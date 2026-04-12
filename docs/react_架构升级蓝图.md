# GIBH-AGENT-V2 语义路由与子智能体架构升级蓝图

**文档性质**：战略级架构设计（Architecture Design Blueprint）  
**版本**：v0.1（草案，供总指挥评审）  
**适用范围**：编排层（Orchestration）全局意图与分发；与《[ARCHITECTURE_CORE_PRINCIPLES.md](./ARCHITECTURE_CORE_PRINCIPLES.md)》协同——**不削弱**资产总线、执行器反射、TaaS 隔离等既有铁律。

---

## 1. 架构演进背景与痛点（Background）

### 1.1 现状概览

当前系统在 `AgentOrchestrator.stream_process` 与 `_classify_global_intent` 中采用 **「多层 if-else + 关键词白名单 + 条件 LLM」** 的混合路由。该模式在快速迭代期具备可实施性，但随着 **Task（组学）**、**HPC（超算 MCP）**、**Chat/ReAct（通用对话与工具）** 三条业务线并行扩张，缺陷已呈结构性暴露。

### 1.2 核心痛点

| 维度 | 问题描述 | 业务影响 |
|------|----------|----------|
| **脆弱的白名单** | 闲聊类关键词（如「你好」）与业务句在同一字符串空间内做子串匹配，易 **误伤或误放**：同一意图的不同措辞可走完全不同的分支。 | 用户体验不可预测；排障需读源码而非读契约。 |
| **意图与文件状态脱节** | 文件嗅探、资产类型与「是否 chat」的规则分散在多处；部分路径依赖 LLM 二次判别，**与首层关键词决策顺序** 交织。 | 「无文件却想做组学」与「有文件却想闲聊」边界模糊；回归测试面指数膨胀。 |
| **首包延迟与感知** | 全局分类前常需 **额外 LLM round-trip**（或多次状态 SSE），且与下游 Planner 串行叠加。 | 用户可见「正在理解您的意图…」停留时间偏长；弱网环境下更明显。 |
| **修改风险与耦合** | 新业务入口往往表现为 **在 Orchestrator 再开一支 if**；与 `target_domain` 快车道、`Skill_Route`、MCP 开关等 **隐式优先级** 共存。 | CR 难以证明「无副作用」；新人上手成本高。 |

### 1.3 战略结论

**废弃「以关键词为中心的脆弱路由树」**，改为 **单一职责的语义路由网关（Semantic Router）** 输出 **结构化 JSON 决策**，再由 **上下文感知的垂直子智能体（Sub-Agents）** 执行。路由层 **不做深度推理、不执行业务工具**；子智能体 **不重复实现资产嗅探与编排契约**（继续委托 `OmicsAssetManager`、既有 Planner/Executor 链）。

---

## 2. 核心架构设计：Semantic Router + Sub-Agents（Core Architecture）

### 2.1 设计原则（与白皮书对齐）

- **资产真相源**：文件与路径状态仍由 **`OmicsAssetManager` + 统一路径规范化** 提供；Router 只消费其摘要，不复制嗅探逻辑。  
- **重型计算不出主进程**：HPC 与组学重型算子仍遵守 **TaaS / Worker / MCP 网关** 边界。  
- **快车道可并存**：`[Skill_Route: …]`、`target_domain` 等可作为 **Router 之前的显式策略层** 或 **Router 输入特征**，但须在文档与代码中 **显式排序**，避免隐式覆盖。

### 2.2 高维度架构流转图（Mermaid）

```mermaid
flowchart TB
    subgraph IN["用户与系统输入"]
        U[("用户自然语言<br/>User_Prompt")]
    end

    subgraph AGG["状态聚合器<br/>Context Aggregator"]
        F[("File_Status<br/>路径列表 / 资产摘要<br/>OmicsAssetManager")]
        M[("MCP_Status<br/>compute_scheduler 等开关")]
        S[("Session_Flags<br/>target_domain / Skill_Route / 历史轮次")]
        AGG_OUT[["聚合上下文<br/>RouterInput JSON"]]
        U --> AGG
        F --> AGG
        M --> AGG
        S --> AGG
        AGG --> AGG_OUT
    end

    subgraph SR["语义路由网关<br/>Semantic Router"]
        direction TB
        R1["仅调用轻量 LLM 或<br/>经校准的小型分类模型"]
        R2["输出标准 JSON<br/>route + confidence + rationale_short"]
        AGG_OUT --> R1 --> R2
    end

    subgraph DEC["决策与策略"]
        D1{"置信度 ≥ τ 且<br/>前提满足?"}
        D2{"最多 1 次重试"}
        D3["Clarify 模式<br/>人机澄清"]
        R2 --> D1
        D1 -->|是| DISPATCH
        D1 -->|否| D2
        D2 -->|仍失败| D3
        D2 -->|重试后成功| DISPATCH
    end

    subgraph DISPATCH["子智能体分发"]
        T["Task Sub-Agent<br/>组学生产线"]
        H["HPC Sub-Agent<br/>MCP 指令链"]
        C["Chat / ReAct Sub-Agent<br/>通用对话与基础工具"]
    end

    DISPATCH --> T
    DISPATCH --> H
    DISPATCH --> C

    subgraph OUT["输出"]
        SSE[("统一 SSE / 状态机<br/>与现有前端契约对齐")]
    end

    T --> SSE
    H --> SSE
    C --> SSE
    D3 --> SSE
```

### 2.3 数据流摘要

1. **输入**：原始 `User_Prompt` + API 层元数据（`enabled_mcps`、`thinking_mode`、`target_domain` 等）。  
2. **聚合**：状态聚合器产出 **稳定 schema** 的 `RouterInput`（见 §3）。  
3. **路由**：Semantic Router **仅** 产出 `RouterDecision` JSON。  
4. **分发**：根据 `route` 字段实例化或委托对应 Sub-Agent；**Clarify** 为独立出口，可挂起会话。  
5. **输出**：各 Sub-Agent 复用现有 SSE 事件类型，避免前端大爆炸式改造。

---

## 3. 语义路由网关详细设计（The Semantic Router）

### 3.1 职责边界（单一职责）

| 做 | 不做 |
|----|------|
| 根据聚合上下文输出 **离散路由标签** + **置信度** + **极短理由**（便于审计与日志） | 不调用组学 Planner、不执行 DAG、不直连 MCP `tools/call` |
| 校验 **路由前提**（例如：`route=task` 且无文件时是否允许「仅规划」） | 不进行多步 ReAct、不生成工作流 JSON |
| 触发 **最多 1 次** 自我纠错式重试（见 §3.5） | 不把「闲聊关键词」与「业务意图」混在同一层子串表中作为主要依据 |

### 3.2 输入层定义（RouterInput）

路由层 **必须** 综合以下字段（缺失项显式为 `null` 或空数组，禁止静默省略）：

| 字段 | 类型（逻辑） | 说明 |
|------|----------------|------|
| `user_prompt` | `string` | 用户当前轮自然语言（已 UTF-8 规范化）。 |
| `file_status` | `object` | **至少包含**：`has_upload: bool`、`resolved_paths: string[]`、`asset_digest: string`（由 `OmicsAssetManager` + 既有 digest 格式化函数产出）。 |
| `mcp_status` | `object` | **至少包含**：`compute_scheduler: bool`、可选 `gateway_connected: bool?`（若健康检查可用）。 |
| `session_flags` | `object` | **可选**：`target_domain`、`skill_route_tool`、`history_turn_count` 等，仅作特征，不替代 Router 决策。 |
| `router_version` | `string` | 便于灰度与回放（如 `semantic-router@2026-04`）。 |

> **说明**：`Skill_Route` 与 `target_domain` 可在 **Router 之前** 短路（保持与现白皮书「技能快车道」一致），或作为 `session_flags` 输入由 Router **显式**输出 `route=skill_fast_lane`（若总指挥选择统一经 Router 审计）。Phase 1 建议 **维持现有短路、Router 不处理 Skill 暗号**，以降低首阶段风险。

### 3.3 输出约束（RouterDecision JSON）

路由网关 **仅允许** 输出下列 **JSON 对象**（单行或可解析块；禁止 Markdown 包裹以外的自然语言）：

```json
{
  "route": "task | hpc | chat | clarify | skill_fast_lane",
  "confidence": 0.0,
  "rationale_short": "string, max 200 chars, zh or en",
  "sub_route": "optional string, e.g. planning_only | execution",
  "requires_clarify": false,
  "clarify_question_id": null
}
```

**字段语义**：

- **`route`**：主分发键；`clarify` 表示不进入任何 Sub-Agent，直接进入人机澄清子状态机。  
- **`confidence`**：`[0,1]` 浮点；由模型 logits 校准或独立评分头给出（Phase 1 可用固定阈值 + LLM 自报置信度 **实验性** 实现）。  
- **`rationale_short`**：仅供日志与合规审计，**不展示给最终用户**（或与 Clarify 问题区分）。  
- **`sub_route`**：Task 线内部细分（如仅规划）；HPC/Chat 可留空。  
- **`requires_clarify`**：与 `route=clarify` 可二选一或组合；具体由状态机统一解释。

**阈值建议（可配置）**：

- `τ_high = 0.75`：高于则直接分发。  
- `τ_low = 0.45`：低于则进入 **一次重试** 或 **Clarify**（见 §3.5）。

### 3.4 路由与三大业务映射（逻辑表）

| `route` | 目标 Sub-Agent | 典型触发语义（示例，非实现硬编码） |
|---------|----------------|--------------------------------------|
| `task` | Task Agent | 多步骤组学分析、工作流执行、含明确 modality 与数据动作且 `file_status` 支持或使用户接受「仅规划」。 |
| `hpc` | HPC Agent | 超算连接、作业提交/查询、集群状态、`execute_hpc_command` 类意图，且 `mcp_status.compute_scheduler == true`。 |
| `chat` | Chat/ReAct Agent | 通用问答、百科、与组学/HPC 无强绑定的工具调用；或未开 HPC 时的兜底对话。 |
| `clarify` | Clarify 子状态 | 前提缺失、意图冲突、置信度不足。 |
| `skill_fast_lane` | （可选）现有 SkillAgent | 若总指挥决定纳入 Router 审计则保留；否则 Phase 1 在 Router 前短路。 |

### 3.5 兜底机制（Fail-safe）

1. **置信度不足**：若 `confidence < τ_high` 且 `≥ τ_low`，允许 **使用「收紧版 system + 更短上下文」重试 1 次**；**禁止**无限重试。  
2. **前提缺失**：例如 `route=task` 且 `has_upload==false` 且用户未选择「仅规划」——**不重试**，直接 `route=clarify` 并生成结构化追问（选项按钮可由现有 `ask_human_for_clarification` 或等价通道承载）。  
3. **HPC 关闭**：若模型倾向 `hpc` 但 `compute_scheduler==false`，**降级为 `chat`** 或 `clarify`（由产品策略二选一，须在配置中写死）。  
4. **解析失败**：JSON 非法 → 计为 **一次失败**；若已消耗重试次数 → **`clarify`**，内容为「系统无法稳定理解请求，请用一句话说明要做：组学分析 / 超算操作 / 普通问答」。  
5. **审计**：每次 Router 调用落结构化日志（输入摘要 hash + 输出 JSON + 是否重试），便于回放「同意图不同措辞」类问题。

---

## 4. 垂直子智能体定义（Sub-Agents Definition）

### 4.1 Task Agent（组学生产线）

- **职责**：沿用现有 **文件检查 → SOPPlanner → DAG/Executor** 链路；**不**在 Task 内嵌套 HPC MCP 或通用闲聊长链。  
- **隔离原则**：Task 入口仅依赖 **RouterDecision.route == task** 与既有 `OmicsAssetManager` 摘要；**禁止** Task 再次解析「是否 chat」的全局意图。  
- **减负目标**：移除/旁路 Task 路径上对「闲聊关键词」「天气」等的重复判断；统一上收至 Router。  
- **稳定性**：作为核心演示与交付链路，**Phase 1 不得破坏默认行为**；仅将「从 Orchestrator 顶部坠入 Task」改为「由 Router 显式派发」。

### 4.2 HPC Agent（超算 MCP）

- **职责**：在 `compute_scheduler` 可用时，将自然语言 **映射为 MCP 工具调用序列**（ReAct / DeepReAct 或经统一 Runner）；严格经 **mcp-gateway**。  
- **边界**：不承载组学 DAG；不替代 `OmicsAssetManager`。  
- **与白皮书关系**：重型命令执行仍在网关/上游；主进程仅编排与参数组装。

### 4.3 Chat / ReAct Agent（通用对话与基础工具）

- **职责**：闲聊、文档型问答、Web/NCBI 等 **非组学重 DAG** 工具；可选 ReAct 循环。  
- **边界**：当 Router 给出 `chat` 时进入；**不**自动升格为 `task`，除非用户在新一轮明确表达且经 Router 重新分类。

### 4.4 Clarify（人机澄清）

- **性质**：可为 **子状态机** 或 **轻量子 Agent**；与 Router 平级或作为 Router 的一等出口。  
- **输出**：结构化追问 + 可点击选项（与现有前端能力对齐）。

---

## 5. 第一阶段施工计划（Phase 1 Execution Plan）

### 5.1 目标

在 **不破坏现有 Task 默认路径** 的前提下，引入 **独立的 Semantic Router 模块**（新包或 `gibh_agent/core/semantic_router.py`），以 **Mock / 契约测试** 验证 JSON 路由正确性；**生产开关** 默认关闭（`FEATURE_SEMANTIC_ROUTER=false`），待灰度开启。

### 5.2 步骤分解

| 步骤 | 交付物 | 验收标准 |
|------|--------|----------|
| **P1-0** | 本文档评审定稿 + `RouterInput` / `RouterDecision` JSON Schema（Pydantic 或 JSON Schema 文件） | 总指挥签字或等效里程碑。 |
| **P1-1** | `SemanticRouter.decide(input: RouterInput) -> RouterDecision` 纯函数 + **单元测试**（≥30 组用例：含「你好+超算」「无文件组学」「仅 HPC」等） | CI 全绿；**无**对 `stream_process` 的默认调用链改动。 |
| **P1-2** | **Mock LLM**：固定返回合法 JSON，测试分发层「根据 route 选择 handler」的表驱动逻辑 | Mock 下覆盖率 ≥ 关键分支 90%。 |
| **P1-3** | 集成测试夹具：从真实 `Orchestrator` **抽取** `RouterInput` 构造器（复用资产嗅探代码路径） | 与现网嗅探结果 **字节级或摘要级** 一致（对同一 fixture）。 |
| **P1-4** | Feature flag 接入：`stream_process` 在 flag 开启时 **先** 调 Router，再分发；关闭时 **零行为差分**（回归快照测试） | 开关 off 时二进制/日志与基线一致。 |
| **P1-5** | 观测性：结构化日志 +（可选）OpenTelemetry span `semantic_router.decide` | 可在测试环境回放一次请求全链路。 |

### 5.3 明确不在 Phase 1 做的事（防范围蔓延）

- 不删除现有 `_classify_global_intent` 实现；**并行共存**，由 flag 选择。  
- 不改 `Executor` 反射规则、不改 `OmicsAssetManager` 核心分类算法。  
- 不大改前端 SSE 协议（仅可增加可选字段，且须向后兼容）。

### 5.4 风险与回滚

- **风险**：Router LLM 延迟或供应商异常。  
- **缓解**：超时后 **降级** 至当前 `_classify_global_intent` 或默认 `chat`（可配置）。  
- **回滚**：关闭 feature flag 即恢复现网行为。

---

## 6. 附录：与现有文档索引

| 文档 | 关系 |
|------|------|
| [ARCHITECTURE_CORE_PRINCIPLES.md](./ARCHITECTURE_CORE_PRINCIPLES.md) | 铁律不变；Router 为编排层增量。 |
| [HPC_MCP_STATUS_BRIEF.md](./HPC_MCP_STATUS_BRIEF.md) | HPC 网关与工具清单运维事实源。 |
| [工具库.md](../工具库.md) | Registry + MCP 附录；Router **不替代**工具注册表。 |

---

**文档维护**：架构组应在每次 Router 阈值或 `route` 枚举变更时更新本文 §3.3–§3.5，并同步 CI fixture 版本号。

*— 全文完 —*
