# 代谢组分析全流程：LLM 调用与耗时归因（精简版）

本文说明用户发起 **代谢组类任务**（含上传表格等）时，系统在**规划—执行—总结**链路上 **何时调用 LLM、调用几次、如何排障**。与 **V2.1** 对齐：编排器已接入 **SemanticRouter**（可关），全局路由与旧版略有差异。

---

## 〇、排障心智（压缩）

1. **串行累加**：任务链上多次 `achat` 墙钟大致相加；体感「十几秒」未必是单次卡死。  
2. **热点**：**数据诊断**、**分析总结**（输入长、`max_tokens` 大，总结可能重试）通常比「意图类单次」更耗时。  
3. **先打日志**：在 `llm_client.py` 的 `achat`/`astream` 记录 **model、耗时**；用 **§五** 关键字对齐阶段。  
4. **工具慢 ≠ LLM 慢**：pandas/绘图/Worker 耗时不计入本文「LLM 次数」。

---

## 一、主链路概览（V2.1）

```
用户输入 + 文件
    ↓
Layer 0 全局路由（见 §1）
    ├ 闲聊关键词 → chat（0 次 LLM）
    ├ target_domain 快车道 → task（0 次 LLM）
    └ 否则
         ├ [默认] SemanticRouter.decide_route → task|chat|hpc|clarify（1～2 次 achat，见下）
         │    clarify → SSE 直接结束，后续 0 LLM
         │    hpc → chat 模式 + 强制 DeepReAct（多轮 achat/工具，次数不定）
         ├ 开关关闭或 SemanticRouter 异常 → _classify_global_intent（0～1 次，含资产嗅探短路）
         └ 无 LLM 客户端 → 回退 _classify_global_intent
    ↓
task：QueryRewriter（0～1 次 achat，澄清合并分支跳过）→ 预检 inspect（无 LLM）
    → 域名意图 _classify_intent（0～1 次 **astream** 流式抽 JSON）→ 步骤意图 _analyze_user_intent（1 次 achat）
    → 无文件：generate_plan（注入 domain/steps 时内部跳过重复意图 LLM）
    → 有文件：inspect → 数据诊断（0～1）→ generate_plan（同上）
    ↓
模板填充（无 LLM）→ 执行器逐步工具（无 LLM，除非工具内嵌）→ 分析总结（1～2）
```

**与旧版差异（全局层）**：未开 `target_domain` 且未命中闲聊关键词时，**默认先走 SemanticRouter**，不再优先执行 `_classify_global_intent` 里「仅工作流原生资产 → 零 LLM 直出 task」的短路。典型「上传代谢表 + 要做全流程」在全局层可能 **多 1 次** 结构化路由 `achat`（`temperature=0.1`，`max_tokens=512`，解析/低置信最多再 **1 次**）。

**旁路（非本代谢任务主链）**  
- **chat**：`_stream_chat_mode` / `_stream_deep_react_chat`（流式 + 可选 tools，次数视轮次而定）。  
- **技能快车道**：`SkillAgent.execute_skill`。  
- **HPC**：除历史「隔离路由」外，SemanticRouter 判 **hpc** 时可在 **本 Orchestrator** 内走 **DeepReAct + compute_scheduler MCP**。

**开关**：环境变量 **`GIBH_ENABLE_SEMANTIC_ROUTER=0|false|no|off`** 关闭语义路由，Layer 0 恢复仅 `_classify_global_intent`。

---

## 二、各阶段 LLM 明细（表格式）

### 1. 全局路由（Layer 0）

| 项目 | 说明 |
|------|------|
| **SemanticRouter（默认）** | `gibh_agent/core/semantic_router.py` → `SemanticRouter.decide_route`；编排器 `stream_process` 在「非关键词、非 target_domain」分支注入 `RouterInput` 并调用。 |
| **调用** | `llm_client.achat`，`temperature=0.1`，`max_tokens=512`；解析失败或 `confidence<0.7` 时 **同轮最多再 1 次**，仍失败则路由为 **clarify**（不再打第三次）。 |
| **输出** | JSON → `route`（task/hpc/chat/clarify/skill_fast_lane）+ `confidence` + `rationale_short`。 |
| **旧路径** | `orchestrator.py` → `_classify_global_intent`：混合资产等 **1 次** `max_tokens=80`；无文件歧义 **1 次** `max_tokens=50`；多条规则短路 **0 次**。 |

### 2. 查询重写（Orchestrator，task 主路径靠前）

| 项目 | 说明 |
|------|------|
| **位置** | `gibh_agent/core/agentic.py` → `QueryRewriter.rewrite`；`orchestrator.stream_process` 在澄清合并分支外调用。 |
| **调用** | `llm_client.achat`，`temperature=0.3`, `max_tokens=256`；失败则回退原始 query，**不计入成功 LLM 轮次**。 |
| **触发** | `self.query_rewriter` 非空且非「澄清续写」合并逻辑时 **1 次**；`query_rewriter` 未初始化则 **0 次**。 |

### 3～4. 域名意图 + 步骤意图（Planner）

| 项目 | 说明 |
|------|------|
| **位置** | `planner.py`：`SOPPlanner._classify_intent`、`_analyze_user_intent` |
| **域名** | **非快车道**：`_classify_intent` 当前实现为 **`llm_client.astream`** + `stream_and_extract_json`（仍计 **1 轮** LLM，非独立 `achat`）。`temperature=0.1`, `max_tokens=512`。 |
| **步骤** | `_analyze_user_intent_core` 内 **1 次** `achat`（`temperature=0.1`, `max_tokens=512`）；JSON 失败走关键词回退，**默认不再二次 achat**。 |
| **快车道** | 有 `target_domain`：**跳过** `_classify_intent`，仍 **`_analyze_user_intent`** 各 **1 轮**。 |

### 5. 文件检查

`FileInspector.inspect_file`：**无 LLM**。

### 6. 数据诊断

`base_agent._perform_data_diagnosis` → DataDiagnostician：**1 次**为主（有缓存可 0）；`temperature`/`max_tokens` 以 `data_diagnostician` 与 agent 实现为准（常见诊断输出较长）。

### 7. generate_plan 内意图

Orchestrator **已传入** `domain_name`、`target_steps` 等时，Planner **跳过**与 §3～4 重复的意图 LLM；未传时内部仍会补跑（兼容旧调用方）。

### 8. 执行器

`WorkflowExecutor`：**无 LLM**（工具内另计）。

### 9. 分析总结

`base_agent._generate_analysis_summary`：**1～2 次**（过短可能重试）。

### 10. MetabolomicsAgent 内其它 `achat`

参数解析、说明文案等分支可能额外调用；**主路径不一定全触发**——以仓库内 **`metabolomics_agent.py` + `achat`** 搜索为准。

---

## 三、调用次数粗算（有文件 + task + generate_plan 已注入 domain/steps）

以下为 **「一轮 LLM 请求」** 计数（`astream` 与 `achat` 均各算 1 次交互，除非失败未发起）。

| 阶段 | 次数 | 备注 |
|------|------|------|
| Layer 0（SemanticRouter 开） | **1～2** | 多数 1 次；低置信/坏 JSON 再 1 次 |
| Layer 0（仅 `_classify_global_intent`） | **0～1** | 原生资产短路 **0**；需 LLM 时单次请求内 **混合资产** 与 **无上传歧义** 两分支 **互斥**，最多 **1 次** `achat`（`max_tokens` 为 80 或 50） |
| QueryRewriter | **0～1** | 编排器已挂载且非澄清续写时通常为 1 |
| 域名意图 `_classify_intent` | **0～1** | 快车道为 0；否则 **astream** 1 轮 |
| 步骤意图 `_analyze_user_intent` | **1** | 快车道与普通路均常见 |
| 数据诊断 | **0～1** | 视缓存 |
| generate_plan 内意图 | **0** | Orchestrator 已传 `domain_name` / `target_steps` 时跳过 |
| 分析总结 | **1～2** | 过短可能重试 |
| **合计（典型：语义路由开 + Rewriter 开 + 有文件执行）** | **约 6～9 次** | 1+1+1+1+1+(0～1 诊断)+1～2 总结；语义路由第二次尝试或总结重试可顶到 **约 10** |
| **Clarifier / Reflector** | **0** | 已在 `AgentOrchestrator` 中构造，**当前 `stream_process` 未调用**，不计入 |

关闭 **`GIBH_ENABLE_SEMANTIC_ROUTER`** 且命中原生资产短路时，全局层可回到 **0**，总次数相应下降。

---

## 四、相关代码索引

| 功能 | 位置 |
|------|------|
| 语义路由 | `gibh_agent/core/semantic_router.py`（`SemanticRouter.decide_route`） |
| Layer 0 集成与开关 | `gibh_agent/core/orchestrator.py`（`ENABLE_SEMANTIC_ROUTER` / `_semantic_router_feature_flag` / `_get_semantic_router_file_status` / `stream_process`） |
| 全局意图（回退） | `orchestrator.py`：`_classify_global_intent` |
| 查询重写 | `gibh_agent/core/agentic.py`：`QueryRewriter.rewrite` |
| Planner 意图 | `planner.py`：`generate_plan`、`_classify_intent`（**astream**）、`_analyze_user_intent`（**achat**） |
| 聊天 / DeepReAct | `orchestrator.py`：`_stream_chat_mode`、`_stream_deep_react_chat` |
| 数据诊断 / 总结 | `base_agent.py`；统计侧见 `data_diagnostician.py` |
| 执行器 | `executor.py`：`WorkflowExecutor` |
| LLM 客户端 | `llm_client.py`：`achat` / `astream` |
| 单测 | `tests/test_semantic_router.py` |

---

## 五、Profiler / 日志关键字

| 关键字 | 含义 |
|--------|------|
| `[Profiler] stream_process 入口` | 请求进入 |
| `🧠 [Semantic Router]` / `[Profiler] 语义路由完成` | SemanticRouter 决策与耗时 |
| `[Profiler] 全局意图分类完成` | `_classify_global_intent` 或 SemanticRouter 异常回退路径 |
| `🔄 [QueryRewriter]` | 查询重写 achat |
| `[Profiler] 意图分类(LLM)完成` / `用户意图分析(LLM)完成` | Planner 意图阶段 |
| `[DataDiagnostician]` / `诊断` 相关 | 数据诊断 LLM |
| `[AnalysisSummary]` | 分析总结 LLM |

---

*文档版本：**V2.1**（已与 `orchestrator.py` / `semantic_router.py` / `planner.py` / `agentic.py` / `base_agent.py` 对照校验；`_classify_intent` 为 **astream** 非 `achat`。）*
