# 代谢组分析全流程：智能体与 LLM 信息收发总结

本文档描述当用户输入 **「我要做代谢组分析全流程」** 并上传 **CSV 文件** 时，智能体后台在**规划**与**执行**过程中，**发送给 LLM 与从 LLM 接收**的信息路径。**首要用途**是：在真实业务场景里，把「用户体感慢」**拆成可观测阶段**，判断瓶颈是否在 LLM、在**哪一次** LLM、以及是**提示词/输出长度**还是**网络与模型**问题。

---

## 〇、排查思路：如何把「慢」归因到 LLM

以下是在生产/联调中排「智能体慢」时的**思考顺序**，与后文各节一一对应；新版架构去掉 `generate_plan` 内重复意图 LLM 后，**总调用次数下降，但单次耗时的主导因素仍在**，不可仅凭「少打了两次」就预期体感线性变快。

### 0.1 延迟从哪里来（心智模型）

1. **串行水线**  
   任务模式下，多数 LLM 调用**前后存在依赖**（例如要先知 domain/steps 再规划，要先跑完工具才有总结），整体耗时近似为各次 LLM **墙钟时间之和**（外加文件 IO、工具计算）。因此：**即使单次都不夸张，5～7 次串行仍可能累加成「十几秒到数十秒」**，这是正常现象，需用日志区分「累加」与「某一次异常拖尾」。

2. **单次调用的三个旋钮**  
   - **输入规模**：长 System、整段 `file_metadata` JSON、多步 `steps_results` 摘要 → 首包延迟（TTFT）与总时延上升。  
   - **输出上限**：`max_tokens` 大（如诊断 1500、总结 2500）→ 生成阶段更长，且与模型速度强相关。  
   - **模型与链路**：同一提示在不同模型、不同 `base_url`、不同并发负载下差异极大；**应先在 `llm_client` 侧打「单次 achat 耗时」再谈 prompt 优化**。

3. **「规划慢」在新版里通常指什么**  
   去掉 `generate_plan` 内重复意图后，**规划前段的主要 LLM 仍往往是 Orchestrator 侧的 `_classify_intent` + `_analyze_user_intent` 各 1 次**（长 System、512 max_tokens）。若在日志里仍看到这两段各出现 **2 次**，再怀疑 `generate_plan` 漏传 `domain_name`/`target_steps`（见 §6、§4.2）。

### 0.2 热点分级（优先看哪里）

| 优先级 | 阶段 | 为何容易慢 |
|--------|------|------------|
| **高** | 数据诊断（§5） | 统计 JSON + 领域指令输入长，`max_tokens=1500`，输出 Markdown 报告。 |
| **高** | 分析总结（§8） | 聚合多步结果，要求长文解读，`max_tokens` 高且可能**重试 1 次** → 最坏接近 2 倍。 |
| **中** | 域名意图 + 步骤意图（§2、§3） | 各 1 次、**通常串行**；System 规则长，User 带 metadata/步骤表。 |
| **低** | 全局 Chat/Task（§1） | 多数场景 0 次或 1 次短调用；仅在走 LLM 路由或混合资产时值得关注。 |

**结论**：用户说「跑到一半卡很久」，在时间轴上多数对应 **诊断** 或 **总结**；说「刚开始就慢」，优先看 **两次意图 LLM** 与是否误进 **全局 LLM 路由**。

### 0.3 推荐归因步骤（可操作）

1. 在 **`llm_client.py` 的 `achat`/`astream` 入口**记录：`model`、消息条数、可选 token 估计、**起止耗时**（或对接口返回的 usage 打日志）。  
2. 用 **§六 Profiler 关键字** 把墙钟时间粗对齐到阶段；再与 `llm_client` 细日志交叉验证「这一段是否真有一次 achat」。  
3. **同文件第二次请求**：若诊断明显变快，说明命中缓存，瓶颈曾主要在诊断 LLM。  
4. 若 **工具阶段很长** 但 LLM 日志稀疏：慢可能不在 LLM（大表 pandas、绘图、远端 Worker）；本文档不覆盖工具内算子，但排障时**避免默认怪 LLM**。  
5. 若走了 **§三 MetabolomicsAgent 内额外 `achat`** 或 **旁路（闲聊/技能/HPC）**：§4.1 的「5～7 次」会低估，需按实际日志计数。

### 0.4 优化判断（与排障的关系）

- **先度量再改 prompt**：无耗时日志时缩短 System 可能得不偿失。  
- **重复意图已治理**：当前主路径的边际收益更多来自 **诊断/总结** 的模型选择、缓存策略、或「是否必须全文生成」。  
- **重试策略**：总结过短重试会直接放大尾延迟，应在日志里统计重试率后再决定是否改阈值或 `max_tokens`。

---

## 一、整体流程概览

```
用户输入 + 上传 CSV
    ↓
① 全局意图 Chat vs Task（多数「工作流原生资产」可走规则短路为 task，不调用 LLM；混合资产等才走 LLM）
    ↓
② 意图阶段：必要时文件预检 (无 LLM) → 域名意图 _classify_intent (LLM)
           → 用户步骤意图 _analyze_user_intent (LLM)
           （快车道 target_domain：跳过 _classify_intent，仍跑 _analyze_user_intent）
    ↓
③ 分支：无文件 → Branch A 预览（generate_plan 带 domain_name/target_steps，内部不再重复意图 LLM）
         有文件 → Path A：正式体检 inspect_file → 数据诊断 (LLM) → generate_plan
    ↓
④ generate_plan：若 Orchestrator 已传入 domain_name、target_steps（及 steps_to_skip），
   则内部跳过 _classify_intent / _analyze_user_intent（不再与②重复）
    ↓
⑤ 生成工作流配置 (无 LLM，代码 DAG / 模板填充)
    ↓
⑥ 执行器按步骤执行工具（编排器本身不调 LLM；个别工具若内嵌 LLM/MCP 另计）
    ↓
⑦ 生成分析总结报告 (LLM，可能重试一次)
```

**旁路（不在本代谢主链路内，但会占用 LLM 或其它后端）**  
- **闲聊模式**：`intent_type == "chat"` 时走 `_stream_chat_mode`（流式 LLM，可挂 MCP tools）。  
- **技能快车道**：命中技能广场等路径时可能走 `SkillAgent.execute_skill`（独立 LLM+工具链）。  
- **compute_scheduler MCP**：`gibh_agent/api/routers/chat.should_use_hpc_isolated_route` 为真时，SSE 走 `hpc_orchestrator`（经 MCP Gateway，不经 `AgentOrchestrator` 本流程）。

---

## 二、各阶段 LLM 调用明细

### 1. 全局意图分类（Chat vs Task）

| 项目 | 说明 |
|------|------|
| **位置** | `gibh_agent/core/orchestrator.py` → `_classify_global_intent()`；入口在 `stream_process` 中「Layer 0」全局路由 |
| **触发** | 闲聊关键词 / BepiPred 等命中可直接 `chat`；`target_domain` 快车道时**整段跳过**本方法，直接 `task` |
| **非 LLM 短路** | 上传资产经 `OmicsAssetManager` 嗅探后，若仅为工作流原生类型（如代谢表格等）→ 直接 `task`；若为非工作流资产且查询模糊 → 直接 `chat`（详见代码内 `ROUTING_*` 分支） |
| **走 LLM 时** | **System**：路由器说明（chat=对话/技能工具链，task=多步骤组学工作流等）。**User**：用户输入 + 「上传资产嗅探摘要」。期望 JSON：`{"type": "chat"}` 或 `{"type": "task"}` |
| **参数** | `temperature=0.1`, `max_tokens=80`（以当前 `llm_client.achat` 调用为准） |
| **耗时特点** | 命中短路时 0 次 LLM；调用时请求仍较短 |

---

### 2. 意图分类（域名 + 模式 + 目标步骤）— 第一次

| 项目 | 说明 |
|------|------|
| **位置** | `gibh_agent/core/planner.py` → `SOPPlanner._classify_intent()` |
| **触发** | 非 `target_domain` 快车道：在文件分支前完成预检后 `planner._classify_intent(refined_query, file_metadata_for_intent)`（见 `orchestrator.py` 意图分析段）；快车道则映射域名后**不调用**本 LLM |
| **发送给 LLM** | **System**: 长段英文规则（Intent Classifier），说明：<br>• 域名：Metabolomics / RNA / Spatial / Radiomics<br>• 模式：EXECUTION / PLANNING<br>• 文件上下文与路由规则（CSV→Metabolomics 等）<br>**User**: `**User Query:**` 用户查询 + `**File Context:** File Uploaded: True` + `**File Metadata:**` 的 JSON（含 file_path、file_type、n_samples、n_features 等） |
| **接收** | JSON：`{"domain_name": "Metabolomics", "mode": "EXECUTION", "target_steps": [...]}` |
| **参数** | `temperature=0.1`, `max_tokens=512` |
| **耗时特点** | System 文本较长，User 带完整 file_metadata，单次调用可能较慢 |

---

### 3. 用户意图分析（选择目标步骤）— 第一次

| 项目 | 说明 |
|------|------|
| **位置** | `gibh_agent/core/planner.py` → `SOPPlanner._analyze_user_intent()` |
| **触发** | Orchestrator 对 `planner._analyze_user_intent` 使用 `async for` 消费（thought / intent_result 等事件）；快车道与传统路均会执行（快车道跳过域名 LLM 但不跳过步骤分析） |
| **发送给 LLM** | **System**: 说明当前域（Metabolomics）及可用步骤 ID 列表（inspect_data, preprocess_data, pca_analysis, ...），要求从用户查询中推断目标步骤与跳过步骤。<br>**User**: 用户查询 + 可用步骤 JSON |
| **接收** | JSON：`{"target_steps": [...], "skip_steps": [...]}` |
| **参数** | `temperature=0.1`, `max_tokens=512` |
| **耗时特点** | 与步骤 2 类似，依赖同一段用户查询和上下文 |

---

### 4. 文件检查（无 LLM）

| 项目 | 说明 |
|------|------|
| **位置** | `gibh_agent/core/file_inspector.py` → `FileInspector.inspect_file()` |
| **说明** | 纯本地逻辑：读取 CSV、解析表头、统计行列数等，**不调用 LLM**。可能被调用两次：一次为意图阶段提供 `file_metadata_for_intent`，一次在 Path A 为诊断与规划提供 `file_metadata`。 |

---

### 5. 数据诊断（生成诊断报告）

| 项目 | 说明 |
|------|------|
| **位置** | `gibh_agent/agents/base_agent.py` → `_perform_data_diagnosis()` → `DataDiagnostician.analyze()` 内部 LLM 调用 |
| **触发** | Path A（有文件、执行分支）中文件检查通过后，若无缓存则调用数据诊断（具体行号随版本变动，搜索 `_perform_data_diagnosis` / `diagnosis` SSE） |
| **发送给 LLM** | **System**: 领域指令（如代谢组学专用 METABO_INSTRUCTION）+ 统计事实（n_samples, n_features, 缺失率等 JSON）。<br>**User**: 由 DataDiagnostician 组装的提示，要求根据统计数据生成「数据诊断报告」和参数建议。 |
| **接收** | 一段 Markdown 格式的诊断报告文本；内部还会解析参数推荐表。 |
| **参数** | `temperature=0.3`, `max_tokens=1500` |
| **耗时特点** | 输入中含完整统计 JSON，输出较长，**单次调用耗时可能明显**。若存在诊断缓存则跳过此次 LLM。 |

---

### 6. 规划器 generate_plan 内部（意图 LLM 是否重复）

| 项目 | 说明 |
|------|------|
| **位置** | `gibh_agent/core/planner.py` → `SOPPlanner.generate_plan()` |
| **触发** | Branch A（无文件预览）与 Path A（有文件执行）均 `async for` 消费 `generate_plan(...)`。 |
| **显式注入（当前实现）** | Orchestrator 传入 `domain_name=`、`target_steps=`（及 `steps_to_skip=`、`target_domain=`、`is_template=` 等）。当 `domain_name` 已给定时 **跳过** `_classify_intent`；当 `target_steps` 已给定时 **跳过** `_analyze_user_intent`，并打日志「避免重复 LLM」。 |
| **仍会出现内部意图 LLM 的情况** | 直接调用 `generate_plan` 且**未**传入 `domain_name` 或 `target_steps` 时，Planner 内部仍会补跑对应 LLM（兼容旧调用方）。 |
| **后续** | 依赖解析、模板生成（`workflow.generate_template()`）、占位符/参数填充等均为代码逻辑，**不调用 LLM**（除非未来扩展）。 |

---

### 7. 执行阶段（无 LLM）

| 项目 | 说明 |
|------|------|
| **位置** | `gibh_agent/core/executor.py` → `WorkflowExecutor.execute_workflow()` → 每步 `execute_step()` |
| **说明** | 按工作流配置依次调用 ToolRegistry 中的工具（如 `metabo_data_validation`, `inspect_data`, `preprocess_data`, `pca_analysis` 等），**不在此阶段调用 LLM**。工具为 Python 函数（pandas、统计、绘图等）。 |

---

### 8. 分析总结报告（执行完成后）

| 项目 | 说明 |
|------|------|
| **位置** | `gibh_agent/agents/base_agent.py` → `_generate_analysis_summary()`；由 Orchestrator 在步骤执行完成后调用（行号随版本变动，可搜索 `_generate_analysis_summary`）。 |
| **触发** | `summary = await target_agent._generate_analysis_summary(steps_results=..., ...)` |
| **发送给 LLM** | **System**: 角色设定（资深生物信息/代谢组解读）。<br>**User**: 各步骤结果摘要（成功/警告/失败）、关键数值与图表路径等，要求生成「深度生物学解释」报告（定量结果、生物学机制、潜在标志物、下一步建议等）。 |
| **接收** | 长段 Markdown 分析报告；若模型返回过短可能触发一次重试。 |
| **参数** | 首次：`temperature=0.3`, `max_tokens=2500`；重试：`max_tokens=2000` |
| **耗时特点** | 输入包含多步骤摘要，输出长，**单次或两次调用都可能较慢**。 |

---

## 三、代谢组 Agent 内部其他 LLM 调用（若被触发）

在数据诊断、参数抽取等路径中，`MetabolomicsAgent` 还可能涉及：

| 场景 | 位置 | 用途 | 参数 |
|------|------|------|------|
| 步骤/参数解析 | `metabolomics_agent.py` 中 `achat` 约 319、965、1058、1251 行附近 | 意图检测、结束步骤、工作流参数等 | temperature 0.1–0.2, max_tokens 64–800 |
| 文件解释 | 约 390 行 | 生成文件说明文本 | temperature=0.3, max_tokens=800 |
| 诊断相关 | 约 1437 行 | 生成/润色诊断类短文本 | temperature=0.2, max_tokens=500 |

这些在「我要做代谢组分析全流程」+ CSV 的**主路径**上不一定都会触发，但与诊断、参数推荐相关的分支可能增加额外 LLM 次数。

---

## 四、LLM 调用次数与耗时排查要点

> **与 §〇 的关系**：§4.1 回答「大概打几次」；**真正把慢落到具体阶段**时，请结合 **§〇 的心智模型、热点分级与归因步骤**。次数变少不等于单次最贵的调用（诊断/总结）消失。

### 4.1 预估调用次数（主路径：代谢 CSV + 任务模式 + Orchestrator 已注入 domain/target）

| 阶段 | 调用次数 | 说明 |
|------|----------|------|
| 全局意图 | 0～1 | 工作流原生资产可短路为 task；否则 1 次 LLM（快车道带 `target_domain` 时跳过全局意图 LLM） |
| 意图分类（Orchestrator 侧） | 0～1 | 快车道 0；普通路 1 |
| 用户意图分析（Orchestrator 侧） | 1 | 快车道与普通路通常均有 |
| 数据诊断 | 0～1 | 无文件分支无诊断；有文件 Path A 且无缓存时 1 |
| 意图分类 / 用户意图（generate_plan 内） | **0** | Orchestrator 已传 `domain_name` + `target_steps` 时跳过（当前主路径） |
| 分析总结 | 1～2 | 正常 1 次，过短时重试 1 次 |
| **合计（典型有文件执行）** | **约 5～7 次** | 视全局是否 LLM、是否命中诊断缓存、总结是否重试而定 |

若某调用方**未**向 `generate_plan` 传入 `domain_name`/`target_steps`，Planner 内部仍可能各多 1 次意图 LLM（与早期文档描述的「重复」现象一致）。

### 4.2 排查运行缓慢时的建议（实操清单）

1. **日志中确认每次 LLM 调用**  
   在 `llm_client.py` 的 `achat`/`astream` 入口打日志（如请求前打 `model`、`messages` 长度或 token 数），并在响应返回后打耗时，便于对应到上述阶段。

2. **重点看规划阶段**  
   • 意图分类、用户意图分析：System + User 较长（尤其带完整 file_metadata），单次可能数秒到十数秒。<br>• **重复调用**：若仍看到同一轮中「意图分类」「用户意图分析」各出现两次，优先检查是否仍有代码路径调用 `generate_plan` 时**漏传** `domain_name` / `target_steps`（当前 Orchestrator 主路径已注入）。

3. **优化方向**  
   • 诊断报告：确认是否命中缓存（同一文件第二次可跳过诊断 LLM）。<br>• 分析总结：若模型经常返回过短导致重试，可考虑调整 prompt 或 max_tokens，减少重试率。<br>• 快车道与全局短路：减少不必要的全局意图 LLM 与域名分类 LLM。

4. **确认 API 与网络**  
   所有 LLM 请求均通过 `gibh_agent/core/llm_client.py` 的 `AsyncOpenAI` 发出，需确认 `base_url`、`api_key` 及网络延迟；模型选择（如 DeepSeek-R1、Qwen 等）也会影响单次耗时。

---

## 五、相关代码索引

| 功能 | 文件: 行号或方法 |
|------|------------------|
| 全局意图 | `orchestrator.py`: `_classify_global_intent()` |
| 意图分类 | `planner.py`: `SOPPlanner._classify_intent()` |
| 用户意图分析 | `planner.py`: `SOPPlanner._analyze_user_intent()` |
| 规划入口 | `planner.py`: `SOPPlanner.generate_plan()`（参数 `domain_name` / `target_steps` / `steps_to_skip` / `is_template` / `target_domain`） |
| 聊天模式流式 | `orchestrator.py`: `_stream_chat_mode()` |
| HPC 隔离聊天 | `hpc_orchestrator.py`: `stream_hpc_direct_chat`（经 MCP，非 `llm_client` 直连时亦不经 Orchestrator DAG） |
| 文件检查 | `file_inspector.py`: `FileInspector.inspect_file()` |
| 数据诊断 | `base_agent.py`: `_perform_data_diagnosis()`；`data_diagnostician.py`（统计）+ LLM 报告 |
| 执行器 | `executor.py`: `WorkflowExecutor.execute_workflow()` / `execute_step()` |
| 分析总结 | `base_agent.py`: `_generate_analysis_summary()` |
| LLM 客户端 | `llm_client.py`: `LLMClient.achat()` / `astream()` |

---

## 六、Profiler 日志关键字（便于对号入座）

后台若已开启 Profiler 探针，可搜索以下日志以对应上述阶段：

| 日志关键字 | 对应阶段 |
|------------|----------|
| `[Profiler] stream_process 入口` | 请求进入 |
| `[Profiler] 全局意图分类完成` | ① 全局意图 LLM 结束 |
| `[Profiler] 意图分类(LLM)完成` | ②/⑥ 意图分类 LLM 结束 |
| `[Profiler] 用户意图分析(LLM)完成` | ③/⑥ 用户意图分析 LLM 结束 |
| `[Profiler] 规划开始` / `规划完成` | ④ 规划阶段 |
| `[Profiler] 工具[xxx] 执行耗时` | ⑦ 各步骤工具执行（非 LLM） |
| `[DataDiagnostician] 调用 LLM` / `LLM调用完成` | ⑤ 数据诊断 LLM |
| `[AnalysisSummary] 开始LLM调用` / `LLM调用完成` | ⑧ 分析总结 LLM |

结合 `monitor-lite.sh` 菜单中的「查看后端性能与耗时日志 (Profiler)」可快速定位哪一阶段耗时最长。

---

*文档版本：已与 `generate_plan` 显式注入 `domain_name`/`target_steps` 后的 Orchestrator 主路径对齐；行号请以仓库内搜索为准。*
