# 代谢组分析全流程：智能体与 LLM 信息收发总结

本文档描述当用户输入 **「我要做代谢组分析全流程」** 并上传 **CSV 文件** 时，智能体后台在**规划**与**执行**过程中，**发送给 LLM 与从 LLM 接收**的信息路径，用于排查运行缓慢原因。

---

## 一、整体流程概览

```
用户输入 + 上传 CSV
    ↓
① 全局意图分类 (LLM)
    ↓
② 文件检查 (无 LLM) + 意图分类 (LLM) + 用户意图分析 (LLM)
    ↓
③ Path A：再次文件检查 → 数据诊断 (LLM) → 规划 generate_plan
    ↓
④ generate_plan 内部：意图分类 (LLM) + 用户意图分析 (LLM)  ← 可能与②重复
    ↓
⑤ 生成工作流配置 (无 LLM，代码 DAG)
    ↓
⑥ 执行器按步骤执行工具 (无 LLM，纯工具调用)
    ↓
⑦ 生成分析总结报告 (LLM，可能重试一次)
```

---

## 二、各阶段 LLM 调用明细

### 1. 全局意图分类（Chat vs Task）

| 项目 | 说明 |
|------|------|
| **位置** | `gibh_agent/core/orchestrator.py` → `_classify_global_intent()` |
| **触发** | `stream_process` 入口，在「工作流复用」判断之后、任务模式分支内 |
| **发送给 LLM** | **System**: `你是一个意图分类助手。只返回JSON格式的意图类型。`<br>**User**: `用户输入: "我要做代谢组分析全流程"` + 判断是「一般聊天」还是「生物信息学分析任务」，只返回 `{"type": "chat"}` 或 `{"type": "task"}` |
| **接收** | 期望 JSON：`{"type": "task"}` |
| **参数** | `temperature=0.1`, `max_tokens=50` |
| **耗时特点** | 请求短、回复短，通常较快 |

---

### 2. 意图分类（域名 + 模式 + 目标步骤）— 第一次

| 项目 | 说明 |
|------|------|
| **位置** | `gibh_agent/core/planner.py` → `SOPPlanner._classify_intent()` |
| **触发** | 有文件时，Orchestrator 在「Path A」之前先做意图分析：`planner._classify_intent(refined_query, file_metadata_for_intent)`（约 1104 行） |
| **发送给 LLM** | **System**: 长段英文规则（Intent Classifier），说明：<br>• 域名：Metabolomics / RNA / Spatial / Radiomics<br>• 模式：EXECUTION / PLANNING<br>• 文件上下文与路由规则（CSV→Metabolomics 等）<br>**User**: `**User Query:**` 用户查询 + `**File Context:** File Uploaded: True` + `**File Metadata:**` 的 JSON（含 file_path、file_type、n_samples、n_features 等） |
| **接收** | JSON：`{"domain_name": "Metabolomics", "mode": "EXECUTION", "target_steps": [...]}` |
| **参数** | `temperature=0.1`, `max_tokens=512` |
| **耗时特点** | System 文本较长，User 带完整 file_metadata，单次调用可能较慢 |

---

### 3. 用户意图分析（选择目标步骤）— 第一次

| 项目 | 说明 |
|------|------|
| **位置** | `gibh_agent/core/planner.py` → `SOPPlanner._analyze_user_intent()` |
| **触发** | 紧接上一步：`planner._analyze_user_intent(refined_query, workflow)`（约 1118 行） |
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
| **触发** | Path A 中文件检查通过后，若无缓存则：`agent_instance._perform_data_diagnosis(file_metadata=..., omics_type="Metabolomics", dataframe=...)`（orchestrator 约 1477 行） |
| **发送给 LLM** | **System**: 领域指令（如代谢组学专用 METABO_INSTRUCTION）+ 统计事实（n_samples, n_features, 缺失率等 JSON）。<br>**User**: 由 DataDiagnostician 组装的提示，要求根据统计数据生成「数据诊断报告」和参数建议。 |
| **接收** | 一段 Markdown 格式的诊断报告文本；内部还会解析参数推荐表。 |
| **参数** | `temperature=0.3`, `max_tokens=1500` |
| **耗时特点** | 输入中含完整统计 JSON，输出较长，**单次调用耗时可能明显**。若存在诊断缓存则跳过此次 LLM。 |

---

### 6. 规划器 generate_plan 内部（可能重复的意图 LLM）

| 项目 | 说明 |
|------|------|
| **位置** | `gibh_agent/core/planner.py` → `SOPPlanner.generate_plan()` |
| **触发** | Orchestrator 在 Path A 中调用：`planner.generate_plan(user_query=refined_query, file_metadata=file_metadata)`（约 1588 行），**未传入已求得的 domain_name / target_steps**。 |
| **内部行为** | • 若 `domain_name` 为空：再次调用 `_classify_intent(user_query, file_metadata)` → **与阶段 2 等价的第二次 LLM**<br>• 若 `target_steps` 为空：再次调用 `_analyze_user_intent(user_query, workflow)` → **与阶段 3 等价的第二次 LLM** |
| **说明** | 当前实现下，Orchestrator 已先做过意图分类与用户意图分析，但未把结果传给 `generate_plan`，导致规划阶段可能**重复两次 LLM 调用**，是排查「规划慢」的重点。 |
| **后续** | 依赖解析、模板生成（`workflow.generate_template()`）等均为代码逻辑，**不调用 LLM**。 |

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
| **位置** | `gibh_agent/agents/base_agent.py` → `_generate_analysis_summary()`；由 Orchestrator 在步骤执行完成后调用（约 676 行）。 |
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
| 步骤/参数解析 | `metabolomics_agent.py` 约 319、962、1055、1248 行 | 从 LLM 返回中解析 JSON（目标步骤、工作流参数等） | temperature 0.1–0.2, max_tokens 64–800 |
| 文件解释 | 约 390 行 | 生成文件说明文本 | temperature=0.3, max_tokens=800 |
| 诊断相关 | 约 1434 行 | 生成/润色诊断类短文本 | temperature=0.2, max_tokens=500 |

这些在「我要做代谢组分析全流程」+ CSV 的**主路径**上不一定都会触发，但与诊断、参数推荐相关的分支可能增加额外 LLM 次数。

---

## 四、LLM 调用次数与耗时排查要点

### 4.1 预估调用次数（当前逻辑）

| 阶段 | 调用次数 | 说明 |
|------|----------|------|
| 全局意图 | 1 | 必现 |
| 意图分类（Orchestrator 侧） | 1 | 有文件时 |
| 用户意图分析（Orchestrator 侧） | 1 | 有文件时 |
| 数据诊断 | 1 | 无缓存时 |
| 意图分类（generate_plan 内） | 1 | domain_name 未传入时 **重复** |
| 用户意图分析（generate_plan 内） | 1 | target_steps 未传入时 **重复** |
| 分析总结 | 1～2 | 正常 1 次，过短时重试 1 次 |
| **合计** | **7～9 次** | 其中 2 次为可优化的重复意图调用 |

### 4.2 排查运行缓慢时的建议

1. **日志中确认每次 LLM 调用**  
   在 `llm_client.py` 的 `achat`/`astream` 入口打日志（如请求前打 `model`、`messages` 长度或 token 数），并在响应返回后打耗时，便于对应到上述阶段。

2. **重点看规划阶段**  
   • 意图分类、用户意图分析：System + User 较长（尤其带完整 file_metadata），单次可能数秒到十数秒。<br>• **重复调用**：若在日志中看到「意图分类」或「用户意图分析」在同一轮请求中出现两次且入参相似，即可确认是 `generate_plan` 未复用 Orchestrator 已有结果导致的重复 LLM。

3. **优化方向（建议）**  
   • 在调用 `planner.generate_plan()` 时传入 Orchestrator 已得到的 `domain_name` 和 `target_steps`，避免在 `generate_plan` 内再次执行 `_classify_intent` 和 `_analyze_user_intent`。<br>• 诊断报告：确认是否命中缓存（同一文件第二次可跳过诊断 LLM）。<br>• 分析总结：若模型经常返回过短导致重试，可考虑调整 prompt 或 max_tokens，减少重试率。

4. **确认 API 与网络**  
   所有 LLM 请求均通过 `gibh_agent/core/llm_client.py` 的 `AsyncOpenAI` 发出，需确认 `base_url`、`api_key` 及网络延迟；模型选择（如 DeepSeek-R1、Qwen 等）也会影响单次耗时。

---

## 五、相关代码索引

| 功能 | 文件: 行号或方法 |
|------|------------------|
| 全局意图 | `orchestrator.py`: `_classify_global_intent()` |
| 意图分类 | `planner.py`: `SOPPlanner._classify_intent()` |
| 用户意图分析 | `planner.py`: `SOPPlanner._analyze_user_intent()` |
| 规划入口 | `planner.py`: `SOPPlanner.generate_plan()` |
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

*文档版本：基于当前代码库梳理，用于排查「代谢组分析全流程」场景下智能体运行缓慢问题。*
