# 方法论指导报告：流式输出与思考呈现 · 外部 HTTP 工具接入

**文档性质**：基于 GIBH-AGENT-V2 已落地实现归纳的扩展指南，供后续同类功能复用与评审。  
**适用读者**：后端（编排 / LLM / 工具）、前端（SSE / 思考区 / 执行记录）。

---

## 一、流式输出与思考过程呈现的方法论

### 1.1 设计原则（防回退）

| 原则 | 说明 |
|------|------|
| **双轨分流** | 「推理/思考」与「对用户可见正文」在协议层分开，避免混进同一段 `message` 污染快照与 UI。 |
| **增量优先** | 流式场景下以 chunk 为单位消费，禁止在收到完整响应前用整包 `innerHTML` 覆盖 DOM。 |
| **跨块安全** | 标签（如 `think`）可能被切在两次 chunk 之间，必须用缓冲区 + 状态机，不能只靠单次 `indexOf`。 |
| **快照一致** | `state_snapshot` 中 `reasoning` 与 `text` 的写入规则须与 SSE `event` 类型一一对应，便于历史恢复。 |

### 1.2 后端分层模型

1. **传输层：`LLMClient.astream`**
   - 统一走 OpenAI 兼容 `chat.completions.create(stream=True)`，产出 `ChatCompletionChunk` 异步迭代器。
   - 参考：`gibh_agent/core/llm_client.py` 中 `astream`。

2. **解析层：`stream_from_llm_chunks`（推荐默认路径）**
   - **轨 A**：若 delta 带 `reasoning_content`（部分通道/模型），直接映射为 `("thought", {"content": ...})`。
   - **轨 B**：对 `delta.content` 做缓冲区 + `in_think_tag`，在 `THINK_OPEN` / `THINK_CLOSE`（与 `stream_utils` 常量一致）之间产出 `thought`，之外产出 `message`。
   - 预留对 `<<<SUGGESTIONS>>>` 等扩展块的分拣，避免打断 think 状态机。
   - 参考：`gibh_agent/core/stream_utils.py` 中 `stream_from_llm_chunks` 及注释中的「尾部保留长度」策略。

3. **编排层：Orchestrator / `_emit_sse`**
   - 将解析层产出的 `(event_type, data)` 转为 `event: …\ndata: JSON\n\n`。
   - **硬约束**：`thought` 只累加到 `state_snapshot["reasoning"]`，`message` 只累加到 `state_snapshot["text"]`，禁止交叉写入。
   - 参考：`gibh_agent/core/orchestrator.py` 中 `_emit_sse` 文档字符串与分支实现。

4. **非流式回退：`achat` + `extract_think_and_content`**
   - 一次性响应时，用正则剥离 think 块，再按需包装成带标签的字符串或拆分字段。
   - 参考：`gibh_agent/core/llm_client.py` 中 `extract_think_and_content`；`gibh_agent/agents/base_agent.py` 中 `chat(..., stream=False)` 路径。

### 1.3 专项场景：技能快车道 + 结构化执行日志（可选）

当需要在**不破坏全局 Planner**的前提下，对单次技能做「填参 → 工具 → 结果」的可视化：

- 使用 **`astream` + `stream_from_llm_chunks`** 抽取思考，将「标签外 JSON」解析为工具参数（与强制 function-call 二选一或作回退）。
- 通过 **`status.content` 前缀 `[AgenticLog] ` + JSON** 向前端传递 `step_start` / `substep_delta` 等动作（具体 schema 与前端 `handleAgenticLog` 对齐）。
- 参考：`gibh_agent/agents/skill_agent.py` 中 `execute_skill`、`SkillAgent._agentic_emit`。

### 1.4 前端配合要点

- **SSE**：按 `event` 类型分发；`thought` 使用 **`appendChild(document.createTextNode)`** 或等价方式增量写入思考容器，避免 XSS 与闪烁。
- **执行记录**：普通任务仍走 `updateProcessLog`；若以 `[AgenticLog] ` 开头则走独立解析分支，**不得**把整段 JSON 当普通步骤文案追加。
- 参考：`services/nginx/html/index.html` 中 `handleServerEvent`（`case 'status'`）、`ensureThoughtBlock`；文档 `docs/FRONTEND_THOUGHT_RENDERING_REPLAY.md` 可作补充阅读。

### 1.5 扩展新流式能力时的检查清单

1. 新路由是否 **`async for`** 消费 LLM 流，并在首包/末包打日志便于排障？  
2. 是否复用 **`stream_from_llm_chunks`**，避免复制一套 think 解析？若模型协议不同，是否单独函数而非改坏公共状态机？  
3. SSE 是否明确 **`thought` / `message` / `status`** 语义，且 `_emit_sse` 快照字段是否同步？  
4. 前端是否避免 **`innerHTML` 全量替换** 思考区？  
5. 若需从流中拆 JSON（如 Planner），是否使用 **`stream_and_extract_json`** 或同等「缓冲 + 闭合检测」模式？（见 `stream_utils.py`）

---

## 二、通过外部 API 文档实现工具调用的方法论

托管侧多为 **HTTP JSON API**（示例：BepiPred3 `POST /predict`）。本仓库的成熟范式是：**文档 → 探针脚本 → 注册工具 → 编排/技能消费**。

### 2.1 从文档到契约

1. **固定端点与动词**：记录 base URL、path、方法（多为 `POST`）、`Content-Type: application/json`。
2. **请求体字段**：与文档逐字段对齐；在代码里用 **TypedDict / Pydantic / 显式 `Dict`** 注释清楚含义与枚举（如 `top_20` / `top_50`）。
3. **响应信封**：区分 **HTTP 层**（status code）与 **业务层**（如顶层 `success`、`data`、`data.status`）。
4. **异步与超时**：长计算任务使用 **`httpx.AsyncClient(timeout=较大值)`**（如 300s），与文档建议一致。

### 2.2 推荐实现步骤（可复用流水线）

| 步骤 | 动作 |
|------|------|
| **A. 纯净探针** | 脱离 Agent，用 `requests`/`httpx` 写最小脚本：`POST` 文档示例 payload，**打印完整 JSON**，确认 `success`、`data.status`、`output_files` 等真实形状。 |
| **B. 环境化配置** | BepiPred3 工具**默认本地子进程**（`BEPIPRED3_ROOT` / `BEPIPRED3_PYTHON` 等）；仅当 `BEPIPRED3_USE_REMOTE_API=1` 时使用 `BEPIPRED3_PREDICT_URL` 等远程覆盖。 |
| **C. 防御性输入** | 文档要求与 LLM 填参常见误差要想在前面（如 FASTA 缺 `>` 导致 422 → 服务端自动补 `>Sequence`）。 |
| **D. 错误与重试** | 对 **5xx、超时、连接错误** 做有限次重试+退避；对 **业务失败 `data.status == failed`** 根据 `error_message` 判断是否 GPU 瞬时繁忙可重试。 |
| **E. 结果归一化** | 将远端返回的 `output_files`（dict 或 list）**统一解析**为智能体易用的字段（如 `html_url` / `csv_url` / `zip_url`），供 Markdown 卡片与前端链接使用。 |
| **F. 注册工具** | `@registry.register(...)` + `args_schema`；描述中写清「何时调用」「参数含义」，便于 Planner / Skill 路由填参。 |
| **G. 用户可见文案** | 将「GPU/算力繁忙」等转译为可操作的中文短句，**避免在面向用户的字符串中出现第三方商业产品名**；避免把原始堆栈甩给用户。 |

**参考实现**：`gibh_agent/tools/protein_tools.py`（`bepipred3_prediction`，本地 CLI 为主）；远程磐石 API 纯净探针见 `tests/test_panshi_api_raw.py`。

### 2.3 典型响应模式（以 BepiPred3 为例）

典型处理顺序：

1. `response.status_code == 200` 且 body 为 JSON。  
2. 顶层 **`success is True`**，否则携带 `message`/`error` 返回 `status: error`。  
3. **`data` 为对象**：若 **`data.status == "completed"`**，从 **`data.output_files`** 抽取下载链接；若为 **`failed`**，读取 **`data.error_message`**，按关键字决定是否重试。  
4. 对 **`data.status` 缺失但已有 `output_files`** 的畸形成功响应，可打日志后**降级按成功**（需评审是否允许）。  
5. 若仍为 **`running`**，当前同步工具可返回明确错误提示「任务仍在运行」（若未来支持轮询，再扩展为 job id 模式）。

解析 `output_files` 时建议同时兼容 **字典键**（`html`/`csv`/`zip`）与 **列表项**（`name`/`type`/`download_url`），参见 `_extract_output_urls`。

### 2.4 扩展新 HTTP 工具时的检查清单

1. 是否已有 **`tests/test_*_api_raw.py`** 一类探针，可在无 API Key 的智能体环境外独立验证？  
2. 是否在工具内打 **`response.text` / 解析后 JSON` 的结构化日志**（注意脱敏与长度）？  
3. 是否所有失败路径返回 **`{"status": "error", "message": "..."}`** 形状，与 `SkillAgent._format_result_markdown` 等上游一致？  
4. 成功路径是否提供 **至少一种可点击结果**（URL）或明确说明「无链接」？  
5. 是否需要 **技能广场 seed**、**路由关键词**、**E2E 测试**（参见 `tests/test_bepipred3_skill.py`）？

---

## 三、两条方法论的结合点

- **流式思考**解决「用户感知到模型在推理」；**外部 HTTP 工具**解决「重计算外包给专用服务」。  
- 技能快车道中：先 **流式抽参 +（可选）[AgenticLog] 可视化**，再 **调用已注册的工具函数**，最后 **Markdown 汇总链接**——三层职责清晰，便于单独替换其中一层而不牵动全局。

---

## 四、文档维护

- 当 `stream_utils`、`_emit_sse` 或外部 API 信封变更时，应同步更新本节并补充「变更日期 + 迁移注意」。  
- 新增模型通道（新 `delta` 字段）时，优先在 **`stream_from_llm_chunks`** 扩展轨，避免各业务复制解析逻辑。

---

*本报告由仓库当前实现抽象而成，不替代具体 API 官方文档；接入新托管服务时务必以最新官方说明为准。*
