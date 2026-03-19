# 严重 BUG 汇报：分析任务路由下「思考过程」(thought) 完全丢失

**汇报日期**: 按架构师三步排查要求完成  
**结论**: 根因在后端——分析任务路径下从未产生并下发 `thought` 事件；前端与 DOM 骨架未见导致丢失或覆盖的代码。

---

## 第一步：后端分析路由的 SSE 推送 (Backend Generator Check)

### 检索结果

- **API 入口**: `server.py` 第 2081 行 `@app.post("/api/chat")`，流式生成器为 `generate_sse()`（约 2245–2290 行），内部直接 `async for event in orchestrator.stream_process(...)` 并 `yield event`，即所有 SSE 事件均来自 `orchestrator.stream_process`。
- **闲聊路由**: `gibh_agent/core/orchestrator.py` 第 299–360 行。当 `intent_type == "chat"` 时，调用：
  ```python
  from .stream_utils import stream_from_llm_chunks
  async for event_type, data in stream_from_llm_chunks(
      llm_client.astream(messages, ...),
      model_name=model_name,
  ):
      yield self._emit_sse(state_snapshot, event_type, data)
  ```
  `stream_from_llm_chunks`（`stream_utils.py` 54–146 行）会将 `reasoning_content` 或 content 中 `<think>...</think>` 内容 yield 为 `("thought", {"content": ...})`，再经 `_emit_sse` 发成 SSE 的 `event: thought`，故闲聊时前端能收到 thought。

- **分析任务路由**: 
  - **直接执行路径**（带 `workflow_data`）：orchestrator 第 381–913 行。仅：初始化执行引擎 → 按步骤 `WorkflowExecutor.execute_workflow()` → 根据结果 yield `status`、`step_result`、`diagnosis`、`result`、`done`。**全程没有调用任何 LLM 流（astream）**，因此**不可能产生 thought**。
  - **规划路径**（无 `workflow_data`）：第 996 行起，经过查询重写、文件规范化、意图分析、Branch A（无文件模板）/ Path A（有文件正式规划）。意图与规划均依赖：
    - `planner._classify_intent` / `planner._analyze_user_intent` / `planner.generate_plan`（`planner.py`）：全部使用 `llm_client.achat(...)`（**非流式**），无 `astream`，无 `stream_from_llm_chunks`。
    - `agentic` 模块（QueryRewriter、Clarifier、Reflector）：同样全部 `achat`，无流式。
  - 因此**在分析任务路径下，后端从未把任何 LLM 输出通过 `stream_from_llm_chunks` 或等价逻辑转成 `event_type == "thought"` 并 yield**。

### 根因结论（第一步）

- **不是**「思考过程被丢弃或拦截」——而是**分析任务分支里根本没有「会 yield thought 的 LLM 流」**。
- 只有 `intent_type == "chat"` 时才会走 `stream_from_llm_chunks(llm_client.astream(...))`；一旦进入任务模式（含直接执行与规划），所有 LLM 调用均为非流式 `achat`，自然没有 thought 事件可推。

### 涉及文件与代码位置

| 文件 | 位置 | 说明 |
|------|------|------|
| `gibh_agent/core/orchestrator.py` | 299–360 行 | 仅 chat 分支调用 `stream_from_llm_chunks` 并 `_emit_sse(..., event_type, data)`，含 thought |
| `gibh_agent/core/orchestrator.py` | 381–913 行 | 直接执行路径：无 LLM 流，仅 status/step_result/diagnosis/result/done |
| `gibh_agent/core/orchestrator.py` | 996 行起 | 规划路径：意图/规划全部走 planner 与 agentic 的 achat |
| `gibh_agent/core/planner.py` | 151, 939, 1284 行等 | 全部 `llm_client.achat`，无 astream |
| `gibh_agent/core/agentic.py` | 82, 210, 288, 429 行 | 全部 `llm_client.achat`，无 astream |
| `gibh_agent/core/stream_utils.py` | 54–146 行 | 唯一将 LLM 流解析为 thought 的生成器，仅在 chat 分支被使用 |

---

## 第二步：前端 SSE 事件分发器 (Frontend Dispatcher Check)

### 检索结果

- **SSE 消费与分发**: `services/nginx/html/index.html` 内，SSE 流在检测到 `text/event-stream` 后进入统一解析逻辑，按 `event` 类型（如 `thought`、`status`、`message` 等）调用 `handleServerEvent(eventType, data)`（约 5587 行起）。
- **thought 处理**: `handleServerEvent` 的 `switch` 中有 `case 'thought':`（约 5615–5644 行）：
  - 从 `data` 中取 `content`/`text`/`reasoning` 得到 `thoughtText`；
  - 若 `window.currentThoughtContainer` 为空则调用 `ensureThoughtBlock()`；
  - 使用 `el.appendChild(document.createTextNode(thoughtText))` 增量追加，**无 innerHTML 覆盖**；
  - 未发现任何 `if (isTask)` / `if (isAnalysis)` 等分支对 `eventType === 'thought'` 做 return 或忽略。

### 结论（第二步）

- 前端**没有**在「分析任务模式」下故意忽略或丢弃 `thought` 事件。
- 分析任务下控制台看不到 thought，是因为**后端在分析路径从未下发 `event: thought`**，而非前端未处理。

---

## 第三步：DOM 暴力覆盖 (DOM Overwrite Check)

### 检索结果

- **AI 气泡骨架**: 
  - 发送消息时通过 `appendMessage('ai', ...)` 创建 AI 行（约 7896–7968 行）。AI 气泡内部为**固定四区**（约 7941–7946 行）：
    - `reasoning-zone reasoning-container`（思考区）
    - `process-log-zone log-container`（执行记录）
    - `text-zone ... text-content-zone`（正文）
    - `cards-zone workflow-container`（卡片）
  - 流式响应开始时（约 5124–5177 行），从当前 AI 行的 `.bubble` 内解析出 `reasoningZone`、`logZone`、`textZone`、`cardZone` 并写入 `window.currentMessageZones`，**没有**用 innerHTML 重写整个 bubble。
- **ProcessLog 渲染**: 约 5278–5306 行。创建 `processLogContainer` 后使用 `logZoneRef.appendChild(processLogContainer)` 追加到 `process-log-zone`，**没有**对 `message-content` 或 `bubble` 做 innerHTML 覆盖。
- **清空逻辑**: 约 5357–5359 行仅清空「文本占位区」：`if (aiContentDiv && !aiContentDiv.closest('.reasoning-zone')) aiContentDiv.innerHTML = ''`，即只清 text-zone，**不动** reasoning-zone。
- **done/error 等收尾**: 使用 `aiContentDiv.innerHTML = ...` 的几处（如 7331、7577、7769、7848）中，`aiContentDiv` 均指向 **text-zone**（约 5169 行 `const aiContentDiv = textZone`），因此只更新正文区，**不会**覆盖 reasoning-zone 或 process-log-zone。

### 结论（第三步）

- **未发现**在分析任务下用 innerHTML 覆盖整个消息容器或整个 bubble 的代码。
- 执行记录通过 appendChild 挂到 `process-log-zone`，与 reasoning-zone、text-zone、cards-zone 隔离；若后端在分析路径也下发 thought，现有 DOM 骨架可以正确容纳「reasoning-zone 在上、process-log-zone 紧随其后」的规范。

---

## 修复思路（待架构师批准后实施）

### 方向 A：在分析任务路径中引入「会 yield thought 的 LLM 流」（推荐）

- **意图/规划前增加「思考阶段」**：在进入文件检查与规划之前，对当前用户 query 做一次**流式** LLM 调用（例如「理解需求并简要说明将如何分析」），通过现有的 `stream_from_llm_chunks(llm_client.astream(...))` 消费，并将得到的 `event_type, data` 用 `_emit_sse(state_snapshot, event_type, data)` 原样转发。这样分析任务也会收到与闲聊一致的 thought 事件，前端无需改逻辑。
- **可选细化**：若希望 thought 仅在「规划前」出现一次，可在 `stream_process` 中在分支「TASK MODE + 无 workflow_data」且尚未调用 planner 的位置，插入上述流式调用并转发；若希望执行阶段也有「正在分析结果」等思考，可在生成 diagnosis 前再增加一段流式 LLM 并同样经 `stream_from_llm_chunks` 转发 thought。

### 方向 B：规划/意图的 achat 结果中解析 <think> 并补发单次 thought

- 在 `planner.generate_plan` 或意图分类返回后，若 LLM 返回内容包含 `<think>...</think>` 或 `reasoning_content`，在 orchestrator 中解析并 **yield 一次** `_emit_sse(state_snapshot, "thought", {"content": ...})`。这样分析路由至少能展示「规划时的推理过程」一段，且不需要把规划改为流式，改动面较小，但 thought 只有一段、无打字机效果。

### 方向 C：前端防御性保留（可选）

- 保持现有 DOM 骨架与事件分发逻辑不变；若后续在分析路径增加 thought 下发，前端已具备正确渲染能力。可在 `ensureThoughtBlock` 中再确认：在分析任务下挂载的 bubble 同样具备 `.reasoning-zone`，且不会被 process-log 的 appendChild 影响。

---

## 汇报小结

| 步骤 | 结论 |
|------|------|
| 第一步 后端 | **根因**：分析任务路径（直接执行 + 规划）从未调用会 yield `thought` 的 LLM 流；仅闲聊分支使用 `stream_from_llm_chunks`，故只有闲聊有 thought。 |
| 第二步 前端 | 前端对 `type === 'thought'` 有完整处理，无按任务模式忽略；看不到 thought 是因为后端未发送。 |
| 第三步 DOM | 未发现用 innerHTML 覆盖整泡或 reasoning-zone；ProcessLog 使用 appendChild 追加到 logZone，符合 reasoning-zone → process-log-zone → text-zone → cards-zone 的规范。 |

**建议**：优先按「方向 A」在分析任务分支增加至少一段流式 LLM 并统一经 `stream_from_llm_chunks` 转发 thought，保证与闲聊路由一致、且不破坏现有 DOM 骨架。若需先小范围验证，可用「方向 B」在规划结果中解析 <think> 并补发一次 thought。
