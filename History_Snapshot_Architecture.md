# 历史记录与高维 JSON 快照系统 — 深度技术报告

> 面向首席系统架构师：数据持久化、API、前端注水与跨设备/跨用户白屏根因分析。

---

## 1. 数据持久化层 (Database Schema & ORM)

### 1.1 核心表与 SQLAlchemy 模型

**文件路径**: `gibh_agent/db/models.py`

| 表名 | 模型类 | 说明 |
|------|--------|------|
| `users` | `User` | 用户表，与历史隔离无直接 FK |
| `sessions` | `Session` | 历史会话，`owner_id` 为 String（username / guest_uuid），无 FK |
| `messages` | `Message` | 消息记录，**高维结构存于 `content`** |
| `assets` | `Asset` | 上传文件元数据 |
| `workflow_templates` | `WorkflowTemplate` | 工作流收藏，`config_json` 为完整配置 |

**Message 表（历史快照核心）**:

```python
class Message(Base):
    __tablename__ = "messages"
    id = Column(Integer, primary_key=True, autoincrement=True)
    session_id = Column(String(64), nullable=False, index=True)
    role = Column(String(32), nullable=False)   # user | agent
    content = Column(JSON, nullable=True)        # 复杂嵌套 JSON，无损存取
    created_at = Column(DateTime, nullable=False, default=datetime.utcnow)
```

- **字段类型**: `content` 使用 SQLAlchemy `JSON`，在 MySQL 侧对应 **JSON 类型**（或部分驱动/版本下为 LONGTEXT + 应用层序列化）。
- **可空性**: `content` 允许 `nullable=True`，即 DB 中可出现 `content = NULL` 的记录。

**WorkflowTemplate 表**:

```python
config_json = Column(JSON, nullable=True)  # 工作流参数配置，复杂嵌套无损
```

- 工作流卡片、步骤、参数等全部塞入单列 JSON，无独立步骤表。

### 1.2 高维结构的序列化方式

| 业务概念 | 存储位置 | 结构概要 | 序列化方式 |
|----------|----------|----------|------------|
| 工作流卡片 (Workflow) | `Message.content.state_snapshot.workflow` | `workflow_config.workflow_data.steps[]`、`file_paths` 等 | 整块 dict 存入 `Message.content` 的 JSON |
| 执行清单 (Checklist/Steps) | `Message.content.state_snapshot.steps` | `[{ step_name, status, duration, step_result, ... }]` | 同上，数组在 JSON 内 |
| 图表/报告 (Report) | `Message.content.state_snapshot.report` | `report_data`、`diagnosis`、`summary`、`steps_details`、`images` 等 | 同上 |
| 纯文本回复 | `Message.content.state_snapshot.text` 或 `content.text` | 字符串 | 同上 |

**后端写入快照的唯一点**（`server.py` 流式结束）：

```python
# server.py L2248-2260
content = {
    "state_snapshot": state_snapshot_for_db if state_snapshot_for_db is not None
    else {"text": "", "workflow": None, "steps": [], "report": None}
}
msg = MessageModel(session_id=req.session_id, role="agent", content=content)
db.add(msg)
db.commit()
```

- `state_snapshot_for_db` 来自 SSE 流中 **最后一次** `event: state_snapshot` 的 `data` 解析结果。
- 若流异常中断、未发出 `state_snapshot` 或解析失败，则 `state_snapshot_for_db` 为 `None`，入库为占位结构 `{text:"", workflow:null, steps:[], report:null}`。
- **无版本号、无 schema 校验**：前端若依赖 `workflow.workflow_config.workflow_data.steps` 等路径，一旦后端写入结构变化或缺失，前端易崩。

### 1.3 潜在存储层风险

- **MySQL JSON 列**: 大文档（如含 base64 图、长步骤列表）可能触及 `max_allowed_packet` 或实际存储/复制性能问题；当前实现**无分片、无截断**。
- **类型与编码**: `gibh_agent/db/connection.py` 使用 `charset=utf8mb4`，JSON 内 emoji/特殊字符可正常存；但若驱动或 ORM 在某环境下将 JSON 反序列化为非标准类型（如 `Decimal`），需依赖 `sanitize_for_json` 在 API 层兜底。
- **content 为 NULL**: 若历史数据或异常导致某条 `Message.content = NULL`，API 返回的 `content` 为 `null`，前端若未显式处理会跳过该条或报错（见第 3 节）。

---

## 2. 后端接口层 (API Endpoints)

### 2.1 历史记录核心链路

| 接口 | 方法 | 用途 | 鉴权 |
|------|------|------|------|
| `GET /api/sessions` | list_sessions | 当前用户会话列表，按 created_at 倒序 | `get_current_owner_id` |
| `GET /api/sessions/{session_id}/messages` | list_messages | 指定会话的**全量消息**（含 content） | 同上 + 校验 session.owner_id |
| `POST /api/chat` | chat_endpoint | 发消息 + SSE 流，结束时写 Session/Message | 同上 |

**鉴权与会话隔离**（`gibh_agent/core/deps.py`）:

```python
def get_current_owner_id(request: Request) -> str:
    # 优先 Bearer JWT -> sub 作为 owner_id
    auth = request.headers.get("Authorization")
    if auth and auth.startswith("Bearer "):
        token = auth[7:].strip()
        if token and token.lower() != "null":
            try:
                payload = decode_access_token(token)
                if payload and payload.get("sub"):
                    return str(payload["sub"]).strip()
            except Exception:
                pass
    # 否则 X-Guest-UUID
    guest = request.headers.get("X-Guest-UUID")
    if guest and guest.strip():
        return guest.strip()
    raise HTTPException(401, ...)
```

- **会话隔离**: `list_messages` 先查 `Session`，仅当 `session.owner_id == owner_id` 时返回消息；否则 403。
- **跨设备/跨用户**: B 设备若未带同一 `owner_id`（如未登录、或 guest_uuid 不同），同一 `session_id` 会 403 或 401，前端只做 `console.warn`，列表/内容为空，表现为“空页面”或白屏。

### 2.2 GET /api/sessions/{session_id}/messages 行为

**文件**: `gibh_agent/api/routers/user_data.py`

- **无分页**: `.order_by(MessageModel.created_at.asc()).all()` 一次拉全量。
- **无超时配置**: 依赖 FastAPI/uvicorn 默认；若单会话消息数百条且每条 content 含大 JSON，响应体与 DB 查询时间都可能很大。
- **序列化容错**: 对每条消息用 `sanitize_for_json(raw)`；若单条序列化失败，捕获后降级为占位 `content: { "text": "[该条消息内容无法解析]", "state_snapshot": null }`，避免整接口 500，但前端会拿到残缺数据。

```python
# user_data.py L100-120
for r in rows:
    raw = {"id": r.id, "session_id": r.session_id, "role": r.role, "content": r.content, ...}
    try:
        out.append(sanitize_for_json(raw))
    except Exception as row_err:
        logger.warning("单条消息序列化失败 message_id=%s: %s", ...)
        out.append(sanitize_for_json({
            "id": ..., "content": {"text": "[该条消息内容无法解析]", "state_snapshot": None},
            ...
        }))
```

- **截断**: 接口层**没有**对 `content` 或 `state_snapshot` 做长度/深度截断；大快照会完整返回。

---

## 3. 前端状态恢复与 DOM 渲染层 (Frontend Hydration)

### 3.1 入口：加载会话并恢复

**文件**: `services/nginx/html/index.html`  
**函数**: `loadSessionAndRestore(sessionId, title)`（约 L3655）

流程概要：

1. 清空 `#chatContent`、`#active-workspace-content`、`chatHistory`、`window.workflowCache` 等。
2. `currentSessionId = sessionId`。
3. `fetch('/api/sessions/' + sessionId + '/messages', { headers: getAuthHeaders() })`。
4. `r.ok ? r.json() : Promise.reject(new Error(r.status))` → 非 2xx 直接进入 catch，仅 `console.warn`，**不渲染任何消息**。
5. 若成功：`for (var i = 0; i < messages.length; i++)` 遍历，按 `msg.role` 分支：
   - **user**: `appendMessage('user', text, ...)`。
   - **agent**: 仅当 `content` 为真值时进入分支，根据 `content.state_snapshot` / `content.events` / 否则 fallback 文本，调用 `createRestoreAIMessageRowWithSnapshot(snapshot, msg.id)`。

关键代码：

```javascript
// index.html L3686-3706
for (var i = 0; i < messages.length; i++) {
    var msg = messages[i];
    var content = msg.content;
    if (typeof content === 'string') {
        try { content = JSON.parse(content); } catch (e) { content = { text: content }; }
    }
    if (msg.role === 'user') {
        var text = (content && content.text) || ...;
        appendMessage('user', text || '', ...);
    } else if (msg.role === 'agent' && content) {
        if (content.state_snapshot) {
            createRestoreAIMessageRowWithSnapshot(content.state_snapshot, msg.id);
        } else if (content.events) { ... }
        else {
            var fallbackText = (content.text != null ? content.text : (content.message || ''));
            createRestoreAIMessageRowWithSnapshot({ text: fallbackText, workflow: null, steps: [], report: null }, msg.id);
        }
    }
}
```

- **致命点 1**: `msg.role === 'agent' && content` — 当 `content === null`（DB 中 NULL）或 `content === undefined` 时，该条 agent 消息**被完全跳过**，不创建任何 DOM，也不报错。
- **致命点 2**: 若 `r.ok` 为 false（401/403/500），catch 仅 warn，界面已被清空，表现为**白屏或空对话**。

### 3.2 反序列化与 content 形态

- 后端返回的 `content` 多为**已解析对象**（ORM 的 JSON 列反序列化后经 `sanitize_for_json` 再 json 序列化返回）；前端 `r.json()` 得到的是对象。
- 仅当某处将 content 存成字符串时，前端才有 `typeof content === 'string'` 的 `JSON.parse` 分支；若解析抛错，降级为 `{ text: content }`，不会抛到外层，但 `state_snapshot` 丢失。

### 3.3 createRestoreAIMessageRowWithSnapshot 与渲染链

**函数**: `createRestoreAIMessageRowWithSnapshot(state_snapshot, messageId)`（约 L3997）

- 入参规整：`var snap = state_snapshot && typeof state_snapshot === 'object' ? state_snapshot : { text: '', workflow: null, steps: [], report: null };`
- 先创建 DOM 行并 **appendChild 到 chatContent**，再根据 `snap` 填文本、工作流卡片、执行清单、报告。

**子调用**:

1. **文本**: `extractAndStripSuggestionTags(snap.text)` → `safeMarkedParse(text)` 写 `textZone.innerHTML`。  
   - **safeMarkedParse**（L2548）：若 `marked` 未定义或无 `parse`，则用简单 replace 做转义与换行；若 `text` 为 `null/undefined`，会变成 `undefined.replace(...)` **抛错**。
2. **工作流卡片**: `renderStaticHistorySnapshot(snap, cardZone, logZone)` → `_renderStaticWorkflowCard(snap.workflow, targetZone)`。  
   - `_renderStaticWorkflowCard` 内有多层可选链式取值（如 `workflowData.workflow_config.workflow_data.steps`），但若 `workflowData` 存在且为畸形对象（如 `steps` 为非数组），后续 `workflowSteps.length`、`forEach` 等仍可能抛错。
3. **报告**: `renderStaticSnapshotReport(reportData, targetContainer, options)`。  
   - 开头有 `if (!targetContainer || !reportData || typeof reportData !== 'object') return;`，相对安全；内部对 `stepsDetails`、`images` 等有 `Array.isArray` 或过滤，但若 `reportData.report_data` 或嵌套字段为畸形类型，仍可能在某处抛错。

**整体 try-catch 情况**:

- `createRestoreAIMessageRowWithSnapshot` 自身**无** try-catch；若其中任一步（如 `safeMarkedParse(undefined)`、`_renderStaticWorkflowCard` 内对畸形 workflow 的遍历）抛错，会**中断调用栈**。
- 在 `loadSessionAndRestore` 的 for 循环中，若某条消息的 `createRestoreAIMessageRowWithSnapshot` 抛错，**后续消息都不会再渲染**，且错误只会在控制台未捕获处暴露，用户看到的是**部分历史 + 白屏或卡住**。
- `_renderStaticWorkflowCard` 在 `renderStaticHistorySnapshot` 中被包在 `try { formId = _renderStaticWorkflowCard(...) } catch (e) { console.warn(...) }` 中，仅工作流卡片单块被保护，**同一消息内的 report/checklist 渲染仍无保护**。

### 3.4 致命弱点小结

| 环节 | 问题 | 后果 |
|------|------|------|
| fetch 失败 (401/403/5xx) | 仅 catch + console.warn，不渲染 | 清空后无内容，白屏 |
| `msg.content === null` | agent 分支被跳过 | 该条消息缺失，若全为这类则空列表 |
| `safeMarkedParse(text)` | 未校验 `text` 类型，可能 `undefined` | 抛错，整条或后续消息不渲染 |
| 单条 `createRestoreAIMessageRowWithSnapshot` 抛错 | 循环无 try-catch | 后续消息全部不渲染，半屏或白屏 |
| `snap.workflow` 结构畸形 | 仅工作流卡片有 try-catch，report/checklist 无 | 仍可能在某子渲染中抛错 |
| 无 snapshot 版本/schema | 后端结构演进或缺失字段 | 前端访问 undefined 属性或类型不符导致运行时错误 |

---

## 4. 架构师视角的漏洞诊断 (Vulnerability Assessment)

### 4.1 为何 A 电脑能恢复、B 电脑空页面？—— 三个最可疑逻辑断层

**（1）鉴权与 identity 不一致（LocalStorage / 设备隔离）**

- 前端所有需鉴权的请求均通过 `getAuthHeaders()` 取 `Authorization: Bearer <token>` 或 `X-Guest-UUID`（见 `index.html` 多处 fetch）。
- Token 与 guest_uuid 来源：登录后写 `localStorage.setItem('access_token', ...)`、`guest_uuid` 在未登录时由前端生成并写入 localStorage。
- **B 设备**：未登录则无同一 `access_token`；若用游客则 `guest_uuid` 为新生成，与 A 设备不同。此时：
  - `GET /api/sessions` 可能返回 B 设备自己的空列表或另一套会话；
  - 若用户通过某种方式拿到 A 的 session_id 在 B 上打开（如分享链接），`GET /api/sessions/{session_id}/messages` 会因 `session.owner_id !== owner_id` 返回 **403**。
- fetch 返回非 ok 时，前端仅 `Promise.reject(new Error(r.status))` → catch 里 `console.warn('[loadSessionAndRestore]', e)`，**不展示任何错误 UI**，且此时 `#chatContent` 已被清空 → **白屏**。
- **结论**: 身份完全依赖前端 LocalStorage，跨设备/跨浏览器必然 identity 不同，同一会话在 B 上要么 403 要么拿到的是 B 的会话列表，是“在 B 电脑上加载出空页面”的**首要嫌疑**。

**（2）单条消息解析或渲染抛错导致整页渲染中断（无 per-message try-catch）**

- `loadSessionAndRestore` 的 for 循环中，对每条消息直接调用 `appendMessage` 或 `createRestoreAIMessageRowWithSnapshot`，**没有** per-message 的 try-catch。
- 若其中一条出现：
  - `content` 为字符串但非合法 JSON，已用 `content = { text: content }` 兜底，一般不会抛；
  - 但若 `content.state_snapshot` 存在且为**畸形结构**（如某字段为不可序列化或前端未预期的类型、或后端曾写入过旧版 schema），在 `createRestoreAIMessageRowWithSnapshot` 内部（如 `safeMarkedParse(snap.text)` 当 `snap.text` 为 undefined、或 report 中某属性访问报错）会抛错。
- 一旦抛错，循环终止，**后续消息全部不渲染**；用户看到的是“只恢复了前几条或空屏”。
- A 电脑若从未加载过该条“问题消息”（例如该消息只在 B 首次拉取时返回），或 A 的缓存/顺序不同，可能不会触发该条，表现为“A 正常、B 白屏或半屏”。

**（3）content 为 NULL 或缺失 state_snapshot 导致 agent 消息整条被跳过**

- DB 中若某条 `Message.content = NULL`（历史脏数据、迁移或写入失败），API 返回 `content: null`。
- 前端逻辑 `msg.role === 'agent' && content` 为 false，该条 agent 消息**直接被跳过**，不调用 `createRestoreAIMessageRowWithSnapshot`。
- 若会话内大量或全部 agent 消息的 content 为 null（例如某次 bug 导致只写了 user 消息没写 agent 的 content），则循环结束后聊天区**几乎没有内容**，仅剩 user 气泡，看起来像“空页面”或数据丢失。
- 不同设备/不同会话下，是否命中“content 为 null”的记录具有随机性，可表现为“A 电脑能恢复、B 电脑空”。

---

## 5. 关键代码索引（便于重构时精确定位）

| 模块 | 文件路径 | 行号/函数 |
|------|----------|-----------|
| Message 模型 | `gibh_agent/db/models.py` | L39-49 `Message` |
| 会话/消息 API | `gibh_agent/api/routers/user_data.py` | L80-126 `list_messages` |
| 鉴权 | `gibh_agent/core/deps.py` | L22-44 `get_current_owner_id` |
| 快照写入 DB | `server.py` | L2248-2262 `content = {"state_snapshot": ...}` |
| 快照聚合 | `gibh_agent/core/orchestrator.py` | L217-221 `state_snapshot = {...}`，L2029 `_sanitize_snapshot_text` |
| 前端加载历史 | `services/nginx/html/index.html` | L3655 `loadSessionAndRestore`，L3682-3711 fetch + 循环 |
| 前端单条恢复 | `services/nginx/html/index.html` | L3997 `createRestoreAIMessageRowWithSnapshot`，L3964 `renderStaticHistorySnapshot`，L3892 `renderStaticSnapshotReport`，L3898 `_renderStaticWorkflowCard` |
| 前端鉴权头 | `services/nginx/html/index.html` | L2439 `getAuthHeaders`，L2422-2455 localStorage token/guest_uuid |
| 前端 Markdown | `services/nginx/html/index.html` | L2548 `safeMarkedParse` |

---

## 6. 重构建议（摘要）

1. **后端**: 为 `state_snapshot` 引入版本或 schema 标识；对 `list_messages` 做分页或按需加载大字段；保证写入失败或流中断时仍写入可解析的占位 content，避免 NULL。
2. **API**: 403/401 时返回明确错误码/文案，前端统一处理为“未授权/无权限”提示，避免静默白屏。
3. **前端**: `loadSessionAndRestore` 内对单条消息用 try-catch 包裹，单条失败时记录并跳过或渲染占位，不中断后续消息；对 `content === null` 的 agent 消息显式渲染占位行；`safeMarkedParse` 对非字符串入参做兜底。
4. **跨设备**: 如需“同一会话多端可读”，需在业务上明确身份与共享策略（账号登录、链接+只读 token 等），而不是依赖设备本地 guest_uuid。

以上为当前历史记录与高维 JSON 快照系统的技术现状与主要风险点，供架构重构决策使用。
