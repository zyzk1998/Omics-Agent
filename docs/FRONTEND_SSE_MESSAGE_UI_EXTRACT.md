# 前端三大模块代码提取：SSE 接收解析、消息数据结构、聊天气泡 UI

本文档从项目 `services/nginx/html/index.html` 中提取并整理以下三部分代码，便于查阅与排查。

**说明**：本项目前端为**纯 HTML + 内联 JavaScript**，无 Vue/React；所有逻辑均在 `index.html` 内。

---

## 一、SSE 接收与解析逻辑

### 1.1 使用方式：fetch API + ReadableStream

**未使用** 原生 `EventSource`，而是使用 **fetch + `res.body.getReader()`** 配合 `TextDecoder` 逐块读取流式响应。

- 请求：`POST /api/chat`，`stream: true`，响应 `Content-Type: text/event-stream`。
- 读取：`const reader = res.body.getReader();`，循环 `await reader.read()`，用 `decoder.decode(value, { stream: true })` 解码后追加到 `buffer`。

### 1.2 核心代码：发起请求与进入 SSE 分支

```javascript
// 位置：index.html 内 sendMessage 逻辑
const res = await fetch('/api/chat', { 
    method: 'POST', 
    headers: { ...getAuthHeaders(), 'Content-Type': 'application/json' }, 
    body: JSON.stringify(payload),
    signal: currentAbortController.signal
});
const contentType = res.headers.get("content-type");

if (contentType && contentType.includes("text/event-stream")) {
    // SSE 流式传输模式
    const reader = res.body.getReader();
    const decoder = new TextDecoder();
    let buffer = '';
    let messageBuffer = '';
    let thoughtBuffer = '';
    // ... 见下方解析循环
}
```

### 1.3 SSE 解析循环：按 \n\n 分割事件，解析 event: / data:

```javascript
// 主循环：读取并解析 SSE 流
while (true) {
    const { done, value } = await reader.read();
    if (done) break;

    const chunk = decoder.decode(value, { stream: true });
    buffer += chunk;

    // 按 SSE 标准格式分割：\n\n 是事件分隔符
    const parts = buffer.split('\n\n');
    buffer = parts.pop() || ''; // 保留最后一个不完整的事件

    for (const part of parts) {
        if (!part.trim()) continue;
        const lines = part.split('\n');
        let eventType = 'message';
        let dataLines = [];
        for (const line of lines) {
            if (line.startsWith('event: ')) {
                eventType = line.substring(7).trim();
            } else if (line.startsWith('data: ')) {
                dataLines.push(line.substring(6));
            } else if (line.trim() && dataLines.length > 0) {
                dataLines.push(line);
            }
        }
        if (dataLines.length > 0) {
            const jsonData = dataLines.join('\n');
            try {
                const data = JSON.parse(jsonData);
                const shouldStop = handleServerEvent(eventType, data);
                if (shouldStop) return;
            } catch (e) {
                console.error('❌ SSE JSON 解析错误:', e);
            }
        }
    }
}
// 剩余 buffer 同样按 event:/data: 解析后调用 handleServerEvent(eventType, data);
```

### 1.4 thought 字段的提取与当前消息状态更新

- **提取**：在 **`handleServerEvent(eventType, data)`** 的 `case 'thought':` 分支中，直接使用服务端下发的 `data.content`，不单独解析“thought 字段”；即 **服务端发 `event: thought` + `data: { "content": "..." }`**，前端用 `data.content` 作为思考内容。
- **更新方式**：
  1. 将 `data.content` 追加到局部变量 `thoughtBuffer`（仅用于本流式会话的累积，非全局状态）。
  2. 调用 **`ensureThoughtBlock()`** 获取或创建当前 AI 气泡内的 `.thought-content` 容器。
  3. 用 **`el.appendChild(document.createTextNode(data.content || ''))`** 把本次内容追加到该容器，并 `el.scrollTop = el.scrollHeight` 保持滚动到底；若外层是 `<details>` 则设置 `block.open = true` 保持展开。

```javascript
// handleServerEvent 内
case 'thought':
    thoughtBuffer += data.content || '';
    try {
        var el = ensureThoughtBlock();
        if (el) {
            el.appendChild(document.createTextNode(data.content || ''));
            el.scrollTop = el.scrollHeight;
            var block = el.closest('.agent-thought-block');
            if (block && block.tagName === 'DETAILS') block.open = true;
        }
    } catch (e) { console.warn('[thought] 渲染失败:', e); }
    break;
```

**结论**：thought 内容来自服务端 SSE 的 `event: thought` 的 `data.content`；前端不解析“thought”字段名，只根据事件类型 `thought` 更新当前 AI 气泡内的 `.thought-content` DOM，无单独的消息状态对象（如 Vue/React 的 data.thought），状态体现在 DOM 上。

---

## 二、消息的数据结构定义 (State Interface)

### 2.1 前端存储一条“对话轮次”的结构（chatHistory）

前端用 **数组 `chatHistory`** 存“一轮对话”的摘要，**每条元素** 为：

```javascript
// 声明（约 1868 行）
let chatHistory = [];

// 在 SSE 流结束、合并完 AI 回复后 push（约 5624 行）
chatHistory.push({ user: text, ai: fullText });
```

即：

- **类型**：`Array<{ user: string, ai: string }>`
- **user**：当轮用户发送的文本。
- **ai**：当轮 AI 的完整回复文本（流式结束后合并得到的 `fullText`）。

该结构主要用于 **发给后端的 `history` 上下文**（如 `payload.history = chatHistory`），**不是** 单条消息的完整结构（无 role、id、state_snapshot 等）。

### 2.2 单条消息在 DOM 上的表现（无独立 JS 对象）

单条消息**没有**单独的“消息状态接口”对象；结构体现在 **DOM** 上：

- **容器**：`#chatContent` 下的 `.message-row`（class 含 `user` 或 `ai`）。
- **用户消息**：`message-row.user` → 内含 `.message-container`、`.bubble`（正文）、可选 `.file-attachment-container`、操作栏等。
- **AI 消息**：`message-row.ai` → 内含 `.message-container`，其下可有：
  - `.process-log-zone`（执行记录）
  - `.text-content-zone` / `.bubble`（正文）
  - `.cards-zone`（工作流卡片等）
  - `.agent-thought-block`（深度思考折叠块）
  - `.message-footer`（耗时、操作按钮）等。
- **持久化相关**：`row.dataset.messageId` 与后端 `message_saved` 事件下发的 `message_id` 绑定。

创建 AI 行时通过 **`appendMessage('ai', ...)`** 或流式过程中创建的 **`window.currentMessageZones`**（含 `logZone`、`textZone`、`cardZone`）指向当前正在更新的那一条 AI 消息的 DOM 区域。

### 2.3 后端历史消息结构（供对比）

从 `loadSessionAndRestore` 可见，**后端返回的一条消息** 大致为：

```javascript
// 来自 GET /api/sessions/:id/messages
{
  id: number,           // 消息 ID
  role: 'user' | 'agent',
  content: {
    state_snapshot: { text, workflow, steps, report }  // 新格式
    // 或 events: [...]   // 旧格式
    // 或 text / message 等
  }
}
```

前端**不**把这份结构存成全局的“消息列表状态”，而是遍历后直接生成 DOM（或注入到 `chatHistory` 的 user/ai 文本）。

---

## 三、聊天气泡的 UI 渲染组件（深度思考 + 执行记录）

本项目为**单文件 HTML+JS**，没有独立 Vue/React 组件文件；“组件”以 **DOM 创建 + 内联样式/类名** 实现。下面给出“深度思考中...”和“⚡ 执行记录 (N)”对应的 **逻辑与样式** 所在位置与代码。

### 3.1 “🧠 深度思考中...” 折叠面板

- **实现**：**`<details>` + `<summary>` + `.thought-content`**，通过 **`ensureThoughtBlock()`** 在**当前 AI 气泡的 `.message-content` 内**动态创建或获取。
- **“深度思考中...” 文案**：**写死的静态文本**，写在 `summary` 的 `innerHTML` 里，**没有** 绑定到变量；仅 `.thought-content` 内为动态追加的思考内容。

```javascript
// ensureThoughtBlock（约 3514–3530 行）
function ensureThoughtBlock() {
    var messageContent = ...; // 当前 AI 的 .message-content
    var block = messageContent.querySelector('.agent-thought-block');
    if (!block) {
        block = document.createElement('details');
        block.className = 'agent-thought-block';
        block.setAttribute('open', '');
        block.innerHTML = '<summary>🧠 深度思考中...</summary><div class="thought-content"></div>';
        messageContent.insertBefore(block, messageContent.firstChild);
    }
    return block.querySelector('.thought-content');
}
```

**CSS（约 1323–1346 行）**：

```css
.agent-thought-block {
    margin-bottom: 16px;
    border-left: 3px solid #e5e7eb;
    padding-left: 12px;
}
.agent-thought-block summary {
    cursor: pointer;
    color: #6b7280;
    font-size: 13px;
    font-weight: 500;
    list-style: none;
}
.agent-thought-block summary::-webkit-details-marker { display: none; }
.thought-content {
    margin-top: 8px;
    font-size: 13px;
    color: #9ca3af;
    font-style: italic;
    line-height: 1.6;
    max-height: 300px;
    overflow-y: auto;
    white-space: pre-wrap;
}
```

### 3.2 “⚡ 执行记录 (N)” 区域

- **实现**：在检测到 SSE 流式响应后，在 **当前 AI 的 `logZone`**（或回退到 `messageContent`）里插入 **`.process-log-container`**，其内为 **`.process-header`**（标题 + 数字） + **`.process-steps.checklist-container`**（步骤列表）。数字 N 来自 **`.step-count`**，由 **`addProcessStep`** 在每次添加/更新步骤时更新。
- **标题**：“⚡ 执行记录” 为**固定文案**，只有括号里的数字是**动态**（绑定到 `processStepCount`，写入 `.step-count`）。

**创建 DOM（约 3359–3368 行）**：

```javascript
processLogContainer = document.createElement('div');
processLogContainer.className = 'process-log-container';
processLogContainer.innerHTML = `
    <div class="process-header" onclick="...">
        <span>⚡ 执行记录 (<span class="step-count">0</span>)</span>
        <span class="look-right-hint">  <span class="arrow-bounce">➔</span></span>
        <span class="icon">▼</span>
    </div>
    <div class="process-steps checklist-container"></div>
`;
processStepsDiv = processLogContainer.querySelector('.process-steps');
// 追加到 logZone 或 messageContent
```

**添加单步（addProcessStep，约 3390–3433 行）**：  
根据 `status`（active/running/done/error）设置图标与副标题，并更新 `processLogContainer.querySelector('.step-count').textContent = processStepCount`。

**CSS（约 840–876 行及 1347–1376 行）**：

```css
.process-log-container {
    margin: 12px 0;
    border: 1px solid #e5e7eb;
    border-radius: 8px;
    background: #f9fafb;
    overflow: hidden;
    animation: fadeIn 0.3s;
}
.process-header {
    padding: 10px 14px;
    background: #f3f4f6;
    border-bottom: 1px solid #e5e7eb;
    display: flex;
    align-items: center;
    justify-content: space-between;
    cursor: pointer;
    ...
}
.checklist-container { display: flex; flex-direction: column; gap: 12px; margin-top: 16px; }
.checklist-item { display: flex; align-items: flex-start; gap: 12px; }
.checklist-icon.icon-pending { color: #d1d5db; }
.checklist-icon.icon-running { color: #3b82f6; animation: spin 1s linear infinite; }
.checklist-icon.icon-success { color: #10b981; }
.checklist-icon.icon-error { color: #ef4444; }
.checklist-title { font-size: 14px; font-weight: 600; color: #1f2329; }
.checklist-subtitle { font-size: 12px; color: #6b7280; margin-top: 4px; }
```

### 3.3 小结

| 区域           | 实现方式           | 是否绑定变量                         |
|----------------|--------------------|--------------------------------------|
| “深度思考中...” | `<details>` + 静态 summary 文案 | 否，静态；仅 `.thought-content` 内容动态追加 |
| “执行记录 (N)”  | `.process-header` + `.step-count` | 是，N 来自 `processStepCount`，步骤列表由 `addProcessStep` 动态生成 |

---

*以上代码与行号均对应 `services/nginx/html/index.html`，行号可能随项目变更略有偏移，以实际文件为准。*
