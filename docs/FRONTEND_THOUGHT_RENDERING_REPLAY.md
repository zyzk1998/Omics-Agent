## 目标
让“纯文本场景”的 `thought` 渲染体验回退到修改前的稳定实现，并把实现经验沉淀为可复刻的“渲染状态机 + DOM 约束”文档。

当前代码区分了两条渲染路径：
1. 纯文本闲聊（`workflowData` 为空）：渲染 `agent-thought-block`，并对 `thought` 事件做增量 `appendChild` 打字机效果。
2. 结构化场景（`workflowData` 非空，如工作流/报告）：渲染 `thinking-box` 到 `reasoning-container`（用以修复层级错位）。

本文件重点记录“纯文本路径”的经验细节，确保后续在非纯文本路径复刻时不会再次丢失关键行为。

---

## 纯文本路径：渲染状态机（SSE）

### 关键全局状态
- `window.currentThoughtContainer`
  - 类型：DOM 元素引用（`.agent-thought-block .thought-content`）
  - 赋值来源：`ensureThoughtBlock()`（创建/定位 thought-content 后写入）
- `thoughtBuffer`
  - 用于累积日志/调试（纯文本路径最终渲染仍依赖增量 `thoughtText`）
- `messageBuffer`
  - 用于累积回答内容，最后通过 `<think>` 清理后整体渲染

### DOM 结构（旧版 agent-thought-block）
当首次收到 `thought` 时，如果当前插入父容器里不存在 `.agent-thought-block`，会创建：
- `details.agent-thought-block`（默认打开）
- `summary` 包含：
  - 开启动画 SVG（`thought-icon-open`，转圈）
  - 关闭箭头 SVG（`thought-icon-closed`）
  - 文本：`<span class="thought-summary-text">思考中...</span>`
- `div.thought-content`：`thought` 增量文本的实际落点

对应的创建逻辑在：
- `ensureThoughtBlock()` 内部：
  - 选择插入容器：`reasoningContainer || messageContent`
  - 若不存在 `.agent-thought-block`，则用 `block.innerHTML = '...'` 创建完整结构
  - 设置 `window.currentThoughtContainer = contentEl`（`.thought-content` 节点）

### thought 事件处理（增量打印的关键）
在 `handleServerEvent(eventType, data)` 的 `case 'thought'` 中：
1. 从 payload 中抽取 `thoughtText`（支持 `string` / `content` / `text` / `reasoning` 等形态）
2. 累积到 `thoughtBuffer`（用于 debug）
3. **纯文本路径的核心动作：**
   - 确保 `window.currentThoughtContainer` 存在（必要时调用 `ensureThoughtBlock()`）
   - **如果本次 `thoughtText` 非空：**
     - `el.appendChild(document.createTextNode(thoughtText))`
     - `el.scrollTop = el.scrollHeight` 保证滚动跟随
     - 找到最近的 `.agent-thought-block` 并强制 `details[open]`：
       - `thoughtDetails.setAttribute('open', '')`

为什么这里必须是增量 `appendChild`？
- `updateChatBubbleWithReasoning()` 会把 thought 渲染为整体覆盖（`textContent = parsed.thinking`），在纯文本体验上会变成“整段替换”，失去旧版打字机的连续感。
- 增量 appendChild 能保证 chunk 对齐与视觉节奏稳定。

### message 事件处理（关闭思考 + 渲染答案）
在 `case 'message'` 中：
1. `finalizeThoughtBlock()`：
   - 找到 `window.currentThoughtContainer` 所在的 `.agent-thought-block`
   - 移除 `open`（并移除 `thought-icon-open`）
   - 把 summary 文本从 `思考中...` 变为 `思考过程`
2. 从 `messageBuffer` 中清理 `<think>` 标签：
   - `messageBuffer.replace(/<\/?think>/g, '')`
3. 使用 `updateChatBubble(cleanedText)`：
   - 只更新 `.textZone` / `aiContentDiv`（不触碰 reasoning-zone 的思考 DOM）

### 新对话历史保留约束（与本次回退的结合点）
为了“历史思考过程永久保留”，当前实现调整为：
- `clearThoughtBlock()` 只重置全局引用，不再对历史 DOM 执行任何 `textContent=''` / `removeAttribute('open')` / 清空 `.thinking-box` 等操作。
- 这样旧消息里的思考节点不会被销毁，用户可以查看历史时光机。

---

## 结构化路径：为什么要避免复用纯文本逻辑
结构化场景（`workflowData` 非空）为了修复层级错位：
- 使用 `updateChatBubbleWithReasoning()` 把思考折叠框挂在 `.reasoning-container` 上方
- 并在函数内移除 `.spinner-border` / “思考中”文本节点

需要注意的差异：
- `updateChatBubbleWithReasoning()` 会创建并覆盖 `thinking-box`（`targetReasoning.innerHTML = ''`）
- 纯文本路径不应该触发它，否则会把体验从“增量打字机”退化为“折叠框整体覆盖”，并且可能破坏历史消息结构的稳定性。

---

## 实施要点清单（用于后续复刻）
当你要在非纯文本场景复刻“纯文本思考成功打印经验”时，请逐项确认：
1. thought 落点必须是独立节点：`.agent-thought-block .thought-content`
2. thought 渲染必须是增量追加：
   - `appendChild(TextNode)` 而不是 `textContent = 全量`
3. 必须滚动跟随：
   - `scrollTop = scrollHeight`
4. thought details 必须保证处于 open 状态：
   - `setAttribute('open','')`
5. 收尾必须在 message 转入答案时执行：
   - `finalizeThoughtBlock()` 修改 summary 与图标
6. 回答渲染必须清理 `<think>` 标签，并且只更新文本区：
   - `updateChatBubble(cleanedText)` 不能碰 reasoning-zone 思考 DOM
7. 历史保留策略：
   - 清空只限全局引用，不允许清空历史节点内容

