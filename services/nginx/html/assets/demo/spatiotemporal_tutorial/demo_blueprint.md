# 沉浸式 Demo：状态机与交互蓝图

本文档为 `services/nginx/html/index.html` 与 `css/main.css` 中「单细胞时空动力学」沉浸式演示的**唯一设计依据**。实现须与此一致。

---

## 1. UI 布局规范

| 项 | 要求 |
|----|------|
| 容器 DOM | `#immersiveDemoEntryWrap`（内联 `display:flex; justify-content:center; width:100%; margin-top/bottom: 10px`） |
| 位置 | **必须**在 `.chat-input-panel` **内部**、作为**最后一个子节点**（在 `.toolbar` 下方）；**禁止**塞进 `.input-area` 或与发送按钮同一 flex 行 |
| 布局 | 整行居中；与 `.chat-input-panel` 的 `flex-direction: column` 自然堆叠 |
| 按钮 | `#btnImmersiveSpatiotemporalDemo`，类名 `btn-immersive-demo`（透明底 + 主题色描边 + `breathing-glow` 周期外发光，见 `main.css`）；初始 `display: none` |
| 可见性 | 仅**未登录**访客展示；已登录隐藏（`window.syncImmersiveDemoEntryVisibility()`：存在 `localStorage.token` 或 `access_token` 时隐藏按钮并 `display:none` 整行 wrap） |

### 按钮视觉（`.btn-immersive-demo`）

- 科研线条风：透明底、`1px solid var(--primary-color)`、文字同色、`border-radius: 6px`、`font-weight: 500`、`font-size: 14px`
- `@keyframes breathing-glow`：`box-shadow` 约 2s 周期弱→强循环（见 `main.css`）；`:hover` 使用 `rgba(var(--primary-color-rgb), 0.1)` 淡底

---

## 2. 状态机枚举 `window.DEMO_STATE`

全局标志：`window._immersiveDemoActive === true` 表示用户已从入口进入沉浸式链路。

| 状态 | 含义 | 进入条件 / 动作 |
|------|------|-----------------|
| `IDLE` | 默认，真实用户模式 | `_immersiveDemoActive === false`，不劫持 `fetch` |
| `ARMED` | 已点沉浸式入口 | 自动填入技能 `prompt_template`（或兜底文案）；`DemoCoach` **指向 `#sendBtn`**，文案引导点击发送 |
| `PLAYING_PLAN` | 规划皮影戏 | 在 `sendMessage` 中由 `ARMED` + 无 `workflowData` 转入；**拦截**真实请求，播放 `demo_planning_stream.txt`（经 `playImmersiveDemoStream`） |
| `PLAN_DONE` | 规划播放结束 | 播放完成后置位；`DemoCoach` **优先指向**「载入官方演示数据」按钮（`#chatContent .btn-use-demo`），若无则回退指向「执行工作流」按钮 |
| `PLAYING_EXEC` | 执行皮影戏 | `PLAN_DONE` 且本次 `sendMessage` 携带 **`workflowData`（执行/载入工作流）** 时转入；**拦截**真实请求，播放 `demo_execution_stream.txt` |
| `COMPLETED` | 全流程结束 | 执行皮影戏结束后；`_immersiveDemoActive = false`，`DEMO_STATE` 可视为结束态；弹出注册/登录引流（`openAuthModal` + Toast） |

**重置**：`clearChatAndReset`（及同类「新会话」）须将 `DEMO_STATE = 'IDLE'`、`_immersiveDemoActive = false`，并 `DemoCoach.hide()`。

**说明**：不再使用中间态别名（如 `GUIDING_EXEC`）；规划结束后的等待统一为 `PLAN_DONE`。

---

## 3. 事件劫持锚点（`sendMessage`）

### 3.1 拦截位置

在 `sendMessage` 内，在构造完 `payload`、打完诊断日志之后，**在任何** `fetch('/api/chat', …)` **之前**，插入沉浸式分支：

```text
// 沉浸式 Demo：皮影戏回放（不发起真实 /api/chat）
if (window._immersiveDemoActive && typeof window.playImmersiveDemoStream === 'function') {
    // ARMED → PLAYING_PLAN → 读 demo_planning_stream.txt
    // PLAN_DONE + workflowData → PLAYING_EXEC → 读 demo_execution_stream.txt
    // … 各分支末尾 return，禁止落入下方 fetch
}
try {
    const res = await fetch('/api/chat', { ... });
```

（行号会随文件变动，以 **`fetch('/api/chat'` 紧上方** 为准。）

### 3.2 如何拦截真实网络请求

1. 满足 `_immersiveDemoActive` 且命中「规划回放」或「执行回放」分支时，**必须 `return`**，不得执行后续 `fetch`。
2. 回放通过 `await playImmersiveDemoStream(...)` 消费静态 SSE 剧本；**解析器须按行状态机**（`event:` / `data:` / 空行分块、**`data:` 后无前缀的续行拼入同一 payload**），再 `JSON.parse` 后调用 `handleServerEvent`。**禁止**仅用 `\n\n` 整块切分（易与 JSON 内换行或长行工具冲突）。
3. 游客登录门控：在 `sendMessage` 靠前位置，仅当 `_immersiveDemoActive && (DEMO_STATE === 'ARMED' || (DEMO_STATE === 'PLAN_DONE' && workflowData))` 时允许跳过 `isLoggedIn()` 检查（与上表一致）。

### 3.3 静态剧本路径

- 规划：`/assets/demo/spatiotemporal_tutorial/demo_planning_stream.txt`
- 执行：`/assets/demo/spatiotemporal_tutorial/demo_execution_stream.txt`

---

## 4. 与 `window.onload` 的衔接

- 调用 `syncImmersiveDemoEntryVisibility()`。
- 绑定 `#btnImmersiveSpatiotemporalDemo` 的 `click`：置 `_immersiveDemoActive = true`、`DEMO_STATE = 'ARMED'`，填 prompt，调 `DemoCoach.show` 指向 `#sendBtn`。

---

*版本：与仓库 HTML/CSS 同步维护。*
