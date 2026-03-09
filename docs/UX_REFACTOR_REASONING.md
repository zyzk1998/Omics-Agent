# 终极 UX 重构 - 深度推理报告

## 任务 1：JWT 密钥重置与 401 优雅降级

**根因**：若 SECRET_KEY 使用 `os.urandom()` 或每次进程启动变化，Docker 重启后旧 JWT 签名密钥变更，已签发的 token 验证失败 → 401。

**思路**：
- 后端：`SECRET_KEY` 从环境变量读取，提供**固定强哈希默认值**，确保重启后密钥不变。
- 前端：在 fetch 响应层面拦截 401（全局封装或各请求 catch），静默清除 `localStorage` 的 `access_token`、`username`，将左下角 UI 恢复为「游客」，并提示「登录已过期，请重新登录」。

**边界**：401 可能来自任意接口，统一清 token 并提示即可；多 tab 下仅当前 tab 即时更新，其他 tab 下次请求再 401 时同样会清 token。

---

## 任务 2：真·时光机（状态快照渲染）

**根因**：通过 for 循环重放 SSE 事件不仅慢，且易导致 DOM 闪烁与逻辑冲突（如多次 append、重复绑定）。

**思路**：
- 从 `content.events` 数组中**解析每条 SSE 行**（格式为 `event: xxx\ndata: {...}`），得到 `eventType` 与 `data`。
- **提取最终态**：顺序遍历 events，对以下字段做“最后一次覆盖”保留：
  - `steps_details`：来自 `data.report_data.steps_details` 或 `data.steps_details`
  - `diagnosis_report`：来自 `data.diagnosis_report` 或 `data.report_data.diagnosis`
  - `report_data`：整段（用于 summary/diagnosis 渲染）
- **瞬间渲染**：将最终 `stepsDetails` 中每项 `status` 归一为 `'success'`（避免任何转圈），一次调用 `renderExecutionSteps(stepsDetails)`；再根据最终 `diagnosis_report` / `report_data` 调用 `renderReportSection`。不调用 `handleServerEvent`，不重放。

**边界**：events 可能为混合顺序；规则为「最后一次出现的有效 payload 获胜」。若整段无 `steps_details`（如纯聊天回复），则只渲染文本区，不渲染步骤。若 `content.events` 缺失或格式异常，可 fallback 为仅渲染文本或保留原重放逻辑。

---

## 任务 3：聊天气泡 CRUD 控制权

**根因**：用户无法对已发送/已生成的消息进行复制、编辑、删除、重新生成，缺乏控制感。

**思路**：
- **UI**：在每个气泡（User / Agent）底部增加 `.message-action-bar`，默认隐藏、Hover 显示，含复制、编辑（仅 User）、重新生成（仅 Agent）、删除（SVG 图标）。
- **后端**：新增 `DELETE /api/messages/{message_id}`，校验该 message 所属 session 的 `owner_id` 与当前用户一致后物理删除。
- **前端**：复制 → `navigator.clipboard.writeText`；删除 → 调 DELETE 成功后 `element.remove()`；重新生成 → 删当前 AI 气泡，取上一条 User 内容，调用 `sendMessage`；编辑 → 内容填回输入框，删除该条及之后所有消息（后端可循环 DELETE 或提供 truncate 接口），前端同步移除对应 DOM。

**边界**：历史恢复时消息带 `id`，新建消息若后端未返回 id 可先不绑定删除或等列表刷新；编辑时「删除该消息及之后」需后端支持（单条 DELETE + 循环删除后续，或一次性 truncate）。

---

## 任务 4：数据资产无感自动刷新

**根因**：规划阶段后端已更新 `Asset.modality`，前端侧栏资产树不刷新，需用户手动刷新页面。

**思路**：在 `handleServerEvent` 的 `case 'done'` 中，在 `handleDoneEvent(data)` 之后**静默调用一次 `fetchAssets()`**，不弹 loading，用户看完报告后侧栏已更新。

**边界**：若 fetchAssets 返回 401，由全局 401 逻辑统一处理（清 token、提示）。

---

## 任务 5：技能广场 QoderWork 级微交互

**根因**：「使用此技能」按钮常驻显示，视觉拥挤。

**思路**：纯 CSS。`.skill-card .skill-use-btn` 默认 `opacity: 0; transform: translateY(10px); transition: all 0.3s ease;`；`.skill-card:hover .skill-use-btn` 设为 `opacity: 1; transform: translateY(0);`。

**边界**：无。
