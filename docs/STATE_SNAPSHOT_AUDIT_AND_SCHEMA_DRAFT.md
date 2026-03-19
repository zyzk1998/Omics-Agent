# 快照遗漏清单与 JSON 扩充结构草案（时光机补全）

**目标**：历史记录 100% 还原当前对话页面（耗时、执行记录文案、思考区、报告建议与质量评估等），**严禁存 DOM 字符串**，仅通过完善 `state_snapshot` JSON 实现。  
**约束**：新结构向前兼容旧历史；读取旧记录不得报错或白屏。

---

## 第一步：快照遗漏字段盘点 (Audit Missing States)

### 1.1 后端写入与前端消费对照

| 数据来源 | 当前 state_snapshot 是否包含 | 前端展示位置 | 恢复时是否可见 |
|----------|------------------------------|--------------|----------------|
| **总耗时 (duration)** | ❌ 否。done 时仅写入 SSE 事件 `data.duration`，未回写 `state_snapshot`；`_start_time` 在快照中但未在入库前转为 duration | 消息 footer 的 `.generation-timer`（如 "12.3s"） | ❌ 历史行 footer 无计时器 |
| **单步耗时 (step.duration)** | ❌ 否。Executor 构建 `step_detail` 时未写入 `duration`/`elapsed`；仅打日志 | 执行记录 checklist 副标题 / 报告手风琴（"耗时 X.Xs"） | ❌ 副标题无“耗时 Xs” |
| **执行记录逐条文案 (process log)** | ❌ 否。status 事件仅推 SSE，未追加到快照 | 左侧 process-log-zone（“正在接收请求…”“正在执行步骤: XXX”“完成 (12.3s)”） | ⚠️ 仅用 steps 渲染静态 checklist，无上述逐条文案 |
| **思考过程 (reasoning)** | ❌ 否。thought 与 message 均追加到 `state_snapshot["text"]`，无独立字段 | reasoning-zone（可折叠“思考过程”） | ❌ 历史不区分思考/正文，无单独思考框 |
| **报告建议 (suggestions)** | ✅ 是。经 diagnosis 的 report_data 合并进 `report.report_data`，含 `suggestions` | 报告下方建议 chips | ⚠️ 数据在快照，但 `renderStaticSnapshotReport` 未渲染 |
| **质量评估 (evaluation)** | ✅ 是。同上，`report.report_data.evaluation` | 报告内质量得分/评语 | ⚠️ 数据在快照，但 `renderStaticSnapshotReport` 未渲染 |
| **工作流卡片** | ✅ 是。`state_snapshot.workflow` | cards-zone | ✅ 可恢复 |
| **步骤结果 (steps / report_data.steps_details)** | ✅ 是。step_result 时写入 `steps` 与 `report.report_data` | 中栏 checklist、右栏报告手风琴+图表 | ✅ 可恢复（缺每步 duration） |
| **诊断/专家报告 (diagnosis/summary)** | ✅ 是。`report.diagnosis`、`report.report_data` | 报告区数据诊断、专家总结 | ✅ 可恢复 |

### 1.2 涉及代码位置（便于落地改）

- **后端**
  - 快照初始与聚合：`gibh_agent/core/orchestrator.py` 约 218–224 行（`state_snapshot` 初始化）、2168–2263 行（`_emit_sse` 内各 event_type 的写入）。
  - done 时 duration 仅写入事件 data、未写入快照：2256–2260 行；快照最终下发 911–912 行，此前未执行 `state_snapshot["duration"] = ...`。
  - 入库：`server.py` 约 2278–2311 行，解析 `event: state_snapshot` 的 data 作为 `content.state_snapshot`。
- **Executor 未写单步耗时**：`gibh_agent/core/executor.py` 约 1414–1416 行（`tool_start = time.time()`、`step_result = self.execute_step(...)`），约 1479–1493 行构建 `step_detail`，**未**设置 `step_detail["duration"]`。
- **前端**
  - 历史恢复入口：`services/nginx/html/index.html` 约 3796–3797、4181–4239 行（`createRestoreAIMessageRowWithSnapshot`）。
  - 静态执行记录/报告：约 3865–3890 行（`renderStaticChecklistToLogZone`，用 `step.duration`/`step.summary`）、3923–3986 行（`renderStaticSnapshotReport`，无 suggestions/evaluation 区块）。
  - 直播时 footer 计时器：约 7205–7281 行（done 时用 `data.duration` 更新 `.generation-timer`）；历史行创建时 footer 无 timer 节点（约 4210 行）。

### 1.3 遗漏清单小结

| 序号 | 遗漏项 | 类型 | 影响 |
|------|--------|------|------|
| 1 | **总耗时 (duration)** | 未写入快照 | 历史消息无“X.Xs”计时 |
| 2 | **单步耗时 (step.duration)** | Executor 未填 | 历史 checklist/报告无“耗时 Xs” |
| 3 | **执行记录逐条 (process_log)** | 未入快照 | 历史无法还原“正在接收请求…”等中间态文案 |
| 4 | **思考区独立文案 (reasoning)** | 与 text 混在一起 | 历史无单独“思考过程”折叠区 |
| 5 | **建议 chips / 质量评估 的渲染** | 数据在 report_data，前端未画 | 历史有数据但看不到建议与评估 |

---

## 第二步：快照 JSON 扩充方案 (Schema Upgrade)

### 2.1 现有结构（兼容基准）

```json
{
  "text": "",
  "workflow": null,
  "steps": [],
  "report": null,
  "_start_time": 1234567890.123
}
```

- `steps`：元素为执行步骤对象，含 `step_id`、`name`、`status`、`summary`、`step_result`、`data`、`plot` 等（无 `duration`）。
- `report`：`{ "diagnosis": "...", "report_data": { "diagnosis", "steps_details", "suggestions", "evaluation", ... } }`。

### 2.2 扩充字段（仅新增键，不删旧键）

| 字段路径 | 类型 | 说明 | 写入时机 |
|----------|------|------|----------|
| **`duration`** | number | 本条消息总耗时（秒），与 done 事件一致 | 在 yield `state_snapshot` 前由 orchestrator 写入：`duration = round(time.time() - _start_time, 1)` |
| **`steps[].duration`** | number | 该步骤执行耗时（秒） | Executor 在每次 `steps_details.append(step_detail)` 前写入 `step_detail["duration"] = round(time.time() - tool_start, 1)` |
| **`process_log`** | array | 执行记录逐条，与直播时 process-log-zone 顺序一致 | 每次 status 推送给前端时，orchestrator 同步 append `{ "content": "...", "state": "running"|"completed"|... }`；done 时追加 `{ "content": "完成 (Xs)", "state": "completed" }` |
| **`reasoning`** | string | 仅思考过程正文（thought 累积），与 `text` 分离 | 在 `_emit_sse` 中，当 `event_type == "thought"` 时除追加到 `text` 外，另做 `state_snapshot["reasoning"] = (state_snapshot.get("reasoning") or "") + (data.get("content") or "")`（可选：thought 只进 reasoning 不进 text，需与产品确认） |

**不存储**：任何 DOM 或 innerHTML；仅存上述结构化数据。

### 2.3 向前兼容

- 读取 `content.state_snapshot` 时，所有新字段按**可选**处理：`duration`/`process_log`/`reasoning` 缺省则不展示对应 UI；`steps[].duration` 缺省则副标题仅用 `summary`。
- 旧记录无 `report.report_data.suggestions`/`evaluation` 时，前端不渲染建议块/评估块即可，不报错。

### 2.4 建议的最终快照结构（示例）

```json
{
  "text": "最终回复正文（仅 message，不含 thought）",
  "reasoning": "思考过程全文（可选，有则恢复 reasoning-zone）",
  "workflow": { "workflow_data": {...}, "template_mode": false },
  "steps": [
    {
      "step_id": "step1",
      "name": "数据验证",
      "status": "success",
      "summary": "通过",
      "duration": 1.2,
      "step_result": {...},
      "data": {...},
      "plot": "/results/..."
    }
  ],
  "process_log": [
    { "content": "正在接收请求...", "state": "start" },
    { "content": "正在执行步骤: 数据验证", "state": "running" },
    { "content": "完成 (12.3s)", "state": "completed" }
  ],
  "duration": 12.3,
  "report": {
    "diagnosis": "...",
    "report_data": {
      "diagnosis": "...",
      "steps_details": [...],
      "suggestions": ["建议1", "建议2"],
      "evaluation": { "score": 85, "comment": "..." }
    }
  },
  "_start_time": 1234567890.123
}
```

**说明**：若希望 `text` 与直播完全一致（含 thought 拼接），可保留当前行为，仅新增 `reasoning` 用于历史时单独渲染思考区；否则可改为 thought 只进 `reasoning`、不进 `text`，由产品决定。

---

## 第三步：前端还原逻辑重构 (History Restoration)

### 3.1 原则

- 保持 **reasoning-zone → process-log-zone → text-zone → cards-zone** 骨架不变。
- 仅根据**扩充后的 state_snapshot** 在对应 zone 内用 **createElement / appendChild** 或已有静态渲染函数填内容，**禁止**对整泡或 reasoning-zone 做 innerHTML 覆盖。

### 3.2 具体改造点

| 位置 | 当前行为 | 改造后 |
|------|----------|--------|
| **createRestoreAIMessageRowWithSnapshot** (约 4181) | footer 仅 action-bar，无计时器 | 若 `snap.duration != null`，在 footer 中插入 `<span class="generation-timer"><i class="bi bi-stopwatch"></i> {duration}s</span>`（与直播 footer 结构一致） |
| **reasoning-zone** | 未填充 | 若 `snap.reasoning` 存在且非空，在 reasoning-zone 内挂载与直播一致的“思考过程”折叠块（如 `.agent-thought-block`），内容为 `snap.reasoning`（纯文本或安全 HTML），不覆盖整泡 |
| **process-log-zone** | 仅调用 `renderStaticChecklistToLogZone(snap.steps)` | 若 `snap.process_log` 存在且长度>0：用 `process_log` 逐条生成与直播一致的 checklist 行（content + state 对应图标），append 到 logZone；否则**回退**为仅用 `snap.steps` 渲染静态 checklist（保持旧记录兼容） |
| **renderStaticChecklistToLogZone** | 副标题用 `step.duration` 或 `step.summary` | 已支持 `step.duration`，无需改逻辑；后端补全 `steps[].duration` 后即可显示“耗时 Xs” |
| **renderStaticSnapshotReport** (约 3893) | 只渲染 diagnosis、summary、steps_details、topImages | 若 `reportData.report_data.suggestions` 存在且为数组，在报告底部增加“后续建议”区块并渲染 chips；若 `reportData.report_data.evaluation` 存在，增加“质量评估”区块（得分+评语） |
| **工作流卡片 / 报告区** | 已按 snap.workflow、snap.report 恢复 | 不变；仅保证 snapshot 中 workflow/report 已含 suggestions/evaluation（当前已有） |

### 3.3 旧记录兼容

- `snap.duration`、`snap.reasoning`、`snap.process_log` 缺省：对应区块不渲染或使用默认（如无 timer、无思考区、执行记录仅用 steps）。
- `steps[].duration` 缺省：`renderStaticChecklistToLogZone` 已有 `step.duration != null ? ... : (step.summary || '')`，无需改。
- `report.report_data.suggestions` / `evaluation` 缺省：不渲染建议/评估区块，不抛错。

---

## 四、实施顺序建议

1. **后端**  
   - Orchestrator：在 yield `state_snapshot` 前写 `state_snapshot["duration"]`；在 status/done 路径中维护并写入 `state_snapshot["process_log"]`；thought 时可选写入 `state_snapshot["reasoning"]`（并视产品决定是否从 `text` 中剥离 thought）。  
   - Executor：在每步 append 前写 `step_detail["duration"] = round(time.time() - tool_start, 1)`。  
2. **前端**  
   - 历史恢复：footer 显示 `snap.duration`；若有 `snap.reasoning` 则恢复 reasoning-zone；若有 `snap.process_log` 则用其渲染 process-log-zone，否则回退 steps。  
   - 静态报告：在 `renderStaticSnapshotReport` 中增加 suggestions、evaluation 的渲染。  
3. **兼容**  
   - 所有新字段读取处做存在性判断，避免旧数据导致报错或白屏。

---

## 五、冗余与边界说明

- **不存 DOM**：process_log 只存 `{ content, state }`，前端用其**重新生成**与直播一致的 DOM 节点，不存 HTML 字符串。  
- **report_data 与 steps**：当前 `steps` 与 `report.report_data.steps_details` 存在重复；保留现状即可（恢复时 checklist 用 steps，报告用手风琴用 report_data.steps_details），不在此次扩充中合并，避免破坏现有逻辑。  
- **_start_time**：入库前可保留（便于调试）或删除；若保留，前端不使用其计算 duration，仅使用显式的 `duration` 字段。

请审查上述遗漏清单与 JSON 扩充方案，确认无冗余、且旧记录兼容策略无误后，再按该草案实施代码修改。
