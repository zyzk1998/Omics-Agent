# SSE 事件类型

## 事件类型列表

| 事件类型 | 说明 | 数据格式 |
|---------|------|---------|
| `status` | 状态更新 | `{ "content": "状态消息", "state": "状态值" }` |
| `message` | 文本消息 | `{ "content": "消息内容" }` |
| `workflow` | 工作流配置 | `{ "workflow_config": {...}, "template_mode": true/false }` |
| `step_result` | 步骤执行结果 | `{ "report_data": {...} }` |
| `diagnosis` | 诊断报告 | `{ "report_data": {...} }` |
| `result` | 最终结果 | `{ "report_data": {...} }` 或 `{ "workflow_config": {...} }` |
| `done` | 完成信号 | `{ "status": "success" }` |
| `error` | 错误信息 | `{ "error": "错误描述", "message": "用户友好的错误消息" }` |

## 状态值 (state)

`status` 事件中的 `state` 字段可能的值：

- `"start"`: 开始处理
- `"analyzing"`: 正在分析
- `"thinking"`: 正在思考
- `"running"`: 正在执行
- `"rendering"`: 正在渲染
- `"generating_report"`: 正在生成报告
- `"completed"`: 执行完成
- `"error"`: 发生错误
- `"async_job_started"`: 异步作业已启动
- `"waiting"`: 等待中

---

**返回**: [SSE 流式响应目录](../README.md) | [API 手册首页](../../README.md)
