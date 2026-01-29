# SSE 流式响应格式

当 `stream: true` 时，`/api/chat` 接口返回 Server-Sent Events (SSE) 格式的流式响应。

## 文档列表

1. [事件类型](event-types.md) - SSE 事件类型详解
2. [事件格式](event-format.md) - SSE 事件格式说明
3. [前端处理](frontend-handling.md) - 前端处理 SSE 流的完整示例

---

## 快速参考

### SSE 事件类型

| 事件类型 | 说明 |
|---------|------|
| `status` | 状态更新 |
| `message` | 文本消息 |
| `workflow` | 工作流配置 |
| `step_result` | 步骤执行结果 |
| `diagnosis` | 诊断报告 |
| `result` | 最终结果 |
| `done` | 完成信号 |
| `error` | 错误信息 |

### SSE 事件格式

```
event: {event_type}
data: {json_data}

```

---

**返回**: [API 手册首页](../README.md)
