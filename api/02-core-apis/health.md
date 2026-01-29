# 健康检查接口

## `GET /api/health`

**说明**: 检查 API 服务状态和组件初始化情况

**请求参数**: 无

**响应示例**:

```json
{
  "status": "ok",
  "service": "GIBH-AGENT-V2",
  "agent_initialized": true,
  "tool_retriever_initialized": true
}
```

**响应字段说明**:
- `status`: 服务状态（"ok" 表示正常）
- `service`: 服务名称
- `agent_initialized`: 智能体是否已初始化
- `tool_retriever_initialized`: 工具检索器是否已初始化

---

**返回**: [核心 API 目录](../README.md) | [API 手册首页](../../README.md)
