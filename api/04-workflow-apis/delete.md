# 删除工作流

## `DELETE /api/workflows/{workflow_id}`

**说明**: 删除工作流

**路径参数**:
- `workflow_id` (int, 必需): 工作流ID

**请求参数**:
- `user_id` (string, 可选): 用户ID，默认 "guest"

**响应示例**:

```json
{
  "status": "success",
  "message": "工作流 123 已删除"
}
```

**错误响应** (404 Not Found):

```json
{
  "detail": "工作流不存在或无权删除"
}
```

---

**返回**: [工作流 API 目录](../README.md) | [API 手册首页](../../README.md)
