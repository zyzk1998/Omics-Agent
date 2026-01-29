# 列出工作流

## `GET /api/workflows/list`

**说明**: 列出用户的所有工作流（书签）

**请求参数**:
- `user_id` (string, 可选): 用户ID，默认 "guest"

**响应示例**:

```json
{
  "status": "success",
  "workflows": [
    {
      "id": 123,
      "name": "My Workflow",
      "workflow_json": { ... },
      "created_at": "2025-01-28T12:00:00"
    }
  ],
  "count": 1
}
```

---

**返回**: [工作流 API 目录](../README.md) | [API 手册首页](../../README.md)
