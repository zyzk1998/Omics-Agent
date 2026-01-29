# 任务历史

## `GET /api/jobs/history`

**说明**: 获取任务执行历史

**请求参数**:
- `user_id` (string, 可选): 用户ID，默认 "guest"
- `status` (string, 可选): 任务状态过滤（如 "success", "failed", "running"）
- `limit` (int, 可选): 返回的任务数量，默认 50

**响应示例**:

```json
{
  "status": "success",
  "jobs": [
    {
      "id": 456,
      "user_id": "guest",
      "workflow_name": "Metabolomics Analysis Pipeline",
      "status": "success",
      "created_at": "2025-01-28T12:00:00",
      "completed_at": "2025-01-28T12:05:00"
    }
  ],
  "count": 1
}
```

---

**返回**: [工作流 API 目录](../README.md) | [API 手册首页](../../README.md)
