# 工作流状态

## `GET /api/workflow/status/{run_id}`

**说明**: 查询工作流状态（兼容旧架构，支持 Celery 异步任务）

**路径参数**:
- `run_id` (string, 必需): 运行ID（Celery 任务ID）

**响应示例**:

```json
{
  "status": "success",
  "completed": true,
  "steps_status": [
    {
      "step_id": "inspect_data",
      "status": "success"
    }
  ],
  "error": null
}
```

**状态值**:
- `"running"`: 正在运行
- `"success"`: 执行成功
- `"failed"`: 执行失败
- `"pending"`: 等待执行

---

**返回**: [工作流 API 目录](../README.md) | [API 手册首页](../../README.md)
