# 规划工作流

## `POST /api/workflows/plan`

**说明**: 规划工作流（plan-first：可以在没有文件的情况下生成工作流）

**请求格式**: `application/json`

**请求体**:

```typescript
interface WorkflowPlanRequest {
  query: string;                      // 用户查询
  file_metadata?: Record<string, any>; // 文件元数据（可选）
  user_id?: string;                    // 用户ID（可选）
}
```

**响应示例**:

```json
{
  "status": "success",
  "workflow": {
    "workflow_name": "Metabolomics Analysis Pipeline",
    "steps": [ ... ]
  },
  "user_id": "guest"
}
```

---

**返回**: [工作流 API 目录](../README.md) | [API 手册首页](../../README.md)
