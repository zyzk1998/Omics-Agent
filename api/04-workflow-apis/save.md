# 保存工作流

## `POST /api/workflows/save`

**说明**: 保存工作流（书签）

**请求格式**: `application/json`

**请求体**:

```typescript
interface WorkflowSaveRequest {
  name: string;                        // 工作流名称
  workflow_json: Record<string, any>;  // 工作流 JSON
  user_id?: string;                    // 用户ID（可选）
}
```

**响应示例**:

```json
{
  "status": "success",
  "workflow_id": 123,
  "message": "工作流 'My Workflow' 已保存"
}
```

---

**返回**: [工作流 API 目录](../README.md) | [API 手册首页](../../README.md)
