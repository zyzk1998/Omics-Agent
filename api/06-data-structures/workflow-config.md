# 工作流配置 (WorkflowConfig)

```typescript
interface WorkflowStep {
  step_id: string;                    // 步骤ID
  tool_id: string;                    // 工具ID
  name: string;                       // 步骤名称
  params: Record<string, any>;        // 步骤参数
}

interface WorkflowConfig {
  workflow_name: string;               // 工作流名称
  steps: WorkflowStep[];               // 步骤列表
  file_paths?: string[];                // 文件路径数组（可选）
}
```

## 使用场景

- 工作流规划接口返回的工作流配置
- 工作流执行接口的请求参数
- SSE 事件中的 `workflow` 事件

---

**返回**: [数据结构目录](../README.md) | [API 手册首页](../../README.md)
