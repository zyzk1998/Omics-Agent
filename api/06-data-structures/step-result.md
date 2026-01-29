# 步骤结果 (StepResult)

```typescript
interface StepResult {
  step_id: string;                     // 步骤ID
  tool_id: string;                     // 工具ID
  name: string;                        // 步骤名称
  summary: string;                     // 步骤摘要
  status: "success" | "failed" | "warning";  // 步骤状态
  plot?: string;                       // 图表路径（可选）
  step_result: {
    step_name: string;
    status: string;
    logs: string;
    data: Record<string, any>;
  };
}
```

## 使用场景

- 工作流执行接口返回的步骤结果
- SSE 事件中的 `step_result` 事件

---

**返回**: [数据结构目录](../README.md) | [API 手册首页](../../README.md)
