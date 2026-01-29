# 执行工作流接口

## `POST /api/execute`

**说明**: 直接执行工作流（不通过聊天接口）

**请求格式**: `application/json`

**请求体**:

```typescript
interface ExecuteRequest {
  workflow_data: {
    workflow_name: string;
    steps: Array<{
      step_id: string;
      tool_id: string;
      name: string;
      params: Record<string, any>;
    }>;
  };
  file_paths: string[];  // 文件路径数组（相对路径或绝对路径）
}
```

**请求示例**:

```json
{
  "workflow_data": {
    "workflow_name": "Metabolomics Analysis Pipeline",
    "steps": [
      {
        "step_id": "inspect_data",
        "tool_id": "inspect_data",
        "name": "数据检查",
        "params": {
          "file_path": "example.csv"
        }
      },
      {
        "step_id": "preprocess_data",
        "tool_id": "preprocess_data",
        "name": "数据预处理",
        "params": {
          "file_path": "example.csv",
          "missing_threshold": "0.5",
          "normalization": "log2",
          "scale": "true"
        }
      }
    ]
  },
  "file_paths": ["guest/20250128_120000/example.csv"]
}
```

**成功响应** (200 OK):

```json
{
  "type": "analysis_report",
  "status": "success",
  "report_data": {
    "status": "success",
    "workflow_name": "Metabolomics Analysis Pipeline",
    "steps_details": [...],
    "steps_results": [...],
    "final_plot": "/results/run_20250128_120000/pca_plot.png",
    "output_dir": "/app/results/run_20250128_120000",
    "diagnosis": "## AI 专家分析报告\n\n..."
  },
  "reply": "✅ 工作流执行完成",
  "thought": "[THOUGHT] 使用 ToolRegistry 动态执行"
}
```

**错误响应** (500 Internal Server Error):

```json
{
  "status": "error",
  "error": "ValueError: No input files provided",
  "message": "工作流执行失败: No input files provided"
}
```

---

**返回**: [核心 API 目录](../README.md) | [API 手册首页](../../README.md)
