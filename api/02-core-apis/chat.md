# 聊天接口

## `POST /api/chat`

**说明**: 处理用户查询，支持多种响应类型（流式/JSON）

**请求格式**: `application/json`

**请求体**:

```typescript
interface ChatRequest {
  message: string;                    // 用户消息（可为空，如果有文件）
  history?: Array<{                   // 对话历史（可选）
    role: "user" | "assistant";
    content: string;
  }>;
  uploaded_files?: Array<{             // 已上传的文件列表（可选）
    name: string;                      // 文件名
    path: string;                     // 文件路径（相对路径或绝对路径）
  }>;
  workflow_data?: {                   // 工作流执行数据（可选）
    workflow_name: string;
    steps: Array<{
      step_id: string;
      tool_id: string;
      name: string;
      params: Record<string, any>;
    }>;
    file_paths: string[];             // 文件路径数组（必需）
  };
  test_dataset_id?: string;           // 测试数据集 ID（可选）
  stream?: boolean;                    // 是否使用流式响应，默认 false
  session_id?: string;                 // 会话ID（可选）
  user_id?: string;                    // 用户ID（可选，默认 "guest"）
}
```

**请求示例**:

```json
{
  "message": "分析这个文件",
  "uploaded_files": [
    {
      "name": "example.csv",
      "path": "guest/20250128_120000/example.csv"
    }
  ],
  "stream": true,
  "user_id": "guest",
  "session_id": "20250128_120000"
}
```

**响应类型**: 根据 `stream` 参数和 `Content-Type` 判断

### JSON 响应（非流式，`stream: false`）

**Content-Type**: `application/json`

**响应类型**: 根据 `type` 字段判断

#### 工作流配置响应

```json
{
  "type": "workflow_config",
  "workflow_data": {
    "workflow_name": "Metabolomics Analysis Pipeline",
    "steps": [...]
  },
  "file_paths": ["guest/20250128_120000/example.csv"],
  "diagnosis_report": "数据质量评估报告...",
  "recommendation": "推荐使用 log2 标准化..."
}
```

#### 分析报告响应

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
  }
}
```

### 流式响应（SSE，`stream: true`）

**Content-Type**: `text/event-stream`

**格式**: Server-Sent Events (SSE)

详细说明请参考 [SSE 流式响应格式](../05-sse-streaming/README.md)

**前端处理示例**:

```javascript
const response = await fetch('/api/chat', {
  method: 'POST',
  headers: { 'Content-Type': 'application/json' },
  body: JSON.stringify({
    message: '分析这个文件',
    uploaded_files: [{ name: 'example.csv', path: 'example.csv' }],
    stream: true
  })
});

const reader = response.body.getReader();
const decoder = new TextDecoder();
let buffer = '';

while (true) {
  const { done, value } = await reader.read();
  if (done) break;
  
  buffer += decoder.decode(value, { stream: true });
  // 解析 SSE 事件...
}
```

---

**返回**: [核心 API 目录](../README.md) | [SSE 流式响应](../05-sse-streaming/README.md) | [API 手册首页](../../README.md)
