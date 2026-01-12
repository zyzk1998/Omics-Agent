# GIBH-AGENT-V2 API 文档

**版本**: v2.0  
**基础 URL**: `http://localhost:8028` (开发环境)  
**协议**: HTTP/1.1  
**数据格式**: JSON (除文件上传外)

---

## 📋 目录

1. [通用说明](#通用说明)
2. [API 端点列表](#api-端点列表)
3. [详细接口文档](#详细接口文档)
4. [数据结构定义](#数据结构定义)
5. [错误处理](#错误处理)
6. [使用示例](#使用示例)

---

## 通用说明

### 请求头

所有 JSON 请求需要设置：
```
Content-Type: application/json
```

文件上传请求使用：
```
Content-Type: multipart/form-data
```

### 响应格式

#### 成功响应

```json
{
  "status": "success",
  "data": { ... }
}
```

#### 错误响应

```json
{
  "status": "error",
  "error": "错误描述",
  "message": "用户友好的错误消息"
}
```

### HTTP 状态码

- `200 OK`: 请求成功
- `400 Bad Request`: 请求参数错误
- `404 Not Found`: 资源不存在
- `413 Payload Too Large`: 文件大小超限
- `500 Internal Server Error`: 服务器内部错误

---

## API 端点列表

| 方法 | 路径 | 说明 |
|------|------|------|
| `GET` | `/` | 返回前端 HTML 页面 |
| `POST` | `/api/upload` | 文件上传（支持多文件） |
| `POST` | `/api/chat` | 聊天接口（支持流式响应） |
| `POST` | `/api/execute` | 执行工作流 |
| `GET` | `/api/logs/stream` | 实时日志流（SSE） |
| `GET` | `/api/logs` | 获取历史日志 |
| `GET` | `/api/workflow/status/{run_id}` | 查询工作流状态 |

---

## 详细接口文档

### 1. 文件上传接口

**端点**: `POST /api/upload`

**说明**: 上传一个或多个文件，支持 10x Genomics 数据（自动识别并分组）

**请求格式**: `multipart/form-data`

**请求参数**:
- `files`: `File[]` - 文件列表（支持多文件上传，最多 20 个）

**支持的文件类型**:
- `.h5ad` - AnnData 格式（单细胞数据）
- `.mtx` - Matrix Market 格式
- `.tsv`, `.csv` - 表格数据
- `.txt` - 文本文件
- `.gz`, `.tar`, `.zip` - 压缩文件

**文件大小限制**: 默认 100MB（可通过环境变量 `MAX_FILE_SIZE` 配置）

**响应格式**:

#### 单个文件上传成功

```json
{
  "status": "success",
  "file_id": "example.csv",
  "file_name": "example.csv",
  "file_path": "/path/to/uploads/example.csv",
  "file_size": 1024,
  "metadata": {
    "file_type": "csv",
    "n_samples": 100,
    "n_features": 50
  },
  "is_10x": false,
  "file_paths": ["example.csv"],
  "file_info": [
    {
      "name": "example.csv",
      "size": 1024,
      "path": "example.csv"
    }
  ],
  "count": 1
}
```

#### 多个文件上传成功

```json
{
  "status": "success",
  "file_paths": ["file1.csv", "file2.csv"],
  "file_info": [
    {
      "name": "file1.csv",
      "size": 1024,
      "path": "file1.csv"
    },
    {
      "name": "file2.csv",
      "size": 2048,
      "path": "file2.csv"
    }
  ],
  "count": 2
}
```

#### 10x Genomics 数据上传成功

```json
{
  "status": "success",
  "is_10x_data": true,
  "group_dir": "10x_data_20241201_120000",
  "files": [
    {
      "file_id": "10x_data_20241201_120000",
      "file_name": "matrix.mtx",
      "file_path": "/path/to/uploads/10x_data_20241201_120000/matrix.mtx",
      "file_size": 1024,
      "metadata": { ... },
      "is_10x": true,
      "group_dir": "10x_data_20241201_120000"
    },
    {
      "file_id": "10x_data_20241201_120000",
      "file_name": "barcodes.tsv",
      "file_path": "/path/to/uploads/10x_data_20241201_120000/barcodes.tsv",
      "file_size": 512,
      "metadata": { ... },
      "is_10x": true,
      "group_dir": "10x_data_20241201_120000"
    }
  ],
  "file_paths": ["10x_data_20241201_120000"],
  "message": "10x数据已保存到: 10x_data_20241201_120000"
}
```

#### 错误响应

```json
{
  "detail": "不允许的文件类型: .exe。允许的类型: .h5ad, .mtx, .tsv, .csv, .txt, .gz, .tar, .zip"
}
```

**状态码**: `400 Bad Request`

```json
{
  "detail": "文件 example.csv 超过最大大小限制 (100MB)"
}
```

**状态码**: `413 Payload Too Large`

---

### 2. 聊天接口

**端点**: `POST /api/chat`

**说明**: 处理用户查询，支持多种响应类型（流式/JSON）

**请求格式**: `application/json`

**请求体**:

```typescript
interface ChatRequest {
  message: string;                    // 用户消息（可为空，如果有文件）
  history?: Array<{                    // 对话历史（可选）
    role: "user" | "assistant";
    content: string;
  }>;
  uploaded_files?: Array<{             // 已上传的文件列表（可选）
    file_name?: string;                 // 文件名（兼容字段）
    file_path?: string;                 // 文件路径（兼容字段）
    name?: string;                      // 文件名（新字段）
    path?: string;                      // 文件路径（新字段）
  }>;
  workflow_data?: {                    // 工作流执行数据（可选）
    workflow_name: string;
    steps: Array<{
      step_id: string;
      tool_id: string;
      name: string;
      params: Record<string, any>;
    }>;
    file_paths: string[];              // 文件路径数组（必需）
  };
  test_dataset_id?: string;            // 测试数据集 ID（可选）
}
```

**请求示例**:

```json
{
  "message": "分析这个文件",
  "history": [],
  "uploaded_files": [
    {
      "name": "example.csv",
      "path": "example.csv"
    }
  ]
}
```

**响应类型**: 根据 `Content-Type` 判断

#### 2.1 JSON 响应（非流式）

**Content-Type**: `application/json`

**响应类型**: 根据 `type` 字段判断

##### 2.1.1 工作流配置响应

```json
{
  "type": "workflow_config",
  "workflow_data": {
    "workflow_name": "Metabolomics Analysis Pipeline",
    "steps": [
      {
        "step_id": "inspect_data",
        "tool_id": "inspect_data",
        "name": "数据检查",
        "step_name": "数据检查",
        "desc": "检查数据文件的基本信息",
        "params": {
          "file_path": "example.csv"
        }
      },
      {
        "step_id": "preprocess_data",
        "tool_id": "preprocess_data",
        "name": "数据预处理",
        "step_name": "数据预处理",
        "desc": "数据预处理：处理缺失值、标准化、缩放",
        "params": {
          "file_path": "example.csv",
          "missing_threshold": "0.5",
          "normalization": "log2",
          "scale": "true"
        }
      }
    ]
  },
  "file_paths": ["example.csv"],
  "recommendation": {
    "summary": "检测到数据包含 77 个样本。数值跨度较大 (0-10000+)。",
    "params": {
      "normalization": {
        "value": "log2",
        "reason": "数值跨度大，建议 Log 变换以符合正态分布"
      },
      "missing_threshold": {
        "value": "0.5",
        "reason": "标准质控阈值"
      },
      "scale": {
        "value": true,
        "reason": "标准化有助于后续分析"
      },
      "n_components": {
        "value": "10",
        "reason": "根据样本数推荐"
      }
    }
  }
}
```

##### 2.1.2 工具配置响应

```json
{
  "type": "tool_config",
  "reply": "请配置以下参数：",
  "tool": {
    "name": "inspect_data",
    "description": "检查数据文件",
    "parameters": [
      {
        "name": "file_path",
        "type": "string",
        "required": true,
        "description": "文件路径"
      }
    ]
  }
}
```

##### 2.1.3 工作流启动响应

```json
{
  "type": "workflow_started",
  "reply": "工作流已启动，正在执行...",
  "run_id": "run_20241201_120000"
}
```

##### 2.1.4 数据选择器响应

```json
{
  "type": "data_selector",
  "reply": "请选择数据集：",
  "datasets": [
    {
      "id": "pbmc_1k_v3",
      "name": "PBMC 1k v3",
      "description": "Peripheral Blood Mononuclear Cells",
      "size": "1.2 GB"
    }
  ]
}
```

##### 2.1.5 工具选择响应

```json
{
  "type": "choice",
  "reply": "请选择要使用的工具：",
  "candidates": [
    {
      "name": "inspect_data",
      "description": "检查数据文件",
      "tool_id": "inspect_data"
    },
    {
      "name": "preprocess_data",
      "description": "预处理数据",
      "tool_id": "preprocess_data"
    }
  ]
}
```

##### 2.1.6 分析报告响应

```json
{
  "type": "analysis_report",
  "status": "success",
  "report_data": {
    "status": "success",
    "workflow_name": "Metabolomics Analysis Pipeline",
    "steps_details": [
      {
        "step_id": "inspect_data",
        "tool_id": "inspect_data",
        "name": "数据检查",
        "summary": "检查完成: 77 个样本, 50 个代谢物",
        "status": "success",
        "step_result": {
          "step_name": "数据检查",
          "status": "success",
          "logs": "检查完成: 77 个样本, 50 个代谢物",
          "data": {
            "summary": {
              "n_samples": 77,
              "n_metabolites": 50,
              "missing_percentage": 2.5
            },
            "preview": [ ... ]
          }
        }
      }
    ],
    "steps_results": [ ... ],
    "final_plot": "/results/run_20241201_120000/pca_plot.png",
    "output_dir": "/path/to/results/run_20241201_120000",
    "diagnosis": "## 数据质量评估\n\n数据质量良好，缺失值比例较低（2.5%）...\n\n## 主要发现\n\n..."
  },
  "diagnosis": "## 数据质量评估\n\n..."
}
```

##### 2.1.7 错误响应

```json
{
  "type": "error",
  "error": "智能体未初始化，请检查配置和日志。",
  "message": "智能体初始化失败，请查看服务器日志获取详细信息"
}
```

**状态码**: `500 Internal Server Error`

#### 2.2 流式响应（SSE）

**Content-Type**: `text/event-stream` 或 `text/plain`

**格式**: Server-Sent Events (SSE) 或纯文本流

**流式内容可能包含**:

1. **思考过程**（可选）:
```
<think>
思考内容...
</think>
```

2. **最终回答**:
```
最终回答内容...
```

3. **数据集 JSON**（测试数据选择时）:
```
<!-- DATASETS_JSON: [{"id":"pbmc_1k_v3","name":"PBMC 1k v3",...}] -->
```

**前端处理示例**:

```javascript
const response = await fetch('/api/chat', {
  method: 'POST',
  headers: { 'Content-Type': 'application/json' },
  body: JSON.stringify({
    message: '分析这个文件',
    uploaded_files: [{ name: 'example.csv', path: 'example.csv' }]
  })
});

const contentType = response.headers.get('content-type');

if (contentType && contentType.includes('application/json')) {
  // JSON 响应
  const data = await response.json();
  handleJsonResponse(data);
} else {
  // 流式响应
  const reader = response.body.getReader();
  const decoder = new TextDecoder();
  
  while (true) {
    const { done, value } = await reader.read();
    if (done) break;
    
    const chunk = decoder.decode(value, { stream: true });
    handleStreamChunk(chunk);
  }
}
```

---

### 3. 执行工作流接口

**端点**: `POST /api/execute`

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
  "file_paths": ["example.csv"]
}
```

**响应格式**:

#### 成功响应

```json
{
  "type": "analysis_report",
  "status": "success",
  "report_data": {
    "status": "success",
    "workflow_name": "Metabolomics Analysis Pipeline",
    "steps_details": [
      {
        "step_id": "inspect_data",
        "tool_id": "inspect_data",
        "name": "数据检查",
        "summary": "检查完成: 77 个样本, 50 个代谢物",
        "status": "success",
        "plot": "/results/run_20241201_120000/inspect_plot.png",
        "step_result": {
          "step_name": "数据检查",
          "status": "success",
          "logs": "检查完成: 77 个样本, 50 个代谢物",
          "data": {
            "summary": {
              "n_samples": 77,
              "n_metabolites": 50,
              "missing_percentage": 2.5
            },
            "preview": [ ... ],
            "images": ["/results/run_20241201_120000/inspect_plot.png"]
          }
        }
      }
    ],
    "steps_results": [ ... ],
    "final_plot": "/results/run_20241201_120000/pca_plot.png",
    "output_dir": "/path/to/results/run_20241201_120000",
    "diagnosis": "## 数据质量评估\n\n..."
  }
}
```

#### 错误响应

```json
{
  "status": "error",
  "error": "ValueError: No input files provided",
  "error_detail": "Traceback (most recent call last):\n...",
  "message": "工作流执行失败: ValueError: No input files provided"
}
```

**状态码**: `500 Internal Server Error`

---

### 4. 实时日志流接口

**端点**: `GET /api/logs/stream`

**说明**: 使用 Server-Sent Events (SSE) 实时推送日志

**请求参数**: 无

**响应格式**: `text/event-stream`

**响应示例**:

```
data: {"timestamp": "2024-12-01T12:00:00", "level": "INFO", "message": "工作流执行开始"}

data: {"timestamp": "2024-12-01T12:00:01", "level": "INFO", "message": "步骤 1/6: 数据检查"}

data: {"timestamp": "2024-12-01T12:00:02", "level": "INFO", "message": "步骤 1 完成"}
```

**前端处理示例**:

```javascript
const eventSource = new EventSource('/api/logs/stream');

eventSource.onmessage = (event) => {
  const log = JSON.parse(event.data);
  console.log(`[${log.level}] ${log.message}`);
};

eventSource.onerror = (error) => {
  console.error('日志流连接错误:', error);
  eventSource.close();
};
```

---

### 5. 获取历史日志接口

**端点**: `GET /api/logs`

**说明**: 获取最近的历史日志

**请求参数**:
- `limit`: `number` (可选，默认 100) - 返回的日志条数

**请求示例**:

```
GET /api/logs?limit=50
```

**响应格式**:

```json
{
  "logs": [
    {
      "timestamp": "2024-12-01T12:00:00",
      "level": "INFO",
      "message": "工作流执行开始"
    },
    {
      "timestamp": "2024-12-01T12:00:01",
      "level": "INFO",
      "message": "步骤 1/6: 数据检查"
    }
  ],
  "total": 1000
}
```

---

### 6. 查询工作流状态接口

**端点**: `GET /api/workflow/status/{run_id}`

**说明**: 查询工作流执行状态（如果使用 Celery 异步执行）

**路径参数**:
- `run_id`: `string` - 工作流运行 ID

**请求示例**:

```
GET /api/workflow/status/run_20241201_120000
```

**响应格式**:

#### 运行中

```json
{
  "status": "running",
  "completed": false,
  "steps_status": [],
  "error": null
}
```

#### 执行中（有进度）

```json
{
  "status": "running",
  "completed": false,
  "steps_status": [
    {
      "step_id": "inspect_data",
      "status": "success",
      "summary": "检查完成"
    },
    {
      "step_id": "preprocess_data",
      "status": "running",
      "summary": "预处理中..."
    }
  ],
  "error": null
}
```

#### 执行成功

```json
{
  "status": "success",
  "completed": true,
  "steps_status": [ ... ],
  "report_data": {
    "status": "success",
    "workflow_name": "Metabolomics Analysis Pipeline",
    "steps_details": [ ... ],
    "diagnosis": "..."
  },
  "error": null
}
```

#### 执行失败

```json
{
  "status": "failed",
  "completed": true,
  "steps_status": [ ... ],
  "error": "ValueError: No input files provided"
}
```

#### 未找到

```json
{
  "status": "not_found",
  "message": "工作流未找到或未使用异步执行"
}
```

**状态码**: `404 Not Found`

---

## 数据结构定义

### 文件信息

```typescript
interface FileInfo {
  file_id?: string;           // 文件 ID（单个文件时）
  file_name: string;           // 文件名
  file_path: string;           // 文件路径（绝对路径）
  file_size: number;          // 文件大小（字节）
  metadata?: {                 // 文件元数据（可选）
    file_type: string;         // 文件类型（csv, h5ad, etc.）
    n_samples?: number;        // 样本数（如果可检测）
    n_features?: number;       // 特征数（如果可检测）
    [key: string]: any;        // 其他元数据字段
  };
  is_10x?: boolean;            // 是否为 10x Genomics 数据
  group_dir?: string;          // 10x 数据组目录（如果是 10x 数据）
}
```

### 工作流步骤

```typescript
interface WorkflowStep {
  step_id: string;             // 步骤 ID（唯一标识）
  tool_id: string;             // 工具 ID
  name: string;                 // 步骤名称（显示用）
  step_name?: string;           // 步骤名称（兼容字段）
  desc?: string;                // 步骤描述
  params: Record<string, any>;  // 步骤参数（键值对）
}
```

### 工作流配置

```typescript
interface WorkflowConfig {
  workflow_name: string;        // 工作流名称
  steps: WorkflowStep[];        // 步骤列表
}
```

### 步骤执行结果

```typescript
interface StepResult {
  step_name: string;            // 步骤名称
  status: "success" | "error" | "running";  // 状态
  logs: string;                 // 日志信息
  data?: {                      // 步骤数据（可选）
    summary?: Record<string, any>;  // 摘要信息
    preview?: any[];            // 预览数据
    images?: string[];           // 图片路径数组
    [key: string]: any;         // 其他数据字段
  };
}
```

### 步骤详情

```typescript
interface StepDetail {
  step_id: string;              // 步骤 ID
  tool_id: string;              // 工具 ID
  name: string;                 // 步骤名称
  summary: string;              // 摘要
  status: "success" | "error" | "running";  // 状态
  plot?: string;                // 图片路径（如果有）
  step_result: StepResult;      // 完整步骤结果
}
```

### 分析报告

```typescript
interface AnalysisReport {
  status: "success" | "error";  // 状态
  workflow_name: string;         // 工作流名称
  steps_details: StepDetail[];   // 步骤详情列表
  steps_results?: StepResult[];  // 步骤结果列表（新格式）
  final_plot?: string;           // 最终图片路径
  output_dir: string;            // 输出目录
  diagnosis?: string;             // AI 诊断报告（Markdown 格式）
}
```

### AI 推荐

```typescript
interface Recommendation {
  summary: string;               // 数据摘要
  params: {                      // 参数推荐
    [paramName: string]: {
      value: string | number | boolean;  // 推荐值
      reason: string;             // 推荐理由
    };
  };
}
```

---

## 错误处理

### 错误响应格式

所有错误响应遵循统一格式：

```json
{
  "status": "error",
  "error": "错误类型: 错误描述",
  "message": "用户友好的错误消息",
  "error_detail": "详细错误信息（可选，开发环境）"
}
```

### 常见错误码

| HTTP 状态码 | 错误类型 | 说明 |
|------------|---------|------|
| `400` | `Bad Request` | 请求参数错误 |
| `404` | `Not Found` | 资源不存在 |
| `413` | `Payload Too Large` | 文件大小超限 |
| `500` | `Internal Server Error` | 服务器内部错误 |

### 错误处理示例

```javascript
try {
  const response = await fetch('/api/chat', {
    method: 'POST',
    headers: { 'Content-Type': 'application/json' },
    body: JSON.stringify(requestData)
  });
  
  if (!response.ok) {
    const errorData = await response.json();
    throw new Error(errorData.message || errorData.error);
  }
  
  const data = await response.json();
  // 处理成功响应
} catch (error) {
  console.error('请求失败:', error);
  // 显示错误消息给用户
}
```

---

## 使用示例

### 完整工作流示例

#### 1. 上传文件

```javascript
// 上传文件
const formData = new FormData();
formData.append('files', fileInput.files[0]);

const uploadResponse = await fetch('/api/upload', {
  method: 'POST',
  body: formData
});

const uploadData = await uploadResponse.json();
console.log('文件上传成功:', uploadData.file_paths);
```

#### 2. 生成工作流配置

```javascript
// 发送聊天请求，生成工作流配置
const chatResponse = await fetch('/api/chat', {
  method: 'POST',
  headers: { 'Content-Type': 'application/json' },
  body: JSON.stringify({
    message: '分析这个文件',
    uploaded_files: [
      {
        name: uploadData.file_name,
        path: uploadData.file_paths[0]
      }
    ]
  })
});

const chatData = await chatResponse.json();

if (chatData.type === 'workflow_config') {
  console.log('工作流配置:', chatData.workflow_data);
  console.log('AI 推荐:', chatData.recommendation);
}
```

#### 3. 执行工作流

```javascript
// 方式1: 通过聊天接口执行（推荐）
const executeResponse = await fetch('/api/chat', {
  method: 'POST',
  headers: { 'Content-Type': 'application/json' },
  body: JSON.stringify({
    message: '执行工作流',
    workflow_data: {
      workflow_name: chatData.workflow_data.workflow_name,
      steps: chatData.workflow_data.steps,
      file_paths: chatData.file_paths
    }
  })
});

// 方式2: 直接调用执行接口
const executeResponse = await fetch('/api/execute', {
  method: 'POST',
  headers: { 'Content-Type': 'application/json' },
  body: JSON.stringify({
    workflow_data: chatData.workflow_data,
    file_paths: chatData.file_paths
  })
});

const executeData = await executeResponse.json();

if (executeData.type === 'analysis_report') {
  console.log('分析报告:', executeData.report_data);
  console.log('AI 诊断:', executeData.report_data.diagnosis);
}
```

#### 4. 处理流式响应

```javascript
const response = await fetch('/api/chat', {
  method: 'POST',
  headers: { 'Content-Type': 'application/json' },
  body: JSON.stringify({
    message: '解释这个文件',
    uploaded_files: [{ name: 'example.csv', path: 'example.csv' }]
  })
});

const reader = response.body.getReader();
const decoder = new TextDecoder();
let fullText = '';
let thinkBuffer = '';
let isThinking = false;

while (true) {
  const { done, value } = await reader.read();
  if (done) break;
  
  const chunk = decoder.decode(value, { stream: true });
  fullText += chunk;
  
  // 检测思考过程
  if (chunk.includes('<think>')) {
    isThinking = true;
    // 显示思考过程 UI
  }
  
  if (chunk.includes('</think>')) {
    isThinking = false;
    const parts = chunk.split('</think>');
    thinkBuffer += parts[0].replace('<think>', '');
    // 更新思考过程 UI
  }
  
  if (isThinking) {
    thinkBuffer += chunk.replace('<think>', '');
    // 更新思考过程 UI
  } else {
    // 更新最终回答 UI
  }
}
```

### TypeScript 类型定义

```typescript
// 请求类型
interface ChatRequest {
  message: string;
  history?: Array<{ role: string; content: string }>;
  uploaded_files?: Array<{ name: string; path: string }>;
  workflow_data?: {
    workflow_name: string;
    steps: WorkflowStep[];
    file_paths: string[];
  };
  test_dataset_id?: string;
}

interface ExecuteRequest {
  workflow_data: {
    workflow_name: string;
    steps: WorkflowStep[];
  };
  file_paths: string[];
}

// 响应类型
type ChatResponse = 
  | WorkflowConfigResponse
  | ToolConfigResponse
  | WorkflowStartedResponse
  | DataSelectorResponse
  | ChoiceResponse
  | AnalysisReportResponse
  | ErrorResponse
  | StreamResponse;

interface WorkflowConfigResponse {
  type: 'workflow_config';
  workflow_data: WorkflowConfig;
  file_paths: string[];
  diagnosis_report?: string;  // Markdown 格式的数据诊断报告（所有 Agent 统一生成）
  recommendation?: Recommendation;  // 参数推荐（Metabolomics Agent 特有）
}

interface AnalysisReportResponse {
  type: 'analysis_report';
  status: 'success';
  report_data: AnalysisReport;
  diagnosis?: string;
}

interface ErrorResponse {
  type: 'error';
  error: string;
  message: string;
}
```

---

## 注意事项

### 1. 文件路径

- **上传接口返回**: 相对路径（相对于 `UPLOAD_DIR`）
- **聊天接口使用**: 相对路径或绝对路径均可
- **结果图片路径**: 以 `/results/` 开头的 URL 路径（前端可直接使用）

### 2. 工作流执行

- **推荐方式**: 通过 `/api/chat` 接口执行（自动处理文件路径）
- **直接执行**: 使用 `/api/execute` 接口（需要确保文件路径正确）

### 3. 流式响应

- **Content-Type**: `text/event-stream` 或 `text/plain`
- **编码**: UTF-8
- **格式**: SSE 或纯文本流

### 4. 错误处理

- **始终检查**: `response.ok` 或 `response.status`
- **解析错误**: 使用 `response.json()` 获取错误详情
- **用户提示**: 使用 `message` 字段显示给用户

### 5. 文件大小限制

- **默认限制**: 100MB
- **配置方式**: 环境变量 `MAX_FILE_SIZE`（字节）
- **错误处理**: 413 状态码 + 错误消息

### 6. 前端调试功能

**调试侧边栏**:
- **触发方式**: 
  - 双击导航栏品牌 logo（"GIBH Qwen Agent"）
  - 或点击导航栏右侧 🐛 图标按钮
- **功能**: 自动捕获并美化显示所有 JSON 响应
- **样式**: 深色主题，固定右侧，可折叠/展开

**使用场景**:
- 调试 API 响应
- 查看完整 JSON 数据结构
- 排查前端/后端数据不一致问题

---

## 更新日志

### v2.0 (2024-12)

- ✅ 新增 AI 推荐系统（`recommendation` 字段）
- ✅ 新增 AI 诊断报告（`diagnosis` 字段）
- ✅ 优化工作流配置生成性能（轻量级预览）
- ✅ 新增调试可见性（前端调试面板 + 后端 JSON 监控）

---

**文档版本**: v2.1  
**最后更新**: 2024年12月（更新前端调试功能说明）  
**维护者**: GIBH-AGENT-V2 开发团队

