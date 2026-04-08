# GIBH Agent API 文档

**版本**: v2.1  
**基础 URL**: `http://localhost:8028` (开发环境)  
**协议**: HTTP/1.1  
**数据格式**: JSON (除文件上传外)

**近期更新**: 账户与资产中台（鉴权、会话/消息、数据资产、工作流收藏）、技能广场（分页/收藏/管理员）、双通道大模型路由（model_name 入参）。

---

## 📋 目录

1. [通用说明](#通用说明)
2. [核心 API 端点](#核心-api-端点)
3. [详细接口文档](#详细接口文档)
4. [鉴权与会话/资产/技能接口](#鉴权与会话资产技能接口)
5. [SSE 流式响应格式](#sse-流式响应格式)
6. [数据结构定义](#数据结构定义)
7. [错误处理](#错误处理)
8. [使用示例](#使用示例)
9. [前端集成指南](#前端集成指南)

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

**鉴权**（需登录或身份时）：
- 登录后：`Authorization: Bearer <access_token>`（来自 `POST /api/auth/login` 返回的 `access_token`）
- 游客：`X-Guest-UUID: <uuid>`（前端生成并持久化，用于会话/资产/工作流收藏隔离）

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
  "message": "用户友好的错误消息",
  "detail": "详细错误信息（开发环境）"
}
```

### HTTP 状态码

- `200 OK`: 请求成功
- `400 Bad Request`: 请求参数错误
- `403 Forbidden`: 权限不足（如文件路径不安全）
- `404 Not Found`: 资源不存在
- `413 Payload Too Large`: 文件大小超限（默认 100MB）
- `500 Internal Server Error`: 服务器内部错误
- `503 Service Unavailable`: 服务不可用（如工具检索器未初始化）

---

## 核心 API 端点

| 方法 | 路径 | 说明 | 响应类型 |
|------|------|------|----------|
| `GET` | `/` | 返回前端 HTML 页面 | HTML |
| `GET` | `/api/health` | API 健康检查 | JSON |
| `GET` | `/api/db/init` | 数据库建表/初始化（部署后可选调用） | JSON |
| `POST` | `/api/upload` | 文件上传（多文件；按 owner_id 隔离并写 Asset） | JSON |
| `POST` | `/api/chat` | 聊天（SSE 流式；支持 model_name 双通道模型路由） | SSE / JSON |
| `POST` | `/api/execute` | 执行工作流 | JSON |
| `GET` | `/api/logs/stream` | 实时日志流（SSE） | SSE |
| `GET` | `/api/logs` | 获取历史日志 | JSON |
| `GET` | `/api/tools/search` | 语义搜索工具 | JSON |
| `GET` | `/api/tools/list` | 列出所有工具 | JSON |
| `GET` | `/api/tools/{tool_name}` | 获取工具 Schema | JSON |
| `POST` | `/api/workflows/plan` | 规划工作流 | JSON |
| `POST` | `/api/workflows/save` | 保存工作流（书签） | JSON |
| `GET` | `/api/workflows/list` | 列出用户工作流 | JSON |
| `DELETE` | `/api/workflows/{workflow_id}` | 删除工作流 | JSON |
| `GET` | `/api/jobs/history` | 获取任务历史 | JSON |
| `GET` | `/api/workflow/status/{run_id}` | 查询工作流状态（可选 Celery） | JSON |
| **鉴权** | | | |
| `POST` | `/api/auth/register` | 用户注册 | JSON |
| `POST` | `/api/auth/login` | 登录（表单 username/password，返回 access_token） | JSON |
| `GET` | `/api/auth/me` | 当前用户信息（需 Bearer） | JSON |
| `POST` | `/api/auth/merge_guest_data` | 游客数据合并到当前用户（需 Bearer） | JSON |
| **会话与消息** | | | |
| `GET` | `/api/sessions` | 当前用户历史会话列表 | JSON |
| `GET` | `/api/sessions/{session_id}/messages` | 指定会话消息列表 | JSON |
| `DELETE` | `/api/messages/{message_id}` | 删除单条消息 | JSON |
| `DELETE` | `/api/sessions/{session_id}` | 删除会话及下属消息 | JSON |
| `PUT` | `/api/sessions/{session_id}` | 重命名会话（title） | JSON |
| **数据资产** | | | |
| `GET` | `/api/assets` | 当前用户数据资产列表 | JSON |
| `POST` | `/api/assets/reclassify` | 重新分类资产（模态等） | JSON |
| `DELETE` | `/api/assets/{asset_id}` | 删除资产记录 | JSON |
| `PUT` | `/api/assets/{asset_id}` | 重命名资产（file_name） | JSON |
| **工作流收藏** | | | |
| `GET` | `/api/workflow_templates` | 当前用户工作流收藏列表 | JSON |
| `POST` | `/api/workflow_templates` | 新建工作流收藏（name + config_json） | JSON |
| `DELETE` | `/api/workflow_templates/{template_id}` | 删除工作流收藏 | JSON |
| `PUT` | `/api/workflow_templates/{template_id}` | 重命名工作流收藏 | JSON |
| **技能广场** | | | |
| `GET` | `/api/skills` | 技能分页列表（main_cat/sub_cat/saved_only/page/size） | JSON |
| `POST` | `/api/skills/{skill_id}/bookmark` | 收藏技能 | JSON |
| `DELETE` | `/api/skills/{skill_id}/bookmark` | 取消收藏 | JSON |
| `POST` | `/api/admin/bootstrap-skills` | 管理员：重新注入系统技能 | JSON |

---

## 详细接口文档

### 1. 健康检查接口

#### `GET /api/health`

**说明**: 检查 API 服务状态和组件初始化情况

**请求参数**: 无

**响应示例**:

```json
{
  "status": "ok",
  "service": "GIBH-AGENT-V2",
  "agent_initialized": true,
  "tool_retriever_initialized": true
}
```

**响应字段说明**:
- `status`: 服务状态（"ok" 表示正常）
- `service`: 服务名称
- `agent_initialized`: 智能体是否已初始化
- `tool_retriever_initialized`: 工具检索器是否已初始化

---

### 2. 文件上传接口

#### `POST /api/upload`

**说明**: 上传一个或多个文件，支持 10x Genomics 数据（自动识别并分组）

**请求格式**: `multipart/form-data`

**请求参数**:
- `files` (File[], 必需): 文件列表（支持多文件上传，最多 20 个）
- `user_id` (string, 可选): 用户ID，默认 "guest"
- `session_id` (string, 可选): 会话ID，未提供时自动生成（格式: `YYYYMMDD_HHMMSS`）

**支持的文件类型**:
- `.h5ad` - AnnData 格式（单细胞数据）
- `.mtx` - Matrix Market 格式
- `.tsv`, `.csv` - 表格数据
- `.txt` - 文本文件
- `.gz`, `.tar`, `.zip` - 压缩文件

**文件大小限制**: 默认 100MB（可通过环境变量 `MAX_FILE_SIZE` 配置）

**10x Genomics 数据自动识别**:
- 如果上传的文件包含 `matrix.mtx`、`barcodes.tsv`、`features.tsv`（或 `genes.tsv`），系统会自动识别为 10x Genomics 数据
- 10x 数据会被保存到独立的子目录中（格式: `10x_data_YYYYMMDD_HHMMSS`）
- 返回的 `file_paths` 将指向该子目录，而不是单个文件

**成功响应** (200 OK):

```json
{
  "status": "success",
  "file_paths": [
    "guest/20250128_120000/example.csv",
    "guest/20250128_120000/10x_data_20250128_120000"
  ],
  "file_info": [
    {
      "name": "example.csv",
      "size": 1024000,
      "path": "guest/20250128_120000/example.csv"
    }
  ],
  "count": 2,
  "user_id": "guest",
  "session_id": "20250128_120000",
  "is_10x_data": true,
  "group_dir": "guest/20250128_120000/10x_data_20250128_120000",
  "files": [
    {
      "file_id": "example.csv",
      "file_name": "example.csv",
      "file_path": "/app/uploads/guest/20250128_120000/example.csv",
      "file_size": 1024000,
      "metadata": {
        "file_type": "csv",
        "n_samples": 100,
        "n_features": 50
      },
      "is_10x": false
    }
  ]
}
```

**响应字段说明**:
- `status`: 操作状态（"success" 表示成功）
- `file_paths`: 文件路径数组（相对路径，相对于 `/app/uploads`）
- `file_info`: 文件信息数组，包含 `name`、`size`、`path`
- `count`: 上传的文件数量
- `user_id`: 用户ID
- `session_id`: 会话ID
- `is_10x_data`: 是否为 10x Genomics 数据（仅当检测到 10x 数据时存在）
- `group_dir`: 10x 数据组目录路径（仅当检测到 10x 数据时存在）
- `files`: 文件详细信息数组（向后兼容字段）

**错误响应** (400 Bad Request):

```json
{
  "detail": "文件 example.csv 超过最大大小限制 (100MB)"
}
```

**错误响应** (413 Payload Too Large):

```json
{
  "detail": "文件 example.csv 超过最大大小限制 (100MB)"
}
```

**错误响应** (403 Forbidden):

```json
{
  "detail": "文件路径不安全：不允许访问基础目录外的文件"
}
```

**前端集成示例**:

```javascript
const formData = new FormData();
for (let file of fileInput.files) {
    formData.append('files', file);
}
formData.append('user_id', 'guest');
formData.append('session_id', '20250128_120000');

const response = await fetch('/api/upload', {
    method: 'POST',
    body: formData
});

const result = await response.json();
if (result.status === 'success') {
    console.log('上传成功:', result.file_paths);
    // 保存 file_paths 用于后续的聊天请求
    uploadedFiles = result.file_paths;
}
```

---

### 3. 聊天接口

#### `POST /api/chat`

**说明**: 处理用户查询，支持多种响应类型（流式/JSON）

**请求格式**: `application/json`

**请求体**:

```typescript
interface ChatRequest {
  message: string;                    // 用户消息（可为空，如果有文件）
  history?: Array<{ role: "user"|"assistant"; content: string }>;
  uploaded_files?: Array<{             // 已上传文件（支持 name/path 或 file_name/file_path）
    name?: string;  file_name?: string;
    path?: string;  file_path?: string;
  }>;
  workflow_data?: { workflow_name: string; steps: WorkflowStep[]; file_paths: string[]; };
  test_dataset_id?: string;
  stream?: boolean;                   // 默认 false；true 时返回 SSE
  session_id?: string;                // 新建会话时不传，由后端生成并落库
  user_id?: string;                   // 可选，默认由鉴权得到 owner_id
  model_name?: string;                // 双通道模型路由：如 "deepseek-ai/DeepSeek-R1"、"Pro/zai-org/GLM-5"、"Pro/moonshotai/Kimi-K2.5"、"Qwen/Qwen3.5-397B-A17B"（无斜杠的 qwen* 走 DashScope）
}
```

**请求示例**:

```json
{
  "message": "分析这个文件",
  "uploaded_files": [
    { "file_name": "example.h5ad", "file_path": "/app/uploads/owner1/batch_id/example.h5ad" }
  ],
  "stream": true,
  "session_id": "uuid-from-backend-or-new",
  "user_id": "owner_id",
  "model_name": "deepseek-ai/DeepSeek-R1"
}
```

**响应类型**: 根据 `stream` 参数和 `Content-Type` 判断

#### 3.1 JSON 响应（非流式，`stream: false`）

**Content-Type**: `application/json`

**响应类型**: 根据 `type` 字段判断

##### 3.1.1 工作流配置响应

当系统生成工作流计划时返回：

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
          "normalization": "log2"
        }
      }
    ]
  },
  "file_paths": ["guest/20250128_120000/example.csv"],
  "diagnosis_report": "数据质量评估报告...",
  "recommendation": "推荐使用 log2 标准化..."
}
```

##### 3.1.2 分析报告响应

当工作流执行完成时返回：

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
        "plot": "/results/run_20250128_120000/inspect_plot.png",
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
            "images": ["/results/run_20250128_120000/inspect_plot.png"]
          }
        }
      }
    ],
    "steps_results": [ ... ],
    "final_plot": "/results/run_20250128_120000/pca_plot.png",
    "output_dir": "/app/results/run_20250128_120000",
    "diagnosis": "## AI 专家分析报告\n\n..."
  }
}
```

##### 3.1.3 错误响应

```json
{
  "type": "error",
  "error": "错误描述",
  "message": "用户友好的错误消息"
}
```

#### 3.2 流式响应（SSE，`stream: true`）

**Content-Type**: `text/event-stream`

**格式**: Server-Sent Events (SSE)

**🔥 CRITICAL: Delta Streaming Protocol**

**重要**: 服务器发送的是 **Delta**（增量）token，而不是累积文本。客户端必须**追加**新内容到缓冲区。

**Delta 流式协议说明**:
- 每个 `message` 事件只包含**新的 token**（增量部分）
- 客户端必须使用 `+=` 操作符将新内容追加到现有缓冲区
- 示例：
  ```
  event: message
  data: {"content": "Ap"}
  
  event: message
  data: {"content": "ple"}  // 只包含新 token，不是 "Apple"
  ```
- 客户端处理方式：
  ```javascript
  let messageBuffer = '';
  
  // 收到 Delta token
  messageBuffer += data.content;  // ✅ 正确：追加
  // messageBuffer = data.content;  // ❌ 错误：替换
  ```

**DeepSeek Chain of Thought (CoT) 支持**:

当使用 DeepSeek-R1 等支持思考过程的模型时，响应可能包含 `<think>...</think>` 标签：

```
event: message
data: {"content": "<think>思考过程内容</think>最终答案"}
```

**客户端处理建议**:
1. 解析 `<think>` 标签，将思考过程和最终答案分离
2. 在 UI 中分别显示思考过程（可折叠）和最终答案
3. 支持标签跨多个 SSE chunk 的情况（标签可能被分割）

**SSE 事件类型**: 详见 [SSE 流式响应格式](#sse-流式响应格式)

**前端处理示例（Delta Streaming）**:

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
let sseBuffer = '';  // SSE 解析缓冲区
let messageBuffer = '';  // 消息内容缓冲区（累积 Delta）
let reasoningBuffer = '';  // 思考过程缓冲区
let isInReasoning = false;  // 是否在思考标签内

while (true) {
  const { done, value } = await reader.read();
  if (done) break;
  
  sseBuffer += decoder.decode(value, { stream: true });
  const lines = sseBuffer.split('\n');
  sseBuffer = lines.pop() || '';  // 保留不完整的行
  
  let currentEventType = null;
  
  for (const line of lines) {
    if (line.startsWith('event: ')) {
      currentEventType = line.substring(7).trim();
    } else if (line.startsWith('data: ')) {
      const dataStr = line.substring(6).trim();
      try {
        const data = JSON.parse(dataStr);
        
        // 🔥 CRITICAL: 处理 Delta token
        if (currentEventType === 'message' && data.content) {
          // 服务器发送的是 Delta，必须追加
          messageBuffer += data.content;
          
          // 解析 Chain of Thought 标签
          const parsed = parseReasoningTags(messageBuffer);
          if (parsed.reasoning) {
            reasoningBuffer = parsed.reasoning;
            messageBuffer = parsed.answer;
          }
          
          // 更新 UI
          updateChatUI(reasoningBuffer, messageBuffer);
        } else {
          handleSSEEvent(currentEventType, data);
        }
      } catch (e) {
        console.error('JSON 解析错误:', e, '数据:', dataStr);
      }
    }
  }
}

// 解析思考标签的辅助函数
function parseReasoningTags(content) {
  const reasoningMatch = content.match(/<think>(.*?)<\/redacted_reasoning>/s);
  if (reasoningMatch) {
    return {
      reasoning: reasoningMatch[1],
      answer: content.replace(reasoningMatch[0], '').trim()
    };
  }
  return { reasoning: null, answer: content };
}
```

---

### 4. 执行工作流接口

#### `POST /api/execute`

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
    "steps_details": [ ... ],
    "steps_results": [ ... ],
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

### 5. 日志流接口

#### `GET /api/logs/stream`

**说明**: 实时日志流（Server-Sent Events）

**请求参数**: 无

**响应格式**: `text/event-stream`

**SSE 事件格式**:

```
data: {"timestamp": "2025-01-28T12:00:00", "level": "INFO", "message": "日志内容", "module": "gibh_agent.core"}\n\n
```

**心跳事件** (保持连接):

```
data: {"type": "heartbeat", "timestamp": "2025-01-28T12:00:00"}\n\n
```

**前端集成示例**:

```javascript
const eventSource = new EventSource('/api/logs/stream');

eventSource.onmessage = function(event) {
  const logEntry = JSON.parse(event.data);
  console.log(`[${logEntry.level}] ${logEntry.message}`);
};

eventSource.onerror = function(error) {
  console.error('日志流错误:', error);
  eventSource.close();
};
```

---

### 6. 获取历史日志接口

#### `GET /api/logs`

**说明**: 获取历史日志

**请求参数**:
- `limit` (int, 可选): 返回的日志条数，默认 100

**响应示例**:

```json
{
  "logs": [
    {
      "timestamp": "2025-01-28T12:00:00",
      "level": "INFO",
      "message": "日志内容",
      "module": "gibh_agent.core"
    }
  ],
  "total": 1000
}
```

---

### 7. 工具检索接口

#### `GET /api/tools/search`

**说明**: 语义搜索工具（基于 ChromaDB + Embeddings）

**请求参数**:
- `query` (string, 必需): 查询文本（自然语言）
- `top_k` (int, 可选): 返回前 k 个最相关的工具，默认 5
- `category` (string, 可选): 类别过滤器（如 "Metabolomics", "scRNA-seq"）

**响应示例**:

```json
{
  "status": "success",
  "query": "数据预处理",
  "count": 3,
  "tools": [
    {
      "name": "preprocess_data",
      "description": "数据预处理工具",
      "category": "Metabolomics",
      "parameters": { ... }
    }
  ]
}
```

**错误响应** (503 Service Unavailable):

```json
{
  "detail": "工具检索器未初始化。请检查 Ollama 服务和依赖是否已安装。"
}
```

---

#### `GET /api/tools/list`

**说明**: 列出所有已注册的工具

**请求参数**: 无

**响应示例**:

```json
{
  "status": "success",
  "count": 15,
  "tools": [
    "inspect_data",
    "preprocess_data",
    "pca_analysis",
    ...
  ]
}
```

---

#### `GET /api/tools/{tool_name}`

**说明**: 获取特定工具的完整 Schema

**路径参数**:
- `tool_name` (string, 必需): 工具名称

**响应示例**:

```json
{
  "status": "success",
  "tool": {
    "name": "preprocess_data",
    "description": "数据预处理工具",
    "category": "Metabolomics",
    "parameters": {
      "file_path": {
        "type": "string",
        "description": "文件路径",
        "required": true
      },
      "missing_threshold": {
        "type": "string",
        "description": "缺失值阈值",
        "default": "0.5"
      }
    }
  }
}
```

**错误响应** (404 Not Found):

```json
{
  "detail": "工具 'unknown_tool' 不存在"
}
```

---

### 8. 工作流管理接口

#### `POST /api/workflows/plan`

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

#### `POST /api/workflows/save`

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

#### `GET /api/workflows/list`

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

#### `DELETE /api/workflows/{workflow_id}`

**说明**: 删除工作流

**路径参数**:
- `workflow_id` (int, 必需): 工作流ID

**请求参数**:
- `user_id` (string, 可选): 用户ID，默认 "guest"

**响应示例**:

```json
{
  "status": "success",
  "message": "工作流 123 已删除"
}
```

**错误响应** (404 Not Found):

```json
{
  "detail": "工作流不存在或无权删除"
}
```

---

### 9. 任务历史接口

#### `GET /api/jobs/history`

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

### 10. 工作流状态查询接口

#### `GET /api/workflow/status/{run_id}`

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

## 鉴权与会话/资产/技能接口

### 11. 数据库初始化

#### `GET /api/db/init`

**说明**: 建表与初始化（部署后若出现表不存在或技能为空可调用）

**请求参数**: 无

**成功响应** (200 OK):

```json
{
  "message": "数据库初始化成功",
  "tables_created": true
}
```

---

### 12. 鉴权接口

#### `POST /api/auth/register`

**说明**: 用户注册

**请求格式**: `application/json`

**请求体**:

```json
{
  "username": "myuser",
  "password": "mypassword"
}
```

**成功响应** (200 OK):

```json
{
  "username": "myuser",
  "message": "注册成功",
  "role": "user"
}
```

**错误响应** (400): `{"detail": "用户名已存在"}` 或 `"用户名和密码不能为空"`

---

#### `POST /api/auth/login`

**说明**: 登录，使用表单 `application/x-www-form-urlencoded`，返回 JWT

**请求体**: `username`, `password`（表单字段）

**成功响应** (200 OK):

```json
{
  "access_token": "eyJ...",
  "token_type": "bearer",
  "username": "myuser",
  "role": "user"
}
```

**错误响应** (401): `{"detail": "用户名或密码错误"}`

---

#### `GET /api/auth/me`

**说明**: 当前用户信息（需 Header: `Authorization: Bearer <access_token>`）

**成功响应** (200 OK):

```json
{
  "username": "myuser",
  "role": "user"
}
```

---

#### `POST /api/auth/merge_guest_data`

**说明**: 登录后将游客身份下的会话/资产/工作流收藏合并到当前用户

**请求体**:

```json
{
  "guest_uuid": "uuid-from-x-guest-uuid"
}
```

**成功响应** (200 OK): 合并结果描述

---

### 13. 会话与消息

#### `GET /api/sessions`

**说明**: 当前用户（owner_id）历史会话，按创建时间倒序；需鉴权或 X-Guest-UUID

**响应示例**:

```json
[
  {
    "id": "session-uuid",
    "owner_id": "user1",
    "title": "新会话",
    "created_at": "2025-03-13T12:00:00"
  }
]
```

---

#### `GET /api/sessions/{session_id}/messages`

**说明**: 指定会话的消息列表；校验会话归属

**响应示例**:

```json
[
  {
    "id": 1,
    "session_id": "session-uuid",
    "role": "user",
    "content": { "text": "..." },
    "created_at": "2025-03-13T12:00:00"
  }
]
```

**错误响应** (403/404): 无权或会话不存在

---

#### `DELETE /api/messages/{message_id}`

**说明**: 删除单条消息；校验消息所属会话归属

**成功响应** (200 OK): `{"status": "success"}` 或类似

---

#### `DELETE /api/sessions/{session_id}` / `PUT /api/sessions/{session_id}`

**说明**: 删除会话（及下属消息）；重命名会话（请求体 `{"title": "新标题"}`）

---

### 14. 数据资产

#### `GET /api/assets`

**说明**: 当前用户数据资产列表；需鉴权或 X-Guest-UUID

**响应示例**:

```json
[
  {
    "id": 1,
    "owner_id": "user1",
    "file_name": "data.h5ad",
    "file_path": "user1/batch/data.h5ad",
    "modality": "RNA",
    "created_at": "2025-03-13T12:00:00"
  }
]
```

---

#### `POST /api/assets/reclassify`

**说明**: 重新分类资产（模态等）

#### `DELETE /api/assets/{asset_id}` / `PUT /api/assets/{asset_id}`

**说明**: 删除资产；重命名（请求体 `{"file_name": "新文件名"}`）

---

### 15. 工作流收藏

#### `GET /api/workflow_templates`

**说明**: 当前用户工作流收藏列表

**响应示例**:

```json
[
  {
    "id": 1,
    "owner_id": "user1",
    "name": "我的流程",
    "config_json": { "workflow_name": "...", "steps": [...] },
    "created_at": "2025-03-13T12:00:00"
  }
]
```

---

#### `POST /api/workflow_templates`

**请求体**:

```json
{
  "name": "收藏名称",
  "config_json": { "workflow_name": "...", "steps": [...] }
}
```

**成功响应** (200 OK): `{"template_id": 1, "status": "success"}`

---

#### `DELETE /api/workflow_templates/{template_id}` / `PUT /api/workflow_templates/{template_id}`

**说明**: 删除收藏；重命名（请求体 `{"name": "新名称"}`）

---

### 16. 技能广场

#### `GET /api/skills`

**说明**: 技能分页列表；仅 `status=approved`；`saved_only=true` 时需鉴权且只返回当前用户收藏

**请求参数** (Query):
- `main_cat` (string, 可选): 大类筛选
- `sub_cat` (string, 可选): 小类筛选
- `saved_only` (bool, 可选): 仅收藏，默认 false
- `page` (int, 可选): 页码，默认 1
- `size` (int, 可选): 每页条数，默认 12

**响应示例**:

```json
{
  "items": [
    {
      "id": 1,
      "name": "技能名称",
      "description": "...",
      "main_category": "多模态组学",
      "sub_category": "转录组学",
      "prompt_template": "...",
      "author_id": "system",
      "created_at": "2025-03-13T12:00:00",
      "saved": false
    }
  ],
  "total": 42
}
```

**错误响应** (401): `saved_only=true` 且未登录时 `{"detail": "查看「我的」收藏需要登录或提供身份"}`

---

#### `POST /api/skills/{skill_id}/bookmark`

**说明**: 收藏技能；防重复；需鉴权

**成功响应** (200 OK):

```json
{
  "status": "success",
  "message": "已添加到我的工具",
  "saved": true
}
```

---

#### `DELETE /api/skills/{skill_id}/bookmark`

**说明**: 取消收藏；需鉴权

**成功响应** (200 OK): `{"status": "success", "message": "已取消收藏", "saved": false}`

---

#### `POST /api/admin/bootstrap-skills`

**说明**: 管理员专用；重新注入系统技能（清理 author_id=system 后补种）

**请求头**: `Authorization: Bearer <access_token>`，且用户 `role=admin`

**成功响应** (200 OK): 补种结果描述

---

## SSE 流式响应格式

当 `stream: true` 时，`/api/chat` 接口返回 Server-Sent Events (SSE) 格式的流式响应。

### SSE 事件类型

| 事件类型 | 说明 | 数据格式 |
|---------|------|---------|
| `status` | 状态更新 | `{ "content": "状态消息", "state": "状态值" }` |
| `message` | 文本消息（**Delta token**） | `{ "content": "增量内容" }` ⚠️ **只包含新 token，客户端必须追加** |
| `workflow` | 工作流配置 | `{ "workflow_config": {...}, "template_mode": true/false }` |
| `step_result` | 步骤执行结果 | `{ "report_data": {...} }` |
| `diagnosis` | 诊断报告 | `{ "report_data": {...} }` |
| `result` | 最终结果 | `{ "report_data": {...} }` 或 `{ "workflow_config": {...} }` |
| `done` | 完成信号 | `{ "status": "success" }` |
| `state_snapshot` | 会话状态快照（入库用） | `{ "text", "workflow", "steps", "report" }` |
| `error` | 错误信息 | `{ "error": "错误描述", "message": "用户友好的错误消息" }` |

**⚠️ 重要提示**: `message` 事件中的 `content` 字段只包含**增量 token**，不是累积文本。客户端必须使用 `buffer += data.content` 来累积内容。

### SSE 事件格式

每个事件遵循标准 SSE 格式：

```
event: {event_type}
data: {json_data}

```

### 状态值 (state)

`status` 事件中的 `state` 字段可能的值：

- `"start"`: 开始处理
- `"analyzing"`: 正在分析
- `"thinking"`: 正在思考
- `"running"`: 正在执行
- `"rendering"`: 正在渲染
- `"generating_report"`: 正在生成报告
- `"completed"`: 执行完成
- `"error"`: 发生错误
- `"async_job_started"`: 异步作业已启动
- `"waiting"`: 等待中

### 前端处理示例（Delta Streaming + CoT）

```javascript
async function handleSSEStream(response) {
  const reader = response.body.getReader();
  const decoder = new TextDecoder();
  let sseBuffer = '';  // SSE 解析缓冲区
  let currentEventType = null;
  
  // 🔥 Delta Streaming: 消息内容缓冲区
  let messageBuffer = '';  // 累积的消息内容
  let reasoningBuffer = '';  // 累积的思考过程
  let isInReasoning = false;  // 是否在思考标签内

  while (true) {
    const { done, value } = await reader.read();
    if (done) break;

    sseBuffer += decoder.decode(value, { stream: true });
    const lines = sseBuffer.split('\n');
    sseBuffer = lines.pop() || '';  // 保留不完整的行

    for (const line of lines) {
      if (line.startsWith('event: ')) {
        currentEventType = line.substring(7).trim();
      } else if (line.startsWith('data: ')) {
        const dataStr = line.substring(6).trim();
        try {
          const data = JSON.parse(dataStr);
          handleSSEEvent(currentEventType, data);
        } catch (e) {
          console.error('JSON 解析错误:', e, '数据:', dataStr);
        }
      }
    }
  }
}

function handleSSEEvent(eventType, data) {
  switch (eventType) {
    case 'status':
      console.log(`[状态] ${data.state}: ${data.content}`);
      updateStatusUI(data.state, data.content);
      break;
      
    case 'message':
      // 🔥 CRITICAL: Delta Streaming - 必须追加，不能替换
      if (data.content) {
        messageBuffer += data.content;  // ✅ 追加 Delta token
        
        // 解析 Chain of Thought 标签（支持跨 chunk）
        const parsed = parseReasoningTagsStream(messageBuffer);
        if (parsed.reasoning) {
          reasoningBuffer = parsed.reasoning;
          messageBuffer = parsed.answer;
        }
        
        // 更新 UI（分别显示思考过程和最终答案）
        updateChatBubbleWithReasoning({
          thinking: reasoningBuffer,
          answer: messageBuffer,
          isComplete: parsed.isComplete
        });
      }
      break;
      
    case 'workflow':
      console.log('[工作流]', data.workflow_config);
      renderWorkflowCard(data.workflow_config);
      break;
    case 'step_result':
      console.log('[步骤结果]', data.report_data);
      renderStepResult(data.report_data);
      break;
    case 'diagnosis':
      console.log('[诊断报告]', data.report_data);
      renderDiagnosis(data.report_data);
      break;
    case 'result':
      console.log('[最终结果]', data);
      renderFinalResult(data);
      break;
    case 'done':
      console.log('[完成]', data.status);
      onStreamComplete(data.status);
      break;
    case 'error':
      console.error('[错误]', data.error);
      showError(data.message || data.error);
      break;
    default:
      console.log(`[未知事件] ${eventType}:`, data);
  }
}

// 🔥 解析思考标签的流式解析器（支持跨 chunk）
function parseReasoningTagsStream(content) {
  const THINK_START = '<think>';
  const THINK_END = '</think>';
  
  const startIndex = content.indexOf(THINK_START);
  const endIndex = content.indexOf(THINK_END);
  
  if (startIndex !== -1 && endIndex !== -1) {
    // 完整的思考标签
    const reasoning = content.substring(
      startIndex + THINK_START.length,
      endIndex
    );
    const answer = content.substring(endIndex + THINK_END.length).trim();
    return {
      reasoning: reasoning,
      answer: answer,
      isComplete: true
    };
  } else if (startIndex !== -1) {
    // 思考标签开始但未结束（跨 chunk）
    const reasoning = content.substring(startIndex + THINK_START.length);
    return {
      reasoning: reasoning,
      answer: '',
      isComplete: false
    };
  } else {
    // 没有思考标签
    return {
      reasoning: null,
      answer: content,
      isComplete: true
    };
  }
}

// 更新聊天气泡（支持思考过程 UI）
function updateChatBubbleWithReasoning(parsed) {
  // 显示思考过程（可折叠）
  if (parsed.thinking) {
    updateThinkingBox(parsed.thinking, parsed.isComplete);
  }
  
  // 显示最终答案
  if (parsed.answer) {
    appendMessage(parsed.answer);
  }
}
```

**Delta Streaming 示例**:

```
event: message
data: {"content": "Ap"}

event: message
data: {"content": "ple"}  // 只包含新 token "ple"，不是 "Apple"

event: message
data: {"content": " is"}

event: message
data: {"content": " a"}

event: message
data: {"content": " fruit"}
```

客户端处理：
```javascript
let buffer = '';
buffer += "Ap";      // buffer = "Ap"
buffer += "ple";     // buffer = "Apple"
buffer += " is";     // buffer = "Apple is"
buffer += " a";      // buffer = "Apple is a"
buffer += " fruit";  // buffer = "Apple is a fruit"
```

---

## 数据结构定义

### 文件信息 (FileInfo)

```typescript
interface FileInfo {
  name: string;        // 文件名
  size: number;        // 文件大小（字节）
  path: string;        // 文件路径（相对路径）
}
```

### 工作流步骤 (WorkflowStep)

```typescript
interface WorkflowStep {
  step_id: string;                    // 步骤ID
  tool_id: string;                    // 工具ID
  name: string;                       // 步骤名称
  params: Record<string, any>;        // 步骤参数
}
```

### 工作流配置 (WorkflowConfig)

```typescript
interface WorkflowConfig {
  workflow_name: string;               // 工作流名称
  steps: WorkflowStep[];               // 步骤列表
  file_paths?: string[];                // 文件路径数组（可选）
}
```

### 步骤结果 (StepResult)

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

---

## 错误处理

### 错误响应格式

所有错误响应遵循以下格式：

```json
{
  "status": "error",
  "error": "错误类型或代码",
  "message": "用户友好的错误消息",
  "detail": "详细错误信息（仅开发环境）"
}
```

### 常见错误码

| HTTP 状态码 | 错误类型 | 说明 | 解决方案 |
|------------|---------|------|---------|
| 400 | `BadRequest` | 请求参数错误 | 检查请求参数格式和必填字段 |
| 403 | `Forbidden` | 文件路径不安全 | 确保文件路径在允许的目录内 |
| 404 | `NotFound` | 资源不存在 | 检查资源ID或路径是否正确 |
| 413 | `PayloadTooLarge` | 文件大小超限 | 减小文件大小或调整 `MAX_FILE_SIZE` 配置 |
| 500 | `InternalServerError` | 服务器内部错误 | 查看服务器日志获取详细信息 |
| 503 | `ServiceUnavailable` | 服务不可用 | 检查服务组件是否已初始化（如工具检索器） |

### 错误处理最佳实践

1. **前端错误处理**:
   ```javascript
   try {
     const response = await fetch('/api/chat', { ... });
     if (!response.ok) {
       const error = await response.json();
       throw new Error(error.message || error.detail || '请求失败');
     }
     const data = await response.json();
     // 处理成功响应
   } catch (error) {
     // 显示用户友好的错误消息
     showError(error.message);
   }
   ```

2. **流式响应错误处理**:
   ```javascript
   // 在 SSE 流中监听 error 事件
   if (eventType === 'error') {
     showError(data.message || data.error);
     // 可以选择继续或中断流
   }
   ```

---

## 使用示例

### 完整工作流示例

```javascript
// 1. 上传文件
const formData = new FormData();
formData.append('files', fileInput.files[0]);
formData.append('user_id', 'guest');

const uploadResponse = await fetch('/api/upload', {
  method: 'POST',
  body: formData
});
const uploadResult = await uploadResponse.json();

// 2. 发送分析请求（流式）
const chatResponse = await fetch('/api/chat', {
  method: 'POST',
  headers: { 'Content-Type': 'application/json' },
  body: JSON.stringify({
    message: '分析这个文件',
    uploaded_files: uploadResult.file_paths.map(path => ({
      name: path.split('/').pop(),
      path: path
    })),
    stream: true,
    user_id: 'guest',
    session_id: uploadResult.session_id
  })
});

// 3. 处理流式响应
const reader = chatResponse.body.getReader();
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

## 前端集成指南

### 1. 文件上传流程

1. 用户选择文件
2. 调用 `/api/upload` 上传文件
3. 保存返回的 `file_paths` 和 `session_id`
4. 在后续请求中使用这些路径

### 2. 聊天流程

1. 构建请求体，包含 `message`、`uploaded_files`、`stream` 等字段
2. 根据 `stream` 参数选择处理方式：
   - `stream: true`: 使用 SSE 流式处理
   - `stream: false`: 使用 JSON 响应
3. 根据响应类型（`type` 字段）处理不同的响应：
   - `workflow_config`: 显示工作流配置卡片
   - `analysis_report`: 显示分析报告
   - `error`: 显示错误消息

### 3. 工作流执行流程

1. 用户确认工作流配置
2. 调用 `/api/execute` 或通过 `/api/chat` 发送 `workflow_data`
3. 监听执行进度（流式响应）或等待完成（JSON 响应）
4. 渲染执行结果和 AI 专家分析报告

### 4. 状态管理建议

- 使用状态管理库（如 Redux、Vuex）管理：
  - 已上传的文件列表
  - 当前工作流配置
  - 执行状态和结果
  - 用户ID 和会话ID

### 5. 错误处理建议

- 实现全局错误处理机制
- 显示用户友好的错误消息
- 记录错误日志用于调试
- 提供重试机制

---

## 附录

### A. 环境变量配置

| 变量名 | 说明 | 默认值 |
|--------|------|--------|
| `UPLOAD_DIR` | 上传文件目录 | `/app/uploads` |
| `RESULTS_DIR` | 结果输出目录 | `/app/results` |
| `MAX_FILE_SIZE` | 最大文件大小（字节） | `104857600` (100MB) |
| `ALLOWED_ORIGINS` | CORS 允许的来源 | `*` |
| `SILICONFLOW_API_KEY` | SiliconFlow API Key | - |
| `SILICONFLOW_MODEL` | SiliconFlow 模型名称 | - |
| `OLLAMA_BASE_URL` | Ollama 服务地址 | `http://localhost:11434` |
| `OLLAMA_EMBEDDING_MODEL` | Ollama Embedding 模型 | `nomic-embed-text` |
| `CHROMA_PERSIST_DIR` | ChromaDB 持久化目录 | `./data/chroma_tools` |

### B. 文件路径说明

- **上传文件路径**: 相对于 `UPLOAD_DIR`，格式: `{user_id}/{session_id}/{filename}`
- **结果文件路径**: 相对于 `RESULTS_DIR`，格式: `run_{timestamp}/{filename}`
- **访问结果文件**: 通过 `/results/{path}` 静态文件服务访问

### C. 多用户与鉴权

- **身份解析**: 优先 `Authorization: Bearer <access_token>`（登录用户），否则 `X-Guest-UUID`（游客）；解析结果为 `owner_id`，用于会话/消息/资产/工作流收藏/技能收藏隔离。
- **目录隔离**: 上传目录按 `owner_id` 隔离: `{UPLOAD_DIR}/{owner_id}/{batch_id}/`；Asset 表记录 `owner_id`。
- **会话/消息**: Session、Message 表含 `owner_id` / `session_id`，新建会话时由后端生成 `session_id` 并落库。
- **工作流收藏**: WorkflowTemplate 表按 `owner_id` 隔离；技能收藏 UserSavedSkill 按 `owner_id` + `skill_id` 防重复。

### D. 10x Genomics 数据特殊处理

- 自动识别 `matrix.mtx`、`barcodes.tsv`、`features.tsv`（或 `genes.tsv`）
- 自动分组保存到独立子目录
- 返回的 `file_paths` 指向组目录，而不是单个文件

---

**文档版本**: v2.1  
**最后更新**: 2025-03-13  
**维护者**: GIBH Agent Team
