# Omics Agent API 文档

**版本**: v2.0 · **基础 URL**: `http://127.0.0.1:8028`（开发环境，请按部署替换）  
**协议**: HTTP/1.1 · **数据格式**: JSON（文件上传为 `multipart/form-data`）

---

## 快速开始

1. 调用 `GET /api/health` 确认服务可用。
2. 使用 `POST /api/upload` 上传组学数据（支持多文件）。
3. 使用 `POST /api/chat` 发起对话；响应为 **SSE** 流式事件（工作流规划、执行状态、结果等）。
4. 确认工作流后，可用 `POST /api/execute` 直接执行，或在前端工作流卡片中分步运行。

---

## 鉴权

| 场景 | 请求头 |
|------|--------|
| 已登录用户 | `Authorization: Bearer <token>` |
| 访客试用 | `X-Guest-UUID: <uuid>`（能力受限） |

所有 JSON 请求需设置：`Content-Type: application/json`。

---

## 核心端点

| 方法 | 路径 | 说明 |
|------|------|------|
| GET | `/api/health` | 健康检查与组件状态 |
| POST | `/api/upload` | 文件上传（默认单文件上限见服务端配置） |
| POST | `/api/chat` | 主对话接口，支持 SSE |
| POST | `/api/execute` | 直接执行工作流 |
| GET | `/api/sessions` | 会话列表 |
| GET | `/api/skills` | 技能广场目录 |
| POST | `/api/workflows/plan` | Plan-first 工作流规划 |

更完整的模块化说明见仓库 `api/README.md`。

---

## 首次调用：聊天接口

### 请求示例

```http
POST /api/chat HTTP/1.1
Content-Type: application/json
Authorization: Bearer <your_token>
```

```json
{
  "message": "请对我的单细胞数据进行质控与聚类分析",
  "session_id": null,
  "uploaded_files": [
    {
      "name": "sample.h5ad",
      "path": "/uploads/<user_id>/sample.h5ad"
    }
  ],
  "mode": "standard",
  "thinking_mode": "fast"
}
```

### 响应

- **Content-Type**: `text/event-stream`（SSE）
- 常见事件类型：`session`、`thought`、`message`、`status`、`workflow`、`step_result`、`done`、`state_snapshot` 等。

前端应使用 `EventSource` 或 `fetch` + 流式读取按行解析 `event:` / `data:` 字段。

---

## 统一响应结构

### 成功（非 SSE JSON 端点）

```json
{
  "status": "success",
  "data": {}
}
```

### 错误

```json
{
  "status": "error",
  "error": "错误代码或类型",
  "message": "用户可读说明",
  "detail": "开发环境附加信息（可选）"
}
```

### 常见 HTTP 状态码

| 状态码 | 含义 |
|--------|------|
| 200 | 成功 |
| 400 | 参数错误 |
| 403 | 权限不足（如路径不安全） |
| 404 | 资源不存在 |
| 413 | 上传体积超限 |
| 500 | 服务器内部错误 |
| 503 | 依赖服务不可用 |

---

## 文件上传

```http
POST /api/upload HTTP/1.1
Content-Type: multipart/form-data
```

表单字段名通常为 `files`（多文件）。响应中包含可用于 `uploaded_files[].path` 的服务器路径。

---

## SSE 事件处理要点

- **`workflow`**：携带 `workflow_config`，前端应渲染可编辑工作流卡片。
- **`status`**：`state` 可为 `running` / `completed`；用于中栏执行记录 checklist。
- **`state_snapshot`**：用于历史回放与「时光机」恢复，勿与 `thought` / `message` 混写。
- **`done`**：流结束信号；应关闭 loading 并固化 UI 状态。

---

## 相关资源

- 结构化手册：`api/README.md`
- 交互式架构图：`/whotowork.html`
- OpenAPI（若已启用）：`/api/docs`、`/api/redoc`
