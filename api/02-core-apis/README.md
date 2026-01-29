# 核心 API

核心 API 是系统最常用的接口，包括文件上传、聊天、工作流执行等。

## 接口列表

1. [健康检查](health.md) - `GET /api/health`
2. [文件上传](upload.md) - `POST /api/upload`
3. [聊天接口](chat.md) - `POST /api/chat`
4. [执行工作流](execute.md) - `POST /api/execute`

---

## 典型使用流程

```
1. 健康检查 → 确认服务可用
2. 文件上传 → 上传数据文件
3. 聊天接口 → 发送分析请求（流式响应）
4. 执行工作流 → 执行确认的工作流
```

详细流程请参考 [使用示例](../08-examples/README.md)。
