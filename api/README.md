# Omics Agent API 手册

**版本**: v2.0  
**基础 URL**: `http://localhost:8028` (开发环境)  
**协议**: HTTP/1.1  
**数据格式**: JSON (除文件上传外)

---

## 📚 目录结构

本 API 手册按功能模块组织，方便快速查找和参考：

### 📖 快速开始
- [01-getting-started/README.md](01-getting-started/README.md) - 通用说明、请求头、响应格式、HTTP 状态码

### 🔌 核心 API
- [02-core-apis/README.md](02-core-apis/README.md) - 核心接口文档
  - [健康检查](02-core-apis/health.md)
  - [文件上传](02-core-apis/upload.md)
  - [聊天接口](02-core-apis/chat.md)
  - [执行工作流](02-core-apis/execute.md)

### 🛠️ 工具 API
- [03-tool-apis/README.md](03-tool-apis/README.md) - 工具检索和管理接口
  - [语义搜索工具](03-tool-apis/search.md)
  - [列出所有工具](03-tool-apis/list.md)
  - [获取工具 Schema](03-tool-apis/get-schema.md)

### 📋 工作流 API
- [04-workflow-apis/README.md](04-workflow-apis/README.md) - 工作流管理接口
  - [规划工作流](04-workflow-apis/plan.md)
  - [保存工作流](04-workflow-apis/save.md)
  - [列出工作流](04-workflow-apis/list.md)
  - [删除工作流](04-workflow-apis/delete.md)
  - [任务历史](04-workflow-apis/jobs-history.md)
  - [工作流状态](04-workflow-apis/status.md)

### 📡 SSE 流式响应
- [05-sse-streaming/README.md](05-sse-streaming/README.md) - Server-Sent Events 流式响应格式
  - [事件类型](05-sse-streaming/event-types.md)
  - [事件格式](05-sse-streaming/event-format.md)
  - [前端处理](05-sse-streaming/frontend-handling.md)

### 📊 数据结构
- [06-data-structures/README.md](06-data-structures/README.md) - 数据结构定义
  - [文件信息](06-data-structures/file-info.md)
  - [工作流配置](06-data-structures/workflow-config.md)
  - [步骤结果](06-data-structures/step-result.md)

### ⚠️ 错误处理
- [07-error-handling/README.md](07-error-handling/README.md) - 错误处理指南
  - [错误响应格式](07-error-handling/error-format.md)
  - [常见错误码](07-error-handling/error-codes.md)
  - [最佳实践](07-error-handling/best-practices.md)

### 💡 使用示例
- [08-examples/README.md](08-examples/README.md) - 完整使用示例
  - [文件上传示例](08-examples/upload-example.md)
  - [聊天接口示例](08-examples/chat-example.md)
  - [工作流执行示例](08-examples/workflow-example.md)
  - [SSE 流式处理示例](08-examples/sse-example.md)

### 📎 附录
- [09-appendix/README.md](09-appendix/README.md) - 附录信息
  - [环境变量配置](09-appendix/environment-variables.md)
  - [文件路径说明](09-appendix/file-paths.md)
  - [多用户支持](09-appendix/multi-user.md)
  - [10x Genomics 数据](09-appendix/10x-genomics.md)

---

## 🚀 快速导航

### 按使用场景查找

**我想上传文件并分析**
1. 阅读 [文件上传](02-core-apis/upload.md)
2. 阅读 [聊天接口](02-core-apis/chat.md)
3. 参考 [完整工作流示例](08-examples/workflow-example.md)

**我想使用流式响应**
1. 阅读 [SSE 流式响应](05-sse-streaming/README.md)
2. 参考 [SSE 流式处理示例](08-examples/sse-example.md)

**我想管理工具**
1. 阅读 [工具 API](03-tool-apis/README.md)
2. 参考 [语义搜索工具](03-tool-apis/search.md)

**我想管理工作流**
1. 阅读 [工作流 API](04-workflow-apis/README.md)
2. 参考 [规划工作流](04-workflow-apis/plan.md)

**我遇到了错误**
1. 阅读 [错误处理](07-error-handling/README.md)
2. 查看 [常见错误码](07-error-handling/error-codes.md)

---

## 📖 完整文档

如需查看完整的单文件文档，请参考根目录下的 [API.md](../API.md)。

---

**文档版本**: v2.0  
**最后更新**: 2025-01-28  
**维护者**: Omics Agent Team
