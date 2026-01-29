# 快速开始

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

## 下一步

- 查看 [核心 API](../02-core-apis/README.md) 开始使用接口
- 查看 [使用示例](../08-examples/README.md) 了解完整工作流
