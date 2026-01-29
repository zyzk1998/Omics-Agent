# 错误响应格式

所有错误响应遵循以下格式：

```json
{
  "status": "error",
  "error": "错误类型或代码",
  "message": "用户友好的错误消息",
  "detail": "详细错误信息（仅开发环境）"
}
```

## 示例

```json
{
  "status": "error",
  "error": "BadRequest",
  "message": "请求参数错误：缺少必需字段 file_paths",
  "detail": "workflow_data.steps[0].params.file_path is required"
}
```

---

**返回**: [错误处理目录](../README.md) | [API 手册首页](../../README.md)
