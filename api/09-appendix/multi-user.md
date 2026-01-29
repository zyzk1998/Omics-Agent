# 多用户支持

系统支持多用户隔离：

## 文件目录隔离

- 每个用户有独立的文件目录: `{UPLOAD_DIR}/{user_id}/`
- 每个会话有独立的子目录: `{UPLOAD_DIR}/{user_id}/{session_id}/`
- 示例: `/app/uploads/guest/20250128_120000/`

## 工作流和任务历史隔离

- 工作流按用户ID隔离存储
- 任务历史按用户ID过滤
- 默认用户ID为 `"guest"`（未登录用户）

## 使用方式

在 API 请求中指定 `user_id` 和 `session_id`：

```javascript
{
  "user_id": "user123",
  "session_id": "20250128_120000",
  ...
}
```

---

**返回**: [附录目录](../README.md) | [API 手册首页](../../README.md)
