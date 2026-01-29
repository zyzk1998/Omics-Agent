# 文件信息 (FileInfo)

```typescript
interface FileInfo {
  name: string;        // 文件名
  size: number;        // 文件大小（字节）
  path: string;        // 文件路径（相对路径）
}
```

## 使用场景

- 文件上传接口返回的文件信息
- 聊天接口中的 `uploaded_files` 字段

---

**返回**: [数据结构目录](../README.md) | [API 手册首页](../../README.md)
