# 文件上传示例

## 完整示例

```javascript
// 1. 用户选择文件
const fileInput = document.getElementById('fileInput');
const files = fileInput.files;

// 2. 构建 FormData
const formData = new FormData();
for (let file of files) {
    formData.append('files', file);
}
formData.append('user_id', 'guest');
formData.append('session_id', '20250128_120000');

// 3. 发送上传请求
const response = await fetch('/api/upload', {
    method: 'POST',
    body: formData
});

// 4. 处理响应
const result = await response.json();
if (result.status === 'success') {
    console.log('上传成功:', result.file_paths);
    console.log('会话ID:', result.session_id);
    
    // 保存文件路径用于后续请求
    uploadedFiles = result.file_paths;
    sessionId = result.session_id;
} else {
    console.error('上传失败:', result.error);
}
```

---

**返回**: [使用示例目录](../README.md) | [API 手册首页](../../README.md)
