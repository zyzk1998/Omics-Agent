# 工作流执行示例

## 完整工作流示例

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

**返回**: [使用示例目录](../README.md) | [API 手册首页](../../README.md)
