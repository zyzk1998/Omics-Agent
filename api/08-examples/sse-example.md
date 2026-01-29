# SSE 流式处理示例

## 完整示例

参考 [SSE 流式响应 - 前端处理](../05-sse-streaming/frontend-handling.md) 获取完整的 SSE 处理代码。

## 简化示例

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
let buffer = '';
let currentEventType = null;

while (true) {
  const { done, value } = await reader.read();
  if (done) break;
  
  buffer += decoder.decode(value, { stream: true });
  const lines = buffer.split('\n');
  buffer = lines.pop() || '';
  
  for (const line of lines) {
    if (line.startsWith('event: ')) {
      currentEventType = line.substring(7).trim();
    } else if (line.startsWith('data: ')) {
      const dataStr = line.substring(6).trim();
      try {
        const data = JSON.parse(dataStr);
        // 根据 currentEventType 处理数据
        handleSSEEvent(currentEventType, data);
      } catch (e) {
        console.error('JSON 解析错误:', e);
      }
    }
  }
}
```

---

**返回**: [使用示例目录](../README.md) | [SSE 流式响应](../05-sse-streaming/README.md) | [API 手册首页](../../README.md)
