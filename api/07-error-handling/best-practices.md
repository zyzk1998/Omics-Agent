# 错误处理最佳实践

## 1. 前端错误处理

```javascript
try {
  const response = await fetch('/api/chat', { ... });
  if (!response.ok) {
    const error = await response.json();
    throw new Error(error.message || error.detail || '请求失败');
  }
  const data = await response.json();
  // 处理成功响应
} catch (error) {
  // 显示用户友好的错误消息
  showError(error.message);
}
```

## 2. 流式响应错误处理

```javascript
// 在 SSE 流中监听 error 事件
if (eventType === 'error') {
  showError(data.message || data.error);
  // 可以选择继续或中断流
}
```

## 3. 重试机制

对于网络错误或临时故障，可以实现重试机制：

```javascript
async function fetchWithRetry(url, options, maxRetries = 3) {
  for (let i = 0; i < maxRetries; i++) {
    try {
      const response = await fetch(url, options);
      if (response.ok) return response;
    } catch (error) {
      if (i === maxRetries - 1) throw error;
      await new Promise(resolve => setTimeout(resolve, 1000 * (i + 1)));
    }
  }
}
```

---

**返回**: [错误处理目录](../README.md) | [API 手册首页](../../README.md)
