# 前端处理 SSE 流

## 完整示例

```javascript
async function handleSSEStream(response) {
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
          handleSSEEvent(currentEventType, data);
        } catch (e) {
          console.error('JSON 解析错误:', e, '数据:', dataStr);
        }
      }
    }
  }
}

function handleSSEEvent(eventType, data) {
  switch (eventType) {
    case 'status':
      console.log(`[状态] ${data.state}: ${data.content}`);
      updateStatusUI(data.state, data.content);
      break;
    case 'message':
      console.log(`[消息] ${data.content}`);
      appendMessage(data.content);
      break;
    case 'workflow':
      console.log('[工作流]', data.workflow_config);
      renderWorkflowCard(data.workflow_config);
      break;
    case 'step_result':
      console.log('[步骤结果]', data.report_data);
      renderStepResult(data.report_data);
      break;
    case 'diagnosis':
      console.log('[诊断报告]', data.report_data);
      renderDiagnosis(data.report_data);
      break;
    case 'result':
      console.log('[最终结果]', data);
      renderFinalResult(data);
      break;
    case 'done':
      console.log('[完成]', data.status);
      onStreamComplete(data.status);
      break;
    case 'error':
      console.error('[错误]', data.error);
      showError(data.message || data.error);
      break;
    default:
      console.log(`[未知事件] ${eventType}:`, data);
  }
}
```

## 使用示例

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

await handleSSEStream(response);
```

---

**返回**: [SSE 流式响应目录](../README.md) | [API 手册首页](../../README.md)
