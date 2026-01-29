# 聊天接口示例

## JSON 响应示例

```javascript
const response = await fetch('/api/chat', {
  method: 'POST',
  headers: { 'Content-Type': 'application/json' },
  body: JSON.stringify({
    message: '分析这个文件',
    uploaded_files: [
      { name: 'example.csv', path: 'guest/20250128_120000/example.csv' }
    ],
    stream: false,
    user_id: 'guest',
    session_id: '20250128_120000'
  })
});

const data = await response.json();

if (data.type === 'workflow_config') {
  // 显示工作流配置
  renderWorkflowCard(data.workflow_data);
} else if (data.type === 'analysis_report') {
  // 显示分析报告
  renderAnalysisReport(data.report_data);
}
```

---

**返回**: [使用示例目录](../README.md) | [SSE 流式处理示例](sse-example.md) | [API 手册首页](../../README.md)
