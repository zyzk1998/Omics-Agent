# SSE 事件格式

## 标准 SSE 格式

每个事件遵循标准 SSE 格式：

```
event: {event_type}
data: {json_data}

```

## 示例

### status 事件

```
event: status
data: {"content": "正在执行步骤: 数据检查...", "state": "running"}

```

### workflow 事件

```
event: workflow
data: {"workflow_config": {"workflow_name": "Metabolomics Analysis Pipeline", "steps": [...]}, "template_mode": false}

```

### step_result 事件

```
event: step_result
data: {"report_data": {"step_id": "inspect_data", "name": "数据检查", "status": "success", ...}}

```

### diagnosis 事件

```
event: diagnosis
data: {"report_data": {"diagnosis": "## AI 专家分析报告\n\n..."}}

```

### done 事件

```
event: done
data: {"status": "success"}

```

---

**返回**: [SSE 流式响应目录](../README.md) | [API 手册首页](../../README.md)
