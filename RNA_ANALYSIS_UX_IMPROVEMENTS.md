# RNA 分析用户体验优化总结

## 修复内容

### 1. 修复 Path 变量错误 ✅

**问题**: `annotation.py` 中 `Path` 变量重复导入，导致 `local variable 'Path' referenced before assignment` 错误。

**修复**: 移除第 73 行的重复导入 `from pathlib import Path`，使用文件顶部已导入的 `Path`。

**文件**: `gibh_agent/tools/rna/annotation.py`

---

### 2. Cell Ranger 步骤用户体验优化 ✅

**问题**: 
- Cell Ranger 步骤耗时很长（30 分钟到数小时），但用户没有收到友好的等待提示
- 原始数据文件通常很大，超出上传限制

**修复**:
1. **检测 Cell Ranger 步骤并发送友好提示**:
   - 在 `orchestrator.py` 中检测到 `cellranger` 步骤时，发送专门的等待提示
   - 提示信息：`"⏳ Cell Ranger 正在后台运行，这可能需要较长时间（通常 30 分钟到数小时），请耐心等待..."`
   - SSE 事件中包含 `show_waiting_bubble: true` 标志，前端可以显示等待气泡

2. **预览模式测试数据询问**:
   - 检测用户意图是否包含 RNA 分析全流程（包括 cellranger）
   - 如果检测到但未上传文件，在预览模式中询问是否使用测试数据
   - 发送 `test_data_question` SSE 事件，包含测试数据路径和信息

**文件**: `gibh_agent/core/orchestrator.py`

**SSE 事件格式**:
```json
{
  "type": "workflow",
  "workflow_config": null,
  "template_mode": true,
  "test_data_question": {
    "type": "test_data_question",
    "message": "检测到您想要进行 RNA 分析全流程，但未上传文件...",
    "question": "是否使用测试数据进行全流程分析？",
    "test_data_path": "/app/test_data/pbmc_1k_v3_fastqs",
    "test_data_info": {
      "name": "PBMC 1k v3",
      "description": "10x Genomics 单细胞 RNA-seq 测试数据",
      "size": "约 100MB"
    }
  }
}
```

---

### 3. 错误信息用户友好化 ✅

**问题**: 执行结果中的错误信息是面向程序员的日志，用户难以理解。

**修复**:
1. **创建错误格式化器** (`error_formatter.py`):
   - 将技术错误信息转换为用户友好的提示
   - 分类错误类型：数据问题、配置问题、网络问题、资源问题等
   - 为每种错误类型提供优化建议
   - 标记哪些错误可以跳过，哪些不能跳过

2. **在 Executor 中应用错误格式化**:
   - 步骤执行失败时，自动格式化错误信息
   - 返回格式化的错误信息，包含：
     - `user_message`: 用户友好的错误消息
     - `error_category`: 错误类别
     - `suggestion`: 优化建议
     - `can_skip`: 是否可以跳过
     - `technical_details`: 技术细节（供调试）

**文件**: 
- `gibh_agent/core/error_formatter.py` (新建)
- `gibh_agent/core/executor.py`

**错误类别**:
- `data_issue`: 数据文件问题（不能跳过）
- `config_issue`: 配置/依赖问题（部分可跳过）
- `network_issue`: 网络/下载问题（可跳过）
- `resource_issue`: 资源不足（可跳过）
- `data_quality`: 数据质量问题（通常可忽略）
- `internal_error`: 内部错误（可跳过）
- `unknown_error`: 未知错误（可跳过）

**示例**:
```json
{
  "status": "error",
  "step_id": "rna_cell_annotation",
  "step_name": "细胞类型注释",
  "user_message": "模型下载失败：细胞类型注释步骤需要下载 CellTypist 模型，但下载失败。",
  "error_category": "network_issue",
  "suggestion": "此步骤可以跳过，不影响其他分析结果。您可以稍后手动下载模型或使用其他注释方法。",
  "can_skip": true,
  "technical_details": "无法下载CellTypist模型: Connection timeout..."
}
```

---

## 前端集成说明

### 1. Cell Ranger 等待气泡

前端需要监听 SSE 事件中的 `show_waiting_bubble` 标志：

```javascript
if (data.show_waiting_bubble) {
  // 显示等待气泡
  showWaitingBubble(data.waiting_message || "正在处理，请稍候...");
}

// 当收到 done 事件时，隐藏气泡
if (eventType === 'done') {
  hideWaitingBubble();
}
```

### 2. 测试数据询问

前端需要处理 `test_data_question` 事件：

```javascript
if (data.test_data_question) {
  const question = data.test_data_question;
  // 显示确认对话框
  const confirmed = await showConfirmDialog({
    title: question.question,
    message: question.message,
    testDataInfo: question.test_data_info
  });
  
  if (confirmed) {
    // 使用测试数据重新发送请求
    sendMessage(null, null, null, null, false, [{
      name: question.test_data_info.name,
      path: question.test_data_path
    }]);
  }
}
```

### 3. 错误信息显示

前端在执行结果卡片中显示格式化的错误信息：

```javascript
if (step.status === 'error') {
  // 显示用户友好的错误消息
  showErrorCard({
    message: step.user_message,
    category: step.error_category,
    suggestion: step.suggestion,
    canSkip: step.can_skip,
    technicalDetails: step.technical_details  // 可折叠显示
  });
}
```

---

## 测试建议

1. **Cell Ranger 等待提示测试**:
   - 触发包含 `rna_cellranger_count` 步骤的工作流
   - 验证前端收到 `show_waiting_bubble: true` 标志
   - 验证等待气泡正确显示和隐藏

2. **测试数据询问测试**:
   - 发送 RNA 分析全流程查询，但不上传文件
   - 验证收到 `test_data_question` 事件
   - 验证确认后使用测试数据执行工作流

3. **错误信息格式化测试**:
   - 触发会失败的步骤（如缺少依赖的步骤）
   - 验证错误信息是用户友好的
   - 验证错误类别、建议和可跳过标志正确

---

## 注意事项

1. **测试数据路径**: 确保 `/app/test_data/pbmc_1k_v3_fastqs` 路径存在且可访问
2. **前端兼容性**: 前端需要更新以处理新的 SSE 事件格式
3. **错误分类**: 错误分类规则可能需要根据实际使用情况调整

---

**修复日期**: 2025-01-28  
**修复人员**: AI Assistant  
**状态**: ✅ 已完成
