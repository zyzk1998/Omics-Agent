# 诊断报告传递流程调试指南

## 问题描述
Metabolomics Agent 的诊断报告没有在前端显示，而 RNA Agent 可以正常显示。

## 调试步骤

### 1. 检查 BaseAgent._perform_data_diagnosis() 返回值
- 查看日志: `✅ [DataDiagnostician] 诊断报告生成成功，长度: XXX`
- 如果看到 `⚠️ [DataDiagnostician] 诊断报告为空`，说明 LLM 返回为空

### 2. 检查 MetabolomicsAgent 接收诊断报告
- 查看日志: `📝 [DEBUG] Diagnosis report generated, length: XXX`
- 如果看到 `⚠️ [DEBUG] Diagnosis report is None or empty`，说明诊断失败

### 3. 检查添加到返回结果
- 查看日志: `📝 [DEBUG] Adding diagnosis_report to result, length: XXX`
- 如果看到 `⚠️ [DEBUG] diagnosis_report is None, NOT adding to result`，说明没有添加到结果

### 4. 检查最终返回结构
- 查看日志: `📤 [DEBUG] MetabolomicsAgent returning result with keys: [...]`
- 确认 `diagnosis_report` 在 keys 列表中
- 查看日志: `📤 [DEBUG] MetabolomicsAgent has diagnosis_report: True/False`

### 5. 检查 API 路由传递
- 查看日志: `📤 返回工作流配置: 包含推荐=..., 包含诊断=...`
- 确认 `包含诊断=True`

## 可能的问题

1. **诊断报告生成失败但未捕获**
   - 检查异常日志: `❌ [CHECKPOINT] Unified diagnosis failed`
   - 确保 `diagnosis_report = None` 在异常处理中设置

2. **诊断报告为空字符串**
   - LLM 可能返回空字符串而不是 None
   - 需要检查 `if diagnosis_report:` 是否能正确判断

3. **返回结构被覆盖**
   - 检查是否有其他地方修改了 `result` 字典
   - 确保 `diagnosis_report` 在最后添加到结果中

## 修复建议

如果诊断报告仍然不显示，检查：
1. 前端控制台日志: `📋 [DEBUG] renderWorkflowForm:`
2. 确认 `hasDiagnosisReport: true`
3. 检查 `diagnosisReport` 的值是否为空字符串
