# 修复诊断报告显示问题

## 问题描述
Metabolomics Agent 的诊断报告没有在前端显示，而 RNA Agent 可以正常显示。

## 根本原因分析

1. **异常处理不完整**：当诊断失败时，`diagnosis_report` 变量可能没有被正确设置为 None
2. **空字符串检查不足**：LLM 可能返回空字符串而不是 None，需要检查 `strip()`
3. **调试信息不足**：无法追踪诊断报告在传递过程中的状态

## 修复内容

### 1. BaseAgent._perform_data_diagnosis()
- ✅ 添加调试日志：打印诊断报告长度和预览
- ✅ 确保返回 None 或有效字符串

### 2. MetabolomicsAgent._generate_workflow_config()
- ✅ 修复异常处理：确保 `diagnosis_report = None` 在异常时设置
- ✅ 添加调试日志：追踪诊断报告生成和添加过程
- ✅ 改进验证逻辑：检查 `diagnosis_report` 是否为有效字符串（非 None、非空、非纯空白）

### 3. RNAAgent._generate_workflow_config()
- ✅ 同步修复：应用相同的验证逻辑和调试日志
- ✅ 确保异常处理正确

## 调试日志位置

### BaseAgent
```
✅ [DataDiagnostician] 诊断报告生成成功，长度: XXX
📝 [DEBUG] Diagnosis report preview: ...
⚠️ [DataDiagnostician] 诊断报告为空
```

### MetabolomicsAgent
```
📝 [DEBUG] Diagnosis report generated, length: XXX
📝 [DEBUG] Diagnosis report preview: ...
⚠️ [DEBUG] Diagnosis report is None or empty
📝 [DEBUG] Adding diagnosis_report to result, length: XXX
⚠️ [DEBUG] diagnosis_report is invalid (None/empty), NOT adding to result
📤 [DEBUG] MetabolomicsAgent returning result with keys: [...]
📤 [DEBUG] MetabolomicsAgent has diagnosis_report: True/False
```

### server.py
```
📤 返回工作流配置: 包含推荐=..., 包含诊断=...
```

## 验证步骤

1. **重启容器**
   ```bash
   sudo docker compose restart api-server worker
   ```

2. **测试 Metabolomics Agent**
   - 上传代谢组学数据文件
   - 查看日志，确认诊断报告生成
   - 检查前端是否显示诊断报告

3. **检查日志**
   - 按照上述日志位置追踪诊断报告传递过程
   - 确认每个步骤都正常执行

## 预期结果

- ✅ Metabolomics Agent 的诊断报告应该在前端显示
- ✅ 诊断报告显示在蓝色卡片中（优先于推荐信息）
- ✅ 诊断报告使用 Markdown 格式渲染

## 如果仍然不显示

1. 检查前端控制台日志：`📋 [DEBUG] renderWorkflowForm:`
2. 确认 `hasDiagnosisReport: true`
3. 检查 `diagnosisReport` 的值是否为空字符串
4. 查看服务器日志，确认诊断报告是否生成
