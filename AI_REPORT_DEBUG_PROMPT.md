# 🔍 AI专家分析报告生成失败问题诊断提示词

## 📋 问题现象

根据控制台日志，AI专家分析报告显示的是后备摘要内容，而不是真正的生信分析报告：

```
diagnosis: "## 分析结果摘要\n\n本次分析完成了 5 个步骤。请查看上方的详细图表和统计结果以获取更深入的生物学解释。\n\n### 关键发现\n- 成功步骤: 5/5\n- 请查看执行结果中的图表和数据表格获取详细分析。"
```

这是 `orchestrator.py` 第508-514行的后备摘要，说明LLM调用可能失败或返回内容过短。

---

## 🎯 根本原因分析

### 可能原因1：`extract_think_and_content` 处理问题

**问题描述**：
- LLM返回的内容可能包含 `<think>` 或 `<think>` 标签
- `extract_think_and_content` 移除标签后，如果主要内容在标签内，`actual_content` 会变空或很短
- 导致 `base_agent.py` 第1325行的检查失败（`len(response.strip()) <= 100`）
- 进入重试逻辑，但重试后可能仍然失败
- 最终返回错误信息，但被 `orchestrator.py` 的后备逻辑替换

**检查方法**：
查看后端日志中是否有：
```
🔍 [AnalysisSummary] 内容提取结果: think_length=XXX, response_length=YYY
⚠️ [AnalysisSummary] LLM 返回内容过短（XX 字符），尝试重新生成...
```

### 可能原因2：LLM调用失败但异常被捕获

**问题描述**：
- LLM调用抛出异常（如超时、网络错误等）
- 异常被 `except Exception as llm_error` 捕获
- 返回了错误信息，但错误信息可能被orchestrator的后备逻辑替换

**检查方法**：
查看后端日志中是否有：
```
❌ [AnalysisSummary] LLM 调用失败
```

### 可能原因3：LLM返回内容格式不符合预期

**问题描述**：
- LLM返回的内容可能不是标准的Markdown格式
- 或者包含了大量格式化字符
- 导致 `extract_think_and_content` 解析失败

---

## 🔧 已实施的修复

### 修复1：改进内容提取后的检查逻辑

**位置**：`gibh_agent/agents/base_agent.py` 第1318-1335行

**改进**：
- 添加详细日志，记录 `original_content`、`think_content`、`response` 的长度
- 如果 `response` 太短但 `original_content` 很长，使用 `original_content`
- 检查是否有思考标签，如果有则返回原始内容（让前端解析）

### 修复2：改进 `extract_think_and_content` 处理逻辑

**位置**：`gibh_agent/core/llm_client.py` 第306-312行

**改进**：
- 如果移除标签后内容过短（<50字符）但原始内容很长（>200字符），保留原始内容
- 添加警告日志，提示可能的主要内容在标签内

### 修复3：增强日志记录

**位置**：
- `gibh_agent/agents/base_agent.py` 第1320-1324行
- `gibh_agent/core/orchestrator.py` 第488行

**改进**：
- 记录原始内容长度
- 记录提取后的内容长度
- 记录内容预览（前300字符）

---

## 📝 调试步骤

### 步骤1：检查后端日志

**查找以下日志**：
```bash
# 检查LLM是否被调用
grep "📞 \[AnalysisSummary\] 开始LLM调用" /path/to/logs

# 检查LLM调用是否完成
grep "✅ \[AnalysisSummary\] LLM调用完成" /path/to/logs

# 检查内容提取结果
grep "🔍 \[AnalysisSummary\] 内容提取结果" /path/to/logs

# 检查是否进入重试逻辑
grep "⚠️ \[AnalysisSummary\] LLM 返回内容过短" /path/to/logs

# 检查最终结果
grep "✅ \[Orchestrator\] AI专家分析报告生成完成" /path/to/logs
grep "⚠️ \[Orchestrator\] 摘要过短" /path/to/logs
```

### 步骤2：分析日志输出

**如果看到**：
```
🔍 [AnalysisSummary] 内容提取结果: think_length=500, response_length=50
⚠️ [AnalysisSummary] LLM 返回内容过短（50 字符），尝试重新生成...
```

**说明**：
- LLM返回了内容，但主要内容在思考标签内
- `extract_think_and_content` 移除标签后，实际内容变短
- 需要检查修复是否生效（应该使用 `original_content`）

**如果看到**：
```
❌ [AnalysisSummary] LLM 调用失败
```

**说明**：
- LLM调用抛出异常
- 需要检查异常类型和错误信息
- 查看 `last_llm_error` 中的详细信息

**如果看到**：
```
⚠️ [Orchestrator] 摘要过短（XX字符），使用结构化后备
```

**说明**：
- LLM返回的内容被orchestrator的后备逻辑替换
- 需要检查为什么 `summary` 长度 < 50字符

### 步骤3：验证修复效果

**预期日志**：
```
📞 [AnalysisSummary] 开始LLM调用，max_tokens=2500...
✅ [AnalysisSummary] LLM调用完成，开始解析响应...
🔍 [AnalysisSummary] 原始内容长度: 1500
🔍 [AnalysisSummary] 内容提取结果: think_length=200, response_length=1300
✅ [AnalysisSummary] 深度生物学解释生成成功，长度: 1300
✅ [Orchestrator] AI专家分析报告生成完成，耗时: 15.23秒，长度: 1300字符
✅ [Orchestrator] summary是有效的生信分析内容，长度: 1300字符
```

**如果仍然失败**：
- 检查 `original_content` 的实际内容
- 检查是否有思考标签
- 检查 `extract_think_and_content` 是否正确处理

---

## 🚀 立即行动项

1. **查看后端日志**，确认LLM是否被调用以及返回的内容长度
2. **检查修复是否生效**，确认是否使用了 `original_content` 当 `response` 太短时
3. **如果问题仍然存在**，需要进一步检查：
   - LLM API是否正常工作
   - API密钥是否有效
   - 网络连接是否正常
   - 请求内容是否过长导致超时

---

## 🔍 关键检查点

| 检查点 | 位置 | 预期结果 | 如果失败 |
|--------|------|----------|----------|
| LLM调用 | base_agent.py:1316 | 有日志"开始LLM调用" | LLM未调用 |
| LLM响应 | base_agent.py:1318 | 有日志"LLM调用完成" | LLM调用失败 |
| 原始内容长度 | base_agent.py:1321 | >500字符 | LLM返回内容不足 |
| 内容提取 | base_agent.py:1322 | response或original_content长度>100 | 进入重试逻辑 |
| 重试结果 | base_agent.py:1354 | response或original_content长度>100 | 返回错误信息 |
| orchestrator检查 | orchestrator.py:506 | summary不为None且>50字符 | 使用后备摘要 |

---

## 💡 临时解决方案

如果问题仍然存在，可以临时禁用内容提取，直接使用原始内容：

```python
# 在 base_agent.py 第1319行后
# 临时禁用内容提取，直接使用原始内容
original_content = completion.choices[0].message.content or ""
if original_content and len(original_content.strip()) > 100:
    logger.info(f"✅ [AnalysisSummary] 使用原始内容（长度: {len(original_content)}）")
    return original_content
```

---

## 📊 预期修复效果

修复后，应该看到：
1. ✅ LLM成功调用并返回长文本（>500字符）
2. ✅ `extract_think_and_content` 正确提取实际内容，或使用原始内容
3. ✅ orchestrator使用LLM返回的内容，而不是后备摘要
4. ✅ 前端显示真正的生信分析报告，包含：
   - 统计概览
   - 关键生物标志物
   - 通路机制解读
   - 结论与建议
