# 🔍 专业调试提示：AI专家分析报告生成失败 + 工具执行错误诊断

## 📋 问题概述

### 问题1：AI专家分析报告降级处理（LLM调用失败）
**现象**：前端显示的是降级处理生成的简单分析报告，而不是LLM生成的深度生物学解释报告。

**目标**：实现分析报告根据（数据诊断报告 + 所有步骤的可视化/非可视化执行结果），给出专业的生信分析以及建议。

### 问题2：工具执行失败
**现象**：
- PLS-DA 分析失败：分组列 'Sample' 不存在于数据中
- 通路富集分析失败：分组列 'Sample' 不存在于数据中
- 数据格式问题：所有列都是数值型，没有分组信息列

### 问题3：控制台日志与前端错误信息不一致
**现象**：后端日志显示的信息与前端页面显示的错误信息不匹配。

---

## 🎯 调试目标

1. **彻底诊断AI专家分析报告LLM调用失败的根本原因**
   - 确认LLM API调用是否真正执行
   - 确认API密钥、网络连接、请求参数是否正确
   - 确认错误类型和具体错误信息
   - 确认错误是否被正确捕获和记录

2. **修复工具执行失败问题**
   - 修复分组列检测逻辑
   - 处理数据格式问题（所有列都是数值型的情况）
   - 改进错误信息传递机制

3. **统一前后端错误信息显示**
   - 确保后端错误信息正确传递到前端
   - 确保前端正确解析和显示错误信息

---

## 🔬 需要收集的调试信息

### 1. AI专家分析报告LLM调用失败诊断

#### A. 后端日志检查点
**位置**：`gibh_agent/agents/base_agent.py` 第1307-1432行

**需要检查的日志**：
```bash
# 检查以下日志是否存在：
1. "📞 [AnalysisSummary] 开始LLM调用，max_tokens=2500..."
2. "✅ [AnalysisSummary] LLM调用完成，开始解析响应..."
3. "⚠️ [AnalysisSummary] self.llm_client 不可用，使用 LLMClientFactory.create_default()"
4. "❌ [AnalysisSummary] LLM 调用失败"
5. "⚠️ [AnalysisSummary] LLM调用失败，生成基于数据的简单分析报告"
```

**关键检查项**：
- [ ] LLM客户端是否正确初始化（`llm_client_to_use.base_url`、`llm_client_to_use.model`）
- [ ] API密钥是否设置（`api_key_set: true/false`）
- [ ] 错误类型（`error_type`）：`AttributeError`、`TimeoutError`、`ConnectionError`、`APIError`等
- [ ] 错误消息（`error_message`）：具体的错误描述
- [ ] 请求参数：`messages`数量、`system message`长度、`user message`长度
- [ ] 执行上下文：成功步骤数、失败步骤数、关键指标是否为空

#### B. 前端SSE错误事件检查
**位置**：`gibh_agent/core/orchestrator.py` 第439-452行

**需要检查**：
- [ ] 是否收到SSE `error`事件
- [ ] `error`事件的内容：`error_type`、`error_message`、`details`、`context`、`possible_causes`
- [ ] 前端是否正确解析和显示错误信息（`services/nginx/html/index.html` 的 `handleErrorEvent` 函数）

#### C. LLM客户端配置检查
**位置**：`.env` 文件、`gibh_agent/core/llm_client.py`

**需要检查**：
- [ ] `SILICONFLOW_API_KEY` 是否设置且有效
- [ ] `SILICONFLOW_MODEL` 是否设置
- [ ] LLM客户端工厂方法 `LLMClientFactory.create_default()` 是否正常工作
- [ ] 网络连接是否正常（能否访问SiliconFlow API）

#### D. 请求内容检查
**位置**：`gibh_agent/agents/base_agent.py` 第1200-1300行（构建messages的部分）

**需要检查**：
- [ ] `key_findings_json` 是否为空或格式错误
- [ ] `messages` 的内容是否过长（超过API限制）
- [ ] `system message` 和 `user message` 的内容是否合理
- [ ] 是否有特殊字符导致JSON序列化失败

### 2. 工具执行失败诊断

#### A. PLS-DA分析失败
**位置**：`gibh_agent/tools/metabolomics/advanced.py` 的 `run_plsda` 函数

**需要检查**：
- [ ] 数据文件的实际列名（使用 `pandas.read_csv` 读取并打印 `df.columns.tolist()`）
- [ ] 分组列检测逻辑（`_detect_group_column_from_file` 函数）
- [ ] 错误信息是否正确传递到前端

**预期修复**：
- 如果数据中没有分组列，应该：
  1. 自动检测可能的分组列（非数值列，唯一值在2-10之间）
  2. 如果找不到，提示用户数据格式要求
  3. 如果所有列都是数值型，提示用户需要添加分组列

#### B. 通路富集分析失败
**位置**：`gibh_agent/tools/metabolomics/pathway_enrichment.py`

**需要检查**：
- [ ] 分组列检测逻辑（同上）
- [ ] 代谢物ID映射是否正确
- [ ] 通路数据库是否可用

#### C. 数据格式验证
**位置**：`gibh_agent/core/executor.py` 的 `_detect_group_column_from_file` 函数

**需要检查**：
- [ ] 数据文件读取是否正确
- [ ] 列类型检测逻辑是否正确
- [ ] 分组列推荐逻辑是否合理

### 3. 前后端错误信息同步

#### A. 后端错误信息格式
**位置**：`gibh_agent/core/executor.py` 的 `execute_workflow` 方法

**需要检查**：
- [ ] 工具执行失败时，`step_result` 中的 `error` 字段是否包含详细错误信息
- [ ] 错误信息格式是否统一（包含错误类型、错误消息、可能原因、建议）

#### B. 前端错误信息显示
**位置**：`services/nginx/html/index.html` 的 `renderExecutionSteps` 函数

**需要检查**：
- [ ] 前端是否正确解析 `step_result.error` 字段
- [ ] 错误信息是否正确显示在UI上
- [ ] 错误信息的格式是否用户友好

---

## 🛠️ 调试步骤

### 步骤1：收集后端日志
```bash
# 1. 查看完整的后端日志（重点关注AI专家分析报告生成部分）
docker logs <container_name> 2>&1 | grep -A 50 "AnalysisSummary\|LLM调用"

# 2. 查看工具执行日志
docker logs <container_name> 2>&1 | grep -A 20 "PLS-DA\|通路富集\|group_column"

# 3. 查看LLM客户端初始化日志
docker logs <container_name> 2>&1 | grep -A 10 "LLMClient\|SILICONFLOW\|api_key"
```

### 步骤2：检查环境变量
```bash
# 检查.env文件
cat .env | grep SILICONFLOW

# 检查容器内的环境变量
docker exec <container_name> env | grep SILICONFLOW
```

### 步骤3：测试LLM API连接
```python
# 创建测试脚本 test_llm_connection.py
import asyncio
import os
from gibh_agent.core.llm_client import LLMClientFactory

async def test_llm():
    try:
        client = LLMClientFactory.create_default()
        print(f"✅ LLM客户端创建成功")
        print(f"   - base_url: {client.base_url}")
        print(f"   - model: {client.model}")
        print(f"   - api_key: {'已设置' if hasattr(client, 'api_key') and client.api_key else '未设置'}")
        
        # 测试简单调用
        messages = [
            {"role": "system", "content": "You are a helpful assistant."},
            {"role": "user", "content": "Say hello in Chinese."}
        ]
        completion = await client.achat(messages, temperature=0.3, max_tokens=100)
        print(f"✅ LLM调用成功: {completion.choices[0].message.content}")
    except Exception as e:
        print(f"❌ LLM调用失败: {type(e).__name__}: {str(e)}")
        import traceback
        traceback.print_exc()

asyncio.run(test_llm())
```

### 步骤4：检查数据文件格式
```python
# 创建测试脚本 check_data_format.py
import pandas as pd
import sys

file_path = sys.argv[1] if len(sys.argv) > 1 else "guest/20260128_112036/cow_diet.csv"

try:
    df = pd.read_csv(file_path, nrows=10)
    print(f"✅ 数据文件读取成功")
    print(f"   - 形状: {df.shape}")
    print(f"   - 列名: {df.columns.tolist()}")
    print(f"   - 列类型:")
    for col in df.columns:
        dtype = df[col].dtype
        n_unique = df[col].nunique()
        print(f"     - {col}: {dtype} (唯一值: {n_unique})")
    
    # 检查分组列
    non_numeric_cols = df.select_dtypes(include=['object']).columns.tolist()
    print(f"\n   - 非数值列: {non_numeric_cols}")
    
    # 检查可能的分组列（唯一值在2-10之间）
    potential_group_cols = []
    for col in df.columns:
        n_unique = df[col].nunique()
        if 2 <= n_unique <= 10:
            potential_group_cols.append((col, n_unique))
    print(f"   - 可能的分组列: {potential_group_cols}")
    
except Exception as e:
    print(f"❌ 数据文件读取失败: {e}")
    import traceback
    traceback.print_exc()
```

### 步骤5：检查前端SSE事件
```javascript
// 在浏览器控制台中运行
// 监听SSE error事件
const eventSource = new EventSource('/api/chat');
eventSource.addEventListener('error', (e) => {
    console.log('SSE Error Event:', e);
    console.log('Error Data:', JSON.parse(e.data));
});
```

---

## 🔧 预期修复方案

### 修复1：AI专家分析报告LLM调用失败

**可能原因和解决方案**：

1. **API密钥问题**
   - **检查**：确认 `.env` 中的 `SILICONFLOW_API_KEY` 是否正确设置
   - **修复**：如果未设置或无效，更新API密钥

2. **网络连接问题**
   - **检查**：确认服务器能否访问SiliconFlow API
   - **修复**：检查网络配置、防火墙设置、代理设置

3. **请求超时**
   - **检查**：确认请求是否超时（查看错误类型是否为 `TimeoutError`）
   - **修复**：增加超时时间或优化请求内容

4. **请求内容过长**
   - **检查**：确认 `messages` 的总长度是否超过API限制
   - **修复**：截断过长的内容，只保留关键信息

5. **LLM客户端初始化失败**
   - **检查**：确认 `LLMClientFactory.create_default()` 是否正常工作
   - **修复**：确保环境变量正确加载，LLM客户端正确初始化

6. **错误信息未正确传递**
   - **检查**：确认 `self.context["last_llm_error"]` 是否正确设置
   - **修复**：确保错误信息通过SSE正确传递到前端

### 修复2：工具执行失败

**解决方案**：

1. **改进分组列检测逻辑**
   - 自动检测可能的分组列（非数值列，唯一值在2-10之间）
   - 如果找不到，提供清晰的错误提示和数据格式要求
   - 如果所有列都是数值型，提示用户需要添加分组列

2. **改进错误信息格式**
   - 统一错误信息格式（包含错误类型、错误消息、可能原因、建议）
   - 确保错误信息正确传递到前端

3. **数据格式验证**
   - 在执行前验证数据格式
   - 提供数据格式示例和修复建议

### 修复3：前后端错误信息同步

**解决方案**：

1. **统一错误信息格式**
   - 后端统一错误信息格式（JSON结构）
   - 前端统一解析和显示逻辑

2. **增强错误信息传递**
   - 确保所有错误信息都通过SSE传递到前端
   - 确保前端正确解析和显示错误信息

---

## 📊 成功标准

### AI专家分析报告
- [ ] LLM API调用成功（后端日志显示 "✅ [AnalysisSummary] LLM调用完成"）
- [ ] 返回的报告长度 > 500字符（不是降级报告）
- [ ] 报告包含深度生物学机制解读
- [ ] 报告包含基于执行结果的建议

### 工具执行
- [ ] PLS-DA分析成功执行或提供清晰的错误提示
- [ ] 通路富集分析成功执行或提供清晰的错误提示
- [ ] 错误信息清晰、用户友好

### 错误信息同步
- [ ] 后端日志和前端显示的错误信息一致
- [ ] 错误信息包含足够的调试信息
- [ ] 错误信息格式统一、用户友好

---

## 🚀 下一步行动

1. **立即执行**：运行上述调试步骤，收集所有必要信息
2. **分析结果**：根据收集的信息，确定根本原因
3. **实施修复**：根据根本原因，实施相应的修复方案
4. **验证修复**：重新测试，确认问题已解决

---

## 📝 注意事项

1. **日志级别**：确保日志级别设置为 `INFO` 或 `DEBUG`，以便收集足够的调试信息
2. **敏感信息**：在日志中不要输出完整的API密钥，只显示是否设置
3. **错误处理**：确保所有错误都被正确捕获和记录
4. **用户体验**：即使LLM调用失败，也要提供有意义的降级报告，而不是技术错误信息

---

**生成时间**：2025-01-28
**问题优先级**：高
**预计修复时间**：2-4小时
