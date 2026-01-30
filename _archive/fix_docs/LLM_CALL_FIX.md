# LLM调用逻辑统一修复方案

## 问题分析

### 规划阶段（正常工作）
- 位置：`gibh_agent/core/planner.py`
- 调用方式：`await self.llm_client.achat(messages, temperature=0.1, max_tokens=2048)`
- LLM客户端来源：通过构造函数传入 `SOPPlanner(tool_retriever, llm_client)`
- 在 `orchestrator.py` 中通过 `_get_llm_client()` 获取，如果失败则使用 `LLMClientFactory.create_default()`

### 数据诊断阶段（可能失败）
- 位置：`gibh_agent/agents/base_agent.py` 的 `_perform_data_diagnosis` 方法
- 调用方式：`await self.llm_client.achat(messages, temperature=0.3, max_tokens=1500)`
- LLM客户端来源：`self.llm_client`（在 Agent 初始化时传入）
- 问题：如果 `self.llm_client` 为 None 或未正确初始化，调用会失败

### 分析报告生成阶段（可能失败）
- 位置：`gibh_agent/agents/base_agent.py` 的 `_generate_analysis_summary` 方法
- 调用方式：`await self.llm_client.achat(messages, temperature=0.3, max_tokens=2500)`
- LLM客户端来源：`self.llm_client`（在 Agent 初始化时传入）
- 问题：如果 `self.llm_client` 为 None 或未正确初始化，调用会失败

## 修复方案

### 方案1：在调用前检查并创建LLM客户端（推荐）

在 `_perform_data_diagnosis` 和 `_generate_analysis_summary` 方法中，如果 `self.llm_client` 为 None，则使用 `LLMClientFactory.create_default()` 创建客户端。

### 方案2：统一使用工厂方法

所有地方都使用 `LLMClientFactory.create_default()` 创建客户端，而不是依赖 agent 的初始化。

## 实施步骤

1. 修改 `base_agent.py` 的 `_perform_data_diagnosis` 方法
2. 修改 `base_agent.py` 的 `_generate_analysis_summary` 方法
3. 添加回退机制：如果 `self.llm_client` 不可用，使用 `LLMClientFactory.create_default()`
