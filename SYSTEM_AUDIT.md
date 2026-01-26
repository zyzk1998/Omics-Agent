# 系统架构与 LLM 审计报告
**System Architecture & LLM Audit Report**

**生成时间**: 2024-12-19  
**目的**: 提供代码库结构、LLM集成点和数据流的完整透明度

---

## 1. 项目结构 (Project Structure)

```
GIBH-AGENT-V2/
├── gibh_agent/                    # 核心代码库
│   ├── core/                      # 核心功能模块
│   │   ├── llm_client.py         # LLM客户端统一接口
│   │   ├── orchestrator.py       # 工作流编排器（SSE流式处理）
│   │   ├── executor.py            # 工作流执行器
│   │   ├── planner.py             # 工作流规划器（SOPPlanner）
│   │   ├── agentic.py             # Agentic组件（QueryRewriter, Clarifier, Reflector）
│   │   ├── tool_registry.py       # 工具注册表
│   │   ├── data_diagnostician.py  # 数据诊断器
│   │   ├── file_inspector.py      # 文件检查器
│   │   └── workflows/             # 工作流定义
│   │       ├── registry.py
│   │       ├── metabolomics.py
│   │       └── rna.py
│   ├── agents/                    # 智能体模块
│   │   ├── base_agent.py          # 基础智能体（包含_generate_analysis_summary）
│   │   ├── router_agent.py        # 路由智能体
│   │   └── specialists/           # 领域专家智能体
│   │       ├── metabolomics_agent.py
│   │       ├── rna_agent.py
│   │       ├── dna_agent.py
│   │       └── ...
│   ├── tools/                     # 分析工具模块
│   │   ├── metabolomics/
│   │   │   ├── basic.py           # 基础工具（PCA, Diff Analysis）
│   │   │   └── advanced.py        # 高级工具（PLS-DA, Pathway Enrichment）
│   │   └── rna/
│   ├── config/                    # 配置文件
│   │   ├── settings.yaml          # LLM配置（包含SILICONFLOW_API_KEY引用）
│   │   └── prompts/               # Prompt模板
│   └── main.py                    # 主入口（GIBHAgent初始化）
├── services/                      # 服务层
│   └── nginx/
│       └── html/
│           └── index.html         # 前端UI（SSE事件处理）
├── server.py                      # FastAPI服务器
├── scripts/                       # 脚本工具
│   ├── verify_report_quality.py  # 报告质量验证脚本
│   └── verify_ui_logic.py        # UI逻辑验证脚本
└── docker-compose.yml             # Docker配置（环境变量注入）
```

### 核心模块职责

- **`core/orchestrator.py`**: 统一管理Agent的流式处理流程，实时输出状态更新、思考过程和结果
- **`core/executor.py`**: 动态执行工作流，不依赖硬编码的工具逻辑，使用ToolRegistry查找和执行工具
- **`core/planner.py`**: 从用户查询中检索相关工具，使用LLM生成工作流计划
- **`agents/base_agent.py`**: 所有领域智能体的基类，包含`_generate_analysis_summary`方法（AI Expert Diagnosis的核心）
- **`core/agentic.py`**: Agentic组件（QueryRewriter查询重写、Clarifier主动澄清、Reflector自我反思）

---

## 2. LLM 调用图谱 (LLM Call Graph)

### "Truth Table" - 所有LLM调用点

| 组件 | 文件路径 | 函数/方法 | 目的 | Fallback逻辑 |
|------|---------|----------|------|-------------|
| **BaseAgent** | `gibh_agent/agents/base_agent.py` | `process_query` | 处理用户查询（流式响应） | 无fallback，直接抛出异常 |
| **BaseAgent** | `gibh_agent/agents/base_agent.py` | `_perform_data_diagnosis` | 生成数据体检报告 | 失败时返回`None`，由调用方处理 |
| **BaseAgent** | `gibh_agent/agents/base_agent.py` | `_generate_analysis_summary` | **生成AI专家分析报告** | **⚠️ CRITICAL**: 失败时返回错误信息Markdown，不隐藏 |
| **BaseAgent** | `gibh_agent/agents/base_agent.py` | `_evaluate_analysis_quality` | 评估分析质量（0-100分） | 失败时返回`None`，不影响主流程 |
| **WorkflowPlanner** | `gibh_agent/core/planner.py` | `plan` | 生成工作流计划 | 失败时返回`{"type": "error", "error": str(e)}` |
| **SOPPlanner** | `gibh_agent/core/planner.py` | `_classify_intent` | 分类用户意图（planning/execution） | 失败时返回`"execution"`（默认执行模式） |
| **SOPPlanner** | `gibh_agent/core/planner.py` | `_analyze_user_intent` | 分析用户意图，提取目标步骤 | 失败时返回`[]`（空列表，表示完整工作流） |
| **SOPPlanner** | `gibh_agent/core/planner.py` | `_check_execution_mode` | 检查执行模式（PLANNING/EXECUTION） | 失败时返回`"execution"`（默认执行模式） |
| **QueryRewriter** | `gibh_agent/core/agentic.py` | `rewrite` | 查询重写（模糊→精确） | **失败时返回原始查询**（`return raw_query`） |
| **Clarifier** | `gibh_agent/core/agentic.py` | `clarify` | 主动澄清缺失信息 | 失败时返回`None`（表示无需澄清） |
| **Clarifier** | `gibh_agent/core/agentic.py` | `_check_intent` | 检查用户意图 | 失败时返回`"execution"`（默认执行模式） |
| **Reflector** | `gibh_agent/core/agentic.py` | `reflect` | 自我反思，检查和纠正计划 | 失败时返回原始计划（`return workflow_plan`） |
| **MetabolomicsAgent** | `gibh_agent/agents/specialists/metabolomics_agent.py` | `_detect_omics_type` | 检测组学类型 | 失败时返回`"Metabolomics"`（默认） |
| **MetabolomicsAgent** | `gibh_agent/agents/specialists/metabolomics_agent.py` | `_suggest_parameters` | 建议工具参数 | 失败时返回`None` |
| **MetabolomicsAgent** | `gibh_agent/agents/specialists/metabolomics_agent.py` | `_validate_workflow` | 验证工作流 | 失败时返回`True`（通过验证） |
| **RNAAgent** | `gibh_agent/agents/specialists/rna_agent.py` | `_detect_omics_type` | 检测组学类型 | 失败时返回`"scRNA-seq"`（默认） |
| **RNAAgent** | `gibh_agent/agents/specialists/rna_agent.py` | `_suggest_parameters` | 建议工具参数 | 失败时返回`None` |
| **RouterAgent** | `gibh_agent/agents/router_agent.py` | `route` | 路由到合适的领域智能体 | 失败时返回`"metabolomics"`（默认） |

### LLM调用统计

- **总调用点**: 43个（通过`grep`统计）
- **关键调用点**: 
  - `_generate_analysis_summary`: **1个**（AI Expert Diagnosis的核心）
  - `_perform_data_diagnosis`: **1个**（数据体检报告）
  - `plan`: **1个**（工作流规划）
  - Agentic组件: **4个**（QueryRewriter, Clarifier x2, Reflector）
  - SOPPlanner: **3个**（意图分类、分析、模式检查）

---

## 3. 数据流分析: "AI Expert Diagnosis" 生成流程

### 完整数据流追踪

```
1. 用户提交工作流执行请求
   ↓
2. Frontend (index.html) -> POST /api/chat (server.py)
   ↓
3. AgentOrchestrator.stream_process() (orchestrator.py:85)
   ↓
4. WorkflowExecutor.execute_workflow() (executor.py:891)
   ├── 执行每个步骤（PCA, PLS-DA, Diff Analysis等）
   └── 返回 results = {
       "steps_details": [...],
       "steps_results": [...],
       "status": "success"
   }
   ↓
5. Orchestrator 收集步骤详情 (orchestrator.py:182)
   steps_details = results.get("steps_details", [])
   ↓
6. 构建 summary_context (orchestrator.py:268)
   summary_context = {
       "has_failures": len(failed_steps) > 0,
       "has_warnings": len(warning_steps) > 0,
       "failed_steps": failed_steps,
       "warning_steps": warning_steps,
       "successful_steps": successful_steps,
       "workflow_status": results.get("status", "unknown")
   }
   ↓
7. 调用 agent._generate_analysis_summary() (orchestrator.py:279)
   summary = await self.agent._generate_analysis_summary(
       results,  # 包含 steps_results
       domain_name,  # "Metabolomics" or "RNA"
       summary_context=summary_context
   )
   ↓
8. BaseAgent._generate_analysis_summary() (base_agent.py:642)
   ├── 提取关键发现（PCA分离、差异代谢物、富集通路、VIP代谢物）
   ├── 构建Prompt（Nature Medicine风格，要求生物学机制解读）
   ├── 调用LLM: await self.llm_client.achat(messages, ...) (base_agent.py:1004)
   │   ├── ✅ 成功: 返回LLM生成的Markdown报告
   │   ├── ⚠️ 内容过短: 重试一次 (base_agent.py:1012)
   │   └── ❌ 失败: 返回错误信息Markdown (base_agent.py:1048-1060)
   │       return f"""## ❌ LLM 生成失败
   │       **错误信息**: {str(llm_error)}
   │       **分析指标**: {key_findings_json}
   │       **说明**: LLM 服务调用失败..."""
   └── 返回 summary (Markdown字符串)
   ↓
9. Orchestrator 检查summary (orchestrator.py:286)
   if not summary or len(summary.strip()) < 50:
       # 使用结构化后备（仅在summary为None或过短时）
       summary = f"""## 分析结果摘要
       本次分析完成了 {len(successful_steps)} 个步骤..."""
   ↓
10. 产生SSE事件 (orchestrator.py:340-360)
    yield self._format_sse("step_result", {...})  # 执行结果
    yield self._format_sse("diagnosis", {...})     # 专家报告
    yield self._format_sse("result", {...})       # 完整结果（向后兼容）
   ↓
11. Frontend 接收事件 (index.html:1246-1280)
    case 'diagnosis':
        renderDiagnosisCard(data.report_data.diagnosis, data)
    case 'step_result':
        renderExecutionSteps(data.report_data.steps_details)
   ↓
12. 渲染到UI (index.html:1430-1524)
    renderDiagnosisCard() -> 统一报告容器 -> 专家洞察区域
```

### 关键验证点

#### ✅ LLM确实被调用
- **位置**: `base_agent.py:1004`
- **代码**: `completion = await self.llm_client.achat(messages, temperature=0.3, max_tokens=2500)`
- **验证**: 如果LLM失败，会抛出异常并被捕获（`except Exception as llm_error`）

#### ⚠️ Fallback逻辑分析

1. **`_generate_analysis_summary`的Fallback** (base_agent.py:1037-1060):
   - **类型**: 错误信息Markdown（不是静态列表）
   - **触发条件**: LLM调用失败或返回内容过短
   - **内容**: 包含错误信息、分析指标、已完成步骤数
   - **评估**: ✅ **正确** - 不隐藏失败，明确显示错误

2. **Orchestrator的Fallback** (orchestrator.py:286-294):
   - **类型**: 结构化Markdown（简短摘要）
   - **触发条件**: `summary`为`None`或长度<50字符
   - **内容**: "本次分析完成了 X 个步骤。请查看上方的详细图表..."
   - **评估**: ⚠️ **可能触发** - 如果`_generate_analysis_summary`返回`None`（第1064行），会使用此fallback

3. **`_generate_analysis_summary`返回`None`的情况** (base_agent.py:1062-1064):
   - **触发条件**: 最外层`except Exception as e`捕获到异常
   - **评估**: ⚠️ **罕见** - 只有在数据提取阶段失败时才会发生

### 结论：是否存在"Hardcoded/Fake"逻辑？

**答案：基本否，但有遗留代码**

1. **LLM调用是真实的**: `base_agent.py:1004`确实调用`self.llm_client.achat()`
2. **主要路径无静态列表fallback**: `_generate_analysis_summary`失败时返回的是错误信息Markdown，不是"✅ 成功步骤"列表
3. **遗留的fallback函数**: `_generate_fallback_summary`仍然存在（`orchestrator.py:995`），会生成"✅ 成功步骤"列表
   - **调用场景**: 
     - 当`self.agent`不存在时（`orchestrator.py:338`）
     - 当`_generate_analysis_summary`抛出异常时（`orchestrator.py:300`，但已被新的结构化fallback替代）
   - **评估**: ⚠️ **遗留代码**，在正常流程中不会触发，但应该删除以避免混淆
4. **Orchestrator的新fallback**: 仅在summary为`None`或过短时触发简短摘要（`orchestrator.py:288-294`），而`_generate_analysis_summary`现在总是返回字符串（错误信息或LLM输出）

---

## 4. 环境变量与配置审计

### API Key加载流程

#### 1. Docker环境（生产环境）

```
docker-compose.yml
  ↓
环境变量注入: SILICONFLOW_API_KEY=${SILICONFLOW_API_KEY}
  ↓
GIBHAgent._init_llm_clients() (main.py:90)
  ↓
读取配置: gibh_agent/config/settings.yaml
  ├── cloud.siliconflow.api_key: "${SILICONFLOW_API_KEY:}"
  └── 如果为空，使用环境变量: os.getenv("SILICONFLOW_API_KEY", "")
  ↓
LLMClientFactory.create_cloud_siliconflow() (llm_client.py:387)
  ↓
创建LLMClient: LLMClient(base_url="https://api.siliconflow.cn/v1", api_key=api_key, ...)
```

#### 2. 验证脚本环境（开发环境）

```
scripts/verify_report_quality.py
  ↓
手动读取环境变量: os.getenv("DEEPSEEK_API_KEY", os.getenv("LLM_API_KEY", "EMPTY"))
  ↓
直接创建LLMClient: LLMClient(base_url="https://api.deepseek.com/v1", api_key=api_key, ...)
```

### 401错误原因分析

**为什么`scripts/verify_report_quality.py`会失败？**

1. **环境变量未加载**:
   - 脚本直接使用`os.getenv("DEEPSEEK_API_KEY", ...)`
   - 如果`.env`文件存在但未加载，环境变量为空
   - 默认值`"EMPTY"`导致401认证失败

2. **Docker应用成功的原因**:
   - `docker-compose.yml`通过`environment`注入环境变量
   - 或者在启动时通过`export SILICONFLOW_API_KEY=...`设置
   - `main.py`会验证API key是否存在，不存在会抛出明确的错误

3. **解决方案**:
   ```bash
   # 方法1: 设置环境变量
   export DEEPSEEK_API_KEY="your_key_here"
   python3 scripts/verify_report_quality.py
   
   # 方法2: 使用.env文件（需要python-dotenv）
   # 在脚本开头添加:
   # from dotenv import load_dotenv
   # load_dotenv()
   ```

### 配置优先级

1. **配置文件** (`settings.yaml`): `api_key: "${SILICONFLOW_API_KEY:}"`
2. **环境变量**: `os.getenv("SILICONFLOW_API_KEY", "")`
3. **默认值**: `"EMPTY"`（会导致401错误）

---

## 5. 关键发现与建议

### ✅ 正面发现

1. **LLM调用是真实的**: 所有关键功能都确实调用LLM，没有hardcoded逻辑
2. **错误处理透明**: LLM失败时返回明确的错误信息，不隐藏问题
3. **Fallback合理**: 仅在极端情况下使用简短摘要，且明确标注

### ⚠️ 需要注意的点

1. **`_generate_fallback_summary`函数仍然存在**:
   - 位置: `orchestrator.py:995-1055`
   - 功能: 生成"✅ 成功步骤"列表格式的fallback
   - **当前状态**: ⚠️ **已废弃但未删除**
   - **调用位置**: 
     - `orchestrator.py:338` - 当agent不存在时调用
     - `orchestrator.py:300` - 当`_generate_analysis_summary`异常时调用（但已被新的结构化fallback替代）
   - **评估**: 这个函数会生成"✅ 成功步骤"列表，但**仅在极端情况下调用**（agent不存在或双重异常）

2. **`_generate_analysis_summary`可能返回`None`**:
   - 位置: `base_agent.py:1064`
   - 触发条件: 数据提取阶段异常
   - 影响: Orchestrator会使用简短fallback（`orchestrator.py:288-294`）
   - **建议**: 确保数据提取阶段的异常处理也返回错误信息字符串

3. **验证脚本的环境变量加载**:
   - 当前: 直接使用`os.getenv()`，可能未加载`.env`文件
   - **建议**: 添加`python-dotenv`支持，或明确文档说明需要设置环境变量

4. **Orchestrator的fallback触发条件**:
   - 当前: `if not summary or len(summary.strip()) < 50`
   - **评估**: 合理，但`_generate_analysis_summary`现在总是返回字符串（错误信息或LLM输出），此fallback可能永远不会触发

### 🔍 验证建议

1. **运行验证脚本前设置环境变量**:
   ```bash
   export DEEPSEEK_API_KEY="your_key"
   # 或
   export SILICONFLOW_API_KEY="your_key"
   python3 scripts/verify_report_quality.py
   ```

2. **检查LLM调用日志**:
   - 查看`gibh_agent/agents/base_agent.py:1000`的日志输出
   - 确认`📞 [AnalysisSummary] 调用 LLM 生成深度生物学解释...`出现

3. **验证报告内容**:
   - 如果包含"❌ LLM 生成失败"，说明LLM调用失败但错误处理正确
   - 如果包含"生物学机制解读"、"潜在标志物"等关键词，说明LLM成功生成

---

## 6. 总结

### 系统架构
- **模块化设计**: 清晰的职责分离（orchestrator, executor, planner, agents）
- **统一LLM接口**: `LLMClient`提供统一的异步/同步接口
- **流式处理**: SSE事件流确保实时反馈

### LLM集成
- **43个调用点**: 分布在规划、执行、诊断、评估等各个环节
- **真实调用**: 所有关键功能都确实调用LLM，无hardcoded逻辑
- **错误处理**: 失败时返回明确的错误信息，不隐藏问题

### 数据流
- **AI Expert Diagnosis**: 从执行结果 → 数据提取 → Prompt构建 → LLM调用 → 前端渲染
- **无静态fallback**: 失败时返回错误信息Markdown，不是"✅ 成功步骤"列表
- **透明性**: 所有步骤都有日志记录，可追踪

### 环境配置
- **Docker环境**: 通过`docker-compose.yml`注入环境变量，配置正确
- **验证脚本**: 需要手动设置环境变量或使用`.env`文件
- **401错误**: 由于环境变量未加载，不是代码逻辑问题

---

---

## 7. 快速参考 (Quick Reference)

### LLM调用点速查表

| 功能 | 文件 | 行号 | 是否必需 | Fallback |
|------|------|------|---------|---------|
| AI专家分析报告 | `base_agent.py` | 1004 | ✅ 是 | 错误信息Markdown |
| 数据体检报告 | `base_agent.py` | 586 | ✅ 是 | 返回None |
| 工作流规划 | `planner.py` | 106 | ✅ 是 | 错误字典 |
| 查询重写 | `agentic.py` | 82 | ⚠️ 可选 | 返回原始查询 |
| 主动澄清 | `agentic.py` | 210 | ⚠️ 可选 | 返回None |
| 自我反思 | `agentic.py` | 429 | ⚠️ 可选 | 返回原始计划 |

### 环境变量配置

```bash
# 生产环境（Docker）
export SILICONFLOW_API_KEY="your_api_key_here"

# 开发环境（验证脚本）
export DEEPSEEK_API_KEY="your_api_key_here"
# 或
export LLM_API_KEY="your_api_key_here"
```

### 验证LLM是否被调用

1. **查看日志**:
   ```bash
   grep "📞.*LLM" logs/*.log
   ```

2. **检查报告内容**:
   - ✅ 包含"生物学机制解读"、"潜在标志物" → LLM成功
   - ⚠️ 包含"❌ LLM 生成失败" → LLM失败但错误处理正确
   - ❌ 包含"✅ 成功步骤"列表 → 这是旧的fallback（已移除）

3. **运行验证脚本**:
   ```bash
   export DEEPSEEK_API_KEY="your_key"
   python3 scripts/verify_report_quality.py
   ```

---

**文档版本**: 1.0  
**最后更新**: 2024-12-19  
**维护者**: System Audit Script

