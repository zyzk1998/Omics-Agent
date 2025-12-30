# GIBH-AGENT 项目继续讨论提示词

## 📋 项目背景

我正在开发 **GIBH-AGENT**，一个基于多模态大模型与微服务架构的单细胞生信分析智能体平台。项目正在从单体脚本重构为**分层多智能体系统**。

**重要**: 新架构代码在独立目录 `/home/ubuntu/GIBH-AGENT-V2/`，原项目 `/home/ubuntu/GIBH-AGENT/` 保持不变。

## 🏗️ 当前项目状态

### 原项目（保持不变）
**位置**: `/home/ubuntu/GIBH-AGENT/services/`
- **架构**: FastAPI + Celery + Redis + vLLM
- **功能**: 单细胞转录组分析（scRNA-seq）
- **双脑设计**: 
  - 逻辑大脑：Qwen3-Coder-30B-AWQ（工作流规划）
  - 视觉大脑：Qwen3-VL-8B（多模态对话）
- **核心文件**:
  - `services/api/src/agent.py`: BioBlendAgent（旧版单体智能体）
  - `services/api/src/celery_app.py`: Celery 任务处理
  - `services/api/src/main.py`: FastAPI 入口
  - `services/api/src/skills/scanpy_local.py`: Scanpy 分析逻辑

### 新架构（独立开发）
**位置**: `/home/ubuntu/GIBH-AGENT-V2/gibh_agent/`
- **核心基础设施**:
  - `core/llm_client.py`: 统一 LLM 客户端（支持本地/云端切换）✅
  - `core/prompt_manager.py`: 提示管理器（Jinja2 模板）✅
  - `core/dispatcher.py`: 任务分发器（本地/Slurm/SSH）✅
- **智能体系统**:
  - `agents/base_agent.py`: 基础智能体抽象类 ✅
  - `agents/router_agent.py`: 路由智能体（识别组学类型）✅
  - `agents/specialists/rna_agent.py`: 转录组智能体（已重构）✅
  - `agents/specialists/dna_agent.py` 等6个：占位符 ⏳
- **工具类**:
  - `tools/cellranger_tool.py`: Cell Ranger 脚本生成器 ✅
  - `tools/scanpy_tool.py`: Scanpy 工作流脚本生成器 ✅
- **配置**:
  - `config/settings.yaml`: 统一配置文件 ✅

## 🎯 重构目标

### 核心改进
1. **支持7种组学模态**: Transcriptomics, Genomics, Epigenomics, Metabolomics, Proteomics, Spatial Omics, Imaging
2. **处理TB级数据**: 控制平面（智能体）vs 数据平面（HPC）分离
3. **灵活LLM切换**: 本地（vLLM）和云端（DeepSeek/SiliconFlow）无缝切换
4. **模块化架构**: 易于扩展和维护

### 关键设计原则
1. **控制平面 vs 数据平面分离**
   - 智能体只处理文件路径（字符串），不处理二进制数据
   - 大数据处理通过 TaskDispatcher 提交到 HPC
2. **统一接口**
   - 所有智能体继承 `BaseAgent`
   - 统一的 LLM 调用接口
3. **配置驱动**
   - LLM 配置可切换（本地/云端）
   - 任务分发方式可配置

## 📊 架构流程

```
用户查询
    ↓
RouterAgent (路由智能体)
    - 分析用户自然语言
    - 识别组学类型（7种模态）
    - 识别用户意图（分析/可视化/解释）
    - 路由到对应的领域智能体
    ↓
Domain Agents (领域智能体)
    ├── RNAAgent ✅ (转录组，已实现)
    ├── DNAAgent ⏳ (基因组，占位符)
    ├── EpigenomicsAgent ⏳ (表观遗传，占位符)
    ├── MetabolomicsAgent ⏳ (代谢组，占位符)
    ├── ProteomicsAgent ⏳ (蛋白质组，占位符)
    ├── SpatialAgent ⏳ (空间组学，占位符)
    └── ImagingAgent ⏳ (影像分析，占位符)
    ↓
Tools (工具类)
    - 只生成脚本，不执行
    - CellRangerTool: 生成 Cell Ranger 脚本
    - ScanpyTool: 生成 Scanpy 工作流脚本
    ↓
TaskDispatcher (任务分发器)
    - 提交脚本到 HPC/服务器
    - 支持：本地执行、Slurm 提交、SSH 远程提交
```

## 🔧 关键技术实现

### LLM 客户端切换
```python
from gibh_agent.core.llm_client import LLMClientFactory

# 本地模型（vLLM）
client = LLMClientFactory.create_local_vllm("qwen3-vl")

# 云端模型（DeepSeek）
client = LLMClientFactory.create_cloud_deepseek()
```

### 任务分发（只处理文件路径）
```python
# 智能体生成脚本（只传路径，不读文件）
script = cellranger_tool.generate_count_script(
    fastq_dir="/data/fastq",  # 只传路径字符串
    sample_id="sample1"
)

# TaskDispatcher 提交执行
task_info = await dispatcher.submit_script(script)
```

### 路由决策
```python
# RouterAgent 自动路由
route_result = await router.process_query(query, files)
# 返回: {
#   "modality": "transcriptomics",
#   "routing": "rna_agent",
#   "confidence": 0.95
# }
```

## 📁 关键文件位置

### 新架构代码（独立目录）
**位置**: `/home/ubuntu/GIBH-AGENT-V2/gibh_agent/`
- `core/llm_client.py`: LLM 客户端
- `core/prompt_manager.py`: 提示管理器
- `core/dispatcher.py`: 任务分发器
- `agents/router_agent.py`: 路由智能体
- `agents/specialists/rna_agent.py`: 转录组智能体
- `tools/cellranger_tool.py`: Cell Ranger 工具
- `tools/scanpy_tool.py`: Scanpy 工具
- `main.py`: 主入口
- `config/settings.yaml`: 配置文件

### 原项目代码（保持不变）
**位置**: `/home/ubuntu/GIBH-AGENT/services/api/src/`
- `agent.py`: 旧版 BioBlendAgent
- `main.py`: FastAPI 入口
- `celery_app.py`: Celery 任务

### 文档
**位置**: `/home/ubuntu/GIBH-AGENT-V2/`
- `CONTINUE_PROJECT_PROMPT.md`: 本文档（继续讨论提示词）
- `PROJECT_SUMMARY.md`: 项目总结
- `REFACTORING_PLAN.md`: 详细重构方案
- `IMPROVEMENT_ANALYSIS.md`: 改进分析
- `SETUP.md`: 设置指南

## ⏳ 当前状态

### ✅ 已完成（Phase 1）
- 核心基础设施（LLMClient, PromptManager, TaskDispatcher）
- 基础智能体框架（BaseAgent, RouterAgent）
- RNAAgent（转录组智能体，重构自现有代码）
- 工具类（CellRangerTool, ScanpyTool）
- 配置文件结构
- 其他6个领域智能体占位符

### ⏳ 待完成
- [ ] 完善 RNAAgent 的完整工作流执行逻辑
- [ ] 与现有 FastAPI 服务集成（可选，保持独立也可以）
- [ ] 实现其他领域智能体（DNAAgent 等）
- [ ] 添加单元测试
- [ ] 性能优化

## 🎯 下一步讨论方向

1. **完善现有功能**: RNAAgent 的完整工作流、错误处理
2. **集成现有服务**: 如何与 FastAPI 服务集成，保持向后兼容
3. **扩展领域智能体**: 实现 DNAAgent 或其他智能体
4. **工具类扩展**: 添加更多生信工具（GATK, MACS2 等）
5. **性能优化**: 路由优化、缓存机制
6. **测试和部署**: 单元测试、集成测试、生产部署

## 💡 关键设计决策

1. **为什么只处理文件路径？**
   - TB级FASTQ文件无法加载到内存
   - 智能体只需要知道文件位置，不需要读取内容
   - 数据处理在HPC集群完成

2. **为什么使用脚本生成？**
   - 脚本可以在HPC集群执行
   - 支持Slurm等任务调度系统
   - 可以重试和监控

3. **为什么使用 OpenAI SDK？**
   - 标准接口，兼容性好
   - 支持本地（vLLM）和云端（DeepSeek等）
   - 社区支持好

## 📝 使用示例

```python
from gibh_agent import create_agent

# 创建智能体
agent = create_agent("gibh_agent/config/settings.yaml")

# 处理查询
result = await agent.process_query(
    query="帮我分析一下这个单细胞数据",
    uploaded_files=[{"name": "sample.h5ad", "path": "/data/sample.h5ad"}]
)

# 结果包含路由信息和处理结果
print(result["routing_info"])  # 路由决策
print(result)  # 智能体响应
```

## ⚠️ 重要说明

- ✅ **新架构独立**: `/home/ubuntu/GIBH-AGENT-V2/` 是独立的新项目
- ✅ **原项目不变**: `/home/ubuntu/GIBH-AGENT/` 完全不受影响
- ✅ **可以独立开发**: 可以单独提交到新的 Git 仓库
- ✅ **建议先理解**: 阅读文档后再开始开发

---

**请基于以上信息继续讨论项目改进和实现细节。**
