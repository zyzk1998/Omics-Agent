# GIBH-AGENT 项目总结

## 📋 项目概述

**GIBH-AGENT** 是一个基于多模态大模型与微服务架构的单细胞生信分析智能体平台。项目正在从单体脚本重构为**分层多智能体系统**，以支持7种组学模态和TB级数据处理。

## 🏗️ 当前架构

### 技术栈
- **前端/网关**: Nginx (反向代理)
- **API 服务**: FastAPI + Gunicorn
- **任务调度**: Celery + Redis
- **推理引擎**: vLLM (Qwen3-VL-8B + Qwen3-Coder-30B-AWQ)
- **数据存储**: 本地文件系统 + ChromaDB
- **计算引擎**: Scanpy (单细胞分析)

### 核心功能
1. **多模态对话**: 支持图文多模态交互
2. **自动化工作流**: 标准单细胞分析流程 (QC -> Normalize -> PCA -> UMAP)
3. **本地化部署**: 模型权重、向量库、数据文件全部本地化
4. **出版级绘图**: 自动生成高分辨率 (300 DPI) 分析图表

## 🔄 重构目标

### 核心改进
1. **支持7种组学模态**: Transcriptomics, Genomics, Epigenomics, Metabolomics, Proteomics, Spatial Omics, Imaging
2. **处理TB级数据**: 通过控制平面/数据平面分离
3. **灵活LLM切换**: 本地（vLLM）和云端（DeepSeek/SiliconFlow）无缝切换
4. **模块化架构**: 易于扩展和维护

### 新架构设计

```
用户查询
    ↓
RouterAgent (路由智能体) - 识别组学类型和意图
    ↓
Domain Agents (领域智能体)
    ├── RNAAgent (转录组) ✅ 已实现
    ├── DNAAgent (基因组) ⏳ 占位符
    ├── EpigenomicsAgent ⏳ 占位符
    ├── MetabolomicsAgent ⏳ 占位符
    ├── ProteomicsAgent ⏳ 占位符
    ├── SpatialAgent ⏳ 占位符
    └── ImagingAgent ⏳ 占位符
    ↓
Tools (工具类) - 生成脚本，不执行
    ├── CellRangerTool ✅
    └── ScanpyTool ✅
    ↓
TaskDispatcher (任务分发器) - 提交到HPC
    ├── 本地执行
    ├── Slurm 提交
    └── SSH 远程提交
```

## 📁 项目结构

### 现有代码（services/）
```
services/
├── api/                    # FastAPI 后端
│   ├── src/
│   │   ├── agent.py        # BioBlendAgent (旧版)
│   │   ├── celery_app.py   # Celery 任务
│   │   ├── main.py         # FastAPI 入口
│   │   └── skills/         # Scanpy 逻辑
├── nginx/                  # 网关配置
└── worker/                 # 异步计算节点
```

### 新架构代码（gibh_agent/）
```
gibh_agent/
├── config/
│   ├── settings.yaml       # 统一配置文件
│   └── prompts/            # 提示词模板
│       └── router.yaml
├── core/
│   ├── llm_client.py       # 统一 LLM 客户端 ✅
│   ├── prompt_manager.py   # 提示管理器 ✅
│   └── dispatcher.py       # 任务分发器 ✅
├── agents/
│   ├── base_agent.py       # 基础智能体类 ✅
│   ├── router_agent.py     # 路由智能体 ✅
│   └── specialists/
│       ├── rna_agent.py   # 转录组智能体 ✅
│       ├── dna_agent.py    # 基因组智能体 ⏳
│       └── ...             # 其他5个智能体 ⏳
├── tools/
│   ├── cellranger_tool.py  # Cell Ranger 工具 ✅
│   └── scanpy_tool.py      # Scanpy 工具 ✅
└── main.py                 # 主入口 ✅
```

## 🎯 核心设计原则

### 1. 控制平面 vs 数据平面分离
- **控制平面**：智能体只处理文件路径（字符串），不处理二进制数据
- **数据平面**：TaskDispatcher 提交脚本到 HPC，处理实际数据

### 2. 统一接口
- 所有智能体继承 `BaseAgent`
- 统一的 LLM 调用接口（LLMClient）
- 统一的工具调用接口

### 3. 配置驱动
- LLM 配置可切换（本地/云端）
- 任务分发方式可配置（本地/Slurm/SSH）
- 提示词模板化管理

## 📊 实施状态

### ✅ 已完成（Phase 1）
- [x] LLMClient 统一客户端
- [x] PromptManager 提示管理器
- [x] TaskDispatcher 任务分发器
- [x] BaseAgent 基础类
- [x] RouterAgent 路由智能体
- [x] RNAAgent 转录组智能体（重构）
- [x] CellRangerTool 和 ScanpyTool
- [x] 配置文件结构
- [x] 其他6个领域智能体占位符

### ⏳ 待完成（Phase 2-4）
- [ ] 完善 RNAAgent 的完整工作流
- [ ] 与现有 FastAPI 服务集成
- [ ] 实现其他领域智能体
- [ ] 添加单元测试
- [ ] 性能优化

## 🔧 关键技术点

### LLM 切换
```python
# 本地模型
client = LLMClientFactory.create_local_vllm("qwen3-vl")

# 云端模型
client = LLMClientFactory.create_cloud_deepseek()
```

### 任务分发
```python
# 智能体生成脚本（只处理路径）
script = cellranger_tool.generate_count_script(
    fastq_dir="/data/fastq",  # 只传路径
    sample_id="sample1"
)

# TaskDispatcher 提交执行
task_info = await dispatcher.submit_script(script)
```

### 路由决策
```python
# RouterAgent 自动路由
route_result = await router.process_query(query, files)
# 返回: {"routing": "rna_agent", "modality": "transcriptomics"}

# 获取对应智能体
target_agent = agents[route_result["routing"]]
result = await target_agent.process_query(query, files)
```

## 📝 关键文件说明

### 配置文件
- `config/settings.yaml`: 统一配置文件，支持环境变量
- `config/prompts/router.yaml`: 路由智能体提示词模板

### 核心组件
- `core/llm_client.py`: LLM 客户端，支持本地/云端切换
- `core/prompt_manager.py`: 提示管理器，使用 Jinja2 模板
- `core/dispatcher.py`: 任务分发器，支持本地/Slurm/SSH

### 智能体
- `agents/router_agent.py`: 路由智能体，识别组学类型
- `agents/specialists/rna_agent.py`: 转录组智能体（已实现）
- `agents/specialists/*.py`: 其他6个领域智能体（占位符）

### 工具类
- `tools/cellranger_tool.py`: Cell Ranger 脚本生成器
- `tools/scanpy_tool.py`: Scanpy 工作流脚本生成器

## 🚀 使用示例

```python
from gibh_agent import create_agent

# 创建智能体
agent = create_agent("config/settings.yaml")

# 处理查询
result = await agent.process_query(
    query="帮我分析一下这个单细胞数据",
    uploaded_files=[{"name": "sample.h5ad", "path": "/data/sample.h5ad"}]
)

# 结果包含路由信息和处理结果
print(result["routing_info"])  # 路由决策
print(result)  # 智能体响应
```

## 📚 文档

- `REFACTORING_PLAN.md`: 详细重构方案
- `IMPLEMENTATION_SUMMARY.md`: 实施总结
- `IMPROVEMENT_ANALYSIS.md`: 改进分析
- `gibh_agent/README.md`: 新架构使用指南

## ⚠️ 注意事项

1. **向后兼容**: 现有 `services/api/src/agent.py` 保持不变，新架构并行运行
2. **渐进式迁移**: 逐步替换旧代码，不一次性重构
3. **配置管理**: 使用 YAML 配置文件，支持环境变量替换
4. **文件路径原则**: 智能体只处理文件路径，不读取二进制数据

## 🎯 下一步计划

1. **完善现有功能**: RNAAgent 完整工作流
2. **集成现有服务**: 与 FastAPI 服务集成
3. **扩展领域智能体**: 按优先级实现其他智能体
4. **添加测试**: 单元测试和集成测试
5. **性能优化**: 路由优化、缓存机制

