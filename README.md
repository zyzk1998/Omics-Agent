# Omics Agent

<div align="center">

![Omics Agent Logo](https://via.placeholder.com/150x150.png?text=Omics+Agent)

**基于多模态大模型与微服务架构的多组学数据分析智能体平台**

**支持 7 种组学模态：转录组、基因组、表观遗传、代谢组、蛋白质组、空间组学、影像分析**

[![Python](https://img.shields.io/badge/Python-3.10%2B-blue?logo=python)](https://www.python.org/)
[![FastAPI](https://img.shields.io/badge/Backend-FastAPI-009688?logo=fastapi)](https://fastapi.tiangolo.com/)
[![Docker](https://img.shields.io/badge/Deploy-Docker-2496ED?logo=docker)](https://www.docker.com/)
[![License](https://img.shields.io/badge/License-Proprietary-red)](LICENSE)

[功能特性](#-功能特性) • [技术架构](#-技术架构) • [快速开始](#-快速开始) • [API文档](#-api-文档) • [项目结构](#-项目结构)

</div>

---

## 📖 项目简介

**Omics Agent** 是一款面向生物信息学研究的**企业级多组学数据分析智能体平台**，采用**多智能体系统（Multi-Agent System, MAS）架构**和**工具增强生成（Tool-Augmented Generation, TAG）范式**，通过自然语言交互实现多组学数据的端到端自动化分析。

系统基于**大语言模型（Large Language Model, LLM）驱动的动态工作流规划**，结合**语义工具检索（Semantic Tool Retrieval）**和**模块化工具执行引擎**，实现了从数据质控、特征提取、统计分析到结果可视化的全流程自动化。核心创新在于将传统的硬编码分析流程转化为**LLM 驱动的智能规划与执行系统**，显著降低了多组学数据分析的技术门槛。

### ✨ 核心亮点

- **🧬 多组学模态支持**：覆盖 7 种主流组学类型（转录组、基因组、表观遗传、代谢组、蛋白质组、空间组学、影像分析），采用**领域特定智能体（Domain-Specific Agents）**架构，实现专业化分析能力
- **🤖 多模态交互能力**：基于 **DeepSeek-V3.2** 多模态大模型，支持自然语言查询、文件上传和图表解读，实现**人机协作式数据分析**
- **⚡ 动态工作流规划**：采用 **Tool-RAG（Tool Retrieval-Augmented Generation）** 架构，通过语义检索和 LLM 规划生成定制化分析流程，无需硬编码模板
- **🔒 数据隐私与安全**：支持**混合部署模式**（本地 LLM + 云端 API），满足科研机构数据安全合规要求
- **📊 出版级可视化**：基于 **Scanpy**、**Matplotlib** 等专业库，自动生成符合学术发表标准的矢量图表（300+ DPI，支持 SVG/PDF 导出）
- **🚀 可扩展架构**：采用**插件化工具系统**和**通用执行引擎**，新增分析工具无需修改核心代码，支持快速迭代和功能扩展

---

## 🏗 技术架构

<div align="center">
  <img src="./System Architecture" width="400" alt="Omics Agent System Architecture">
  <br>
  <em>Omics Agent 微服务交互架构图</em>
</div>

<br>

系统采用**前后端分离的微服务架构**，遵循**领域驱动设计（DDD）**和**关注点分离（SoC）**原则，各组件通过 Docker Compose 进行容器化编排：

| 架构层次 | 技术选型 | 设计说明 |
| :--- | :--- | :--- |
| **表现层（Presentation Layer）** | HTML5 + Bootstrap 5 + Marked.js | 响应式 Web 界面，支持 Markdown 渲染、代码高亮和实时流式更新（SSE） |
| **应用服务层（Application Layer）** | FastAPI + Gunicorn + Uvicorn | 异步 HTTP 服务，基于 ASGI 协议实现高并发请求处理，支持流式响应（Server-Sent Events） |
| **业务逻辑层（Business Logic Layer）** | Python 3.10+ + Pydantic v2 | 领域模型封装、参数验证、工作流编排和智能体协调 |
| **计算执行层（Execution Layer）** | Celery + Redis + Python Workers | 分布式任务队列，支持异步任务调度、结果持久化和状态监控 |
| **推理服务层（Inference Layer）** | SiliconFlow API / Local LLM | 大语言模型推理服务，支持多模态输入（文本、图像、文件）和流式输出 |
| **数据持久层（Persistence Layer）** | 本地文件系统 + SQLite（可选） | 用户数据隔离、会话管理、结果存储和元数据管理 |

### 架构演进

**当前版本（V2）**：已实现**动态 Tool-RAG 架构**，支持工具自动发现、语义检索、动态规划和通用执行

```
用户查询
    ↓
WorkflowPlanner (工作流规划器) - LLM 驱动的动态规划
    ├── ToolRetriever.retrieve() → 语义检索相关工具
    ├── 注入工具 schema 到提示词
    └── 生成 JSON 工作流计划
    ↓
Domain Agents (领域智能体)
    ├── RNAAgent (转录组) ✅ 已实现
    ├── DNAAgent (基因组) ⏳ 占位符
    ├── EpigenomicsAgent (表观遗传) ⏳ 占位符
    ├── MetabolomicsAgent (代谢组) ✅ 已实现
    ├── ProteomicsAgent (蛋白质组) ⏳ 占位符
    ├── SpatialAgent (空间组学) ⏳ 占位符
    └── ImagingAgent (影像分析) ⏳ 占位符
    ↓
ToolRegistry (工具注册表) - 模块化插件系统
    ├── 自动发现和注册工具
    ├── Pydantic 参数验证
    └── 单例模式管理
    ↓
WorkflowExecutor (工作流执行器) - 通用执行引擎
    ├── 动态工具查找
    ├── 参数验证和执行
    └── 数据流处理
    ↓
Tools (工具类) - 按领域组织
    ├── general/file_inspector.py ✅
    ├── metabolomics/preprocessing.py ✅
    ├── metabolomics/statistics.py ✅
    ├── metabolomics/plotting.py ✅
    ├── CellRangerTool ✅
    └── ScanpyTool ✅
```

**核心特性**:
- 🔍 **自动工具发现（Auto Tool Discovery）**: 基于 Python 装饰器和反射机制，递归遍历工具目录，自动注册工具到全局注册表
- 🔎 **语义工具检索（Semantic Tool Retrieval）**: 采用 **ChromaDB** 向量数据库和 **Ollama Embeddings**，实现基于语义相似度的工具推荐，支持自然语言查询
- 🧠 **LLM 驱动的工作流规划（LLM-Driven Workflow Planning）**: 将工具 Schema 注入到 LLM 提示词中，通过结构化输出（JSON）生成可执行的工作流计划
- ⚙️ **通用执行引擎（Universal Execution Engine）**: 基于**策略模式**和**依赖注入**，实现工具无关的执行框架，支持动态工具查找、参数验证和结果聚合
- 📦 **模块化工具系统（Modular Tool System）**: 按领域（Domain）组织工具，遵循**单一职责原则（SRP）**，支持插件化扩展和版本管理

---

## 系统图

<div align="center">
  <img src="./Data flow" width="800" alt="Data Flow Pipeline">
  <br>
  <em>智能体数据处理与逻辑流向图</em>
</div>

## ⚡ 快速开始

### 1. 环境准备

- **操作系统**: Linux (Ubuntu 20.04+ 推荐)
- **硬件资源**: 
  - CPU: 8 cores+
  - RAM: 32GB+ (生信分析内存消耗大)
  - GPU: 可选（如需本地 LLM 推理）
- **软件依赖**: Docker, Docker Compose

### 2. 部署步骤

```bash
# 1. 克隆仓库
git clone https://github.com/zyzk1998/GIBH-AGENT-V2.git
cd GIBH-AGENT-V2
# 注意：项目文件夹名称保持为 GIBH-AGENT-V2（避免导入错误），但产品名称为 Omics Agent

# 2. 配置环境变量（可选）
# 编辑 docker-compose.yml 中的 SILICONFLOW_API_KEY 和 SILICONFLOW_MODEL

# 3. 启动服务
# 首次启动会自动构建镜像，耗时较长请耐心等待
docker compose up -d --build

# 4. 验证状态
docker compose logs -f api-server
# 等待出现 "Application startup complete" 字样
```

### 3. 访问服务

- **Web 界面**: `http://localhost:8028`
- **API 文档**: `http://localhost:8028/api/docs` (Swagger UI)

---

## 📂 项目结构

```text
GIBH-AGENT-V2/
├── docker-compose.yml              # 容器编排配置
├── server.py                       # FastAPI 服务器入口
├── requirements.txt                # Python 依赖
│
├── services/                       # 服务配置
│   ├── api/                        # API 服务 Dockerfile
│   │   └── Dockerfile
│   └── nginx/                      # 前端静态文件
│       └── html/
│           └── index.html          # 前端页面
│
├── gibh_agent/                     # 核心智能体系统
│   ├── core/                       # 核心基础设施
│   │   ├── llm_client.py          # 统一 LLM 客户端（支持本地/云端切换）
│   │   ├── prompt_manager.py      # 提示管理器（Jinja2 模板）
│   │   ├── dispatcher.py          # 任务分发器（本地/Slurm/SSH）
│   │   ├── file_inspector.py      # 文件检测和元数据生成
│   │   ├── celery_app.py          # Celery 异步任务配置
│   │   ├── tool_registry.py       # 🔥 工具注册系统（ToolRegistry）
│   │   ├── tool_retriever.py       # 🔥 工具检索系统（ChromaDB + Embeddings）
│   │   ├── planner.py             # 🔥 工作流规划器（LLM 驱动）
│   │   └── executor.py             # 🔥 工作流执行器（通用执行引擎）
│   │
│   ├── agents/                     # 智能体系统
│   │   ├── base_agent.py          # 基础智能体抽象类
│   │   ├── router_agent.py         # 路由智能体（识别组学类型）
│   │   └── specialists/            # 领域智能体
│   │       ├── rna_agent.py       # 转录组智能体 ✅
│   │       ├── dna_agent.py       # 基因组智能体 ⏳
│   │       ├── epigenomics_agent.py # 表观遗传智能体 ⏳
│   │       ├── metabolomics_agent.py # 代谢组智能体 ⏳
│   │       ├── proteomics_agent.py # 蛋白质组智能体 ⏳
│   │       ├── spatial_agent.py   # 空间组学智能体 ⏳
│   │       └── imaging_agent.py    # 影像分析智能体 ⏳
│   │
│   ├── tools/                      # 🔥 模块化工具系统（自动发现）
│   │   ├── __init__.py            # 自动发现和加载系统
│   │   ├── general/               # 通用工具
│   │   │   └── file_inspector.py  # 文件检查工具 ✅
│   │   ├── metabolomics/           # 代谢组学工具
│   │   │   ├── preprocessing.py   # 数据预处理 ✅
│   │   │   ├── statistics.py      # 统计分析（PCA, 差异分析）✅
│   │   │   └── plotting.py        # 可视化（火山图, 热图）✅
│   │   ├── cellranger_tool.py     # Cell Ranger 脚本生成器 ✅
│   │   ├── scanpy_tool.py         # Scanpy 工作流脚本生成器 ✅
│   │   └── ...
│   │
│   └── config/                     # 配置文件
│       ├── settings.yaml          # 统一配置文件
│       └── prompts/                # 提示词模板
│
├── data/                           # 数据目录
│   ├── uploads/                    # 用户上传的原始数据
│   ├── results/                    # 分析结果产出
│   └── redis/                      # Redis 持久化数据
│
└── results/                        # 分析结果（挂载点）
```

---

## 🛠 维护与排查

### 常用命令

```bash
# 服务重启
docker compose restart api-server worker

# 查看日志
docker compose logs -f api-server
docker compose logs -f worker

# 查看服务状态
docker compose ps

# 清理并重建
docker compose down -v
docker compose up -d --build

# 进入容器调试
docker compose exec api-server bash
```

### 常见问题

1. **502 Bad Gateway**: 检查 API 服务是否正常启动
```bash
   docker compose logs api-server
   ```

2. **文件上传失败**: 检查 `data/uploads` 目录权限
   ```bash
   sudo chmod -R 777 data/uploads
   ```

3. **Redis 连接失败**: 检查 Redis 服务状态
   ```bash
   docker compose ps redis
   docker compose logs redis
   ```

---

## 📋 功能特性

### 1. 文件上传与管理

- ✅ **多文件批量上传**：支持同时选择多个文件
- ✅ **文件暂存机制**：选择文件后暂存，点击发送时统一上传
- ✅ **10x Genomics 数据检测**：自动识别并分组 `matrix.mtx`、`barcodes.tsv`、`features.tsv`
- ✅ **文件类型检测**：自动识别单细胞数据、FASTQ、CSV 等格式
- ✅ **会话级文件注册表**：支持多文件上下文管理，解决上下文停滞问题

### 2. 智能对话分析

- ✅ **自然语言交互**：通过对话描述分析需求
- ✅ **动态工作流规划**：🔥 LLM 驱动的智能工作流生成（无需硬编码模板）
- ✅ **语义工具检索**：🔥 基于 ChromaDB 的工具语义检索
- ✅ **工具选择**：智能推荐合适的分析工具（Cell Ranger、Scanpy 等）
- ✅ **流式响应**：实时显示分析进度和结果

### 3. 多智能体路由

- ✅ **意图识别**：自动识别用户查询的组学类型
- ✅ **智能路由**：将查询路由到对应的领域智能体
- ✅ **策略模式提示词**：领域特定指令隔离，消除领域幻觉
- ✅ **扩展性**：易于添加新的组学模态支持

### 4. 模块化工具系统 🔥

- ✅ **自动发现**：递归遍历目录，自动发现和注册工具
- ✅ **装饰器注册**：使用 `@registry.register` 装饰器定义工具
- ✅ **参数验证**：Pydantic v2 自动生成和验证参数 schema
- ✅ **按领域组织**：工具按领域（general, metabolomics, etc.）组织
- ✅ **可扩展性**：新增工具只需创建文件，无需修改其他代码

### 5. 通用工作流执行 🔥

- ✅ **动态执行**：工具无关的执行引擎，支持任意注册工具
- ✅ **参数验证**：自动验证工具参数（Pydantic）
- ✅ **数据流处理**：处理步骤间的数据流传递
- ✅ **错误处理**：优雅降级，结构化错误消息

### 6. 异步任务处理

- ✅ **Celery 任务队列**：处理耗时的生信分析任务
- ✅ **任务状态监控**：实时查看任务执行状态
- ✅ **结果持久化**：分析结果自动保存到 `results/` 目录

---

## 🔧 API 文档

### 交互式 API 文档

启动服务后，访问以下地址查看完整的 API 文档：

- **Swagger UI**: `http://localhost:8028/api/docs` - 交互式 API 测试界面
- **ReDoc**: `http://localhost:8028/api/redoc` - 可读性更强的 API 文档

### 完整 API 文档

详细的 API 接口文档请参考：[**API.md**](API.md)

该文档包含：
- 所有 API 端点的详细说明（请求/响应格式、参数说明、错误处理）
- SSE 流式响应格式和事件类型
- 数据结构定义（TypeScript 接口）
- 前端集成指南和最佳实践
- 错误处理策略和常见问题解决方案

### 主要 API 端点概览

| 端点 | 方法 | 说明 | 响应类型 |
|------|------|------|----------|
| `/api/upload` | POST | 文件上传（支持多文件、10x Genomics 自动识别） | JSON |
| `/api/chat` | POST | 聊天接口（支持流式响应 SSE） | SSE / JSON |
| `/api/execute` | POST | 直接执行工作流 | JSON |
| `/api/tools/search` | GET | 语义搜索工具 | JSON |
| `/api/workflows/plan` | POST | 规划工作流（plan-first 模式） | JSON |
| `/api/health` | GET | 健康检查和组件状态 | JSON |

---

## 📄 版权说明

Copyright © 2025 Omics Agent Team. All Rights Reserved.
本项目为商业版代码，未经授权禁止商用分发。

---

## 🔗 相关文档

### 核心文档
- [**API.md**](API.md) - 📚 **完整 API 接口文档**（前后端交接必备，包含所有接口的输入输出、错误处理和使用示例）
- [架构重构总结](ARCHITECTURE_REFACTORING.md) - 🔥 Tool-RAG 动态工作流系统详细说明
- [项目总结](PROJECT_SUMMARY.md) - 项目概述和架构说明

### 开发文档
- [重构方案](REFACTORING_PLAN.md) - 详细的重构计划
- [快速参考](QUICK_REFERENCE.md) - 常用命令和配置
- [Docker 部署](DOCKER_DEPLOYMENT.md) - Docker 部署详细指南

### 前端集成
- [**lite.html**](services/nginx/html/lite.html) - 瘦前端演示页面，展示如何与后端 API 集成

---

<div align="center">

**Made with ❤️ by Omics Agent Team**

</div>
