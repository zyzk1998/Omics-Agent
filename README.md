# GIBH-AGENT-V2 (重构版本)

<div align="center">

![GIBH-AGENT Logo](https://via.placeholder.com/150x150.png?text=GIBH-AGENT)

**基于多模态大模型与微服务架构的生物信息学分析智能体平台**

[![Python](https://img.shields.io/badge/Python-3.10%2B-blue?logo=python)](https://www.python.org/)
[![FastAPI](https://img.shields.io/badge/Backend-FastAPI-009688?logo=fastapi)](https://fastapi.tiangolo.com/)
[![Docker](https://img.shields.io/badge/Deploy-Docker-2496ED?logo=docker)](https://www.docker.com/)
[![License](https://img.shields.io/badge/License-Proprietary-red)](LICENSE)

[功能特性](#-功能特性) • [技术架构](#-技术架构) • [快速开始](#-快速开始) • [API文档](#-api-文档) • [项目结构](#-项目结构)

</div>

---

## 📖 项目简介

**GIBH-AGENT-V2** 是一款企业级生物信息学分析平台，旨在通过自然语言交互（Chat）实现多组学数据的全流程自动化分析。

系统采用 **DeepSeek-V3.2** 多模态大模型作为核心大脑，结合 **Scanpy**、**Cell Ranger** 等强大的计算引擎，让科研人员可以通过对话完成从数据质控（QC）、降维聚类到细胞注释的复杂分析任务。

### ✨ 核心亮点

- **🤖 多模态交互**：支持图文对话，不仅能听懂"帮我分析这个数据"，还能识别并解读生信图表。
- **⚡ 自动化工作流**：内置标准单细胞分析 Pipeline (QC -> Normalize -> PCA -> Neighbors -> UMAP -> Clustering)。
- **🔒 数据隐私安全**：支持本地化部署（Local LLM）和云端 API（SiliconFlow）灵活切换，保障科研数据安全。
- **📊 出版级绘图**：自动生成符合 SCI 发表标准的矢量图表（300 DPI+）。
- **🚀 多智能体架构**：从单体脚本重构为分层多智能体系统，支持7种组学模态扩展。

---

## 🏗 技术架构

系统采用前后端分离的微服务架构，各组件通过 Docker Compose 编排：

| 组件 | 技术选型 | 说明 |
| :--- | :--- | :--- |
| **前端层** | HTML + Bootstrap + Marked.js | 响应式 Web 界面，支持 Markdown 渲染和代码高亮 |
| **应用层** | FastAPI + Gunicorn | 高并发异步 API 服务，处理业务逻辑 |
| **计算层** | Celery + Redis | 分布式任务队列，处理耗时的生信分析任务 |
| **推理层** | SiliconFlow API (DeepSeek-V3.2) | 云端大模型推理服务，支持多模态对话 |
| **存储层** | 本地文件系统 | 用户上传数据、分析结果存储 |

### 架构演进

**当前版本（V2）**：正在从单体架构重构为**分层多智能体系统**

```
用户查询
    ↓
RouterAgent (路由智能体) - 识别组学类型和意图
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
Tools (工具类) - 生成分析脚本
    ├── CellRangerTool ✅ 已实现
    └── ScanpyTool ✅ 已实现
    ↓
TaskDispatcher (任务分发器) - 提交到 HPC 集群
    ├── 本地执行
    ├── Slurm 提交
    └── SSH 远程提交
```

---

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
│   │   └── celery_app.py          # Celery 异步任务配置
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
│   ├── tools/                      # 工具类
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

### 2. 智能对话分析

- ✅ **自然语言交互**：通过对话描述分析需求
- ✅ **工作流配置**：自动生成标准分析流程配置
- ✅ **工具选择**：智能推荐合适的分析工具（Cell Ranger、Scanpy 等）
- ✅ **流式响应**：实时显示分析进度和结果

### 3. 多智能体路由

- ✅ **意图识别**：自动识别用户查询的组学类型
- ✅ **智能路由**：将查询路由到对应的领域智能体
- ✅ **扩展性**：易于添加新的组学模态支持

### 4. 异步任务处理

- ✅ **Celery 任务队列**：处理耗时的生信分析任务
- ✅ **任务状态监控**：实时查看任务执行状态
- ✅ **结果持久化**：分析结果自动保存到 `results/` 目录

---

## 🔧 API 文档

启动服务后，访问以下地址查看完整的 API 文档：

- **Swagger UI**: `http://localhost:8028/api/docs`
- **ReDoc**: `http://localhost:8028/api/redoc`

### 主要 API 端点

- `POST /api/chat` - 发送聊天消息，获取分析结果
- `POST /api/upload` - 上传文件（支持多文件）
- `GET /api/health` - 健康检查

---

## 📄 版权说明

Copyright © 2025 GIBH-AGENT Team. All Rights Reserved.
本项目为商业版代码，未经授权禁止商用分发。

---

## 🔗 相关文档

- [项目总结](PROJECT_SUMMARY.md) - 项目概述和架构说明
- [重构方案](REFACTORING_PLAN.md) - 详细的重构计划
- [快速参考](QUICK_REFERENCE.md) - 常用命令和配置
- [Docker 部署](DOCKER_DEPLOYMENT.md) - Docker 部署详细指南

---

<div align="center">

**Made with ❤️ by GIBH-AGENT Team**

</div>
