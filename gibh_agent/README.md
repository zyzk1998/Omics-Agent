# GIBH-AGENT 重构版本

## 📖 概述

这是 GIBH-AGENT 的重构版本，从单体脚本升级为**分层多智能体系统**，支持7种组学模态和大数据处理。

## 🏗️ 架构设计

### 核心原则

1. **控制平面 vs 数据平面分离**
   - 智能体只处理文件路径（字符串），不处理二进制数据
   - 大数据处理通过 TaskDispatcher 提交到 HPC 集群

2. **统一接口**
   - 所有智能体继承 `BaseAgent`
   - 统一的 LLM 调用接口
   - 统一的工具调用接口

3. **可扩展性**
   - 新增组学模态只需添加新的 Specialist Agent
   - 新增工具只需实现对应的 Tool 类

### 架构层次

```
用户查询
    ↓
RouterAgent (路由智能体)
    ↓
Domain Agents (领域智能体)
    ├── RNAAgent (转录组)
    ├── DNAAgent (基因组)
    ├── EpigenomicsAgent (表观遗传)
    ├── MetabolomicsAgent (代谢组)
    ├── ProteomicsAgent (蛋白质组)
    ├── SpatialAgent (空间组学)
    └── ImagingAgent (影像分析)
    ↓
Tools (工具类)
    ├── CellRangerTool (生成脚本)
    ├── ScanpyTool (生成脚本)
    └── ...
    ↓
TaskDispatcher (任务分发器)
    ├── 本地执行
    ├── Slurm 提交
    └── SSH 远程提交
```

## 🚀 快速开始

### 1. 安装依赖

```bash
pip install openai jinja2 pyyaml paramiko
```

### 2. 配置

编辑 `config/settings.yaml`：

```yaml
llm:
  default: "local"  # 或 "cloud"
  local:
    logic:
      base_url: "http://localhost:8001/v1"
      model: "qwen3-coder-awq"
    vision:
      base_url: "http://localhost:8000/v1"
      model: "qwen3-vl"

dispatcher:
  type: "local"  # 或 "slurm", "ssh"
```

### 3. 使用示例

```python
from gibh_agent import create_agent
import asyncio

async def main():
    # 创建智能体
    agent = create_agent("config/settings.yaml")
    
    # 处理查询
    result = await agent.process_query(
        query="帮我分析一下这个单细胞数据",
        uploaded_files=[
            {"name": "sample.h5ad", "path": "/data/sample.h5ad"}
        ]
    )
    
    print("路由信息:", result.get("routing_info"))
    print("处理结果:", result)

asyncio.run(main())
```

## 📁 目录结构

```
gibh_agent/
├── config/
│   ├── settings.yaml          # 配置文件
│   └── prompts/               # 提示词模板
│       └── router.yaml
├── core/
│   ├── llm_client.py          # LLM 客户端
│   ├── prompt_manager.py      # 提示管理器
│   └── dispatcher.py          # 任务分发器
├── agents/
│   ├── base_agent.py          # 基础智能体
│   ├── router_agent.py        # 路由智能体
│   └── specialists/           # 领域智能体
│       ├── rna_agent.py
│       ├── dna_agent.py
│       └── ...
├── tools/
│   ├── cellranger_tool.py     # Cell Ranger 工具
│   └── scanpy_tool.py         # Scanpy 工具
└── main.py                    # 主入口
```

## 🔧 核心组件

### LLMClient

统一 LLM 客户端，支持本地和云端切换：

```python
from gibh_agent.core.llm_client import LLMClientFactory

# 本地模型
client = LLMClientFactory.create_local_vllm("qwen3-vl")

# 云端模型
client = LLMClientFactory.create_cloud_deepseek()
```

### PromptManager

提示管理器，使用模板系统：

```python
from gibh_agent.core.prompt_manager import PromptManager

manager = PromptManager("config/prompts")
prompt = manager.get_prompt("rna_expert", {
    "file_path": "/data/sample.h5ad",
    "user_intent": "单细胞分析"
})
```

### TaskDispatcher

任务分发器，处理大数据：

```python
from gibh_agent.core.dispatcher import TaskDispatcher

dispatcher = TaskDispatcher(config)
task_info = await dispatcher.submit_script(
    script_content=script,
    script_name="job.sh"
)
```

## 📝 扩展指南

### 添加新的领域智能体

1. 在 `agents/specialists/` 创建新文件
2. 继承 `BaseAgent`
3. 实现 `process_query` 方法
4. 在 `main.py` 中注册

### 添加新的工具

1. 在 `tools/` 创建新工具类
2. 实现脚本生成方法（只生成脚本，不执行）
3. 在对应的智能体中使用

## 🔄 迁移路径

1. **保持向后兼容**：现有 API 保持不变
2. **渐进式迁移**：逐步替换旧代码
3. **配置驱动**：通过配置文件切换功能

## 📚 文档

- [重构方案](./REFACTORING_PLAN.md)
- [实施总结](./IMPLEMENTATION_SUMMARY.md)

## ⚠️ 注意事项

1. **文件路径 vs 二进制数据**：智能体只处理文件路径，不读取文件内容
2. **脚本生成 vs 执行**：工具类只生成脚本，通过 TaskDispatcher 执行
3. **配置管理**：使用 YAML 配置文件，支持环境变量

## 🎯 下一步

- [ ] 完善各领域智能体的实现
- [ ] 添加更多工具类
- [ ] 与现有 FastAPI 服务集成
- [ ] 添加单元测试
- [ ] 性能优化

