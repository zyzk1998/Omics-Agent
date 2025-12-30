# GIBH-AGENT-V2 设置指南

## 📋 项目说明

这是 GIBH-AGENT 的重构版本，**完全独立于原项目**。

- ✅ 原项目 (`/home/ubuntu/GIBH-AGENT/`) 保持不变
- ✅ 新架构代码在本目录独立开发
- ✅ 可以单独提交到新的 Git 仓库

## 🔧 环境设置

### 1. 安装依赖

```bash
cd /home/ubuntu/GIBH-AGENT-V2
pip install openai jinja2 pyyaml paramiko
```

### 2. 配置环境变量

创建 `.env` 文件（可选）：

```bash
# LLM 配置
export VLLM_LOGIC_URL="http://localhost:8001/v1"
export LLM_LOGIC_MODEL="qwen3-coder-awq"
export VLLM_VL_URL="http://localhost:8000/v1"
export LLM_VL_MODEL="qwen3-vl"

# 云端 LLM（可选）
export DEEPSEEK_API_KEY="your_key_here"
export SILICONFLOW_API_KEY="your_key_here"

# 路径配置
export UPLOAD_DIR="/app/uploads"
export RESULTS_DIR="/app/results"
```

### 3. 编辑配置文件

编辑 `gibh_agent/config/settings.yaml`，根据你的环境调整配置。

## 🚀 快速测试

```python
from gibh_agent import create_agent
import asyncio

async def test():
    agent = create_agent("gibh_agent/config/settings.yaml")
    result = await agent.process_query(
        query="你好",
        uploaded_files=[]
    )
    print(result)

asyncio.run(test())
```

## 📁 目录说明

- `gibh_agent/`: 新架构代码
- `CONTINUE_PROJECT_PROMPT.md`: 继续讨论提示词
- `PROJECT_SUMMARY.md`: 项目总结
- `REFACTORING_PLAN.md`: 重构方案

## ⚠️ 注意事项

1. **独立开发**: 本目录与原项目完全独立，可以自由修改
2. **先理解再开发**: 建议先阅读文档，理解架构
3. **Git 管理**: 建议创建独立的 Git 仓库管理

## 🔄 Git 初始化（可选）

```bash
cd /home/ubuntu/GIBH-AGENT-V2
git init
git add .
git commit -m "Initial commit: GIBH-AGENT V2 refactored architecture"

# 推送到新仓库
git remote add origin <your-new-repo-url>
git push -u origin main
```

