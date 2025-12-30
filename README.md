# GIBH-AGENT-V2 (重构版本)

## 📋 说明

这是 **GIBH-AGENT** 的重构版本，**完全独立于原项目**。

- ✅ **原项目保持不变**: `/home/ubuntu/GIBH-AGENT/` 目录完全不受影响
- ✅ **新架构独立**: 本目录包含重构后的新架构代码
- ✅ **可以独立开发**: 可以单独提交到新的 Git 仓库

## 📁 项目位置

### 原项目（保持不变）
**位置**: `/home/ubuntu/GIBH-AGENT/`
- `services/` - 现有代码（FastAPI + Celery）
- `benchmark.py` - 基准测试
- 其他原有文件

### 新架构（独立开发）
**位置**: `/home/ubuntu/GIBH-AGENT-V2/` （本目录）
- `gibh_agent/` - 新架构代码
- 文档文件

## 🏗️ 项目结构

```
GIBH-AGENT-V2/
├── README.md                      # 本文件
├── CONTINUE_PROJECT_PROMPT.md     # 继续讨论提示词 ⭐
├── PROJECT_SUMMARY.md             # 项目总结
├── REFACTORING_PLAN.md            # 重构方案
├── IMPROVEMENT_ANALYSIS.md        # 改进分析
├── SETUP.md                       # 设置指南
├── PROJECT_LOCATION.md            # 项目位置说明
├── CLEANUP_NOTICE.md             # 清理说明
├── SAVE_TO_LOCAL.md              # 保存到本地指南
│
└── gibh_agent/                    # 新架构代码
    ├── README.md                  # 使用指南
    ├── main.py                    # 主入口
    ├── core/                      # 核心基础设施
    │   ├── llm_client.py         # LLM 客户端
    │   ├── prompt_manager.py     # 提示管理器
    │   └── dispatcher.py         # 任务分发器
    ├── agents/                    # 智能体系统
    │   ├── base_agent.py         # 基础智能体
    │   ├── router_agent.py       # 路由智能体
    │   └── specialists/          # 领域智能体
    │       ├── rna_agent.py     # 转录组智能体 ✅
    │       └── ...               # 其他6个智能体 ⏳
    ├── tools/                     # 工具类
    │   ├── cellranger_tool.py    # Cell Ranger 工具 ✅
    │   └── scanpy_tool.py        # Scanpy 工具 ✅
    └── config/                    # 配置文件
        ├── settings.yaml         # 统一配置文件
        └── prompts/              # 提示词模板
```

## 🚀 快速开始

### 在新对话框中继续讨论

1. 打开 `CONTINUE_PROJECT_PROMPT.md`
2. 复制全部内容
3. 粘贴到新对话框
4. 开始讨论

### 查看项目

- **项目概述**: 阅读 `PROJECT_SUMMARY.md`
- **重构方案**: 阅读 `REFACTORING_PLAN.md`
- **代码使用**: 参考 `gibh_agent/README.md`

## 📝 与原项目的关系

- **原项目**: `/home/ubuntu/GIBH-AGENT/` （保持不变）
- **新架构**: `/home/ubuntu/GIBH-AGENT-V2/` （本目录，独立开发）

两个项目可以并行存在，互不影响。

## 🔄 Git 仓库建议

### 选项1: 新建独立仓库（推荐）
```bash
cd /home/ubuntu/GIBH-AGENT-V2
git init
git add .
git commit -m "Initial commit: GIBH-AGENT V2 refactored architecture"
# 然后推送到新的远程仓库
```

### 选项2: 作为原项目的分支
```bash
cd /home/ubuntu/GIBH-AGENT
git checkout -b refactored-architecture
# 然后将 GIBH-AGENT-V2 的内容合并进来
```

## ⚠️ 重要说明

### 关于原项目目录下的 gibh_agent

如果原项目目录 `/home/ubuntu/GIBH-AGENT/` 下有 `gibh_agent/` 目录，那是之前创建的。你可以：

1. **删除它**（推荐，保持原项目完全不变）：
   ```bash
   rm -rf /home/ubuntu/GIBH-AGENT/gibh_agent/
   ```

2. **保留它**（如果你想在原项目中也保留一份）：
   - 两个目录互不影响
   - 新项目代码在 `/home/ubuntu/GIBH-AGENT-V2/gibh_agent/`

### 开发建议

- ✅ 本目录是**独立的新项目**，不会影响原项目
- ✅ 可以自由修改和开发
- ✅ 建议先阅读文档，理解架构后再开始开发
- ✅ 原项目代码在 `/home/ubuntu/GIBH-AGENT/services/` 保持不变

## 📦 保存到本地

查看 `SAVE_TO_LOCAL.md` 了解如何保存到本地。

## 🎯 下一步

1. 阅读 `CONTINUE_PROJECT_PROMPT.md` 了解完整项目信息
2. 阅读 `SETUP.md` 了解如何设置环境
3. 开始开发新架构功能
