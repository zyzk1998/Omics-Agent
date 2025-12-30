# 快速参考

## 📍 项目位置

- **原项目**: `/home/ubuntu/GIBH-AGENT/` （保持不变）
- **新架构**: `/home/ubuntu/GIBH-AGENT-V2/` （独立开发）

## 📚 关键文档

1. **CONTINUE_PROJECT_PROMPT.md** ⭐ - 继续讨论提示词（最重要）
2. **PROJECT_SUMMARY.md** - 项目总结
3. **REFACTORING_PLAN.md** - 重构方案
4. **SETUP.md** - 设置指南

## 🔧 核心代码

- `gibh_agent/core/llm_client.py` - LLM 客户端
- `gibh_agent/core/prompt_manager.py` - 提示管理器
- `gibh_agent/core/dispatcher.py` - 任务分发器
- `gibh_agent/agents/router_agent.py` - 路由智能体
- `gibh_agent/agents/specialists/rna_agent.py` - 转录组智能体

## 🚀 快速使用

```python
from gibh_agent import create_agent

agent = create_agent("gibh_agent/config/settings.yaml")
result = await agent.process_query("你好")
```

## ⚠️ 注意事项

- ✅ 原项目完全不受影响
- ✅ 新架构独立开发
- ✅ 可以自由修改
