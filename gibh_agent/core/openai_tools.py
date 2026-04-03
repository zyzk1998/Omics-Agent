# -*- coding: utf-8 -*-
"""将 ToolRegistry 中的工具转为 OpenAI Chat Completions `tools` 格式。"""
from typing import Any, Dict, List, Optional

from gibh_agent.core.tool_registry import registry

HPC_MCP_TOOL_PREFIX = "hpc_mcp_"
COMPUTE_SCHEDULER_MCP_KEY = "compute_scheduler"

# 仅追加，不替换核心 system（由调用方拼在已有 chat_system 末尾）
HPC_CHAT_SYSTEM_APPEND = (
    "\n\n【超算调度授权】你现在已通过 MCP 协议连接至国家级超级计算机中心。你可以直接使用**自然语言**与超算进行交互。"
    "当用户提出重型计算需求时，请调用 `hpc_mcp_` 开头的工具，并将用户的意图转化为清晰的自然语言指令作为参数传递给超算。超算会自动理解并执行。"
    "\n\n【异步任务操作规范（SOP）】超算任务通常是异步的。当你调用「提交任务」类工具并获得 Job ID（或等价作业标识）后，"
    "**绝对不要**在同一轮对话里盲目、反复循环调用状态查询工具来刷进度。"
    "你必须立刻用中文回复用户，明确告知：「任务已提交，Job ID 为 XXX。您可以随时让我查询该任务的最新进度。」"
    "只有当用户**明确要求**查询进度、状态或结果时，你才调用对应的状态/队列查询类工具。"
)


def _enabled_mcp_set(enabled_mcps: Optional[List[str]]) -> set:
    return {str(x).strip() for x in (enabled_mcps or []) if x is not None and str(x).strip()}


def apply_hpc_mcp_tool_policy(names: List[str], enabled_mcps: Optional[List[str]] = None) -> List[str]:
    """
    仅当 enabled_mcps 含 compute_scheduler 时，向工具列表注入 registry 内全部 hpc_mcp_*；
    否则从列表中剔除所有 hpc_mcp_*（防 LLM 侧漏挂）。
    """
    allow = COMPUTE_SCHEDULER_MCP_KEY in _enabled_mcp_set(enabled_mcps)
    out = [n for n in names if not str(n).startswith(HPC_MCP_TOOL_PREFIX)]
    if allow:
        seen = set(out)
        for n in registry.list_tools():
            if n.startswith(HPC_MCP_TOOL_PREFIX) and n not in seen:
                seen.add(n)
                out.append(n)
    return out


def hpc_chat_system_suffix(enabled_mcps: Optional[List[str]]) -> str:
    if COMPUTE_SCHEDULER_MCP_KEY not in _enabled_mcp_set(enabled_mcps):
        return ""
    return HPC_CHAT_SYSTEM_APPEND


def tool_names_to_openai_tools(names: List[str]) -> List[Dict[str, Any]]:
    out: List[Dict[str, Any]] = []
    for name in names:
        meta = registry.get_metadata(name)
        if not meta:
            continue
        try:
            schema = meta.args_schema.model_json_schema()
        except Exception:
            schema = {"type": "object", "properties": {}}
        out.append(
            {
                "type": "function",
                "function": {
                    "name": meta.name,
                    "description": (meta.description or "")[:4096],
                    "parameters": schema,
                },
            }
        )
    return out
