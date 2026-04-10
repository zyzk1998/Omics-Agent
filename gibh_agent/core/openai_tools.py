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

# 与 docs/hpc_mcp_tools_catalog.json 及网关 tools/list 对齐；帮助模型直接命中 function.name，减少「先想再搜」
HPC_MCP_TOOL_ROUTING_CHEAT_SHEET = (
    "\n\n【hpc_mcp_* 工具速查 — 直接选用，勿用 execute 代替已有专用工具】"
    "以下为当前 OpenAI `tools` 中**实际出现的 function.name**（均以 `hpc_mcp_` 开头）。"
    "用户话里出现同义说法或英文原名时，应**立即**调用对应工具，无需先在对话中复述整份工具表。\n"
    "**HPC（超算）** "
    "探活/连通性 → `hpc_mcp_test_hpc_connection`；"
    "系统与集群概况、超算状态、资源信息（无具体作业 ID）→ **优先** `hpc_mcp_get_hpc_system_info`；"
    "作业列表、队列、当前跑哪些任务 → `hpc_mcp_list_jobs`；"
    "单个作业状态 → `hpc_mcp_get_job_status`（需作业 ID）；"
    "作业输出/日志 → `hpc_mcp_get_job_output`（需作业 ID）；"
    "取消作业 → `hpc_mcp_cancel_job`（需作业 ID）；"
    "提交批处理脚本 → `hpc_mcp_submit_hpc_job`；"
    "测序 FastQC 类检测提交 → `hpc_mcp_submit_data_conversion_job`；"
    "上传/下载/列目录/读远程文件 → `hpc_mcp_upload_file_to_hpc`、`hpc_mcp_download_file_from_hpc`、"
    "`hpc_mcp_list_hpc_directory`、`hpc_mcp_read_hpc_file`；"
    "sinfo、squeue、pestat、节点占用、自定义 Shell（**仅当**上述专用工具无法覆盖时）→ `hpc_mcp_execute_hpc_command`（参数为完整 Linux 命令字符串）。\n"
    "**工作站** 名称规律为 `hpc_mcp_*_workstation*` 或 `hpc_mcp_run_*_on_workstation`："
    "探活 `hpc_mcp_test_workstation_connection`；系统信息 `hpc_mcp_get_workstation_system_info`；"
    "传文件与目录 `hpc_mcp_upload_file_to_workstation`、`hpc_mcp_download_file_from_workstation`、"
    "`hpc_mcp_list_workstation_directory`、`hpc_mcp_read_workstation_file`；"
    "远程命令 `hpc_mcp_execute_workstation_command`；FastQC `hpc_mcp_run_fastqc_on_workstation`；"
    "通用容器作业 `hpc_mcp_run_container_job_on_workstation`。\n"
    "**名称映射**：用户若说远端短名（如 `get_hpc_system_info`、`list_jobs`），实际调用须加前缀：`hpc_mcp_get_hpc_system_info`、`hpc_mcp_list_jobs`。\n"
    "**禁止**：在已有 `list_jobs`/`get_job_status` 等专用工具可满足时，用 `execute_hpc_command` 代替；无作业 ID 时不要编造 ID。"
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
    return HPC_CHAT_SYSTEM_APPEND + HPC_MCP_TOOL_ROUTING_CHEAT_SHEET


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
