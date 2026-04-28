# -*- coding: utf-8 -*-
"""HPC MCP 配置 API：重连与状态查询。"""
from __future__ import annotations

import json
import logging
import shlex
from typing import Any, Dict, List, Optional

from fastapi import APIRouter
from pydantic import BaseModel, Field

from gibh_agent.core.hpc_scanner import HpcFileScanner
from gibh_agent.mcp_client import hpc_mcp_manager

logger = logging.getLogger(__name__)

router = APIRouter(tags=["config"])
_hpc_scanner = HpcFileScanner()


class MCPConfigBody(BaseModel):
    url: str = Field(default="http://127.0.0.1:8001/mcp", description="MCP Server URL")
    transport: str = Field(default="streamableHttp", description="streamableHttp | sse")


class HPCCommandBody(BaseModel):
    command: str = Field(..., description="需要在 HPC 侧执行的 shell 命令")


class HPCPathBody(BaseModel):
    file_path: str = Field(..., description="需要校验是否存在的绝对路径")


def _flatten_content_to_text(content: Any) -> str:
    if isinstance(content, str):
        return content
    if isinstance(content, list):
        parts: List[str] = []
        for item in content:
            if isinstance(item, str):
                parts.append(item)
            elif isinstance(item, dict):
                txt = item.get("text")
                if isinstance(txt, str):
                    parts.append(txt)
                else:
                    try:
                        parts.append(json.dumps(item, ensure_ascii=False))
                    except Exception:  # noqa: BLE001
                        parts.append(str(item))
            else:
                parts.append(str(item))
        return "\n".join([p for p in parts if p])
    return str(content or "")


async def _run_hpc_shell_command(command: str) -> Dict[str, Any]:
    status = await hpc_mcp_manager.status_dict_async()
    tool_names = status.get("tools") or []
    if not isinstance(tool_names, list):
        tool_names = []
    tool_names = [str(x).strip() for x in tool_names if str(x).strip()]

    preferred = [t for t in tool_names if any(k in t.lower() for k in ("shell", "bash", "command", "exec", "terminal", "run"))]
    candidates = preferred or tool_names
    if not candidates:
        return {"status": "error", "message": "HPC MCP 未连接或无可用命令工具。", "stdout": ""}

    arg_variants = [
        {"command": command},
        {"cmd": command},
        {"script": command},
        {"input": command},
        {"query": command},
        {"bash_command": command},
        {"shell_command": command},
    ]
    last_error: Optional[str] = None
    for tool_name in candidates:
        for args in arg_variants:
            result = await hpc_mcp_manager.call_tool(tool_name, args)
            if result.get("status") == "success":
                text = _flatten_content_to_text(result.get("content") or [])
                return {"status": "success", "tool_name": tool_name, "stdout": text}
            msg = str(result.get("message") or "").strip()
            if msg:
                last_error = f"{tool_name}: {msg}"
    return {"status": "error", "message": last_error or "MCP 命令执行失败。", "stdout": ""}


@router.get("/mcp/status")
async def get_mcp_status():
    return await hpc_mcp_manager.status_dict_async()


@router.post("/mcp")
async def post_mcp_config(body: MCPConfigBody):
    try:
        await hpc_mcp_manager.connect(body.url.strip(), body.transport.strip())
        return {
            "status": "success",
            "message": "已调度后台断开并重连（不阻塞主线程）",
            "config": {"url": body.url.strip().rstrip("/"), "transport": body.transport},
        }
    except Exception as e:  # noqa: BLE001
        logger.warning("POST /api/config/mcp 失败（仍返回 200 侧信息）: %s", e, exc_info=True)
        return {
            "status": "error",
            "message": str(e),
        }


@router.post("/mcp/hpc-files/sync")
async def sync_hpc_files(body: HPCCommandBody):
    cmd = (body.command or "").strip() or "ls -1"
    result = await _run_hpc_shell_command(cmd)
    stdout = result.get("stdout") or ""
    raw_files = [str(x).strip() for x in str(stdout).splitlines() if str(x).strip()]
    cleaned = await _hpc_scanner.filter_files(raw_files)
    return {
        "status": result.get("status"),
        "tool_name": result.get("tool_name"),
        "message": result.get("message"),
        "stdout": stdout,
        "files": cleaned,
    }


@router.post("/mcp/hpc-files/check")
async def check_hpc_file_exists(body: HPCPathBody):
    file_path = (body.file_path or "").strip()
    if not file_path:
        return {"status": "error", "exists": False, "message": "file_path 不能为空。"}
    quoted = shlex.quote(file_path)
    result = await _run_hpc_shell_command(f"test -e {quoted} && echo __EXISTS__ || echo __MISSING__")
    stdout = (result.get("stdout") or "").strip()
    exists = "__EXISTS__" in stdout and "__MISSING__" not in stdout
    return {
        "status": result.get("status"),
        "exists": exists,
        "stdout": stdout,
        "message": result.get("message"),
    }
