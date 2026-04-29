# -*- coding: utf-8 -*-
"""HPC MCP 配置 API：重连与状态查询。"""
from __future__ import annotations

import json
import logging
import shlex
from typing import Any, Dict, List, Optional

from fastapi import APIRouter
from pydantic import BaseModel, Field

from gibh_agent.mcp_client import hpc_mcp_manager

logger = logging.getLogger(__name__)

router = APIRouter(tags=["config"])
# _hpc_scanner = HpcFileScanner()
# 业务新要求：HPC 资产必须“原汁原味”全量展示，不做任何智能过滤。


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


def _extract_mcp_output_text(raw_text: Any) -> str:
    """
    MCP 可能返回 JSON wrapper（如 {"success": true, "output": "..."}）。
    这里强制先提取 output/stderr/stdout 等字段，再做 splitlines，避免把 wrapper 当文件名解析。
    """
    txt = str(raw_text or "").strip()
    if not txt:
        return ""
    try:
        obj = json.loads(txt)
        if isinstance(obj, dict):
            for key in ("output", "stdout", "stderr", "message"):
                v = obj.get(key)
                if isinstance(v, str) and v.strip():
                    return v
            return txt
    except json.JSONDecodeError:
        # MCP 直接返回纯文本 stdout 时，按原文透传，绝不抛异常
        return txt
    except Exception:
        return txt
    return txt


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


def _normalize_hpc_rel_path(raw: str) -> str:
    s = str(raw or "").strip()
    if not s:
        return ""
    if s.startswith("./"):
        s = s[2:]
    if s in (".", "./"):
        return ""
    s = s.replace("\\", "/")
    while "//" in s:
        s = s.replace("//", "/")
    return s.strip("/")


def _build_hpc_tree_from_find_output(stdout: str) -> List[Dict[str, Any]]:
    """
    解析 find 输出（带类型前缀：d/f/l）为树结构。
    行格式：<type>\t<path>，例如：d\t./dir, f\t./dir/a.txt
    """
    rows: List[Dict[str, str]] = []
    for line in str(stdout or "").splitlines():
        ln = line.strip()
        if not ln:
            continue
        if "\t" in ln:
            t, p = ln.split("\t", 1)
        else:
            t, p = "f", ln
        t = str(t or "f").strip().lower()[:1]
        p = _normalize_hpc_rel_path(p)
        if not p:
            continue
        rows.append({"type": t, "path": p})

    root: Dict[str, Any] = {"children": {}}

    def ensure_dir_node(parts: List[str]) -> Dict[str, Any]:
        cur = root
        built: List[str] = []
        for part in parts:
            built.append(part)
            children = cur.setdefault("children", {})
            nxt = children.get(part)
            if not nxt:
                rel = "/".join(built)
                nxt = {
                    "name": part,
                    "type": "directory",
                    "path": f"/{rel}",
                    "children": {},
                }
                children[part] = nxt
            cur = nxt
        return cur

    for row in rows:
        rel_path = row["path"]
        parts = [p for p in rel_path.split("/") if p]
        if not parts:
            continue
        is_dir = row["type"] == "d"
        if is_dir:
            ensure_dir_node(parts)
            continue
        parent = ensure_dir_node(parts[:-1]) if len(parts) > 1 else root
        children = parent.setdefault("children", {})
        name = parts[-1]
        abs_path = "/" + rel_path
        children[name] = {
            "name": name,
            "type": "file",
            "path": abs_path,
        }

    def to_list(node: Dict[str, Any]) -> List[Dict[str, Any]]:
        out: List[Dict[str, Any]] = []
        children = node.get("children", {})
        keys = sorted(children.keys(), key=lambda x: x.lower())
        for k in keys:
            child = children[k]
            if child.get("type") == "directory":
                out.append(
                    {
                        "name": child.get("name"),
                        "type": "directory",
                        "path": child.get("path"),
                        "children": to_list(child),
                    }
                )
            else:
                out.append(
                    {
                        "name": child.get("name"),
                        "type": "file",
                        "path": child.get("path"),
                    }
                )
        return out

    return to_list(root)


def _normalize_hpc_base_path(path_value: Optional[str]) -> str:
    p = str(path_value or "").strip()
    if not p:
        return "."
    return p


def _join_hpc_path(base_path: str, name: str) -> str:
    b = str(base_path or ".").strip()
    n = str(name or "").strip()
    if not n:
        return b or "."
    if b in (".", ""):
        return n
    if b == "/":
        return "/" + n.lstrip("/")
    return b.rstrip("/") + "/" + n.lstrip("/")


def _build_hpc_items_from_ls(stdout: str, base_path: str) -> List[Dict[str, str]]:
    rows: List[Dict[str, str]] = []
    text = _extract_mcp_output_text(stdout)
    for line in str(text or "").splitlines():
        s = line.strip()
        if not s:
            continue
        is_dir = s.endswith("/")
        name = s[:-1] if is_dir else s
        if not name:
            continue
        rows.append(
            {
                "name": name,
                "type": "directory" if is_dir else "file",
                "path": _join_hpc_path(base_path, name),
            }
        )
    rows.sort(key=lambda x: (x.get("type") != "directory", x.get("name", "").lower()))
    return rows


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
    # cleaned = await _hpc_scanner.filter_files(raw_files)
    # 业务要求：禁用 HpcFileScanner 过滤，直接返回全量原始文件列表。
    cleaned = raw_files
    return {
        "status": result.get("status"),
        "tool_name": result.get("tool_name"),
        "message": result.get("message"),
        "stdout": stdout,
        "files": cleaned,
    }


@router.get("/mcp/hpc-files/tree")
async def get_hpc_files_tree(path: Optional[str] = "."):
    # 兼容旧路由名，行为改为单层 list（FileZilla 风格逐级下钻）
    base_path = _normalize_hpc_base_path(path)
    cmd = f"ls -1pa {shlex.quote(base_path)}"
    logger.info("[HPC MCP] 准备执行命令: %s", cmd)
    result = await _run_hpc_shell_command(cmd)
    stdout = result.get("stdout") or ""
    logger.debug("[HPC MCP] 原始返回结果: %s", stdout)
    if result.get("status") != "success":
        return {
            "status": "error",
            "message": result.get("message") or "HPC 文件列表拉取失败",
            "path": base_path,
            "items": [],
            "raw": stdout,
            "command": cmd,
        }
    items = _build_hpc_items_from_ls(stdout, base_path)
    return {
        "status": "success",
        "path": base_path,
        "items": items,
        "raw": stdout,
        "command": cmd,
    }


@router.get("/mcp/hpc-files/list")
async def list_hpc_files(path: Optional[str] = "."):
    # 新路由：语义明确的单层目录列表接口
    return await get_hpc_files_tree(path=path)


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
