from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, List

from gibh_agent.core.tool_registry import registry
from gibh_agent.core.utils import safe_tool_execution
from gibh_agent.core.workspace_context import get_workspace_context


def _resolve_under_workspace(relative_or_file: str) -> Path:
    ctx = get_workspace_context()
    workspace_root_raw = str(ctx.get("workspace_path") or "").strip()
    if not workspace_root_raw:
        raise ValueError("当前未挂载本地项目工作区")

    workspace_root = Path(workspace_root_raw).expanduser().resolve()
    candidate_raw = str(relative_or_file or "").strip()
    candidate = Path(candidate_raw) if candidate_raw else Path(".")
    if not candidate.is_absolute():
        candidate = workspace_root / candidate
    resolved = candidate.expanduser().resolve()
    try:
        resolved.relative_to(workspace_root)
    except ValueError as exc:
        raise ValueError("非法路径：仅允许访问当前工作区根目录内的文件") from exc
    return resolved


@registry.register(
    name="list_local_directory",
    description=(
        "列出当前本地项目工作区目录内容。input_dir 为空表示工作区根目录；"
        "也可传工作区内子目录相对路径。仅允许访问已挂载工作区范围内路径。"
    ),
    category="Workspace",
    output_type="json",
)
@safe_tool_execution
def list_local_directory(input_dir: str = "") -> Dict[str, Any]:
    target = _resolve_under_workspace(input_dir)
    if not target.exists():
        return {"status": "error", "message": "目录不存在"}
    if not target.is_dir():
        return {"status": "error", "message": "目标不是目录"}

    items: List[Dict[str, Any]] = []
    for child in sorted(target.iterdir(), key=lambda p: (not p.is_dir(), p.name.lower())):
        items.append(
            {
                "name": child.name,
                "path": str(child),
                "type": "dir" if child.is_dir() else "file",
            }
        )
    return {
        "status": "success",
        "message": f"已列出目录: {target}",
        "directory": str(target),
        "items": items,
    }


@registry.register(
    name="read_local_file",
    description=(
        "读取当前本地项目工作区内文本文件内容。file_path 支持工作区相对路径或绝对路径（仅限工作区内）；"
        "max_lines 控制返回行数上限，超出将截断。"
    ),
    category="Workspace",
    output_type="json",
)
@safe_tool_execution
def read_local_file(file_path: str, max_lines: int = 1000) -> Dict[str, Any]:
    target = _resolve_under_workspace(file_path)
    if not target.exists():
        return {"status": "error", "message": "文件不存在"}
    if not target.is_file():
        return {"status": "error", "message": "目标不是文件"}

    if max_lines <= 0:
        max_lines = 1000
    max_lines = min(max_lines, 5000)

    with target.open("r", encoding="utf-8", errors="replace") as f:
        lines = f.readlines()
    total_lines = len(lines)
    truncated = total_lines > max_lines
    shown_lines = lines[:max_lines]
    content = "".join(shown_lines)

    msg = f"读取成功: {target}"
    if truncated:
        msg += f"（文件共 {total_lines} 行，仅返回前 {max_lines} 行）"

    return {
        "status": "success",
        "message": msg,
        "file_path": str(target),
        "total_lines": total_lines,
        "returned_lines": len(shown_lines),
        "truncated": truncated,
        "content": content,
    }
