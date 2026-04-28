from __future__ import annotations

import threading
from pathlib import Path
from typing import Any, Dict, Optional

_LOCK = threading.Lock()
_WORKSPACE_CONTEXT: Dict[str, Any] = {}


def set_workspace_context(workspace_path: str) -> Dict[str, Any]:
    resolved = Path(workspace_path).expanduser().resolve()
    result_dir = (resolved / "result").resolve()
    payload = {
        "workspace_path": str(resolved),
        "workspace_name": resolved.name,
        "result_path": str(result_dir),
    }
    with _LOCK:
        _WORKSPACE_CONTEXT.clear()
        _WORKSPACE_CONTEXT.update(payload)
    return dict(payload)


def get_workspace_context() -> Dict[str, Any]:
    with _LOCK:
        return dict(_WORKSPACE_CONTEXT)


def build_workspace_system_instruction(context: Optional[Dict[str, Any]] = None) -> str:
    ctx = context or get_workspace_context()
    workspace_path = str(ctx.get("workspace_path") or "").strip()
    result_path = str(ctx.get("result_path") or "").strip()
    if not workspace_path or not result_path:
        return ""
    return (
        "\n\n[系统强制指令]: 用户当前已挂载本地项目。项目根目录绝对路径为: "
        f"{workspace_path}。你可以使用提供的本地文件工具查看该目录。"
        "你的所有数据分析产出、生成的图表、保存的中间文件，"
        f"**绝对必须**写入到 {result_path} 目录下。绝不允许写在其他位置。"
    )
