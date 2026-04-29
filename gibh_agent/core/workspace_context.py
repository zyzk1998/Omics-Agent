from __future__ import annotations

import threading
from pathlib import Path
from typing import Any, Dict, List, Optional

# 文件树 API：不遍历这些目录名（与本地 Explorer 常见忽略一致）
WORKSPACE_TREE_IGNORE_DIRS = frozenset({".git", "__pycache__", "node_modules", ".venv", "venv"})

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
    if not workspace_path:
        return ""
    return (
        "\n\n用户当前已挂载本地项目，根目录绝对路径为: "
        f"{workspace_path}。请根据用户的具体要求决定文件的保存路径。"
        f"如果用户没有指定，请默认保存在 {workspace_path}/result/ 目录下。"
        "你生成或修改的任何文件，必须在回复中明确告知用户其绝对路径。"
    )


def build_workspace_file_tree(root: Path) -> List[Dict[str, Any]]:
    """
    递归构建工作区根目录下的文件树 JSON。
    忽略 WORKSPACE_TREE_IGNORE_DIRS 中的目录名；忽略以 ``.`` 开头的文件。
    """
    if not root.exists() or not root.is_dir():
        return []

    nodes: List[Dict[str, Any]] = []
    try:
        entries = sorted(root.iterdir(), key=lambda p: (not p.is_dir(), p.name.lower()))
    except (OSError, PermissionError):
        return []

    for entry in entries:
        name = entry.name
        try:
            if entry.is_dir():
                if name in WORKSPACE_TREE_IGNORE_DIRS:
                    continue
                nodes.append(
                    {
                        "name": name,
                        "type": "directory",
                        "children": build_workspace_file_tree(entry),
                    }
                )
            else:
                if name.startswith("."):
                    continue
                nodes.append(
                    {
                        "name": name,
                        "type": "file",
                        "path": str(entry.resolve()),
                    }
                )
        except (OSError, PermissionError):
            continue
    return nodes
