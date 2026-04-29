import logging
import os
import platform
from pathlib import Path
from datetime import datetime
from typing import Any, Dict, List

import uvicorn
from fastapi import FastAPI, HTTPException
from fastapi.middleware.cors import CORSMiddleware
from pydantic import BaseModel


logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
)
logger = logging.getLogger("omics-local-sidecar")

app = FastAPI(title="Omics Agent Local Sidecar", version="1.0.0")
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)


class WorkspaceInitRequest(BaseModel):
    workspace_path: str


class ReadFileRequest(BaseModel):
    file_path: str
    max_lines: int = 1000


class WriteFileRequest(BaseModel):
    file_path: str
    content: str


_workspace_context: Dict[str, Any] = {}
_SKIP_DIRS = {".git", ".venv", "node_modules", "__pycache__", ".mypy_cache", ".pytest_cache"}
_MAX_DEPTH = 6
_MAX_ENTRIES_PER_DIR = 500


def _workspace_name_from_path(path_obj: Path) -> str:
    name = path_obj.name.strip()
    if name:
        return name
    return str(path_obj)


def _set_workspace_context(workspace_path: str) -> Dict[str, Any]:
    root = Path(workspace_path).expanduser().resolve()
    result_dir = root / "result"
    info = {
        "workspace_name": _workspace_name_from_path(root),
        "workspace_path": str(root),
        "result_path": str(result_dir),
        "updated_at": datetime.now().isoformat(),
        "source": "local-sidecar",
    }
    _workspace_context.clear()
    _workspace_context.update(info)
    return dict(_workspace_context)


def _safe_sort_key(path_obj: Path) -> tuple:
    return (0 if path_obj.is_dir() else 1, path_obj.name.lower())


def _build_tree(root: Path, depth: int = 0) -> List[Dict[str, Any]]:
    if depth > _MAX_DEPTH:
        return []

    nodes: List[Dict[str, Any]] = []
    try:
        children = sorted(list(root.iterdir()), key=_safe_sort_key)[:_MAX_ENTRIES_PER_DIR]
    except Exception as exc:
        logger.warning("[Workspace Tree] 读取目录失败: %s (%s)", root, exc)
        return []

    for child in children:
        name = child.name
        if name in _SKIP_DIRS:
            continue
        if child.is_dir():
            nodes.append(
                {
                    "type": "directory",
                    "name": name,
                    "path": str(child),
                    "children": _build_tree(child, depth + 1),
                }
            )
        else:
            nodes.append(
                {
                    "type": "file",
                    "name": name,
                    "path": str(child),
                }
            )
    return nodes


@app.get("/health")
async def health() -> Dict[str, Any]:
    return {"status": "ok", "service": "omics-local-sidecar", "platform": platform.system()}


@app.post("/api/workspace/init")
async def init_workspace(payload: WorkspaceInitRequest) -> Dict[str, Any]:
    raw_path = str(payload.workspace_path or "").strip()
    logger.info("[Workspace Init] 收到原始路径: %s", payload.workspace_path)
    logger.info(
        "[Workspace Init] 当前后端操作系统: %s | 后端工作目录: %s",
        platform.system(),
        os.getcwd(),
    )
    if not raw_path:
        raise HTTPException(status_code=400, detail="workspace_path 不能为空")

    resolved = Path(raw_path).expanduser().resolve()
    if not resolved.is_absolute():
        raise HTTPException(status_code=400, detail="workspace_path 必须是绝对路径")
    if not resolved.exists() or not resolved.is_dir():
        raise HTTPException(
            status_code=400,
            detail=(
                "初始化工作区失败: 路径不存在或不是目录。\n"
                f"(后端尝试访问的路径为: {resolved}，当前后端系统为: {platform.system()})"
            ),
        )

    result_dir = resolved / "result"
    result_dir.mkdir(parents=True, exist_ok=True)

    workspace_info = _set_workspace_context(str(resolved))
    logger.info("✅ [Workspace] initialized: %s", workspace_info["workspace_path"])
    return {"status": "success", "workspace": workspace_info}


@app.get("/api/workspace/tree")
async def workspace_tree() -> Dict[str, Any]:
    root_raw = str(_workspace_context.get("workspace_path") or "").strip()
    if not root_raw:
        raise HTTPException(status_code=400, detail="未挂载本地工作区")

    root = Path(root_raw).expanduser().resolve()
    if not root.exists() or not root.is_dir():
        raise HTTPException(status_code=400, detail="工作区路径不存在或不是目录")

    tree = _build_tree(root)
    return {"status": "success", "workspace_path": str(root), "tree": tree}


@app.get("/api/workspace/list")
async def workspace_list() -> Dict[str, Any]:
    return await workspace_tree()


@app.post("/api/tools/read_file")
async def read_file_tool(payload: ReadFileRequest) -> Dict[str, Any]:
    raw_path = str(payload.file_path or "").strip()
    if not raw_path:
        raise HTTPException(status_code=400, detail="file_path 不能为空")
    target = Path(raw_path).expanduser().resolve()
    if not target.exists() or not target.is_file():
        raise HTTPException(status_code=400, detail=f"文件不存在或不可读: {target}")

    max_lines = int(payload.max_lines or 1000)
    if max_lines <= 0:
        max_lines = 1000
    max_lines = min(max_lines, 5000)
    with target.open("r", encoding="utf-8", errors="replace") as f:
        lines = f.readlines()
    total = len(lines)
    shown = lines[:max_lines]
    return {
        "status": "success",
        "tool_name": "read_local_file",
        "file_path": str(target),
        "total_lines": total,
        "returned_lines": len(shown),
        "truncated": total > max_lines,
        "content": "".join(shown),
    }


@app.post("/api/tools/write_file")
async def write_file_tool(payload: WriteFileRequest) -> Dict[str, Any]:
    raw_path = str(payload.file_path or "").strip()
    if not raw_path:
        raise HTTPException(status_code=400, detail="file_path 不能为空")
    target = Path(raw_path).expanduser().resolve()
    target.parent.mkdir(parents=True, exist_ok=True)
    text = str(payload.content or "")
    target.write_text(text, encoding="utf-8")
    return {
        "status": "success",
        "tool_name": "write_local_file",
        "file_path": str(target),
        "bytes_written": len(text.encode("utf-8")),
        "message": f"已写入本地文件: {target}",
    }


if __name__ == "__main__":
    uvicorn.run(app, host="127.0.0.1", port=8019, log_level="info")
