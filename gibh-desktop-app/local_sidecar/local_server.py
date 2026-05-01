import asyncio
import logging
import os
import platform
from pathlib import Path
from datetime import datetime
from typing import Any, Dict, List, Optional

import httpx
import uvicorn
from fastapi import FastAPI, HTTPException, Query
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


class UploadToCloudRequest(BaseModel):
    """由 Sidecar 从用户磁盘读取并 multipart 转发至云端 /api/upload。"""

    local_file_path: str
    upload_base_url: str
    authorization: Optional[str] = None
    x_guest_uuid: Optional[str] = None


class CheckFileRequest(BaseModel):
    """POST /api/tools/check_file 使用 JSON 传递路径，避免 GET query 对 Windows 反斜杠等字符的编码坑。"""

    path: str


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


def _check_file_result(raw_path: str) -> Dict[str, Any]:
    """供服务端 Docker / 前端称重：本地文件/目录是否存在（勿与 UPLOAD_DIR 拼接）。"""
    path_in = str(raw_path or "").strip()
    if not path_in:
        raise HTTPException(status_code=400, detail="path 不能为空")
    target = Path(path_in).expanduser()
    try:
        target = target.resolve()
    except (OSError, ValueError):
        pass
    exists = target.exists()
    size_bytes: Optional[int] = None
    if exists and target.is_file():
        try:
            size_bytes = int(target.stat().st_size)
        except OSError:
            size_bytes = None
    return {
        "status": "success",
        "exists": bool(exists),
        "is_file": target.is_file() if exists else False,
        "is_dir": target.is_dir() if exists else False,
        "path": str(target),
        "size_bytes": size_bytes,
    }


@app.get("/api/tools/check_file")
async def check_file_tool_get(path: str = Query(..., description="主机侧绝对路径")) -> Dict[str, Any]:
    return _check_file_result(path)


@app.post("/api/tools/check_file")
async def check_file_tool_post(payload: CheckFileRequest) -> Dict[str, Any]:
    """推荐：JSON body 传 path，避免 GET 查询串破坏含 \\、#、& 的 Windows 路径。"""
    return _check_file_result(payload.path)


def _primary_path_from_upload_json(data: Dict[str, Any]) -> str:
    fps = data.get("file_paths")
    if isinstance(fps, list) and fps:
        return str(fps[0])
    files = data.get("files")
    if isinstance(files, list) and files:
        fp = files[0].get("file_path") if isinstance(files[0], dict) else None
        if fp:
            return str(fp)
    return ""


@app.post("/api/tools/upload_to_cloud")
async def upload_to_cloud(payload: UploadToCloudRequest) -> Dict[str, Any]:
    """
    将本地文件以 multipart/form-data 转发到宿主云端 POST /api/upload（字段名 files），
    绕过浏览器对大文件与路径的限制。
    """
    raw = str(payload.local_file_path or "").strip()
    if not raw:
        raise HTTPException(status_code=400, detail="local_file_path 不能为空")
    target = Path(raw).expanduser()
    try:
        target = target.resolve()
    except (OSError, ValueError):
        pass
    if not target.exists() or not target.is_file():
        raise HTTPException(status_code=400, detail=f"本地文件不存在或不是文件: {target}")

    base = str(payload.upload_base_url or "").strip().rstrip("/")
    if not base.startswith("http://") and not base.startswith("https://"):
        raise HTTPException(status_code=400, detail="upload_base_url 无效")

    upload_url = f"{base}/api/upload"
    headers: Dict[str, str] = {}
    if payload.authorization:
        headers["Authorization"] = payload.authorization.strip()
    if payload.x_guest_uuid:
        headers["X-Guest-UUID"] = payload.x_guest_uuid.strip()

    timeout = httpx.Timeout(connect=30.0, read=3600.0, write=3600.0, pool=30.0)

    async with httpx.AsyncClient(timeout=timeout, follow_redirects=True) as client:
        with target.open("rb") as fp:
            files = {"files": (target.name, fp, "application/octet-stream")}
            resp = await client.post(upload_url, files=files, headers=headers)

    if resp.status_code >= 400:
        raise HTTPException(
            status_code=502,
            detail=f"云端上传失败 HTTP {resp.status_code}: {resp.text[:800]}",
        )

    try:
        data = resp.json()
    except Exception:
        raise HTTPException(status_code=502, detail="云端返回非 JSON")

    cloud_path = _primary_path_from_upload_json(data if isinstance(data, dict) else {})
    if not cloud_path:
        raise HTTPException(status_code=502, detail="无法从上传响应解析 file_path")

    return {
        "status": "success",
        "cloud_path": cloud_path,
        "raw": data,
    }


@app.get("/api/tools/check_env")
async def check_env() -> Dict[str, Any]:
    """占位：后续可检测本机 Python/conda。当前固定无环境以触发沙盒流程。"""
    return {"has_env": False}


@app.post("/api/tools/install_sandbox")
async def install_sandbox() -> Dict[str, Any]:
    """本地沙盒安装存根。"""
    await asyncio.sleep(1.2)
    return {"status": "success", "message": "轻量级计算沙盒已配置完成（模拟）"}


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
