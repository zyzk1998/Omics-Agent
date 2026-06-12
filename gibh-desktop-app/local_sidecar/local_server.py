import asyncio
import logging
import os
import platform
import re
import subprocess
import threading
import time
from pathlib import Path
from datetime import datetime
from typing import Any, Dict, List, Optional

import httpx
import uvicorn
from fastapi import FastAPI, HTTPException, Query, Request
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import FileResponse
from pydantic import BaseModel, Field
from starlette.responses import Response


logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
)
logger = logging.getLogger("omics-local-sidecar")

app = FastAPI(title="Omics Agent Local Sidecar", version="1.0.0")

# 本地 Sidecar 仅监听 127.0.0.1；前端可能来自 8018 Web、Electron 或远程 nginx 反代页面。
# allow_credentials=False + 显式 Origin 白名单（禁止裸 * 配合 credentials，避免 ACAO 缺失）。
def _load_sidecar_cors_explicit_origins() -> List[str]:
    base = [
        "http://127.0.0.1:8018",
        "http://localhost:8018",
        "http://127.0.0.1:8019",
        "http://localhost:8019",
        "http://127.0.0.1",
        "http://localhost",
        "null",
    ]
    extra = os.getenv("OMICS_SIDECAR_EXTRA_ORIGINS", "").strip()
    if extra:
        for item in extra.split(","):
            o = item.strip().rstrip("/")
            if o and o not in base:
                base.append(o)
    return base


_SIDECAR_CORS_EXPLICIT_ORIGINS = _load_sidecar_cors_explicit_origins()
_SIDECAR_CORS_ORIGIN_REGEX = (
    r"https?://(?:127\.0\.0\.1|localhost)(?::\d+)?"
    r"|https?://192\.168\.\d{1,3}\.\d{1,3}(?::\d+)?"
    r"|https?://10\.\d{1,3}\.\d{1,3}\.\d{1,3}(?::\d+)?"
    r"|https?://172\.(?:1[6-9]|2\d|3[01])\.\d{1,3}\.\d{1,3}(?::\d+)?"
)

app.add_middleware(
    CORSMiddleware,
    allow_origins=_SIDECAR_CORS_EXPLICIT_ORIGINS,
    allow_origin_regex=_SIDECAR_CORS_ORIGIN_REGEX,
    allow_credentials=False,
    allow_methods=["GET", "POST", "PUT", "PATCH", "DELETE", "OPTIONS", "HEAD"],
    allow_headers=[
        "Authorization",
        "Content-Type",
        "X-Guest-UUID",
        "X-Requested-With",
        "Accept",
        "Origin",
        "Access-Control-Request-Method",
        "Access-Control-Request-Headers",
    ],
    expose_headers=["Content-Disposition", "Content-Length", "Content-Type"],
    max_age=86400,
)


def _sidecar_cors_allow_origin(request: Request) -> str:
    origin = (request.headers.get("origin") or "").strip()
    if not origin:
        return "*"
    if origin == "null":
        return "null"
    if origin in _SIDECAR_CORS_EXPLICIT_ORIGINS:
        return origin
    if re.fullmatch(_SIDECAR_CORS_ORIGIN_REGEX, origin):
        return origin
    return "*"


@app.middleware("http")
async def sidecar_cors_fallback_middleware(request: Request, call_next):
    """
    兜底 CORS：确保 FileResponse、HTTPException 与 OPTIONS 预检均返回跨域头。
    修复历史配置 allow_credentials+\\* 导致 CORSMiddleware 不写 ACAO 的问题。
    """
    allow_origin = _sidecar_cors_allow_origin(request)
    cors_base = {
        "Access-Control-Allow-Origin": allow_origin,
        "Access-Control-Allow-Methods": "GET, POST, PUT, PATCH, DELETE, OPTIONS, HEAD",
        "Access-Control-Allow-Headers": (
            request.headers.get("access-control-request-headers")
            or "Authorization, Content-Type, X-Guest-UUID, X-Requested-With, Accept, Origin"
        ),
        "Access-Control-Max-Age": "86400",
    }
    if allow_origin != "*":
        cors_base["Vary"] = "Origin"

    if request.method == "OPTIONS":
        return Response(status_code=204, headers=cors_base)

    response = await call_next(request)
    for key, value in cors_base.items():
        if key.lower() == "access-control-allow-headers" and value == "*":
            response.headers.setdefault(key, "*")
        else:
            response.headers.setdefault(key, value)
    return response


class WorkspaceInitRequest(BaseModel):
    workspace_path: str


class ReadFileRequest(BaseModel):
    file_path: str
    max_lines: int = 1000


class WriteFileRequest(BaseModel):
    file_path: str
    content: str = ""
    content_b64: Optional[str] = None


class UploadToCloudRequest(BaseModel):
    """由 Sidecar 从用户磁盘读取并 multipart 转发至云端 /api/upload。"""

    local_file_path: str
    upload_base_url: str
    authorization: Optional[str] = None
    x_guest_uuid: Optional[str] = None


class CheckFileRequest(BaseModel):
    """POST /api/tools/check_file 使用 JSON 传递路径，避免 GET query 对 Windows 反斜杠等字符的编码坑。"""

    path: str


class WorkspaceFilePathRequest(BaseModel):
    """POST /api/files/* 使用 JSON 传递路径，与 check_file 一致，避免 Windows 反斜杠在 GET query 中被破坏。"""

    path: str


class SilentDeployFileSpec(BaseModel):
    dest_name: str
    content_b64: str


class SilentDeployRequest(BaseModel):
    session_id: str
    session_title: Optional[str] = None
    session_folder: Optional[str] = None
    folder_timestamp: Optional[str] = None
    host_mount_path: Optional[str] = None
    files: List[SilentDeployFileSpec] = Field(default_factory=list)


class DeployUrlItem(BaseModel):
    dest_name: str
    download_url: str
    authorization: Optional[str] = None
    x_guest_uuid: Optional[str] = None


class SilentDeployFromUrlRequest(BaseModel):
    session_id: str
    session_title: Optional[str] = None
    session_folder: Optional[str] = None
    folder_timestamp: Optional[str] = None
    host_mount_path: Optional[str] = None
    items: List[DeployUrlItem] = Field(default_factory=list)
    files: List[SilentDeployFileSpec] = Field(default_factory=list)


class LocalResultCheckRequest(BaseModel):
    host_mount_path: Optional[str] = None
    session_id: str
    session_title: Optional[str] = None
    session_folder: Optional[str] = None
    folder_timestamp: Optional[str] = None


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
    return {
        "status": "ok",
        "service": "omics-local-sidecar",
        "platform": platform.system(),
        "version": "1.1.0",
        "capabilities": [
            "workspace_tree",
            "check_file",
            "files_download",
            "files_read_text",
            "upload_to_cloud",
            "write_file",
            "write_bytes",
            "silent_deploy",
            "silent_deploy_from_url",
            "local_result_check",
            "fs_open_file",
            "fs_show_in_explorer",
            "fs_delete_file",
            "fs_copy_file",
        ],
    }


@app.post("/api/shutdown")
async def shutdown_sidecar() -> Dict[str, Any]:
    """Electron 退出前调用：短暂返回响应后 os._exit，释放 8019，避免残留僵尸 Sidecar。"""

    def _exit_soon() -> None:
        time.sleep(0.08)
        os._exit(0)

    threading.Thread(target=_exit_soon, daemon=True).start()
    return {"status": "ok", "message": "terminating"}


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


def _resolve_workspace_file(path: str) -> Path:
    raw = str(path or "").strip().strip('"').strip("'")
    if not raw:
        raise HTTPException(status_code=400, detail="path 不能为空")
    target = Path(raw).expanduser().resolve()
    if not target.is_file():
        raise HTTPException(status_code=400, detail=f"文件不存在或不可读: {target}")
    root_raw = (_workspace_context.get("workspace_path") or "").strip()
    if root_raw:
        root = Path(root_raw).expanduser().resolve()
        try:
            target.relative_to(root)
            return target
        except ValueError:
            pass
    return target


def _resolve_workspace_entry(path: str) -> Path:
    """解析工作区内的文件或目录路径（供 show_in_explorer 等）。"""
    raw = str(path or "").strip().strip('"').strip("'")
    if not raw:
        raise HTTPException(status_code=400, detail="path 不能为空")
    target = Path(raw).expanduser().resolve()
    if not target.exists():
        raise HTTPException(status_code=400, detail=f"路径不存在或不可读: {target}")
    root_raw = (_workspace_context.get("workspace_path") or "").strip()
    if root_raw:
        root = Path(root_raw).expanduser().resolve()
        try:
            target.relative_to(root)
            return target
        except ValueError:
            pass
    return target


def _workspace_file_media_type(target: Path) -> str:
    ext = target.suffix.lower().lstrip(".")
    return {
        "png": "image/png",
        "jpg": "image/jpeg",
        "jpeg": "image/jpeg",
        "jfif": "image/jpeg",
        "jpe": "image/jpeg",
        "pjpeg": "image/jpeg",
        "gif": "image/gif",
        "webp": "image/webp",
        "svg": "image/svg+xml",
        "svgz": "image/svg+xml",
        "bmp": "image/bmp",
        "dib": "image/bmp",
        "ico": "image/x-icon",
        "cur": "image/x-icon",
        "tif": "image/tiff",
        "tiff": "image/tiff",
        "avif": "image/avif",
        "heic": "image/heic",
        "heif": "image/heif",
        "apng": "image/apng",
        "jxl": "image/jxl",
        "jp2": "image/jp2",
        "j2k": "image/jp2",
        "jpx": "image/jpx",
        "psd": "image/vnd.adobe.photoshop",
        "eps": "application/postscript",
        "emf": "image/emf",
        "wmf": "image/wmf",
        "pdf": "application/pdf",
    }.get(ext, "application/octet-stream")


def _workspace_file_download_response(target: Path) -> FileResponse:
    return FileResponse(
        path=str(target),
        media_type=_workspace_file_media_type(target),
        filename=target.name,
    )


@app.get("/api/files/download")
async def download_workspace_file(path: str = Query(...)):
    """本地工作区文件内联下载/预览（图片、PDF 等）。推荐 POST + JSON body 传 path（Windows 路径更可靠）。"""
    target = _resolve_workspace_file(path)
    return _workspace_file_download_response(target)


@app.post("/api/files/download")
async def download_workspace_file_post(payload: WorkspaceFilePathRequest):
    """JSON body 传 path，避免 GET 查询串破坏含 \\、#、& 的 Windows 路径（与 check_file 一致）。"""
    target = _resolve_workspace_file(payload.path)
    return _workspace_file_download_response(target)


@app.post("/api/tools/download_file")
async def download_file_tool_alias(payload: WorkspaceFilePathRequest):
    """与 POST /api/files/download 等价；兼容旧客户端或工具链别名。"""
    target = _resolve_workspace_file(payload.path)
    return _workspace_file_download_response(target)


@app.post("/api/files/read_text")
async def read_workspace_file_text_post(
    payload: WorkspaceFilePathRequest,
    max_bytes: int = Query(512_000, ge=1024, le=2_000_000),
):
    """JSON body 传 path，避免 Windows GET query 编码问题。"""
    target = _resolve_workspace_file(payload.path)
    return _read_workspace_file_text_result(target, max_bytes)


@app.get("/api/files/read_text")
async def read_workspace_file_text(
    path: str = Query(...),
    max_bytes: int = Query(512_000, ge=1024, le=2_000_000),
):
    target = _resolve_workspace_file(path)
    return _read_workspace_file_text_result(target, max_bytes)


def _read_workspace_file_text_result(target: Path, max_bytes: int) -> Dict[str, Any]:
    ext = target.suffix.lower().lstrip(".")
    allowed = {
        "txt", "md", "markdown", "py", "csv", "tsv", "json", "yaml", "yml",
        "html", "htm", "xml", "log", "ini", "sh", "js", "ts", "sql", "toml", "env",
    }
    if ext not in allowed:
        raise HTTPException(status_code=400, detail="该扩展名不支持文本预览")
    size = target.stat().st_size
    truncated = size > max_bytes
    content = target.read_bytes()[:max_bytes].decode("utf-8", errors="replace")
    return {
        "status": "success",
        "file_path": str(target),
        "file_name": target.name,
        "content": content,
        "truncated": truncated,
        "size_bytes": size,
    }


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
    if payload.content_b64:
        import base64

        try:
            data = base64.b64decode(str(payload.content_b64))
        except (ValueError, TypeError) as exc:
            raise HTTPException(status_code=400, detail=f"content_b64 无效: {exc}") from exc
        target.write_bytes(data)
        return {
            "status": "success",
            "tool_name": "write_local_file",
            "file_path": str(target),
            "bytes_written": len(data),
            "message": f"已写入本地文件: {target}",
        }
    text = str(payload.content or "")
    target.write_text(text, encoding="utf-8")
    return {
        "status": "success",
        "tool_name": "write_local_file",
        "file_path": str(target),
        "bytes_written": len(text.encode("utf-8")),
        "message": f"已写入本地文件: {target}",
    }


@app.post("/api/tools/write_bytes")
async def write_bytes_tool(payload: WriteFileRequest) -> Dict[str, Any]:
    """二进制落盘（content_b64）；与 write_file 带 content_b64 等价，便于能力探测。"""
    if not payload.content_b64:
        raise HTTPException(status_code=400, detail="content_b64 不能为空")
    return await write_file_tool(payload)


def _assert_loopback_client(request: Request) -> None:
    """原生 FS 操作仅允许本机访问（Sidecar 绑定 127.0.0.1）。"""
    client = request.client
    host = (client.host if client else "") or ""
    if host not in ("127.0.0.1", "::1", "localhost"):
        raise HTTPException(status_code=403, detail="原生文件操作仅允许本机访问")


def _assert_deploy_client(request: Request) -> None:
    """
    静默落盘/检测：允许本机 + Docker 网桥（API 容器经 host.docker.internal 回调 Sidecar）。
    """
    client = request.client
    host = (client.host if client else "") or ""
    if host in ("127.0.0.1", "::1", "localhost"):
        return
    if os.getenv("OMICS_SIDECAR_ALLOW_DOCKER", "1") == "1":
        if host.startswith("172.") or host.startswith("192.168.") or host.startswith("10."):
            return
    _assert_loopback_client(request)


def _open_file_native(target: Path) -> None:
    path_str = str(target)
    system = platform.system()
    if system == "Windows":
        os.startfile(path_str)  # type: ignore[attr-defined]
    elif system == "Darwin":
        subprocess.run(["open", path_str], check=False)
    else:
        subprocess.run(["xdg-open", path_str], check=False)


def _show_in_explorer_native(target: Path) -> None:
    path_str = str(target.resolve())
    system = platform.system()
    if target.is_dir():
        if system == "Windows":
            subprocess.run(["explorer", os.path.normpath(path_str)], check=False)
        elif system == "Darwin":
            subprocess.run(["open", path_str], check=False)
        else:
            subprocess.run(["xdg-open", path_str], check=False)
        return
    if system == "Windows":
        subprocess.run(["explorer", "/select,", os.path.normpath(path_str)], check=False)
    elif system == "Darwin":
        subprocess.run(["open", "-R", path_str], check=False)
    else:
        parent = str(target.parent.resolve())
        subprocess.run(["xdg-open", parent], check=False)


@app.post("/api/fs/open_file")
async def fs_open_file(request: Request, payload: WorkspaceFilePathRequest) -> Dict[str, Any]:
    """使用系统默认程序打开本地文件。"""
    _assert_loopback_client(request)
    target = _resolve_workspace_file(payload.path)
    try:
        _open_file_native(target)
    except OSError as exc:
        logger.warning("[FS] open_file failed: %s (%s)", target, exc)
        raise HTTPException(status_code=500, detail=f"无法打开文件: {exc}") from exc
    return {"status": "success", "path": str(target), "action": "open_file"}


@app.post("/api/fs/show_in_explorer")
async def fs_show_in_explorer(request: Request, payload: WorkspaceFilePathRequest) -> Dict[str, Any]:
    """在文件资源管理器中定位并选中文件或打开目录。"""
    _assert_loopback_client(request)
    target = _resolve_workspace_entry(payload.path)
    try:
        _show_in_explorer_native(target)
    except OSError as exc:
        logger.warning("[FS] show_in_explorer failed: %s (%s)", target, exc)
        raise HTTPException(status_code=500, detail=f"无法打开文件所在目录: {exc}") from exc
    return {"status": "success", "path": str(target), "action": "show_in_explorer"}


@app.post("/api/fs/delete_file")
async def fs_delete_file(request: Request, payload: WorkspaceFilePathRequest) -> Dict[str, Any]:
    """删除本地文件（仅限普通文件，非目录）。"""
    _assert_loopback_client(request)
    target = _resolve_workspace_file(payload.path)
    if target.is_dir():
        raise HTTPException(status_code=400, detail="不支持删除目录，请仅指定文件路径")
    try:
        target.unlink()
    except OSError as exc:
        logger.warning("[FS] delete_file failed: %s (%s)", target, exc)
        raise HTTPException(status_code=500, detail=f"删除失败: {exc}") from exc
    return {"status": "success", "path": str(target), "action": "delete_file"}


@app.post("/api/fs/copy_file")
async def fs_copy_file(request: Request, payload: WorkspaceFilePathRequest) -> Dict[str, Any]:
    """返回可复制到剪贴板的文本（路径；小文本文件则返回内容）。"""
    _assert_loopback_client(request)
    target = _resolve_workspace_file(payload.path)
    clipboard_text = str(target)
    text_ext = {
        "txt", "md", "markdown", "json", "csv", "tsv", "yaml", "yml", "log", "py", "html", "htm", "xml",
    }
    try:
        size = target.stat().st_size
        if target.suffix.lower().lstrip(".") in text_ext and size <= 512_000:
            clipboard_text = target.read_text(encoding="utf-8", errors="replace")
    except OSError:
        pass
    return {
        "status": "success",
        "path": str(target),
        "clipboard_text": clipboard_text,
        "action": "copy_file",
    }


def _resolve_host_mount_path(raw: Optional[str]) -> Path:
    host = str(raw or _workspace_context.get("workspace_path") or "").strip()
    if not host:
        raise HTTPException(status_code=400, detail="未挂载本地工作区")
    root = Path(host).expanduser().resolve()
    if not root.is_dir():
        raise HTTPException(status_code=400, detail=f"工作区路径不存在: {root}")
    return root


def _session_result_dir_on_host(
    host_root: Path,
    *,
    session_id: str,
    session_title: Optional[str],
    session_folder: Optional[str] = None,
    folder_timestamp: Optional[str] = None,
) -> Dict[str, Any]:
    folder = _safe_session_folder_name(
        session_title,
        session_id,
        session_folder=session_folder,
        folder_timestamp=folder_timestamp,
    )
    session_root = host_root / folder
    upload_dir = session_root / "upload"
    result_dir = session_root / "result"
    result_dir.mkdir(parents=True, exist_ok=True)
    upload_dir.mkdir(parents=True, exist_ok=True)
    return {
        "session_folder": folder,
        "session_root": str(session_root.resolve()),
        "upload_dir": str(upload_dir.resolve()),
        "result_dir": str(result_dir.resolve()),
        "result_dir_path": result_dir,
    }


def _safe_session_folder_name(
    session_title: Optional[str],
    session_id: str,
    *,
    session_folder: Optional[str] = None,
    folder_timestamp: Optional[str] = None,
) -> str:
    explicit = str(session_folder or "").strip()
    if explicit:
        return explicit
    sid = str(session_id or "").strip() or "unknown_session"
    raw = str(session_title or "").strip()
    if not raw:
        return sid
    cleaned = re.sub(r'[<>:"/\\|?*\x00-\x1f]', "_", raw)
    cleaned = re.sub(r"\s+", " ", cleaned).strip(" .")
    if not cleaned:
        return sid
    ts = str(folder_timestamp or "").strip() or datetime.now().strftime("%Y%m%d%H%M")
    suffix = f"-{ts}"
    max_base = max(1, 80 - len(suffix))
    if len(cleaned) > max_base:
        cleaned = cleaned[:max_base].rstrip(" .")
    return f"{cleaned}{suffix}" if cleaned else sid


@app.post("/api/workspace/silent_deploy")
async def workspace_silent_deploy(request: Request, payload: SilentDeployRequest) -> Dict[str, Any]:
    """宿主机静默落盘：将 inline Base64 小文件写入 {mount}/{session}/result/。"""
    _assert_deploy_client(request)
    host_root = _resolve_host_mount_path(payload.host_mount_path)
    tree = _session_result_dir_on_host(
        host_root,
        session_id=payload.session_id,
        session_title=payload.session_title,
        session_folder=payload.session_folder,
        folder_timestamp=payload.folder_timestamp,
    )
    result_dir: Path = tree["result_dir_path"]
    changed = _write_silent_deploy_inline_files(result_dir, payload.files or [])
    if not changed:
        raise HTTPException(status_code=400, detail="无有效文件可落盘")
    return _silent_deploy_success_payload(host_root, tree, payload.session_id, changed)


def _write_silent_deploy_inline_files(result_dir: Path, files: List[SilentDeployFileSpec]) -> List[str]:
    import base64

    changed: List[str] = []
    for spec in files or []:
        name = str(spec.dest_name or "").strip().replace("\\", "/").lstrip("/")
        if not name or ".." in name.split("/"):
            continue
        dest = result_dir / name
        dest.parent.mkdir(parents=True, exist_ok=True)
        try:
            dest.write_bytes(base64.b64decode(spec.content_b64))
        except (ValueError, OSError) as exc:
            raise HTTPException(status_code=400, detail=f"写入失败 {name}: {exc}") from exc
        changed.append(str(dest.resolve()))
    return changed


async def _stream_download_to_path(
    *,
    url: str,
    dest: Path,
    authorization: Optional[str] = None,
    x_guest_uuid: Optional[str] = None,
) -> None:
    headers: Dict[str, str] = {}
    if authorization:
        headers["Authorization"] = str(authorization).strip()
    if x_guest_uuid:
        headers["X-Guest-UUID"] = str(x_guest_uuid).strip()
    timeout = httpx.Timeout(connect=30.0, read=3600.0, write=3600.0, pool=30.0)
    dest.parent.mkdir(parents=True, exist_ok=True)
    async with httpx.AsyncClient(timeout=timeout, follow_redirects=True) as client:
        async with client.stream("GET", url, headers=headers) as resp:
            if resp.status_code >= 400:
                body_preview = (await resp.aread())[:800]
                raise HTTPException(
                    status_code=502,
                    detail=f"云端下载失败 HTTP {resp.status_code}: {body_preview!r}",
                )
            with dest.open("wb") as fp:
                async for chunk in resp.aiter_bytes(256 * 1024):
                    fp.write(chunk)


def _silent_deploy_success_payload(
    host_root: Path,
    tree: Dict[str, Any],
    session_id: str,
    changed: List[str],
) -> Dict[str, Any]:
    result_dir = tree["result_dir_path"]
    logger.info("[SilentDeploy/Sidecar] session=%s copied=%d -> %s", session_id, len(changed), result_dir)
    return {
        "status": "success",
        "changed_paths": changed,
        "copied_paths": changed,
        "result_dir": str(result_dir.resolve()),
        "destination_dir": str(result_dir.resolve()),
        "session_root": tree["session_root"],
        "session_folder": tree["session_folder"],
        "mount_tree": {
            "session_folder": tree["session_folder"],
            "session_root": tree["session_root"],
            "upload_dir": tree["upload_dir"],
            "result_dir": tree["result_dir"],
        },
        "mount_path": str(host_root),
    }


@app.post("/api/workspace/silent_deploy_from_url")
async def workspace_silent_deploy_from_url(
    request: Request,
    payload: SilentDeployFromUrlRequest,
) -> Dict[str, Any]:
    """大文件中继：Sidecar 从云端 staging URL 流式下载并写入本地 result/（浏览器不经手字节）。"""
    _assert_deploy_client(request)
    host_root = _resolve_host_mount_path(payload.host_mount_path)
    tree = _session_result_dir_on_host(
        host_root,
        session_id=payload.session_id,
        session_title=payload.session_title,
        session_folder=payload.session_folder,
        folder_timestamp=payload.folder_timestamp,
    )
    result_dir: Path = tree["result_dir_path"]
    changed: List[str] = []

    for item in payload.items or []:
        name = str(item.dest_name or "").strip().replace("\\", "/").lstrip("/")
        url = str(item.download_url or "").strip()
        if not name or ".." in name.split("/") or not url.startswith(("http://", "https://")):
            continue
        dest = result_dir / name
        await _stream_download_to_path(
            url=url,
            dest=dest,
            authorization=item.authorization,
            x_guest_uuid=item.x_guest_uuid,
        )
        changed.append(str(dest.resolve()))

    if payload.files:
        changed.extend(_write_silent_deploy_inline_files(result_dir, payload.files))

    if not changed:
        raise HTTPException(status_code=400, detail="无有效文件可落盘")

    return _silent_deploy_success_payload(host_root, tree, payload.session_id, changed)


@app.post("/api/workspace/local_result_check")
async def workspace_local_result_check(request: Request, payload: LocalResultCheckRequest) -> Dict[str, Any]:
    """检测宿主机挂载目录下是否已有 sft_corpus 语料。"""
    _assert_deploy_client(request)
    host_root = _resolve_host_mount_path(payload.host_mount_path)
    tree = _session_result_dir_on_host(
        host_root,
        session_id=payload.session_id,
        session_title=payload.session_title,
        session_folder=payload.session_folder,
        folder_timestamp=payload.folder_timestamp,
    )
    result_dir: Path = tree["result_dir_path"]
    corpus_glob = sorted(result_dir.glob("sft_corpus*.json"))
    return {
        "status": "success",
        "has_local_result": bool(corpus_glob),
        "mount_path": str(host_root),
        "host_mount_path": str(host_root),
        "result_dir": str(result_dir.resolve()),
        "session_root": tree["session_root"],
        "session_folder": tree["session_folder"],
        "sft_corpus_paths": [str(p.resolve()) for p in corpus_glob],
        "mount_tree": {
            "session_folder": tree["session_folder"],
            "session_root": tree["session_root"],
            "upload_dir": tree["upload_dir"],
            "result_dir": tree["result_dir"],
        },
        "deploy_mode": "sidecar",
    }


if __name__ == "__main__":
    _port = int(os.getenv("OMICS_SIDECAR_PORT", "8019"))
    uvicorn.run(app, host="127.0.0.1", port=_port, log_level="info")
