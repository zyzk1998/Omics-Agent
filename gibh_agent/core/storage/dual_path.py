# -*- coding: utf-8 -*-
"""
双路存储：
- 永久路：用户挂载卷 {session}/upload|result（永不主动删除）
- 临时路：RESULTS_DIR / UPLOAD_DIR 运行时缓存（由 StorageManager 定时清理）
"""
from __future__ import annotations

import logging
import os
import shutil
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence

from gibh_agent.core.storage.session_paths import resolve_session_mount_tree

logger = logging.getLogger(__name__)

_LIGHT_PREVIEW_EXT = frozenset({".png", ".jpg", ".jpeg", ".gif", ".webp", ".svg", ".pdf", ".csv", ".json", ".md"})


def get_temp_results_root() -> Path:
    return Path(os.getenv("RESULTS_DIR", "/app/results")).expanduser().resolve()


def get_temp_upload_root() -> Path:
    return Path(os.getenv("UPLOAD_DIR", "/app/uploads")).expanduser().resolve()


def _file_row(path: Path, *, zone: str) -> Optional[Dict[str, Any]]:
    try:
        if not path.is_file():
            return None
        st = path.stat()
    except OSError:
        return None
    ext = path.suffix.lower()
    return {
        "name": path.name,
        "path": str(path.resolve()),
        "size_bytes": st.st_size,
        "mtime": st.st_mtime,
        "zone": zone,
        "kind": "upload" if zone == "upload" else "result",
        "previewable": ext in _LIGHT_PREVIEW_EXT,
    }


def _resolve_path(raw: str, *, results_root: Path, upload_root: Path) -> Optional[Path]:
    s = str(raw or "").strip()
    if not s:
        return None
    p = Path(s).expanduser()
    if not p.is_absolute():
        if s.startswith("uploads/") or s.startswith("/uploads/"):
            p = (upload_root / s.lstrip("/")).resolve()
        elif s.startswith("results/") or s.startswith("/results/"):
            rel = s.lstrip("/")
            if rel.startswith("results/"):
                rel = rel[len("results/") :]
            p = (results_root / rel).resolve()
        else:
            p = (results_root / s.lstrip("/")).resolve()
    try:
        return p.resolve()
    except OSError:
        return p


def _path_under_session_runtime(path: Path, session_id: str, *, results_root: Path, upload_root: Path) -> bool:
    """防穿透：路径必须落在当前 session 专属目录或显式引用路径下。"""
    sid = str(session_id or "").strip()
    if not sid:
        return False
    try:
        resolved = path.resolve()
    except OSError:
        return False
    res_s = str(resolved)
    roots = (
        results_root / sid,
        results_root / "corpus_archive" / sid,
        results_root / "corpus_hitl" / sid,
    )
    for base in roots:
        try:
            if base.is_dir() and resolved.is_relative_to(base):
                return True
        except (OSError, ValueError):
            if res_s.startswith(str(base.resolve())):
                return True
    if upload_root.is_dir():
        try:
            if resolved.is_relative_to(upload_root.resolve()) and sid in res_s:
                return True
        except (OSError, ValueError):
            if res_s.startswith(str(upload_root.resolve())) and sid in res_s:
                return True
    return False


def _collect_run_dirs_from_references(
    referenced_paths: Sequence[str],
    *,
    results_root: Path,
    upload_root: Path,
) -> set[Path]:
    """仅从快照/消息显式引用的路径反推可扫描的 run_* 目录（禁止全量 glob）。"""
    allowed: set[Path] = set()
    res_root = results_root.resolve()
    for raw in referenced_paths:
        p = _resolve_path(raw, results_root=results_root, upload_root=upload_root)
        if not p:
            continue
        try:
            resolved = p.resolve()
        except OSError:
            continue
        parts = resolved.parts
        for i, part in enumerate(parts):
            if part.startswith("run_"):
                candidate = Path(*parts[: i + 1])
                try:
                    if candidate.resolve().parent == res_root:
                        allowed.add(candidate.resolve())
                except OSError:
                    pass
                break
    return allowed


def list_session_runtime_files(
    session_id: str,
    *,
    referenced_paths: Optional[Sequence[str]] = None,
) -> Dict[str, List[Dict[str, Any]]]:
    """扫描临时缓存区中与 session 关联的上传/结果文件（严格 session 隔离）。"""
    sid = str(session_id or "").strip()
    uploads: List[Dict[str, Any]] = []
    results: List[Dict[str, Any]] = []
    if not sid:
        return {"uploads": uploads, "results": results}

    refs = [str(p).strip() for p in (referenced_paths or []) if str(p).strip()]
    results_root = get_temp_results_root()
    upload_root = get_temp_upload_root()
    seen_res: set[str] = set()
    seen_up: set[str] = set()

    def _append_result(row: Optional[Dict[str, Any]]) -> None:
        if not row:
            return
        key = row.get("path") or ""
        if key and key not in seen_res:
            seen_res.add(key)
            results.append(row)

    def _append_upload(row: Optional[Dict[str, Any]]) -> None:
        if not row:
            return
        key = row.get("path") or ""
        if key and key not in seen_up:
            seen_up.add(key)
            uploads.append(row)

    for sub in ("corpus_archive", "corpus_hitl"):
        d = results_root / sub / sid
        if d.is_dir():
            for p in sorted(d.rglob("*")):
                _append_result(_file_row(p, zone="result"))

    sess_dir = results_root / sid
    if sess_dir.is_dir():
        for p in sorted(sess_dir.rglob("*")):
            _append_result(_file_row(p, zone="result"))

    allowed_runs = _collect_run_dirs_from_references(refs, results_root=results_root, upload_root=upload_root)

    for run_dir in sorted(allowed_runs):
        if not run_dir.is_dir():
            continue
        for p in sorted(run_dir.rglob("*")):
            _append_result(_file_row(p, zone="result"))

    for raw in refs:
        p = _resolve_path(raw, results_root=results_root, upload_root=upload_root)
        if not p or not p.is_file():
            continue
        if _path_under_session_runtime(p, sid, results_root=results_root, upload_root=upload_root):
            row = _file_row(p, zone="result")
            _append_result(row)
            continue
        try:
            resolved = p.resolve()
            in_allowed_run = False
            for run_dir in allowed_runs:
                try:
                    if resolved.is_relative_to(run_dir):
                        in_allowed_run = True
                        break
                except (OSError, ValueError):
                    if str(resolved).startswith(str(run_dir)):
                        in_allowed_run = True
                        break
            if in_allowed_run:
                _append_result(_file_row(p, zone="result"))
        except OSError:
            continue

    if upload_root.is_dir():
        for p in sorted(upload_root.rglob("*")):
            if sid not in str(p):
                continue
            _append_upload(_file_row(p, zone="upload"))

    return {"uploads": uploads[:80], "results": results[:120]}


def mirror_uploads_to_session_mount(
    *,
    mount_path: str,
    session_title: Optional[str],
    session_id: str,
    upload_paths: Sequence[str],
    folder_timestamp: Optional[str] = None,
    session_created_at: Optional[Any] = None,
) -> List[str]:
    """将用户上传文件复制到挂载卷 upload/ 子目录，返回新增/更新路径列表。"""
    tree = resolve_session_mount_tree(
        mount_path,
        session_title=session_title,
        session_id=session_id,
        folder_timestamp=folder_timestamp,
        session_created_at=session_created_at,
    )
    dest_root = Path(tree["upload_dir"])
    dest_root.mkdir(parents=True, exist_ok=True)
    changed: List[str] = []
    seen_names: set[str] = set()
    for raw in upload_paths:
        src = Path(str(raw).expanduser())
        if not src.is_file():
            continue
        name = src.name
        if name in seen_names:
            stem, suf = src.stem, src.suffix
            name = f"{stem}_{len(seen_names)}{suf}"
        seen_names.add(name)
        dest = dest_root / name
        try:
            shutil.copy2(str(src), str(dest))
            changed.append(str(dest.resolve()))
        except OSError as exc:
            logger.warning("[DualPath] skip upload mirror %s -> %s: %s", src, dest, exc)
    return changed


def notify_session_files_changed(
    *,
    session_id: str,
    changed_paths: Sequence[str],
    mount_tree: Optional[Dict[str, str]] = None,
    file_statuses: Optional[Dict[str, str]] = None,
) -> Dict[str, Any]:
    """构造供 API/SSE 下发的前端文件高亮通知载荷。"""
    paths = [str(p) for p in changed_paths if p]
    statuses = file_statuses or {}
    changed_files: List[Dict[str, str]] = []
    for raw in paths:
        resolved_key = raw
        try:
            resolved_key = str(Path(raw).expanduser().resolve())
        except OSError:
            pass
        st = statuses.get(raw) or statuses.get(resolved_key) or "added"
        if st not in ("added", "modified"):
            st = "added"
        changed_files.append({"path": raw, "status": st})
    payload: Dict[str, Any] = {
        "event": "session_files_changed",
        "session_id": session_id,
        "changed_paths": paths,
        "changed_files": changed_files,
        "highlight": True,
    }
    if mount_tree:
        payload["mount_tree"] = mount_tree
    return payload
