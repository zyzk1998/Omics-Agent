# -*- coding: utf-8 -*-
"""从会话消息与快照汇总 upload/result 文件清单（供 API 与前端手风琴）。"""
from __future__ import annotations

import json
import os
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence

from gibh_agent.core.storage.dual_path import (
    get_temp_results_root,
    list_session_runtime_files,
)
from gibh_agent.core.storage.session_paths import resolve_session_mount_tree


def _collect_paths_from_obj(obj: Any, out: set[str]) -> None:
    if obj is None:
        return
    if isinstance(obj, str):
        s = obj.strip()
        if s.startswith(("/", "uploads/", "results/")) or os.path.isabs(s):
            if any(s.lower().endswith(ext) for ext in (".png", ".jpg", ".jpeg", ".csv", ".tsv", ".fastq", ".mtx", ".json", ".md", ".pdf", ".html", ".h5ad", ".nii.gz")):
                out.add(s)
            elif os.path.isfile(s):
                out.add(s)
        return
    if isinstance(obj, dict):
        for k, v in obj.items():
            if k in ("path", "file_path", "data_path", "image_path", "mask_path", "sequence_or_path"):
                if isinstance(v, str) and v.strip():
                    out.add(v.strip())
            _collect_paths_from_obj(v, out)
        return
    if isinstance(obj, (list, tuple)):
        for item in obj:
            _collect_paths_from_obj(item, out)


def _message_content_dict(msg: Any) -> Optional[Dict[str, Any]]:
    if isinstance(msg, dict):
        content = msg.get("content", msg)
        if isinstance(content, str):
            try:
                content = json.loads(content)
            except (json.JSONDecodeError, TypeError):
                return None
        return content if isinstance(content, dict) else None
    content = getattr(msg, "content", None)
    if isinstance(content, str):
        try:
            content = json.loads(content)
        except (json.JSONDecodeError, TypeError):
            return None
    return content if isinstance(content, dict) else None


def extract_all_paths_from_messages(messages: Sequence[Any]) -> List[str]:
    """从会话全部消息（含 agent 快照）收集路径引用，供运行时文件防穿透过滤。"""
    paths: set[str] = set()
    for msg in messages:
        content = _message_content_dict(msg)
        if content:
            _collect_paths_from_obj(content, paths)
            for key in ("uploaded_files", "attachments", "input_draft_attachments", "files", "images"):
                _collect_paths_from_obj(content.get(key), paths)
        tool_calls = None
        if isinstance(msg, dict):
            tool_calls = msg.get("tool_calls")
        else:
            tool_calls = getattr(msg, "tool_calls", None)
        for tc in tool_calls or []:
            fn = tc.get("function") if isinstance(tc, dict) else getattr(tc, "function", None)
            if fn is None:
                continue
            args_raw = fn.get("arguments") if isinstance(fn, dict) else getattr(fn, "arguments", None)
            if not args_raw:
                continue
            try:
                args_obj = json.loads(args_raw) if isinstance(args_raw, str) else args_raw
            except (json.JSONDecodeError, TypeError):
                continue
            _collect_paths_from_obj(args_obj, paths)
        metadata = None
        if isinstance(msg, dict):
            metadata = msg.get("metadata")
        else:
            metadata = getattr(msg, "metadata", None) or getattr(msg, "message_metadata", None)
        if metadata:
            _collect_paths_from_obj(metadata, paths)
    return sorted(paths)


def extract_upload_paths_from_messages(messages: Sequence[Any]) -> List[str]:
    paths: set[str] = set()
    for msg in messages:
        content = getattr(msg, "content", msg)
        if isinstance(content, str):
            try:
                content = json.loads(content)
            except (json.JSONDecodeError, TypeError):
                continue
        if not isinstance(content, dict):
            continue
        role = getattr(msg, "role", None)
        if role is None and isinstance(msg, dict):
            role = msg.get("role")
        if role and role != "user":
            continue
        _collect_paths_from_obj(content, paths)
        for key in ("uploaded_files", "attachments", "input_draft_attachments", "files"):
            _collect_paths_from_obj(content.get(key), paths)
    return sorted(paths)


def _scan_mount_zone(root: Path, zone: str, limit: int = 100) -> List[Dict[str, Any]]:
    rows: List[Dict[str, Any]] = []
    if not root.is_dir():
        return rows
    for p in sorted(root.rglob("*")):
        if not p.is_file():
            continue
        try:
            st = p.stat()
        except OSError:
            continue
        rows.append(
            {
                "name": p.name,
                "path": str(p.resolve()),
                "size_bytes": st.st_size,
                "mtime": st.st_mtime,
                "zone": zone,
                "kind": zone,
                "storage_tier": "permanent",
                "previewable": p.suffix.lower() in {".png", ".jpg", ".jpeg", ".gif", ".webp", ".svg", ".pdf", ".csv", ".json", ".md"},
            }
        )
        if len(rows) >= limit:
            break
    return rows


def build_session_files_inventory(
    *,
    session_id: str,
    session_title: Optional[str],
    mount_path: Optional[str],
    messages: Sequence[Any],
    include_uploads: bool = False,
) -> Dict[str, Any]:
    """汇总临时缓存 + 挂载永久目录下的会话输出文件（默认不含 upload）。"""
    sid = str(session_id or "").strip()
    referenced_paths = extract_all_paths_from_messages(messages)
    runtime = list_session_runtime_files(sid, referenced_paths=referenced_paths)

    uploads: List[Dict[str, Any]] = []
    results: List[Dict[str, Any]] = list(runtime.get("results") or [])

    if include_uploads:
        uploads = list(runtime.get("uploads") or [])
        upload_paths = extract_upload_paths_from_messages(messages)
        seen_upload = {u.get("path") for u in uploads}
        upload_root = Path(os.getenv("UPLOAD_DIR", "/app/uploads")).expanduser()
        for raw in upload_paths:
            p = Path(raw)
            if not p.is_absolute():
                p = (upload_root / raw.lstrip("/")).resolve()
            if not p.is_file() or str(p) in seen_upload:
                continue
            try:
                st = p.stat()
            except OSError:
                continue
            uploads.append(
                {
                    "name": p.name,
                    "path": str(p.resolve()),
                    "size_bytes": st.st_size,
                    "mtime": st.st_mtime,
                    "zone": "upload",
                    "kind": "upload",
                    "storage_tier": "runtime",
                    "previewable": p.suffix.lower() in {".png", ".jpg", ".jpeg", ".gif", ".webp", ".svg", ".pdf", ".csv", ".json", ".md"},
                }
            )
            seen_upload.add(str(p))

    mount_tree = None
    if mount_path:
        mount_tree = resolve_session_mount_tree(
            mount_path,
            session_title=session_title,
            session_id=sid,
        )
        res_rows = _scan_mount_zone(Path(mount_tree["result_dir"]), "result")
        seen_res = {r.get("path") for r in results}
        for row in res_rows:
            if row["path"] not in seen_res:
                results.append(row)
        if include_uploads:
            up_rows = _scan_mount_zone(Path(mount_tree["upload_dir"]), "upload")
            seen_up = {u.get("path") for u in uploads}
            for row in up_rows:
                if row["path"] not in seen_up:
                    uploads.append(row)

    results_root = get_temp_results_root()
    return {
        "session_id": sid,
        "session_title": session_title or sid,
        "uploads": uploads if include_uploads else [],
        "results": results,
        "outputs": results,
        "mount_tree": mount_tree,
        "temp_cache_root": str(results_root),
    }
