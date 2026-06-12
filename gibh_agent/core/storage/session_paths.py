# -*- coding: utf-8 -*-
"""会话维度挂载目录树：{mount}/{session_title|session_id}/upload|result。"""
from __future__ import annotations

import re
import unicodedata
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, Optional

_INVALID_FS_CHARS = re.compile(r'[<>:"/\\|?*\x00-\x1f]')
_MULTI_SPACE = re.compile(r"\s+")


def compact_folder_timestamp(dt: Optional[datetime] = None) -> str:
    """紧凑时间戳，用于防同名会话目录冲突：%Y%m%d%H%M（如 202606111628）。"""
    return (dt or datetime.now()).strftime("%Y%m%d%H%M")


def resolve_session_folder_name(
    session_title: Optional[str],
    session_id: str,
    *,
    folder_timestamp: Optional[str] = None,
    session_created_at: Optional[datetime] = None,
) -> str:
    """
    将会话 AI 标题或 session_id 转为安全文件夹名。
    标题非空时在末尾追加 ``-{YYYYMMDDHHMM}`` 以防同名会话覆盖。
    标题为空或清洗后为空时回退 session_id（不追加时间戳）。
    """
    sid = str(session_id or "").strip() or "unknown_session"
    raw = str(session_title or "").strip()
    if not raw:
        return sid
    normalized = unicodedata.normalize("NFKC", raw)
    normalized = _INVALID_FS_CHARS.sub("_", normalized)
    normalized = _MULTI_SPACE.sub(" ", normalized).strip(" .")
    if not normalized:
        return sid
    ts = str(folder_timestamp or "").strip() or compact_folder_timestamp(session_created_at)
    suffix = f"-{ts}"
    max_base = max(1, 80 - len(suffix))
    if len(normalized) > max_base:
        normalized = normalized[:max_base].rstrip(" .")
    return f"{normalized}{suffix}" if normalized else sid


def resolve_session_mount_tree(
    mount_path: str,
    *,
    session_title: Optional[str] = None,
    session_id: str = "",
    folder_timestamp: Optional[str] = None,
    session_created_at: Optional[datetime] = None,
    session_folder: Optional[str] = None,
) -> Dict[str, Any]:
    """返回会话在挂载根下的 upload/ 与 result/ 绝对路径。"""
    base = Path(str(mount_path)).expanduser().resolve()
    folder = str(session_folder or "").strip() or resolve_session_folder_name(
        session_title,
        session_id,
        folder_timestamp=folder_timestamp,
        session_created_at=session_created_at,
    )
    session_root = base / folder
    upload_dir = session_root / "upload"
    result_dir = session_root / "result"
    return {
        "session_folder": folder,
        "session_root": str(session_root),
        "upload_dir": str(upload_dir),
        "result_dir": str(result_dir),
    }


def session_upload_dir(mount_path: str, *, session_title: Optional[str], session_id: str) -> Path:
    tree = resolve_session_mount_tree(mount_path, session_title=session_title, session_id=session_id)
    return Path(tree["upload_dir"])


def session_result_dir(mount_path: str, *, session_title: Optional[str], session_id: str) -> Path:
    tree = resolve_session_mount_tree(mount_path, session_title=session_title, session_id=session_id)
    return Path(tree["result_dir"])
