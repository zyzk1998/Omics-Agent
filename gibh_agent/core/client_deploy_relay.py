# -*- coding: utf-8 -*-
"""浏览器中继落盘：小文件 inline Base64，大文件 staging + URL 流式下载。"""
from __future__ import annotations

import base64
import logging
import os
import secrets
import shutil
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

from gibh_agent.core.local_sidecar_client import read_file_as_b64
from gibh_agent.core.storage.mount_path_resolver import normalize_workspace_mount_root
from gibh_agent.core.storage.session_paths import (
    compact_folder_timestamp,
    resolve_session_folder_name,
)

logger = logging.getLogger(__name__)

INLINE_MAX_BYTES = int(os.getenv("CLIENT_DEPLOY_INLINE_MAX_BYTES", str(1024 * 1024)))
STAGING_TTL_SECONDS = int(os.getenv("CLIENT_DEPLOY_STAGING_TTL_SECONDS", "7200"))


def staging_root_dir() -> Path:
    return Path(os.getenv("CLIENT_DEPLOY_STAGING_DIR", "/app/data/client_deploy_staging"))


@dataclass(frozen=True)
class DeployFileEntry:
    dest_name: str
    source_path: Path


@dataclass
class StagingRecord:
    owner_id: str
    session_id: str
    root: Path
    expires_at: float
    file_map: Dict[str, Path]


_staging_registry: Dict[str, StagingRecord] = {}


def _sanitize_stored_name(dest_name: str) -> str:
    name = str(dest_name or "").strip().replace("\\", "/").lstrip("/")
    if not name or ".." in name.split("/"):
        raise ValueError(f"非法落盘文件名: {dest_name}")
    return name.replace("/", "__")


def _purge_expired_staging() -> None:
    now = time.time()
    expired = [tok for tok, rec in _staging_registry.items() if rec.expires_at <= now]
    for tok in expired:
        rec = _staging_registry.pop(tok, None)
        if rec and rec.root.is_dir():
            try:
                shutil.rmtree(rec.root, ignore_errors=True)
            except OSError:
                pass


def register_staging_record(
    *,
    token: str,
    owner_id: str,
    session_id: str,
    root: Path,
    file_map: Dict[str, Path],
) -> None:
    _purge_expired_staging()
    _staging_registry[token] = StagingRecord(
        owner_id=owner_id,
        session_id=session_id,
        root=root,
        expires_at=time.time() + STAGING_TTL_SECONDS,
        file_map=file_map,
    )


def resolve_staged_file(token: str, stored_name: str, *, owner_id: str) -> Path:
    _purge_expired_staging()
    rec = _staging_registry.get(token)
    if not rec or rec.expires_at <= time.time():
        raise FileNotFoundError("staging token 无效或已过期")
    if rec.owner_id != owner_id:
        raise PermissionError("无权访问该 staging token")
    key = str(stored_name or "").strip()
    if not key or ".." in key or "/" in key or "\\" in key:
        raise FileNotFoundError("staging 文件名非法")
    path = rec.file_map.get(key)
    if not path or not path.is_file():
        raise FileNotFoundError(f"staging 文件不存在: {stored_name}")
    resolved = path.resolve()
    if not str(resolved).startswith(str(rec.root.resolve())):
        raise PermissionError("staging 路径越界")
    return resolved


def collect_ingestion_deploy_entries(
    *,
    archive_path: str,
    bundle_dir: Optional[str] = None,
) -> List[DeployFileEntry]:
    entries: List[DeployFileEntry] = []
    archive = Path(str(archive_path)).expanduser()
    if archive.is_file():
        entries.append(DeployFileEntry(dest_name=archive.name, source_path=archive))
    if bundle_dir:
        src_bundle = Path(str(bundle_dir)).expanduser()
        if src_bundle.is_dir():
            for p in sorted(src_bundle.rglob("*")):
                if p.is_file():
                    rel = p.relative_to(src_bundle).as_posix()
                    entries.append(
                        DeployFileEntry(
                            dest_name=f"{src_bundle.name}/{rel}",
                            source_path=p,
                        )
                    )
    return entries


def collect_corpus_bundle_entries(bundle: Dict[str, Any]) -> List[DeployFileEntry]:
    entries: List[DeployFileEntry] = []

    hitl_path = str(bundle.get("corpus_hitl_path") or "").strip()
    if hitl_path:
        src = Path(hitl_path).expanduser()
        if src.is_file():
            entries.append(DeployFileEntry(dest_name=src.name, source_path=src))
            if src.name != "sft_corpus.json":
                entries.append(DeployFileEntry(dest_name="sft_corpus.json", source_path=src))

    dataset_path = str(bundle.get("dataset_path") or "").strip()
    if dataset_path:
        src_ds = Path(dataset_path).expanduser()
        if src_ds.is_file():
            entries.append(DeployFileEntry(dest_name=src_ds.name, source_path=src_ds))

    archive_dir = str(bundle.get("archive_dir") or "").strip()
    if archive_dir:
        src_dir = Path(archive_dir).expanduser()
        if src_dir.is_dir():
            for p in sorted(src_dir.rglob("*")):
                if p.is_file():
                    rel = p.relative_to(src_dir).as_posix()
                    entries.append(
                        DeployFileEntry(
                            dest_name=f"{src_dir.name}/{rel}",
                            source_path=p,
                        )
                    )
    return entries


def build_client_relay_deploy_package(
    *,
    owner_id: str,
    session_id: str,
    session_title: Optional[str],
    host_mount_path: str,
    entries: List[DeployFileEntry],
    folder_timestamp: Optional[str] = None,
    session_created_at: Optional[Any] = None,
) -> Dict[str, Any]:
    """
    构建浏览器 → Sidecar 中继载荷。
    - 单文件 <= INLINE_MAX_BYTES：inline Base64（files）
    - 更大文件：复制至 staging，返回 download_url（Sidecar 流式拉取，浏览器不经手字节）
    """
    sid = str(session_id or "").strip() or "anonymous"
    host = normalize_workspace_mount_root(host_mount_path)
    oid = str(owner_id or "").strip() or "anonymous"
    ts = str(folder_timestamp or "").strip() or None
    if not ts and session_created_at is not None:
        try:
            ts = compact_folder_timestamp(session_created_at)
        except (TypeError, ValueError, AttributeError):
            ts = None
    if not ts:
        ts = compact_folder_timestamp()
    session_folder = resolve_session_folder_name(
        session_title,
        sid,
        folder_timestamp=ts,
    )

    inline_files: List[Dict[str, str]] = []
    download_items: List[Dict[str, Any]] = []
    token: Optional[str] = None
    staging_root: Optional[Path] = None
    file_map: Dict[str, Path] = {}

    for entry in entries:
        src = entry.source_path.expanduser()
        if not src.is_file():
            continue
        size = int(src.stat().st_size)
        if size <= INLINE_MAX_BYTES:
            inline_files.append(
                {
                    "dest_name": entry.dest_name,
                    "content_b64": read_file_as_b64(src),
                    "size_bytes": size,
                }
            )
            continue

        if token is None:
            token = secrets.token_urlsafe(24)
            staging_root = staging_root_dir() / oid[:16] / sid[:16] / token
            staging_root.mkdir(parents=True, exist_ok=True)

        stored = _sanitize_stored_name(entry.dest_name)
        staged_path = staging_root / stored
        staged_path.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(str(src), str(staged_path))
        file_map[stored] = staged_path
        download_items.append(
            {
                "dest_name": entry.dest_name,
                "download_url": f"/api/ingestion/deploy-staging/{token}/{stored}",
                "size_bytes": size,
            }
        )

    if token and staging_root is not None:
        register_staging_record(
            token=token,
            owner_id=oid,
            session_id=sid,
            root=staging_root,
            file_map=file_map,
        )

    transport = "inline"
    if download_items and inline_files:
        transport = "mixed"
    elif download_items:
        transport = "url"

    package: Dict[str, Any] = {
        "transport": transport,
        "session_id": sid,
        "session_title": session_title or sid,
        "session_folder": session_folder,
        "folder_timestamp": ts,
        "host_mount_path": host,
        "files": inline_files,
        "download_items": download_items,
        "inline_max_bytes": INLINE_MAX_BYTES,
    }
    if token:
        package["staging_token"] = token
    return package
