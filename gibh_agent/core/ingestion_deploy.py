# -*- coding: utf-8 -*-
"""入库物理落盘：会话维度 upload/result 树 + tar.gz 归档。"""
from __future__ import annotations

import logging
import shutil
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence

from gibh_agent.core.ingestion_mount_discovery import validate_ingestion_mount_path
from gibh_agent.core.client_deploy_relay import (
    build_client_relay_deploy_package,
    collect_ingestion_deploy_entries,
)
from gibh_agent.core.local_sidecar_client import sidecar_silent_deploy
from gibh_agent.core.storage.dual_path import mirror_uploads_to_session_mount
from gibh_agent.core.storage.session_paths import resolve_session_mount_tree

logger = logging.getLogger(__name__)


def deploy_artifacts_to_mount(
    *,
    archive_path: str,
    mount_path: str,
    session_id: str = "",
    owner_id: str = "",
    bundle_dir: Optional[str] = None,
    session_title: Optional[str] = None,
    upload_paths: Optional[Sequence[str]] = None,
    folder_timestamp: Optional[str] = None,
    session_created_at: Optional[Any] = None,
) -> Dict[str, Any]:
    """
    将 tar.gz 归档与打包工作目录复制到挂载目录的会话 result/ 子树；
    可选同步 upload_paths 至 upload/。

    目标结构::

        {mount_path}/{session_title|session_id}/
            upload/    # 用户上传
            result/    # bundle_*.tar.gz + bundle_*/...
    """
    validated = validate_ingestion_mount_path(mount_path)
    if not validated:
        raise ValueError(f"挂载路径不可用或不可写: {mount_path}")

    archive = Path(str(archive_path)).expanduser()
    if not archive.is_file():
        raise ValueError(f"归档文件不存在: {archive}")

    tree = resolve_session_mount_tree(
        validated,
        session_title=session_title,
        session_id=session_id or (f"owner-{owner_id[:8]}" if owner_id else "anonymous"),
        folder_timestamp=folder_timestamp,
        session_created_at=session_created_at,
    )
    upload_dir = Path(tree["upload_dir"])
    result_dir = Path(tree["result_dir"])
    upload_dir.mkdir(parents=True, exist_ok=True)
    result_dir.mkdir(parents=True, exist_ok=True)

    changed_paths: List[str] = []

    if upload_paths:
        mirrored = mirror_uploads_to_session_mount(
            mount_path=validated,
            session_title=session_title,
            session_id=session_id,
            upload_paths=upload_paths,
            folder_timestamp=folder_timestamp,
            session_created_at=session_created_at,
        )
        changed_paths.extend(mirrored)

    dest_file = result_dir / archive.name
    shutil.copy2(str(archive), str(dest_file))
    changed_paths.append(str(dest_file.resolve()))

    bundle_dest: Optional[Path] = None
    if bundle_dir:
        src_bundle = Path(str(bundle_dir)).expanduser()
        if src_bundle.is_dir():
            bundle_dest = result_dir / src_bundle.name
            if bundle_dest.exists():
                shutil.rmtree(bundle_dest)
            shutil.copytree(str(src_bundle), str(bundle_dest))
            changed_paths.append(str(bundle_dest.resolve()))
            logger.info(
                "[IngestionDeploy] bundle copied session=%s dest=%s",
                session_id or "na",
                bundle_dest,
            )

    logger.info(
        "[IngestionDeploy] archive copied session=%s archive=%s dest=%s",
        session_id or "na",
        archive.name,
        dest_file,
    )
    return {
        "strategy": "local_volume",
        "destination": str(dest_file.resolve()),
        "destination_dir": str(result_dir.resolve()),
        "session_root": tree["session_root"],
        "upload_dir": tree["upload_dir"],
        "result_dir": tree["result_dir"],
        "session_folder": tree["session_folder"],
        "bundle_destination": str(bundle_dest.resolve()) if bundle_dest else None,
        "copied_paths": changed_paths,
        "changed_paths": changed_paths,
        "mount_path": validated,
        "mount_tree": tree,
    }


def build_sidecar_deploy_package(
    *,
    archive_path: str,
    host_mount_path: str,
    session_id: str = "",
    owner_id: str = "",
    bundle_dir: Optional[str] = None,
    session_title: Optional[str] = None,
    folder_timestamp: Optional[str] = None,
    session_created_at: Optional[Any] = None,
) -> Dict[str, Any]:
    """构建供浏览器中继至 Local Sidecar 的落盘载荷（小文件 inline，大文件 staging URL）。"""
    archive = Path(str(archive_path)).expanduser()
    if not archive.is_file():
        raise ValueError(f"归档文件不存在: {archive}")

    entries = collect_ingestion_deploy_entries(
        archive_path=str(archive_path),
        bundle_dir=bundle_dir,
    )
    if not entries:
        raise ValueError("无有效文件可落盘")

    return build_client_relay_deploy_package(
        owner_id=owner_id,
        session_id=session_id or (f"owner-{owner_id[:8]}" if owner_id else "anonymous"),
        session_title=session_title,
        host_mount_path=host_mount_path,
        entries=entries,
        folder_timestamp=folder_timestamp,
        session_created_at=session_created_at,
    )


def delivery_from_sidecar_body(body: Dict[str, Any], *, host_mount_path: str) -> Dict[str, Any]:
    changed = list(body.get("changed_paths") or body.get("copied_paths") or [])
    mount_tree = body.get("mount_tree") or {}
    return {
        "strategy": "local_sidecar",
        "deploy_mode": "client_relay",
        "destination": changed[0] if changed else "",
        "destination_dir": body.get("result_dir") or body.get("destination_dir"),
        "session_root": body.get("session_root"),
        "result_dir": body.get("result_dir"),
        "session_folder": body.get("session_folder"),
        "bundle_destination": None,
        "copied_paths": changed,
        "changed_paths": changed,
        "mount_path": host_mount_path,
        "host_mount_path": host_mount_path,
        "mount_tree": mount_tree,
    }


def deploy_artifacts_via_sidecar(
    *,
    archive_path: str,
    host_mount_path: str,
    session_id: str = "",
    owner_id: str = "",
    bundle_dir: Optional[str] = None,
    session_title: Optional[str] = None,
    folder_timestamp: Optional[str] = None,
    session_created_at: Optional[Any] = None,
) -> Dict[str, Any]:
    """经 Local Sidecar 将归档写入宿主机挂载目录（仅同机 Docker 可达时可用；默认改走浏览器中继）。"""
    package = build_sidecar_deploy_package(
        archive_path=archive_path,
        host_mount_path=host_mount_path,
        session_id=session_id,
        owner_id=owner_id,
        bundle_dir=bundle_dir,
        session_title=session_title,
        folder_timestamp=folder_timestamp,
        session_created_at=session_created_at,
    )
    body = sidecar_silent_deploy(
        session_id=package["session_id"],
        session_title=package["session_title"],
        session_folder=package.get("session_folder"),
        folder_timestamp=package.get("folder_timestamp"),
        host_mount_path=host_mount_path,
        file_specs=package["files"],
    )
    return delivery_from_sidecar_body(body, host_mount_path=host_mount_path)
