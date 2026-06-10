# -*- coding: utf-8 -*-
"""入库物理落盘：tar.gz + 解包目录树复制到用户挂载卷。"""
from __future__ import annotations

import logging
import shutil
from pathlib import Path
from typing import Any, Dict, List, Optional

from gibh_agent.core.ingestion_mount_discovery import validate_ingestion_mount_path

logger = logging.getLogger(__name__)


def deploy_artifacts_to_mount(
    *,
    archive_path: str,
    mount_path: str,
    session_id: str = "",
    owner_id: str = "",
    bundle_dir: Optional[str] = None,
) -> Dict[str, Any]:
    """
    将 tar.gz 归档与打包工作目录（含 corpus_archive 等）复制到挂载目录。

    目标结构::
        {mount_path}/gibh_ingest/{session_id}/bundle_*.tar.gz
        {mount_path}/gibh_ingest/{session_id}/bundle_*/...
    """
    validated = validate_ingestion_mount_path(mount_path)
    if not validated:
        raise ValueError(f"挂载路径不可用或不可写: {mount_path}")

    archive = Path(str(archive_path)).expanduser()
    if not archive.is_file():
        raise ValueError(f"归档文件不存在: {archive}")

    base = Path(validated)
    sub = str(session_id or "").strip() or (f"owner-{owner_id[:8]}" if owner_id else "anonymous")
    dest_dir = base / "gibh_ingest" / sub
    dest_dir.mkdir(parents=True, exist_ok=True)

    dest_file = dest_dir / archive.name
    shutil.copy2(str(archive), str(dest_file))

    copied: List[str] = [str(dest_file.resolve())]
    bundle_dest: Optional[Path] = None

    if bundle_dir:
        src_bundle = Path(str(bundle_dir)).expanduser()
        if src_bundle.is_dir():
            bundle_dest = dest_dir / src_bundle.name
            if bundle_dest.exists():
                shutil.rmtree(bundle_dest)
            shutil.copytree(str(src_bundle), str(bundle_dest))
            copied.append(str(bundle_dest.resolve()))
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
        "destination_dir": str(dest_dir.resolve()),
        "bundle_destination": str(bundle_dest.resolve()) if bundle_dest else None,
        "copied_paths": copied,
        "mount_path": validated,
    }
