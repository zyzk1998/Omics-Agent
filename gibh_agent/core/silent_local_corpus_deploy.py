# -*- coding: utf-8 -*-
"""本地项目挂载静默落盘：语料/报告复制至 {mount}/{session}/result/。"""
from __future__ import annotations

import logging
import shutil
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

from gibh_agent.core.client_deploy_relay import (
    build_client_relay_deploy_package,
    collect_corpus_bundle_entries,
)
from gibh_agent.core.local_sidecar_client import (
    sidecar_check_local_result,
)
from gibh_agent.core.storage.mount_path_resolver import (
    is_ephemeral_container_path,
    is_host_side_mount_path,
    resolve_container_writable_mount,
    resolve_host_workspace_path,
)
from gibh_agent.core.storage.session_paths import resolve_session_mount_tree
from gibh_agent.core.workspace_context import get_workspace_context
from gibh_agent.utils.path_resolver import is_windows_abs_path

logger = logging.getLogger(__name__)

SFT_CORPUS_GLOB = "sft_corpus*.json"


def resolve_local_project_mount_path(
    workspace_context: Optional[Dict[str, Any]] = None,
) -> Optional[str]:
    """返回 API 容器内可写的本地项目挂载根（非 /app/results 临时区）。"""
    return resolve_container_writable_mount(None, workspace_context=workspace_context)


def resolve_deploy_strategy(
    workspace_context: Optional[Dict[str, Any]] = None,
    *,
    local_workspace_mounted: bool = False,
) -> Tuple[Optional[str], Optional[str]]:
    """
    返回 (mode, path)：mode 为 ``container`` 或 ``sidecar``。
    Windows/不可容器直达的路径走 Sidecar；同进程可写目录走 container。
    """
    ctx = workspace_context if workspace_context is not None else get_workspace_context()
    host = resolve_host_workspace_path(workspace_context=ctx)
    host_is_sidecar = bool(host and is_windows_abs_path(host))

    if host_is_sidecar and (local_workspace_mounted or bool(ctx.get("host_workspace_path"))):
        return "sidecar", host

    container = resolve_container_writable_mount(None, workspace_context=ctx, allow_default=not local_workspace_mounted)
    if container:
        return "container", container

    if host:
        return "sidecar", host
    return None, None


def find_sft_corpus_in_result_dir(result_dir: Path) -> List[Path]:
    if not result_dir.is_dir():
        return []
    return sorted(result_dir.glob(SFT_CORPUS_GLOB))


def check_local_session_result_exists(
    *,
    mount_path: str,
    session_id: str,
    session_title: Optional[str] = None,
    session_created_at: Optional[Any] = None,
    folder_timestamp: Optional[str] = None,
) -> Dict[str, Any]:
    """检测挂载树下 result/ 是否已有 sft_corpus 语料文件。"""
    container = resolve_container_writable_mount(mount_path)
    if container:
        tree = resolve_session_mount_tree(
            container,
            session_title=session_title,
            session_id=session_id,
            session_created_at=session_created_at,
            folder_timestamp=folder_timestamp,
        )
        result_dir = Path(tree["result_dir"])
        corpus_files = find_sft_corpus_in_result_dir(result_dir)
        return {
            "has_local_result": bool(corpus_files),
            "mount_path": container,
            "host_mount_path": resolve_host_workspace_path(mount_path) or mount_path,
            "result_dir": tree["result_dir"],
            "session_root": tree["session_root"],
            "session_folder": tree["session_folder"],
            "sft_corpus_paths": [str(p.resolve()) for p in corpus_files],
            "mount_tree": tree,
            "deploy_mode": "container",
        }

    host = resolve_host_workspace_path(mount_path)
    if host:
        tree = resolve_session_mount_tree(
            host,
            session_title=session_title,
            session_id=session_id,
            session_created_at=session_created_at,
            folder_timestamp=folder_timestamp,
        )
        check = sidecar_check_local_result(
            host_mount_path=host,
            session_id=session_id,
            session_title=session_title,
            session_folder=tree["session_folder"],
        )
        check.setdefault("host_mount_path", host)
        check.setdefault("deploy_mode", "sidecar")
        return check

    return {
        "has_local_result": False,
        "mount_path": mount_path,
        "message": "挂载路径不可用或不可写",
    }


def _collect_bundle_file_specs(bundle: Dict[str, Any]) -> List[Dict[str, str]]:
    from gibh_agent.core.local_sidecar_client import read_file_as_b64

    specs: List[Dict[str, str]] = []
    for entry in collect_corpus_bundle_entries(bundle):
        src = entry.source_path.expanduser()
        if src.is_file():
            specs.append({"dest_name": entry.dest_name, "content_b64": read_file_as_b64(src)})
    return specs


def _build_client_relay_delivery(
    *,
    host_mount_path: str,
    session_id: str,
    session_title: Optional[str],
    bundle: Dict[str, Any],
    owner_id: str = "",
    session_created_at: Optional[Any] = None,
    folder_timestamp: Optional[str] = None,
) -> Optional[Dict[str, Any]]:
    """构建浏览器 → Local Sidecar 中继落盘载荷（小文件 inline，大文件 URL）。"""
    entries = collect_corpus_bundle_entries(bundle)
    if not entries:
        return None
    sid = str(session_id or "").strip() or "unknown"
    package = build_client_relay_deploy_package(
        owner_id=owner_id or sid,
        session_id=sid,
        session_title=session_title,
        host_mount_path=host_mount_path,
        entries=entries,
        folder_timestamp=folder_timestamp,
        session_created_at=session_created_at,
    )
    return {
        "strategy": "client_sidecar_relay",
        "deploy_mode": "client_relay",
        "client_deploy_package": package,
        "mount_path": package.get("host_mount_path"),
        "host_mount_path": package.get("host_mount_path"),
        "changed_paths": [],
        "copied_paths": [],
    }


def _deploy_to_container(
    *,
    mount_path: str,
    session_id: str,
    session_title: Optional[str],
    bundle: Dict[str, Any],
    session_created_at: Optional[Any] = None,
    folder_timestamp: Optional[str] = None,
) -> Optional[Dict[str, Any]]:
    if is_ephemeral_container_path(mount_path):
        logger.warning("[SilentDeploy] reject ephemeral mount: %s", mount_path)
        return None

    sid = str(session_id or "").strip() or "unknown"
    tree = resolve_session_mount_tree(
        mount_path,
        session_title=session_title,
        session_id=sid,
        session_created_at=session_created_at,
        folder_timestamp=folder_timestamp,
    )
    result_dir = Path(tree["result_dir"])
    result_dir.mkdir(parents=True, exist_ok=True)

    changed: List[str] = []

    hitl_path = str(bundle.get("corpus_hitl_path") or "").strip()
    if hitl_path:
        src = Path(hitl_path).expanduser()
        if src.is_file():
            dest = result_dir / src.name
            shutil.copy2(str(src), str(dest))
            changed.append(str(dest.resolve()))
            canonical = result_dir / "sft_corpus.json"
            if canonical.name != dest.name:
                shutil.copy2(str(src), str(canonical))
                changed.append(str(canonical.resolve()))

    archive_dir = str(bundle.get("archive_dir") or "").strip()
    if archive_dir:
        src_dir = Path(archive_dir).expanduser()
        if src_dir.is_dir():
            dest_dir = result_dir / src_dir.name
            if dest_dir.exists():
                shutil.rmtree(dest_dir)
            shutil.copytree(str(src_dir), str(dest_dir))
            changed.append(str(dest_dir.resolve()))

    dataset_path = str(bundle.get("dataset_path") or "").strip()
    if dataset_path:
        src_ds = Path(dataset_path).expanduser()
        if src_ds.is_file():
            dest_ds = result_dir / src_ds.name
            shutil.copy2(str(src_ds), str(dest_ds))
            changed.append(str(dest_ds.resolve()))

    if not changed:
        return None

    return {
        "strategy": "local_silent",
        "deploy_mode": "container",
        "destination_dir": str(result_dir.resolve()),
        "session_root": tree["session_root"],
        "session_folder": tree["session_folder"],
        "mount_path": mount_path,
        "mount_tree": tree,
        "changed_paths": changed,
        "copied_paths": changed,
    }



def silent_deploy_corpus_bundle_to_local_mount(
    *,
    session_id: str,
    session_title: Optional[str],
    bundle: Dict[str, Any],
    workspace_context: Optional[Dict[str, Any]] = None,
    local_workspace_mounted: bool = False,
    owner_id: str = "",
    session_created_at: Optional[Any] = None,
    folder_timestamp: Optional[str] = None,
) -> Optional[Dict[str, Any]]:
    """
    当用户已挂载本地项目目录时，将语料 JSON 与归档目录静默复制至
    ``{mount}/{session_folder}/result/``（容器卷或经 Sidecar 写入宿主机）。
    """
    ctx = workspace_context if workspace_context is not None else get_workspace_context()
    if not local_workspace_mounted and not ctx.get("host_workspace_path") and not ctx.get("workspace_path"):
        return None

    host = resolve_host_workspace_path(workspace_context=ctx)
    if local_workspace_mounted and not host and not resolve_container_writable_mount(None, workspace_context=ctx, allow_default=False):
        logger.warning("[SilentDeploy] local workspace mounted but no resolvable host/container path")
        return None

    mode, target = resolve_deploy_strategy(workspace_context=ctx, local_workspace_mounted=local_workspace_mounted)
    if mode == "container" and target:
        return _deploy_to_container(
            mount_path=target,
            session_id=session_id,
            session_title=session_title,
            bundle=bundle,
            session_created_at=session_created_at,
            folder_timestamp=folder_timestamp,
        )
    if mode == "sidecar" and target:
        return _build_client_relay_delivery(
            host_mount_path=target,
            session_id=session_id,
            session_title=session_title,
            bundle=bundle,
            owner_id=owner_id,
            session_created_at=session_created_at,
            folder_timestamp=folder_timestamp,
        )

    raw = str(ctx.get("host_workspace_path") or ctx.get("workspace_path") or "").strip()
    if raw and (is_windows_abs_path(raw) or is_host_side_mount_path(raw)):
        return _build_client_relay_delivery(
            host_mount_path=raw,
            session_id=session_id,
            session_title=session_title,
            bundle=bundle,
            owner_id=owner_id,
            session_created_at=session_created_at,
            folder_timestamp=folder_timestamp,
        )
    logger.debug("[SilentDeploy] skip: no writable mount or sidecar path")
    return None
