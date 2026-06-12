# -*- coding: utf-8 -*-
"""三轨合一数据入库 API（语料守门人 + 打包 + 物理路由推送）。"""
from __future__ import annotations

import logging
import traceback
from typing import Any, Dict, List, Optional

from fastapi import APIRouter, Depends, HTTPException
from fastapi.responses import FileResponse
from pydantic import BaseModel, Field
from sqlalchemy.orm import Session as OrmSession

from gibh_agent.core.client_deploy_relay import resolve_staged_file
from gibh_agent.core.data_packager import build_artifacts_archive
from gibh_agent.core.deps import get_current_owner_id
from gibh_agent.core.hitl_resume import _extract_steps_details, _load_latest_agent_snapshot
from gibh_agent.core.hitl_session_registry import get_hitl_session
from gibh_agent.core.ingestion_deploy import (
    build_sidecar_deploy_package,
    deploy_artifacts_to_mount,
    deploy_artifacts_via_sidecar,
)
from gibh_agent.core.ingestion_mount_discovery import (
    discover_ingestion_mount_paths,
    get_default_ingestion_mount_path,
    validate_ingestion_mount_path,
)
from gibh_agent.core.ingestion_router import IngestionRouterError, deliver_archive
from gibh_agent.core.silent_local_corpus_deploy import (
    check_local_session_result_exists,
    resolve_local_project_mount_path,
)
from gibh_agent.core.storage.dual_path import notify_session_files_changed
from gibh_agent.core.storage.session_file_index import (
    build_session_files_inventory,
    extract_upload_paths_from_messages,
)
from gibh_agent.core.storage.mount_path_resolver import (
    resolve_container_writable_mount,
    resolve_host_workspace_path,
    resolve_ingestion_mount_target,
)
from gibh_agent.core.user_database_settings_store import get_database_mount_config, save_database_mount_config
from gibh_agent.db.connection import get_db_session
from gibh_agent.db.models import Message as MessageModel
from gibh_agent.db.models import Session as SessionModel
from gibh_agent.utils.ls_client import LabelStudioClient, LabelStudioClientError

logger = logging.getLogger(__name__)

router = APIRouter(prefix="/api/ingestion", tags=["Ingestion"])


class IngestionTriggerBody(BaseModel):
    session_id: Optional[str] = Field(default=None, description="关联会话 ID")
    skip_hitl: bool = Field(default=False, description="跳过 Label Studio 专家标注直接入库")
    mount_path: Optional[str] = Field(
        default=None,
        description="入库落盘路径（宿主机 Windows 路径或容器内绝对路径）",
    )
    workspace_path: Optional[str] = Field(
        default=None,
        description="本地项目工作区根路径（Sidecar 落盘；与 mount_path 二选一或互补）",
    )
    client_track: Optional[str] = Field(
        default=None,
        description="web | local_sidecar；local_sidecar 时由浏览器中继 Sidecar 落盘",
    )
    client_deploy_via_browser: bool = Field(
        default=True,
        description="local_sidecar 时是否由浏览器调用本机 Sidecar（远程 API 场景必选）",
    )
    persist_mount_path: bool = Field(
        default=False,
        description="是否将 mount_path 写入用户数据库挂载配置",
    )


def _resolve_expert_report(snap: Dict[str, Any]) -> str:
    ex = snap.get("execution_snapshot") if isinstance(snap.get("execution_snapshot"), dict) else {}
    md = (ex.get("expert_report_markdown") or snap.get("expert_report_markdown") or "").strip()
    if md:
        return md
    rep = snap.get("report") if isinstance(snap.get("report"), dict) else {}
    rd = rep.get("report_data") if isinstance(rep.get("report_data"), dict) else {}
    return str(rd.get("report") or rep.get("report") or "").strip()


def _resolve_corpus_bundle(pack: Dict[str, Any]) -> Dict[str, Any]:
    manifest = pack.get("manifest") if isinstance(pack.get("manifest"), dict) else {}
    return {
        "modality": manifest.get("corpus_modality"),
        "archive_dir": manifest.get("corpus_archive_dir"),
        "file_count": len(manifest.get("files") or []),
    }


def _resolve_session_title(db: OrmSession, session_id: str, owner_id: str) -> Optional[str]:
    meta = _resolve_session_meta(db, session_id, owner_id)
    return meta.get("title") if meta else None


def _resolve_session_meta(db: OrmSession, session_id: str, owner_id: str) -> Optional[Dict[str, Any]]:
    if not session_id:
        return None
    row = (
        db.query(SessionModel)
        .filter(SessionModel.id == session_id, SessionModel.owner_id == owner_id)
        .first()
    )
    if not row:
        return None
    title = str(getattr(row, "title", "") or "").strip()
    return {
        "title": title or None,
        "created_at": getattr(row, "created_at", None),
    }


def _ingestion_error_payload(
    *,
    exc: BaseException,
    mount_type: str,
    archive_path: Optional[str] = None,
    manifest: Optional[Dict[str, Any]] = None,
    corpus_bundle: Optional[Dict[str, Any]] = None,
) -> Dict[str, Any]:
    tb = traceback.format_exc()
    err_msg = getattr(exc, "message", None) or str(exc)
    err_details = getattr(exc, "details", None)
    return {
        "status": "error",
        "message": err_msg,
        "error_type": type(exc).__name__,
        "traceback": tb,
        "stack_trace": tb,
        "mount_type": mount_type,
        "archive_path": archive_path,
        "manifest": manifest,
        "corpus": corpus_bundle or {},
        "error_details": err_details,
    }


def _resolve_ingestion_mount_config(
    owner_id: str,
    *,
    mount_path_override: Optional[str] = None,
    workspace_path: Optional[str] = None,
    local_workspace_mounted: bool = False,
    persist: bool = False,
) -> tuple[Dict[str, Any], Optional[str]]:
    """
    合并用户设置、前端传入路径与系统探测默认路径。
    返回 (cfg, error_message)。
    """
    from gibh_agent.core.workspace_context import register_session_workspace

    cfg = get_database_mount_config(owner_id)

    if workspace_path:
        register_session_workspace(str(workspace_path).strip())

    effective_mount = str(mount_path_override or workspace_path or "").strip() or None
    mount_type, target_path, err = resolve_ingestion_mount_target(
        effective_mount,
        local_workspace_mounted=local_workspace_mounted or bool(workspace_path),
    )

    if err and not effective_mount:
        if not ((cfg.get("local_volume") or {}).get("mount_path")):
            default_path = get_default_ingestion_mount_path()
            if default_path:
                cfg = dict(cfg)
                cfg["mount_type"] = "local_volume"
                cfg["local_volume"] = {"mount_path": default_path}
                mount_type, target_path, err = "local_volume", default_path, None

    if err or not mount_type or not target_path:
        if err:
            return cfg, err
        return cfg, "尚未配置业务数据库挂载，且未能自动探测到可用路径。"

    cfg = dict(cfg)
    cfg["mount_type"] = mount_type
    if mount_type == "local_sidecar":
        cfg["local_volume"] = {"mount_path": target_path, "host_path": target_path}
    else:
        cfg["local_volume"] = {"mount_path": target_path}

    if persist and mount_type == "local_volume":
        try:
            save_database_mount_config(owner_id, cfg)
        except ValueError as exc:
            logger.warning("[Ingestion] persist mount_path failed: %s", exc)

    return cfg, None


@router.get("/discover-mount")
async def discover_mount(
    owner_id: str = Depends(get_current_owner_id),
) -> Dict[str, Any]:
    """
    探测 Docker 容器内可写挂载目录，供前端首次一键入库确认弹窗预填。
    """
    paths = discover_ingestion_mount_paths()
    default_path = paths[0]["path"] if paths else None
    return {
        "status": "success",
        "default_path": default_path,
        "paths": paths,
        "hint": "此为容器内路径，对应宿主机 docker-compose 挂载卷；请勿填写 Windows 盘符如 D:\\project",
        "owner_id_prefix": owner_id[:8] if owner_id else "",
    }


@router.get("/sessions/{session_id}/files")
async def list_session_files_for_ingestion(
    session_id: str,
    owner_id: str = Depends(get_current_owner_id),
    db: OrmSession = Depends(get_db_session),
) -> Dict[str, Any]:
    """返回当前会话 upload/result 文件清单（临时缓存 + 挂载永久目录）。"""
    session = db.query(SessionModel).filter(SessionModel.id == session_id).first()
    if not session:
        return {"status": "error", "message": "会话不存在"}
    if session.owner_id != owner_id:
        return {"status": "error", "message": "无权访问该会话"}

    cfg = get_database_mount_config(owner_id)
    mount_path = str((cfg.get("local_volume") or {}).get("mount_path") or "")
    messages = db.query(MessageModel).filter(MessageModel.session_id == session_id).all()
    inventory = build_session_files_inventory(
        session_id=session_id,
        session_title=getattr(session, "title", None),
        mount_path=mount_path or None,
        messages=messages,
    )
    return {"status": "success", **inventory}


@router.get("/sessions/{session_id}/local-result-check")
async def check_session_local_result(
    session_id: str,
    mount_path: Optional[str] = None,
    owner_id: str = Depends(get_current_owner_id),
    db: OrmSession = Depends(get_db_session),
) -> Dict[str, Any]:
    """检测会话结果是否已静默落盘至本地项目挂载目录的 result/ 下。"""
    session = db.query(SessionModel).filter(SessionModel.id == session_id).first()
    if not session:
        return {"status": "error", "message": "会话不存在"}
    if session.owner_id != owner_id:
        return {"status": "error", "message": "无权访问该会话"}

    resolved_mount = str(mount_path or "").strip() or resolve_local_project_mount_path() or ""
    if not resolved_mount:
        cfg = get_database_mount_config(owner_id)
        resolved_mount = str((cfg.get("local_volume") or {}).get("mount_path") or "").strip()
    if not resolved_mount:
        resolved_mount = resolve_host_workspace_path() or ""

    if not resolved_mount:
        return {
            "status": "success",
            "has_local_result": False,
            "message": "未配置本地项目挂载路径",
        }

    check = check_local_session_result_exists(
        mount_path=resolved_mount,
        session_id=session_id,
        session_title=getattr(session, "title", None),
        session_created_at=getattr(session, "created_at", None),
    )
    return {"status": "success", **check}


@router.post("/trigger")
async def trigger_ingestion(
    body: IngestionTriggerBody,
    owner_id: str = Depends(get_current_owner_id),
    db: OrmSession = Depends(get_db_session),
) -> Dict[str, Any]:
    """
    一键入库：Corpus Gatekeeper 盘点/代偿生成语料 → Data Packager 打归档包 → 物理落盘到挂载卷。
    """
    mount_type = "local_volume"
    try:
        local_mounted = (body.client_track or "").strip().lower() == "local_sidecar" or bool(body.workspace_path)
        cfg, cfg_err = _resolve_ingestion_mount_config(
            owner_id,
            mount_path_override=body.mount_path,
            workspace_path=body.workspace_path,
            local_workspace_mounted=local_mounted,
            persist=body.persist_mount_path,
        )
        if cfg_err:
            discovered = discover_ingestion_mount_paths()
            return {
                "status": "error",
                "message": cfg_err,
                "needs_settings": not body.mount_path and not discovered,
                "needs_mount_confirm": bool(discovered) and not body.mount_path,
                "discovered_paths": discovered,
                "default_path": discovered[0]["path"] if discovered else None,
            }

        mount_type = cfg.get("mount_type") or "local_volume"

        sid = str(body.session_id or "").strip()
        session_meta = _resolve_session_meta(db, sid, owner_id) if sid else None
        session_title = session_meta.get("title") if session_meta else None
        session_created_at = session_meta.get("created_at") if session_meta else None
        expert_md = ""
        steps_details: List[Dict[str, Any]] = []
        hitl_annotations = None
        hitl_meta: Dict[str, Any] = {}
        output_dir = None
        snap: Dict[str, Any] = {}
        upload_paths: List[str] = []

        if sid:
            messages = db.query(MessageModel).filter(MessageModel.session_id == sid).all()
            upload_paths = extract_upload_paths_from_messages(messages)
            _msg, snap = _load_latest_agent_snapshot(db, sid, owner_id)
            if snap:
                expert_md = _resolve_expert_report(snap)
                steps_details = _extract_steps_details(snap)
                hitl = snap.get("hitl") if isinstance(snap.get("hitl"), dict) else {}
                reg = get_hitl_session(sid) or {}
                hitl_meta = {**reg, **hitl}
                output_dir = reg.get("output_dir") or hitl.get("output_dir")
                pid = hitl.get("project_id") or reg.get("project_id")
                if pid and not body.skip_hitl:
                    try:
                        hitl_annotations = LabelStudioClient().export_annotations(int(pid))
                    except LabelStudioClientError as exc:
                        logger.warning("[Ingestion] LS export skipped: %s", exc)

        pack = build_artifacts_archive(
            session_id=sid or f"anonymous-{owner_id[:8]}",
            owner_id=owner_id,
            expert_report_markdown=expert_md,
            steps_details=steps_details,
            hitl_annotations=hitl_annotations,
            hitl_meta=hitl_meta,
            output_dir=output_dir,
            skip_hitl=body.skip_hitl,
        )
        corpus_bundle = _resolve_corpus_bundle(pack)

        archive_path = pack.get("archive_path")
        if pack.get("status") != "success" or not archive_path:
            return {
                "status": "error",
                "message": pack.get("message") or "归档打包失败",
                "mount_type": mount_type,
                "corpus": corpus_bundle,
            }

        if mount_type == "local_sidecar":
            host_mount = str(
                (cfg.get("local_volume") or {}).get("host_path")
                or (cfg.get("local_volume") or {}).get("mount_path")
                or ""
            )
            use_browser_relay = body.client_deploy_via_browser or (
                (body.client_track or "").strip().lower() == "local_sidecar"
            )
            if use_browser_relay:
                deploy_package = build_sidecar_deploy_package(
                    archive_path=str(archive_path),
                    host_mount_path=host_mount,
                    session_id=sid,
                    owner_id=owner_id,
                    bundle_dir=pack.get("bundle_dir"),
                    session_title=session_title,
                    session_created_at=session_created_at,
                )
                files_notify = notify_session_files_changed(
                    session_id=sid,
                    changed_paths=[],
                )
                files_notify["client_deploy_package"] = deploy_package
                return {
                    "status": "client_deploy",
                    "message": "归档已就绪，请由本机 Sidecar 落盘",
                    "job_id": f"ingest-{owner_id[:8]}-{sid[:8] if sid else 'na'}",
                    "mount_type": mount_type,
                    "auto_ingestion_enabled": bool(cfg.get("is_auto_ingestion_enabled")),
                    "skip_hitl": body.skip_hitl,
                    "archive_path": archive_path,
                    "manifest": pack.get("manifest"),
                    "corpus": corpus_bundle,
                    "deploy_package": deploy_package,
                    "session_title": session_title,
                    "files_notification": files_notify,
                }
            delivery = deploy_artifacts_via_sidecar(
                archive_path=str(archive_path),
                host_mount_path=host_mount,
                session_id=sid,
                owner_id=owner_id,
                bundle_dir=pack.get("bundle_dir"),
                session_title=session_title,
                session_created_at=session_created_at,
            )
        elif mount_type == "local_volume":
            mount_path = str((cfg.get("local_volume") or {}).get("mount_path") or "")
            delivery = deploy_artifacts_to_mount(
                archive_path=str(archive_path),
                mount_path=mount_path,
                session_id=sid,
                owner_id=owner_id,
                bundle_dir=pack.get("bundle_dir"),
                session_title=session_title,
                upload_paths=upload_paths,
                session_created_at=session_created_at,
            )
        else:
            delivery = deliver_archive(
                str(archive_path),
                cfg,
                session_id=sid,
                owner_id=owner_id,
            )

        files_notify = notify_session_files_changed(
            session_id=sid,
            changed_paths=delivery.get("changed_paths") or delivery.get("copied_paths") or [],
            mount_tree=delivery.get("mount_tree"),
        )

        logger.info(
            "[Ingestion] delivered owner=%s session=%s dest=%s skip_hitl=%s modality=%s",
            owner_id,
            sid,
            delivery.get("destination"),
            body.skip_hitl,
            corpus_bundle.get("modality"),
        )
        dest_display = (
            delivery.get("destination_dir")
            or delivery.get("destination")
            or delivery.get("host_mount_path")
            or archive_path
        )
        return {
            "status": "success",
            "message": f"入库成功：{dest_display}",
            "job_id": f"ingest-{owner_id[:8]}-{sid[:8] if sid else 'na'}",
            "mount_type": mount_type,
            "auto_ingestion_enabled": bool(cfg.get("is_auto_ingestion_enabled")),
            "skip_hitl": body.skip_hitl,
            "archive_path": archive_path,
            "manifest": pack.get("manifest"),
            "corpus": corpus_bundle,
            "delivery": delivery,
            "session_title": session_title,
            "session_folder": delivery.get("session_folder"),
            "files_notification": files_notify,
        }
    except (IngestionRouterError, ValueError) as exc:
        logger.warning(
            "[Ingestion] deliver failed owner=%s session=%s mount=%s: %s",
            owner_id,
            body.session_id,
            mount_type,
            exc,
        )
        return _ingestion_error_payload(exc=exc, mount_type=mount_type)
    except Exception as exc:
        logger.exception(
            "[Ingestion] fatal owner=%s session=%s: %s",
            owner_id,
            body.session_id,
            exc,
        )
        return _ingestion_error_payload(exc=exc, mount_type=mount_type)


@router.get("/deploy-staging/{token}/{stored_name}")
async def download_deploy_staging_file(
    token: str,
    stored_name: str,
    owner_id: str = Depends(get_current_owner_id),
) -> FileResponse:
    """
    浏览器中继落盘：Sidecar 流式拉取大文件（token 绑定 owner，TTL 过期自动失效）。
    """
    try:
        path = resolve_staged_file(token, stored_name, owner_id=owner_id)
    except PermissionError as exc:
        raise HTTPException(status_code=403, detail=str(exc)) from exc
    except FileNotFoundError as exc:
        raise HTTPException(status_code=404, detail=str(exc)) from exc
    return FileResponse(
        path=str(path),
        media_type="application/octet-stream",
        filename=path.name,
    )
