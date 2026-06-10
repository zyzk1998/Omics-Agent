# -*- coding: utf-8 -*-
"""三轨合一数据入库 API（语料守门人 + 打包 + 物理路由推送）。"""
from __future__ import annotations

import logging
from typing import Any, Dict, List, Optional

from fastapi import APIRouter, Depends
from pydantic import BaseModel, Field
from sqlalchemy.orm import Session as OrmSession

from gibh_agent.core.corpus_gatekeeper import ensure_corpus_ready_for_ingestion
from gibh_agent.core.data_packager import build_artifacts_archive
from gibh_agent.core.ingestion_deploy import deploy_artifacts_to_mount
from gibh_agent.core.ingestion_mount_discovery import (
    discover_ingestion_mount_paths,
    get_default_ingestion_mount_path,
    validate_ingestion_mount_path,
)
from gibh_agent.core.ingestion_router import IngestionRouterError, deliver_archive
from gibh_agent.core.deps import get_current_owner_id
from gibh_agent.core.hitl_resume import _extract_steps_details, _load_latest_agent_snapshot
from gibh_agent.core.hitl_session_registry import get_hitl_session
from gibh_agent.core.user_database_settings_store import get_database_mount_config, save_database_mount_config
from gibh_agent.db.connection import get_db_session
from gibh_agent.utils.ls_client import LabelStudioClient, LabelStudioClientError

logger = logging.getLogger(__name__)

router = APIRouter(prefix="/api/ingestion", tags=["Ingestion"])


class IngestionTriggerBody(BaseModel):
    session_id: Optional[str] = Field(default=None, description="关联会话 ID")
    skip_hitl: bool = Field(default=False, description="跳过 Label Studio 专家标注直接入库")
    mount_path: Optional[str] = Field(
        default=None,
        description="容器内挂载绝对路径（前端 localStorage 记忆；优先于设置面板）",
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


def _resolve_ingestion_mount_config(
    owner_id: str,
    *,
    mount_path_override: Optional[str] = None,
    persist: bool = False,
) -> tuple[Dict[str, Any], Optional[str]]:
    """
    合并用户设置、前端传入路径与系统探测默认路径。
    返回 (cfg, error_message)。
    """
    cfg = get_database_mount_config(owner_id)
    validated: Optional[str] = None

    if mount_path_override:
        validated = validate_ingestion_mount_path(mount_path_override)
        if not validated:
            return cfg, f"挂载路径不可用（须为容器内存在且可写的目录）: {mount_path_override}"
        cfg = dict(cfg)
        cfg["mount_type"] = "local_volume"
        cfg["local_volume"] = {"mount_path": validated}
        if persist:
            try:
                save_database_mount_config(owner_id, cfg)
            except ValueError as exc:
                logger.warning("[Ingestion] persist mount_path failed: %s", exc)
    elif not ((cfg.get("local_volume") or {}).get("mount_path")):
        default_path = get_default_ingestion_mount_path()
        if default_path:
            cfg = dict(cfg)
            cfg["mount_type"] = "local_volume"
            cfg["local_volume"] = {"mount_path": default_path}

    mount_type = cfg.get("mount_type") or "local_volume"
    has_config = False
    if mount_type == "local_volume":
        has_config = bool((cfg.get("local_volume") or {}).get("mount_path"))
    elif mount_type == "hpc_slurm":
        has_config = bool((cfg.get("hpc_slurm") or {}).get("host"))
    elif mount_type == "api_url":
        has_config = bool((cfg.get("api_url") or {}).get("endpoint"))

    if not has_config:
        return cfg, "尚未配置业务数据库挂载，且未能自动探测到可用容器路径。"
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


@router.post("/trigger")
async def trigger_ingestion(
    body: IngestionTriggerBody,
    owner_id: str = Depends(get_current_owner_id),
    db: OrmSession = Depends(get_db_session),
) -> Dict[str, Any]:
    """
    一键入库：Corpus Gatekeeper 盘点/代偿生成语料 → Data Packager 打归档包 → 物理落盘到挂载卷。
    """
    cfg, cfg_err = _resolve_ingestion_mount_config(
        owner_id,
        mount_path_override=body.mount_path,
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
    expert_md = ""
    steps_details: List[Dict[str, Any]] = []
    hitl_annotations = None
    hitl_meta: Dict[str, Any] = {}
    output_dir = None
    snap: Dict[str, Any] = {}
    corpus_bundle: Dict[str, Any] = {}

    if sid:
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

    if sid and snap:
        from server import agent

        corpus_bundle = await ensure_corpus_ready_for_ingestion(
            session_id=sid,
            snap=snap,
            agent_root=agent,
            expert_report_md=expert_md,
            steps_details=steps_details,
            hitl_annotations=hitl_annotations,
            hitl_meta=hitl_meta,
            output_dir=output_dir,
            skip_hitl=body.skip_hitl,
        )
        if corpus_bundle.get("status") != "success":
            logger.warning(
                "[Ingestion] corpus gatekeeper partial session=%s: %s",
                sid,
                corpus_bundle.get("message"),
            )

    pack = build_artifacts_archive(
        session_id=sid or f"anonymous-{owner_id[:8]}",
        owner_id=owner_id,
        expert_report_markdown=expert_md,
        steps_details=steps_details,
        hitl_annotations=hitl_annotations,
        hitl_meta=hitl_meta,
        output_dir=output_dir,
        skip_hitl=body.skip_hitl,
        corpus_archive_dir=corpus_bundle.get("archive_dir"),
        corpus_modality=corpus_bundle.get("modality"),
    )

    archive_path = pack.get("archive_path")
    if pack.get("status") != "success" or not archive_path:
        return {
            "status": "error",
            "message": pack.get("message") or "归档打包失败",
            "mount_type": mount_type,
        }

    try:
        if mount_type == "local_volume":
            mount_path = str((cfg.get("local_volume") or {}).get("mount_path") or "")
            delivery = deploy_artifacts_to_mount(
                archive_path=str(archive_path),
                mount_path=mount_path,
                session_id=sid,
                owner_id=owner_id,
                bundle_dir=pack.get("bundle_dir"),
            )
        else:
            delivery = deliver_archive(
                str(archive_path),
                cfg,
                session_id=sid,
                owner_id=owner_id,
            )
    except (IngestionRouterError, ValueError) as exc:
        err_msg = getattr(exc, "message", None) or str(exc)
        err_details = getattr(exc, "details", None)
        logger.warning(
            "[Ingestion] deliver failed owner=%s session=%s mount=%s: %s",
            owner_id,
            sid,
            mount_type,
            err_msg,
        )
        return {
            "status": "error",
            "message": err_msg,
            "mount_type": mount_type,
            "archive_path": archive_path,
            "manifest": pack.get("manifest"),
            "corpus": corpus_bundle,
            "error_details": err_details,
        }

    logger.info(
        "[Ingestion] delivered owner=%s session=%s dest=%s skip_hitl=%s modality=%s",
        owner_id,
        sid,
        delivery.get("destination"),
        body.skip_hitl,
        corpus_bundle.get("modality"),
    )
    return {
        "status": "success",
        "message": f"入库成功：{delivery.get('destination', archive_path)}",
        "job_id": f"ingest-{owner_id[:8]}-{sid[:8] if sid else 'na'}",
        "mount_type": mount_type,
        "auto_ingestion_enabled": bool(cfg.get("is_auto_ingestion_enabled")),
        "skip_hitl": body.skip_hitl,
        "archive_path": archive_path,
        "manifest": pack.get("manifest"),
        "corpus": corpus_bundle,
        "delivery": delivery,
    }
