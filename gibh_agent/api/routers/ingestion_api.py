# -*- coding: utf-8 -*-
"""三轨合一数据入库 API（打包 + 物理路由推送）。"""
from __future__ import annotations

import logging
from typing import Any, Dict, List, Optional

from fastapi import APIRouter, Depends
from pydantic import BaseModel, Field
from sqlalchemy.orm import Session as OrmSession

from gibh_agent.core.data_packager import build_artifacts_archive
from gibh_agent.core.ingestion_router import IngestionRouterError, deliver_archive
from gibh_agent.core.deps import get_current_owner_id
from gibh_agent.core.hitl_resume import _extract_steps_details, _load_latest_agent_snapshot
from gibh_agent.core.hitl_session_registry import get_hitl_session
from gibh_agent.core.user_database_settings_store import get_database_mount_config
from gibh_agent.db.connection import get_db_session
from gibh_agent.utils.ls_client import LabelStudioClient, LabelStudioClientError

logger = logging.getLogger(__name__)

router = APIRouter(prefix="/api/ingestion", tags=["Ingestion"])


class IngestionTriggerBody(BaseModel):
    session_id: Optional[str] = Field(default=None, description="关联会话 ID")
    skip_hitl: bool = Field(default=False, description="跳过 Label Studio 专家标注直接入库")


def _resolve_expert_report(snap: Dict[str, Any]) -> str:
    ex = snap.get("execution_snapshot") if isinstance(snap.get("execution_snapshot"), dict) else {}
    md = (ex.get("expert_report_markdown") or snap.get("expert_report_markdown") or "").strip()
    if md:
        return md
    rep = snap.get("report") if isinstance(snap.get("report"), dict) else {}
    rd = rep.get("report_data") if isinstance(rep.get("report_data"), dict) else {}
    return str(rd.get("report") or rep.get("report") or "").strip()


@router.post("/trigger")
def trigger_ingestion(
    body: IngestionTriggerBody,
    owner_id: str = Depends(get_current_owner_id),
    db: OrmSession = Depends(get_db_session),
) -> Dict[str, Any]:
    """
    一键入库：先经 Data Packager 打归档包，再经 IngestionRouter 推送到用户配置目标。
    """
    cfg = get_database_mount_config(owner_id)
    mount_type = cfg.get("mount_type") or "local_volume"
    has_config = False
    if mount_type == "local_volume":
        has_config = bool((cfg.get("local_volume") or {}).get("mount_path"))
    elif mount_type == "hpc_slurm":
        has_config = bool((cfg.get("hpc_slurm") or {}).get("host"))
    elif mount_type == "api_url":
        has_config = bool((cfg.get("api_url") or {}).get("endpoint"))

    if not has_config:
        return {
            "status": "error",
            "message": "尚未配置业务数据库挂载，请先在设置中完成关联。",
            "needs_settings": True,
        }

    sid = str(body.session_id or "").strip()
    expert_md = ""
    steps_details: List[Dict[str, Any]] = []
    hitl_annotations = None
    hitl_meta: Dict[str, Any] = {}
    output_dir = None

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

    archive_path = pack.get("archive_path")
    if pack.get("status") != "success" or not archive_path:
        return {
            "status": "error",
            "message": pack.get("message") or "归档打包失败",
            "mount_type": mount_type,
        }

    try:
        delivery = deliver_archive(
            str(archive_path),
            cfg,
            session_id=sid,
            owner_id=owner_id,
        )
    except IngestionRouterError as exc:
        logger.warning(
            "[Ingestion] deliver failed owner=%s session=%s mount=%s: %s",
            owner_id,
            sid,
            exc.mount_type,
            exc.message,
        )
        return {
            "status": "error",
            "message": exc.message,
            "mount_type": mount_type,
            "archive_path": archive_path,
            "manifest": pack.get("manifest"),
            "error_details": exc.details,
        }

    logger.info(
        "[Ingestion] delivered owner=%s session=%s dest=%s skip_hitl=%s",
        owner_id,
        sid,
        delivery.get("destination"),
        body.skip_hitl,
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
        "delivery": delivery,
    }
