# -*- coding: utf-8 -*-
"""HITL Webhook 与 resume_from_hitl API。"""
from __future__ import annotations

import logging
import os
from typing import Any, Dict, Optional

from fastapi import APIRouter, Depends, HTTPException, Request
from fastapi.responses import JSONResponse, StreamingResponse
from pydantic import BaseModel, Field
from sqlalchemy.orm import Session as OrmSession

from gibh_agent.core.deps import get_current_owner_id
from gibh_agent.core.hitl_session_registry import resolve_session_by_project
from gibh_agent.core.session_runtime import (
    SESSION_COMPLETED,
    SESSION_RUNNING,
    SESSION_WAITING_FOR_HITL,
    normalize_session_status,
)
from gibh_agent.db.connection import get_db_session
from gibh_agent.db.models import Session as SessionModel

logger = logging.getLogger(__name__)

router = APIRouter(tags=["HITL"])

_WEBHOOK_SECRET = os.getenv("LABEL_STUDIO_WEBHOOK_SECRET", "").strip()


class ResumeFromHitlBody(BaseModel):
    project_id: Optional[int] = Field(default=None, description="Label Studio 项目 ID（可选，默认从快照读取）")
    trigger: str = Field(default="frontend", description="唤醒来源：frontend | webhook | manual_confirm")
    skip: bool = Field(default=False, description="用户主动跳过专家复核（仅关闭入口，不生成最终版）")


def _snapshot_allows_soft_hitl_resume(db: OrmSession, session_id: str, owner_id: str) -> bool:
    """软 HITL：初稿报告已生成但 hitl_pending 仍为真时允许 resume。"""
    from gibh_agent.core.hitl_resume import _load_latest_agent_snapshot

    _msg, snap = _load_latest_agent_snapshot(db, session_id, owner_id)
    if not snap:
        return False
    if snap.get("hitl_resumed"):
        return False
    return bool(snap.get("hitl_pending") or snap.get("hitl"))


def _session_allows_hitl_resume(
    db: OrmSession,
    session: SessionModel,
    owner_id: str,
) -> bool:
    st = normalize_session_status(getattr(session, "status", None))
    if st == SESSION_WAITING_FOR_HITL:
        return True
    if st in (SESSION_COMPLETED, SESSION_RUNNING):
        return _snapshot_allows_soft_hitl_resume(db, session.id, owner_id)
    return False


def _verify_webhook(request: Request) -> bool:
    if not _WEBHOOK_SECRET:
        return True
    token = request.headers.get("X-Label-Studio-Token") or request.headers.get("Authorization", "")
    token = token.replace("Bearer ", "").replace("Token ", "").strip()
    return token == _WEBHOOK_SECRET


def _extract_project_id(payload: Dict[str, Any]) -> Optional[int]:
    for key in ("project", "project_id"):
        val = payload.get(key)
        if isinstance(val, dict) and val.get("id") is not None:
            return int(val["id"])
        if val is not None and not isinstance(val, dict):
            try:
                return int(val)
            except (TypeError, ValueError):
                pass
    task = payload.get("task")
    if isinstance(task, dict) and task.get("project") is not None:
        try:
            return int(task["project"])
        except (TypeError, ValueError):
            pass
    return None


def _is_annotation_submit_event(payload: Dict[str, Any]) -> bool:
    action = str(payload.get("action") or payload.get("event") or "").upper()
    if not action:
        return False
    triggers = (
        "ANNOTATION",
        "ANNOTATIONS",
        "TASK_COMPLETED",
        "TASK_FINISHED",
        "SUBMIT",
        "COMPLETED",
        "EXPORT",
    )
    return any(t in action for t in triggers)


@router.post("/api/hitl/ls-session-bridge")
async def bridge_label_studio_session():
    """
    将 api-server 侧已 bootstrap 的 Label Studio Django session 写入浏览器 Cookie，
    使同源 /label-studio/ iframe 免二次登录即可打开项目页。
    """
    from gibh_agent.utils.ls_client import DEFAULT_LS_BASE_URL, LabelStudioClientError, _auto_bootstrap_session

    base_url = os.getenv("LABEL_STUDIO_URL", DEFAULT_LS_BASE_URL).rstrip("/")
    try:
        session = _auto_bootstrap_session(base_url, timeout_sec=15.0)
    except LabelStudioClientError as exc:
        raise HTTPException(status_code=503, detail=str(exc)) from exc

    resp = JSONResponse(
        {
            "status": "success",
            "message": "Label Studio 会话已桥接到浏览器",
        }
    )
    for key in ("sessionid", "csrftoken"):
        val = session.cookies.get(key)
        if not val:
            continue
        resp.set_cookie(
            key=key,
            value=val,
            path="/",
            httponly=(key == "sessionid"),
            samesite="lax",
        )
    if not session.cookies.get("sessionid"):
        raise HTTPException(status_code=503, detail="Label Studio sessionid 未获取，无法桥接登录态")
    return resp


@router.post("/api/hitl/webhook")
async def label_studio_webhook(
    request: Request,
    db: OrmSession = Depends(get_db_session),
):
    """
    Label Studio Webhook 入口。

    在 LS 项目 Settings → Webhooks 中配置：
    URL: http://<api-host>/api/hitl/webhook
    事件：ANNOTATION_CREATED / ANNOTATIONS_CREATED / TASK_COMPLETED 等
    可选 Header: X-Label-Studio-Token: <LABEL_STUDIO_WEBHOOK_SECRET>
    """
    if not _verify_webhook(request):
        raise HTTPException(status_code=403, detail="Webhook 鉴权失败")

    try:
        payload = await request.json()
    except Exception:
        payload = {}

    if not isinstance(payload, dict):
        payload = {}

    if not _is_annotation_submit_event(payload):
        return {"status": "ignored", "message": "非标注提交类事件，已忽略", "action": payload.get("action")}

    project_id = _extract_project_id(payload)
    if project_id is None:
        return {"status": "error", "message": "无法解析 project_id"}

    session_id = resolve_session_by_project(project_id)
    if not session_id:
        return {"status": "error", "message": f"未找到 project_id={project_id} 关联的 session"}

    session = db.query(SessionModel).filter(SessionModel.id == session_id).first()
    if not session:
        return {"status": "error", "message": "会话不存在"}

    st = normalize_session_status(getattr(session, "status", None))
    if not _session_allows_hitl_resume(db, session, session.owner_id or ""):
        return {
            "status": "ignored",
            "message": f"会话状态为 {st}，且快照无 hitl_pending，无法 Webhook 唤醒",
            "session_id": session_id,
        }

    logger.info("[HITL Webhook] project=%s session=%s action=%s", project_id, session_id, payload.get("action"))

    from gibh_agent.core.hitl_resume import stream_resume_from_hitl
    from server import agent
    from gibh_agent.core.orchestrator import AgentOrchestrator

    upload_dir = os.getenv("UPLOAD_DIR", "/app/uploads")
    orchestrator = AgentOrchestrator(agent, upload_dir=upload_dir)

    async def _gen():
        async for chunk in stream_resume_from_hitl(
            orchestrator=orchestrator,
            session_id=session_id,
            owner_id=session.owner_id,
            db=db,
            project_id=project_id,
            trigger="webhook",
        ):
            yield chunk

    return StreamingResponse(_gen(), media_type="text/event-stream")


@router.post("/api/sessions/{session_id}/resume_from_hitl")
async def resume_from_hitl(
    session_id: str,
    body: ResumeFromHitlBody,
    owner_id: str = Depends(get_current_owner_id),
    db: OrmSession = Depends(get_db_session),
):
    """
    前端 / iframe postMessage 回调：专家标注完成后唤醒 pipeline，SSE 流式返回最终专家报告。
    """
    sid = str(session_id or "").strip()
    session = db.query(SessionModel).filter(SessionModel.id == sid).first()
    if not session:
        raise HTTPException(status_code=404, detail="会话不存在")
    if session.owner_id != owner_id:
        raise HTTPException(status_code=403, detail="无权操作该会话")

    st = normalize_session_status(getattr(session, "status", None))
    if not _session_allows_hitl_resume(db, session, owner_id):
        raise HTTPException(
            status_code=409,
            detail=f"会话状态为 {st}，且快照无 hitl_pending，无法唤醒（可能已完成 HITL 或未进入复核）",
        )

    if body.skip:
        from gibh_agent.core.session_runtime import SESSION_COMPLETED, set_session_status

        set_session_status(db, sid, SESSION_COMPLETED, owner_id=owner_id)
        return {"status": "skipped", "message": "已跳过专家复核"}

    from gibh_agent.core.hitl_resume import stream_resume_from_hitl
    from server import agent
    from gibh_agent.core.orchestrator import AgentOrchestrator

    upload_dir = os.getenv("UPLOAD_DIR", "/app/uploads")
    orchestrator = AgentOrchestrator(agent, upload_dir=upload_dir)

    async def _gen():
        async for chunk in stream_resume_from_hitl(
            orchestrator=orchestrator,
            session_id=sid,
            owner_id=owner_id,
            db=db,
            project_id=body.project_id,
            trigger=body.trigger or "frontend",
        ):
            yield chunk

    return StreamingResponse(
        _gen(),
        media_type="text/event-stream",
        headers={"Cache-Control": "no-cache", "X-Accel-Buffering": "no"},
    )
