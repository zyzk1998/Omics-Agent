"""
用户反馈：提交与管理员查询
POST /api/feedbacks — 需 JWT 或 X-Guest-UUID（与业务接口一致）
GET  /api/admin/feedbacks — 管理员列表
"""
from __future__ import annotations

import logging
from datetime import datetime
from typing import Any, Dict, List, Optional

from fastapi import APIRouter, Depends, HTTPException, status
from pydantic import BaseModel, Field
from sqlalchemy.orm import Session

from gibh_agent.core.deps import get_current_admin_user, get_current_owner_id
from gibh_agent.core.user_notifications import create_user_notification, notify_all_admins
from gibh_agent.db.connection import get_db_session
from gibh_agent.db.models import User, UserFeedback

logger = logging.getLogger(__name__)

router = APIRouter(prefix="/api", tags=["Feedback"])


class FeedbackCreate(BaseModel):
    type: str = Field(..., max_length=64, description="issue | suggestion | error | other")
    content: str = Field(..., max_length=16000)
    error_context: Optional[str] = Field(None, max_length=32000)
    timestamp: Optional[str] = Field(None, max_length=64)


class FeedbackOut(BaseModel):
    model_config = {"from_attributes": True}

    id: int
    owner_id: str
    feedback_type: str
    content: str
    error_context: Optional[str]
    client_timestamp: Optional[str]
    status: Optional[str] = "open"
    created_at: Optional[datetime]


class FeedbackResolveAction(BaseModel):
    action: str = Field(..., pattern="^(acknowledge|resolve)$")


@router.post("/feedbacks", response_model=Dict[str, Any])
def submit_feedback(
    body: FeedbackCreate,
    owner_id: str = Depends(get_current_owner_id),
    db: Session = Depends(get_db_session),
):
    row = UserFeedback(
        owner_id=owner_id,
        feedback_type=(body.type or "other").strip()[:64],
        content=body.content.strip(),
        error_context=(body.error_context or "").strip() or None,
        client_timestamp=(body.timestamp or "").strip()[:64] or None,
        status="open",
    )
    db.add(row)
    db.flush()
    preview = (body.content or "").strip()[:120]
    if len((body.content or "").strip()) > 120:
        preview += "…"
    try:
        notify_all_admins(
            db,
            ntype="admin_feedback_new",
            title="新用户反馈",
            content=(
                f"来自「{owner_id}」的{row.feedback_type}反馈（ID={row.id}）：{preview or '(无正文)'}"
            ),
            commit=False,
        )
    except Exception as e:
        logger.warning("通知管理员新反馈失败 id=%s: %s", row.id, e)
    db.commit()
    db.refresh(row)
    logger.info("用户反馈已入库: id=%s owner=%s type=%s", row.id, owner_id, row.feedback_type)
    return {"status": "success", "id": row.id}


@router.get("/admin/feedbacks", response_model=List[FeedbackOut])
def list_feedbacks_for_admin(
    db: Session = Depends(get_db_session),
    _admin: User = Depends(get_current_admin_user),
):
    rows = db.query(UserFeedback).order_by(UserFeedback.created_at.desc()).limit(500).all()
    return rows


@router.post("/admin/feedback/{feedback_id}/resolve")
def resolve_feedback(
    feedback_id: int,
    body: FeedbackResolveAction,
    db: Session = Depends(get_db_session),
    _admin: User = Depends(get_current_admin_user),
) -> Dict[str, Any]:
    """管理员处理反馈：acknowledge（收到）| resolve（已落实），并写入用户通知。"""
    row = db.query(UserFeedback).filter(UserFeedback.id == feedback_id).first()
    if not row:
        raise HTTPException(status_code=status.HTTP_404_NOT_FOUND, detail="反馈不存在")

    act = (body.action or "").strip().lower()
    cur = (row.status or "open").strip().lower()
    owner = (row.owner_id or "").strip()

    if act == "acknowledge":
        if cur == "acknowledged":
            return {"ok": True, "feedback_id": feedback_id, "status": "acknowledged", "message": "已是收到状态"}
        if cur == "resolved":
            raise HTTPException(status_code=400, detail="反馈已落实，无法回退为收到")
        row.status = "acknowledged"
        ntype = "feedback_acknowledged"
        title = "反馈已收到"
        content = "💌 您的反馈已收到，我们正在认真评估。"
    elif act == "resolve":
        if cur == "resolved":
            return {"ok": True, "feedback_id": feedback_id, "status": "resolved", "message": "已是落实状态"}
        row.status = "resolved"
        ntype = "feedback_resolved"
        title = "反馈已落实"
        content = "🎉 感谢您的建议！您反馈的问题现已落实/修复，快去体验吧！"
    else:
        raise HTTPException(status_code=400, detail="action 必须是 acknowledge 或 resolve")

    if owner:
        try:
            create_user_notification(
                db,
                user_id=owner,
                ntype=ntype,
                title=title,
                content=content,
                commit=False,
            )
        except Exception as e:
            logger.warning(
                "写入反馈通知失败 feedback_id=%s owner=%s action=%s: %s",
                feedback_id,
                owner,
                act,
                e,
            )

    db.commit()
    db.refresh(row)
    logger.info(
        "管理员处理反馈: id=%s action=%s status=%s owner=%s by=%s",
        feedback_id,
        act,
        row.status,
        owner,
        _admin.username,
    )
    return {
        "ok": True,
        "feedback_id": feedback_id,
        "status": row.status,
        "action": act,
        "message": "已更新并通知用户",
    }
