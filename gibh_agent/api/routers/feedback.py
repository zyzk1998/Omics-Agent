"""
用户反馈：提交与管理员查询
POST /api/feedbacks — 需 JWT 或 X-Guest-UUID（与业务接口一致）
GET  /api/admin/feedbacks — 管理员列表
"""
from __future__ import annotations

import logging
from datetime import datetime
from typing import Any, Dict, List, Optional

from fastapi import APIRouter, Depends
from pydantic import BaseModel, Field
from sqlalchemy.orm import Session

from gibh_agent.core.deps import get_current_admin_user, get_current_owner_id
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
    created_at: Optional[datetime]


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
    )
    db.add(row)
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
