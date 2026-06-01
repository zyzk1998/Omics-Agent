# -*- coding: utf-8 -*-
"""
用户站内通知 API

- GET  /api/user/notifications/unread — 未读数量 + 最近未读列表
- GET  /api/user/notifications — 最近消息（含已读）
- POST /api/user/notifications/{id}/read — 单条已读
- POST /api/user/notifications/read-all — 全部已读
"""
from __future__ import annotations

import logging
from typing import Any, Dict, List, Optional

from fastapi import APIRouter, Depends, HTTPException, Query
from pydantic import BaseModel, Field
from sqlalchemy.orm import Session

from gibh_agent.core.deps import get_current_owner_id
from gibh_agent.core.user_notifications import (
    count_unread_notifications,
    list_recent_notifications,
    notification_to_dict,
)
from gibh_agent.db.connection import get_db_session
from gibh_agent.db.models import UserNotification

logger = logging.getLogger(__name__)

router = APIRouter(prefix="/api/user", tags=["UserNotifications"])


class MarkReadBody(BaseModel):
    mark_all: bool = Field(False, description="为 true 时忽略 path id，全部标已读")


@router.get("/notifications/unread")
def get_unread_notifications(
    owner_id: str = Depends(get_current_owner_id),
    db: Session = Depends(get_db_session),
) -> Dict[str, Any]:
    unread_count = count_unread_notifications(db, owner_id)
    rows = (
        db.query(UserNotification)
        .filter(UserNotification.user_id == owner_id, UserNotification.is_read == 0)
        .order_by(UserNotification.created_at.desc())
        .limit(20)
        .all()
    )
    return {
        "unread_count": unread_count,
        "items": [notification_to_dict(r) for r in rows],
    }


@router.get("/notifications")
def list_notifications(
    owner_id: str = Depends(get_current_owner_id),
    db: Session = Depends(get_db_session),
    limit: int = Query(30, ge=1, le=100),
) -> Dict[str, Any]:
    items = list_recent_notifications(db, owner_id, limit=limit)
    return {
        "unread_count": count_unread_notifications(db, owner_id),
        "items": items,
    }


@router.post("/notifications/{notification_id}/read")
def mark_notification_read(
    notification_id: int,
    owner_id: str = Depends(get_current_owner_id),
    db: Session = Depends(get_db_session),
) -> Dict[str, Any]:
    row = (
        db.query(UserNotification)
        .filter(UserNotification.id == notification_id, UserNotification.user_id == owner_id)
        .first()
    )
    if not row:
        raise HTTPException(status_code=404, detail="通知不存在")
    row.is_read = 1
    db.commit()
    return {"status": "success", "id": row.id, "is_read": True}


@router.post("/notifications/read-all")
def mark_all_notifications_read(
    owner_id: str = Depends(get_current_owner_id),
    db: Session = Depends(get_db_session),
) -> Dict[str, Any]:
    updated = (
        db.query(UserNotification)
        .filter(UserNotification.user_id == owner_id, UserNotification.is_read == 0)
        .update({UserNotification.is_read: 1}, synchronize_session=False)
    )
    db.commit()
    return {"status": "success", "updated": updated}
