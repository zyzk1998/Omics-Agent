# -*- coding: utf-8 -*-
"""用户站内通知：写入与查询辅助。"""
from __future__ import annotations

import logging
from typing import Any, Dict, List, Optional

from sqlalchemy.orm import Session

from gibh_agent.db.models import UserNotification

logger = logging.getLogger(__name__)


def create_user_notification(
    db: Session,
    *,
    user_id: str,
    ntype: str,
    title: str,
    content: str,
    commit: bool = True,
) -> UserNotification:
    uid = (user_id or "").strip()
    if not uid:
        raise ValueError("user_id 不能为空")
    row = UserNotification(
        user_id=uid,
        type=(ntype or "system_notice").strip()[:64],
        title=(title or "").strip()[:512],
        content=(content or "").strip(),
        is_read=0,
    )
    db.add(row)
    if commit:
        db.commit()
        db.refresh(row)
    logger.info("user_notification: to=%s type=%s id=%s", uid, row.type, row.id)
    return row


def notification_to_dict(row: UserNotification) -> Dict[str, Any]:
    return {
        "id": row.id,
        "user_id": row.user_id,
        "type": row.type,
        "title": row.title,
        "content": row.content,
        "is_read": bool(row.is_read),
        "created_at": row.created_at.isoformat() if row.created_at else None,
    }


def list_recent_notifications(
    db: Session, user_id: str, *, limit: int = 30
) -> List[Dict[str, Any]]:
    rows = (
        db.query(UserNotification)
        .filter(UserNotification.user_id == user_id)
        .order_by(UserNotification.created_at.desc())
        .limit(max(1, min(limit, 100)))
        .all()
    )
    return [notification_to_dict(r) for r in rows]


def count_unread_notifications(db: Session, user_id: str) -> int:
    return (
        db.query(UserNotification)
        .filter(UserNotification.user_id == user_id, UserNotification.is_read == 0)
        .count()
    )


def list_admin_usernames(db: Session) -> List[str]:
    """已审核通过的管理员 username 列表（用于待办通知）。"""
    from gibh_agent.db.models import User

    rows = db.query(User).filter(User.role == "admin").all()
    out: List[str] = []
    for u in rows:
        raw = getattr(u, "approval_status", None)
        st = (str(raw).strip().lower() if raw is not None else "")
        if st not in ("", "approved"):
            continue
        un = (u.username or "").strip()
        if un:
            out.append(un)
    return out


def notify_all_admins(
    db: Session,
    *,
    ntype: str,
    title: str,
    content: str,
    commit: bool = False,
) -> int:
    """向所有管理员写入同一条待办通知；返回写入条数。"""
    admins = list_admin_usernames(db)
    if not admins:
        logger.info("notify_all_admins: 无可用管理员，跳过 type=%s", ntype)
        return 0
    for un in admins:
        create_user_notification(
            db,
            user_id=un,
            ntype=ntype,
            title=title,
            content=content,
            commit=False,
        )
    if commit:
        db.commit()
    return len(admins)
