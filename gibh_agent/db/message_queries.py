# -*- coding: utf-8 -*-
"""messages 表查询优化：避免 ORDER BY 时加载超大 JSON content 导致 sort buffer 溢出。"""
from __future__ import annotations

import logging
from typing import Any, Dict, Optional, Tuple

from sqlalchemy.orm import Session as OrmSession

from gibh_agent.db.models import Message as MessageModel, Session as SessionModel

logger = logging.getLogger(__name__)


def get_latest_agent_message_id(db: OrmSession, session_id: str) -> Optional[int]:
    """仅按 id 排序取最新 agent 消息，不 SELECT content（防 MySQL 1038）。"""
    sid = str(session_id or "").strip()
    if not sid:
        return None
    try:
        row_id = (
            db.query(MessageModel.id)
            .filter(MessageModel.session_id == sid, MessageModel.role == "agent")
            .order_by(MessageModel.created_at.desc(), MessageModel.id.desc())
            .limit(1)
            .scalar()
        )
        return int(row_id) if row_id is not None else None
    except (TypeError, ValueError):
        return None


def get_latest_agent_message(db: OrmSession, session_id: str) -> Optional[MessageModel]:
    """两步查询：先取 id，再按主键加载完整行。"""
    mid = get_latest_agent_message_id(db, session_id)
    if mid is None:
        return None
    return db.query(MessageModel).filter(MessageModel.id == mid).first()


def load_latest_agent_snapshot(
    db: OrmSession,
    session_id: str,
    owner_id: str,
) -> Tuple[Optional[MessageModel], Dict[str, Any]]:
    """加载会话最新 agent 消息的 state_snapshot（含 ownership 校验）。"""
    sid = str(session_id or "").strip()
    if not sid:
        return None, {}
    session = db.query(SessionModel).filter(SessionModel.id == sid).first()
    if not session or session.owner_id != owner_id:
        return None, {}
    msg = get_latest_agent_message(db, sid)
    if not msg:
        return None, {}
    content = msg.content if isinstance(msg.content, dict) else {}
    snap = content.get("state_snapshot") if isinstance(content.get("state_snapshot"), dict) else content
    return msg, snap if isinstance(snap, dict) else {}
