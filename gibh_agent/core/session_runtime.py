"""
会话运行时：status 状态机（idle / running / completed / failed）与 DB 写入。
"""
from __future__ import annotations

import logging
from typing import Optional

from sqlalchemy.orm import Session as OrmSession
from sqlalchemy.orm.attributes import flag_modified

from gibh_agent.db.models import Session as SessionModel

logger = logging.getLogger(__name__)

SESSION_IDLE = "idle"
SESSION_RUNNING = "running"
SESSION_COMPLETED = "completed"
SESSION_FAILED = "failed"

_VALID = frozenset({SESSION_IDLE, SESSION_RUNNING, SESSION_COMPLETED, SESSION_FAILED})


def normalize_session_status(raw: Optional[str], *, default: str = SESSION_IDLE) -> str:
    s = (raw or "").strip().lower()
    if s in _VALID:
        return s
    return default


def set_session_status(
    db: OrmSession,
    session_id: str,
    status: str,
    *,
    owner_id: Optional[str] = None,
) -> bool:
    sid = (session_id or "").strip()
    if not sid:
        return False
    st = normalize_session_status(status)
    q = db.query(SessionModel).filter(SessionModel.id == sid)
    if owner_id:
        q = q.filter(SessionModel.owner_id == owner_id)
    row = q.first()
    if not row:
        return False
    row.status = st
    flag_modified(row, "status")
    try:
        db.commit()
        logger.info("[SessionRuntime] session=%s status=%s", sid, st)
        return True
    except Exception as exc:
        db.rollback()
        logger.warning("[SessionRuntime] set status failed session=%s: %s", sid, exc)
        return False
