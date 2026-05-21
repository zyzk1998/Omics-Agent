"""
sessions.status：会话级执行状态（idle / running / completed / failed）幂等补齐。
"""
from __future__ import annotations

import logging
from typing import Optional

from sqlalchemy import inspect as sa_inspect, text
from sqlalchemy.engine import Engine
from sqlalchemy.exc import OperationalError, ProgrammingError

logger = logging.getLogger(__name__)


def _is_duplicate_column_error(exc: BaseException) -> bool:
    s = str(getattr(exc, "orig", exc) or exc).lower()
    if "duplicate column" in s:
        return True
    if "1060" in s:
        return True
    if "already exists" in s and "column" in s:
        return True
    return False


def ensure_sessions_status_column(engine: Optional[Engine]) -> None:
    if engine is None:
        return
    insp = sa_inspect(engine)
    if not insp.has_table("sessions"):
        return
    col_names = {c["name"] for c in insp.get_columns("sessions")}
    if "status" not in col_names:
        try:
            with engine.begin() as conn:
                conn.execute(
                    text(
                        "ALTER TABLE sessions ADD COLUMN status VARCHAR(32) NOT NULL DEFAULT 'idle'"
                    )
                )
            logger.info("[DB] 已添加列 sessions.status")
        except (OperationalError, ProgrammingError) as e:
            if not _is_duplicate_column_error(e):
                logger.error("[DB] 添加 sessions.status 失败: %s", e, exc_info=True)
                raise
            logger.debug("[DB] sessions.status 已存在（并发/重复执行）")
    with engine.begin() as conn:
        conn.execute(text("UPDATE sessions SET status = 'idle' WHERE status IS NULL OR status = ''"))
