"""
sessions.composer_draft：主页面输入框草稿（文本 + 附件药丸元数据）幂等补齐。
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


def ensure_sessions_composer_draft_column(engine: Optional[Engine]) -> None:
    if engine is None:
        return
    insp = sa_inspect(engine)
    if not insp.has_table("sessions"):
        return
    col_names = {c["name"] for c in insp.get_columns("sessions")}
    if "composer_draft" in col_names:
        return
    try:
        with engine.begin() as conn:
            conn.execute(text("ALTER TABLE sessions ADD COLUMN composer_draft JSON NULL"))
        logger.info("[DB] 已添加列 sessions.composer_draft")
    except (OperationalError, ProgrammingError) as e:
        if not _is_duplicate_column_error(e):
            logger.error("[DB] 添加 sessions.composer_draft 失败: %s", e, exc_info=True)
            raise
        logger.debug("[DB] sessions.composer_draft 已存在（并发/重复执行）")
