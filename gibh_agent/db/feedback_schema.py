"""
用户反馈表字段：幂等补齐 user_feedbacks.status。

- 供 server 启动建表后调用，避免 ORM 与表结构不一致。
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


def ensure_user_feedbacks_status_column(engine: Optional[Engine]) -> None:
    if engine is None:
        return
    insp = sa_inspect(engine)
    if not insp.has_table("user_feedbacks"):
        return
    col_names = {c["name"] for c in insp.get_columns("user_feedbacks")}

    if "status" not in col_names:
        try:
            with engine.begin() as conn:
                conn.execute(
                    text(
                        "ALTER TABLE user_feedbacks ADD COLUMN status VARCHAR(32) NULL"
                    )
                )
            logger.info("[DB] 已添加列 user_feedbacks.status")
        except (OperationalError, ProgrammingError) as e:
            if not _is_duplicate_column_error(e):
                logger.error(
                    "[DB] 添加 user_feedbacks.status 失败: %s", e, exc_info=True
                )
                raise
            logger.debug("[DB] user_feedbacks.status 已存在（并发/重复执行）")

    with engine.begin() as conn:
        conn.execute(
            text(
                "UPDATE user_feedbacks SET status = 'open' "
                "WHERE status IS NULL"
            )
        )
    logger.info("[DB] user_feedbacks.status 迁移已就绪（NULL→open 回填已执行）")
