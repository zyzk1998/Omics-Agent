# -*- coding: utf-8 -*-
"""messages 表索引：加速 session_id + role 下按 created_at 取最新一条。"""
from __future__ import annotations

import logging

from sqlalchemy import text
from sqlalchemy.engine import Engine

logger = logging.getLogger(__name__)

_INDEX_NAME = "idx_messages_session_role_created"


def ensure_messages_session_role_created_index(engine: Engine) -> None:
    """幂等创建 (session_id, role, created_at) 复合索引。"""
    try:
        with engine.connect() as conn:
            exists = conn.execute(
                text(
                    "SELECT COUNT(1) FROM information_schema.statistics "
                    "WHERE table_schema = DATABASE() "
                    "AND table_name = 'messages' "
                    "AND index_name = :idx"
                ),
                {"idx": _INDEX_NAME},
            ).scalar()
            if exists:
                return
            conn.execute(
                text(
                    f"CREATE INDEX {_INDEX_NAME} "
                    "ON messages (session_id, role, created_at)"
                )
            )
            conn.commit()
            logger.info("[DB] created index %s on messages", _INDEX_NAME)
    except Exception as exc:
        logger.warning("[DB] ensure messages index skipped: %s", exc)
