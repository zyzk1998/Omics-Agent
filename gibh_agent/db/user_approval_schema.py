"""
用户审核字段：幂等补齐 users.approval_status / users.email，并回填老数据。

- 供 server 启动建表后、以及 get_db_session 首次打开连接时调用，避免 ORM 与表结构不一致导致注册 500。
- 不使用 ADD COLUMN ... COMMENT，以兼容 SQLite 与部分 MySQL/MariaDB 配置。
"""
from __future__ import annotations

import logging

from sqlalchemy import inspect as sa_inspect, text
from typing import Optional

from sqlalchemy.engine import Engine
from sqlalchemy.exc import OperationalError, ProgrammingError

logger = logging.getLogger(__name__)


def _is_duplicate_column_error(exc: BaseException) -> bool:
    s = str(getattr(exc, "orig", exc) or exc).lower()
    if "duplicate column" in s:
        return True
    if "1060" in s:  # MySQL ER_DUP_FIELDNAME
        return True
    if "already exists" in s and "column" in s:
        return True
    return False


def ensure_users_approval_columns(engine: Optional[Engine]) -> None:
    if engine is None:
        return
    insp = sa_inspect(engine)
    if not insp.has_table("users"):
        return
    col_names = {c["name"] for c in insp.get_columns("users")}

    if "approval_status" not in col_names:
        try:
            with engine.begin() as conn:
                conn.execute(
                    text(
                        "ALTER TABLE users ADD COLUMN approval_status VARCHAR(32) NULL"
                    )
                )
            logger.info("[DB] 已添加列 users.approval_status")
        except (OperationalError, ProgrammingError) as e:
            if not _is_duplicate_column_error(e):
                logger.error(
                    "[DB] 添加 users.approval_status 失败: %s", e, exc_info=True
                )
                raise
            logger.debug("[DB] users.approval_status 已存在（并发/重复执行）")

    if "email" not in col_names:
        try:
            with engine.begin() as conn:
                conn.execute(
                    text("ALTER TABLE users ADD COLUMN email VARCHAR(512) NULL")
                )
            logger.info("[DB] 已添加列 users.email")
        except (OperationalError, ProgrammingError) as e:
            if not _is_duplicate_column_error(e):
                logger.error("[DB] 添加 users.email 失败: %s", e, exc_info=True)
                raise
            logger.debug("[DB] users.email 已存在（并发/重复执行）")

    with engine.begin() as conn:
        conn.execute(
            text(
                "UPDATE users SET approval_status = 'approved' "
                "WHERE approval_status IS NULL"
            )
        )
    logger.info("[DB] users 审核字段迁移已就绪（NULL→approved 回填已执行）")
