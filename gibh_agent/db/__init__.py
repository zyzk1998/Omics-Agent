"""
数据库模块

- database: 工作流/任务历史（SQLite）
- connection: Phase 1 用户与资产中台 MySQL 连接（SQLAlchemy，容错）
"""

from .database import WorkflowDB, get_db

# Phase 1: MySQL 连接（容错；导入失败时不影响 WorkflowDB）
try:
    from .connection import (
        engine as mysql_engine,
        SessionLocal as MySQLSessionLocal,
        Base as MySQLBase,
        get_db_session,
        is_available as mysql_is_available,
    )
    _mysql_ok = True
except Exception:
    mysql_engine = None
    MySQLSessionLocal = None
    MySQLBase = None
    get_db_session = None
    mysql_is_available = lambda: False
    _mysql_ok = False

__all__ = [
    "WorkflowDB",
    "get_db",
    "mysql_engine",
    "MySQLSessionLocal",
    "MySQLBase",
    "get_db_session",
    "mysql_is_available",
]

