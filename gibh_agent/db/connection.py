"""
Phase 1 - 用户与资产中台：MySQL 连接与 SQLAlchemy 引擎

- create_engine 为惰性连接，不在模块加载时做连通性探测，不将 engine 置为 None。
- 使用 pymysql 驱动，utf8mb4；pool_pre_ping 防止长连接失效。
- 真实就绪探测与建表由 server.py 启动时智能等待完成。
"""
import os
import threading
from typing import Generator

from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker, declarative_base

# 从环境变量读取，兼容 Docker 与本地
MYSQL_HOST = os.environ.get("MYSQL_HOST", "localhost")
MYSQL_PORT = os.environ.get("MYSQL_PORT", "3306")
MYSQL_DATABASE = os.environ.get("MYSQL_DATABASE", "omics_agent")
MYSQL_USER = os.environ.get("MYSQL_USER", "omics")
MYSQL_PASSWORD = os.environ.get("MYSQL_PASSWORD", "omics_agent_pwd")

DATABASE_URL = (
    f"mysql+pymysql://{MYSQL_USER}:{MYSQL_PASSWORD}@{MYSQL_HOST}:{MYSQL_PORT}/{MYSQL_DATABASE}"
    "?charset=utf8mb4"
)

Base = declarative_base()
engine = create_engine(
    DATABASE_URL,
    pool_pre_ping=True,
    pool_recycle=300,
    pool_size=5,
    max_overflow=10,
    echo=False,
)
SessionLocal = sessionmaker(autocommit=False, autoflush=False, bind=engine)

_users_approval_schema_lock = threading.Lock()
_users_approval_schema_done = False


def is_available() -> bool:
    """兼容旧接口；engine/SessionLocal 始终存在，由 server 启动时等待 MySQL 就绪。"""
    return True


def get_db_session() -> Generator:
    global _users_approval_schema_done
    db = SessionLocal()
    try:
        if not _users_approval_schema_done:
            with _users_approval_schema_lock:
                if not _users_approval_schema_done:
                    from gibh_agent.db.user_approval_schema import (
                        ensure_users_approval_columns,
                    )

                    ensure_users_approval_columns(engine)
                    _users_approval_schema_done = True
        yield db
    finally:
        db.close()
