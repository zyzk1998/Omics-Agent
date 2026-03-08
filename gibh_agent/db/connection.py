"""
Phase 1 - 用户与资产中台：MySQL 连接与 SQLAlchemy 引擎（防阻塞容错设计）

- 若 MySQL 未启动或不可达，仅打 Warning 日志，不抛异常，保证 FastAPI 正常启动。
- 业务在使用 get_db_session() 时需自行判断 engine 是否可用（is_available()）。
- 使用 pymysql 驱动，支持 utf8mb4；pool_pre_ping 防止长连接失效。
"""
import os
import logging
from typing import Generator, Optional

logger = logging.getLogger(__name__)

# 从环境变量读取，兼容 Docker 与本地
MYSQL_HOST = os.environ.get("MYSQL_HOST", "localhost")
MYSQL_PORT = os.environ.get("MYSQL_PORT", "3306")
MYSQL_DATABASE = os.environ.get("MYSQL_DATABASE", "omics_agent")
MYSQL_USER = os.environ.get("MYSQL_USER", "omics")
MYSQL_PASSWORD = os.environ.get("MYSQL_PASSWORD", "omics_agent_pwd")

# SQLAlchemy 对象：连接失败时为 None，不阻塞应用启动
engine = None
SessionLocal = None
Base = None

try:
    from sqlalchemy import create_engine, text
    from sqlalchemy.orm import sessionmaker, declarative_base

    Base = declarative_base()
    # pymysql 驱动，utf8mb4
    url = (
        f"mysql+pymysql://{MYSQL_USER}:{MYSQL_PASSWORD}@{MYSQL_HOST}:{MYSQL_PORT}/{MYSQL_DATABASE}"
        "?charset=utf8mb4"
    )
    engine = create_engine(
        url,
        pool_pre_ping=True,
        pool_recycle=300,
        pool_size=5,
        max_overflow=10,
        echo=False,
    )
    # 仅做一次连通性探测，失败则置空 engine，绝不抛异常
    try:
        with engine.connect() as conn:
            conn.execute(text("SELECT 1"))
    except Exception as e:
        logger.warning(
            "⚠️ [DB] MySQL 连通性检查失败，将禁用数据库依赖（应用仍可正常启动）: %s",
            e,
        )
        engine = None

    if engine is not None:
        SessionLocal = sessionmaker(autocommit=False, autoflush=False, bind=engine)
        logger.info("✅ [DB] MySQL 连接就绪 (%s)", MYSQL_DATABASE)
except ImportError as e:
    logger.warning(
        "⚠️ [DB] 未安装 SQLAlchemy/pymysql，数据库功能不可用（应用仍可正常启动）: %s",
        e,
    )
except Exception as e:
    logger.warning(
        "⚠️ [DB] MySQL 连接初始化失败，应用将继续以无 DB 模式运行: %s",
        e,
    )
    engine = None
    SessionLocal = None
    if Base is None:
        try:
            from sqlalchemy.orm import declarative_base
            Base = declarative_base()
        except ImportError:
            Base = None


def is_available() -> bool:
    """是否已成功连接 MySQL（供业务层判断是否使用 DB）。"""
    return engine is not None and SessionLocal is not None


# 兼容旧调用
is_db_available = is_available


def get_db_session() -> Generator:
    """
    供 FastAPI Depends 使用的 Session 生成器。
    MySQL 不可用时抛出 RuntimeError，业务层应通过 is_available() 判断后再依赖。
    """
    if SessionLocal is None:
        raise RuntimeError("数据库不可用（MySQL 未连接或未安装依赖）")
    session = SessionLocal()
    try:
        yield session
    finally:
        session.close()
