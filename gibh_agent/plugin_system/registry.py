"""
动态技能持久化：ORM 定义与注册逻辑。
【重要】在 Base.metadata.create_all 之前执行: import gibh_agent.plugin_system.registry  # noqa: F401
"""
from __future__ import annotations

import json
import logging
import os
from datetime import datetime
from typing import Any, Dict, Optional

from sqlalchemy import Column, DateTime, Integer, String, Text, UniqueConstraint, text
from sqlalchemy.engine import Engine
from sqlalchemy.orm import Session
from sqlalchemy.types import JSON

from gibh_agent.db.connection import Base

logger = logging.getLogger(__name__)

_DEFAULT_WORKER_ROUTE = "/api/dynamic/run"
MAX_PROMPT_CHUNK = 12000


class DynamicSkillPlugin(Base):
    """用户上传的动态技能：prompt=仅提示词；script=含 main.py，由 worker 执行。"""

    __tablename__ = "dynamic_skill_plugins"
    __table_args__ = (UniqueConstraint("name", name="uq_dynamic_skill_plugins_name"),)

    id = Column(Integer, primary_key=True, autoincrement=True)
    name = Column(String(255), nullable=False, index=True)
    display_name = Column(String(512), nullable=True)
    description = Column(Text, nullable=True)
    parameters_schema = Column(JSON, nullable=True)
    skill_type = Column(String(16), nullable=False, default="script")  # prompt | script
    script_path = Column(String(2048), nullable=True)
    extract_dir = Column(String(2048), nullable=True)
    worker_route = Column(String(512), nullable=False, default=_DEFAULT_WORKER_ROUTE)
    author_id = Column(String(255), nullable=False, index=True)
    status = Column(String(32), nullable=False, default="approved")
    created_at = Column(DateTime, nullable=False, default=datetime.utcnow)


def migrate_dynamic_skill_plugins_schema(engine: Engine) -> None:
    """MySQL：为已存在表补齐 skill_type，并将 script_path 改为可空（幂等、忽略已存在错误）。"""
    stmts = [
        "ALTER TABLE dynamic_skill_plugins ADD COLUMN skill_type VARCHAR(16) NOT NULL DEFAULT 'script'",
        "ALTER TABLE dynamic_skill_plugins MODIFY COLUMN script_path VARCHAR(2048) NULL",
    ]
    with engine.begin() as conn:
        for sql in stmts:
            try:
                conn.execute(text(sql))
                logger.info("plugin_system migrate: %s", sql[:60])
            except Exception as e:
                logger.debug("plugin_system migrate skip (expected if applied): %s", e)


def default_worker_base_url() -> str:
    return (os.environ.get("PYSKILLS_BASE_URL") or "http://worker-pyskills:8000").rstrip("/")


def persist_plugin(
    db: Session,
    *,
    name: str,
    display_name: Optional[str],
    description: Optional[str],
    parameters_schema: Dict[str, Any],
    skill_type: str,
    script_path: Optional[str],
    extract_dir: Optional[str],
    author_id: str,
    worker_route: Optional[str] = None,
    status: str = "approved",
) -> DynamicSkillPlugin:
    """写入 dynamic_skill_plugins；name 全局唯一；script 类必须提供 script_path。"""
    st = (skill_type or "script").strip().lower()
    if st not in ("prompt", "script"):
        raise ValueError("skill_type 必须是 prompt 或 script")
    if st == "script" and not (script_path and str(script_path).strip()):
        raise ValueError("script 类型技能必须提供 main.py 的 script_path")
    route = (worker_route or _DEFAULT_WORKER_ROUTE).strip() or _DEFAULT_WORKER_ROUTE
    sp = os.path.abspath(script_path) if script_path else None

    row = db.query(DynamicSkillPlugin).filter(DynamicSkillPlugin.name == name).first()
    if row:
        if row.author_id != author_id:
            raise ValueError(f"技能标识 name={name!r} 已被其他用户占用")
        row.display_name = display_name
        row.description = description
        row.parameters_schema = parameters_schema
        row.skill_type = st
        row.script_path = sp
        row.extract_dir = extract_dir
        row.worker_route = route
        row.status = status
        db.commit()
        db.refresh(row)
        logger.info("plugin_system: 更新动态技能 name=%s id=%s type=%s", name, row.id, st)
        return row

    row = DynamicSkillPlugin(
        name=name,
        display_name=display_name,
        description=description,
        parameters_schema=parameters_schema,
        skill_type=st,
        script_path=sp,
        extract_dir=extract_dir,
        worker_route=route,
        author_id=author_id,
        status=status,
    )
    db.add(row)
    db.commit()
    db.refresh(row)
    logger.info("plugin_system: 注册动态技能 name=%s id=%s type=%s", name, row.id, st)
    return row


def plugin_to_skill_plaza_payload(row: DynamicSkillPlugin) -> Dict[str, Any]:
    """供技能广场列表合并展示（与现有 Skill JSON 形状对齐）。"""
    schema_json = json.dumps(row.parameters_schema or {}, ensure_ascii=False)[:MAX_PROMPT_CHUNK]
    prompt = (
        f"[Skill_Route: execute_dynamic_skill]\n"
        f"[Dynamic_Plugin]\n"
        f"plugin_id={row.id}\n"
        f"name={row.name}\n"
        f"skill_type={row.skill_type}\n"
        f"parameters_schema={schema_json}\n"
    )
    return {
        "id": f"dyn_{row.id}",
        "name": row.display_name or row.name,
        "description": row.description or "",
        "main_category": "动态插件",
        "sub_category": "用户上传",
        "prompt_template": prompt,
        "author_id": row.author_id,
        "created_at": row.created_at.isoformat() if row.created_at else None,
        "saved": False,
        "is_dynamic_plugin": True,
        "plugin_db_id": row.id,
        "skill_type": row.skill_type,
    }
