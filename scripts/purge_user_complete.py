#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
按 username 物理抹除用户相关数据（不可逆）：
- MySQL: messages → sessions → assets(含文件) → workflow_templates → user_saved_skills →
         他人对该用户 skills 的收藏 → skills → dynamic_skill_plugins(含目录) → users
- SQLite workflows.db: saved_workflows、job_history（user_id 与 username 一致时）

用法:
  PYTHONPATH=. python scripts/purge_user_complete.py --username mars --yes
"""
from __future__ import annotations

import argparse
import logging
import os
import shutil
import sqlite3
import sys
from pathlib import Path

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
logger = logging.getLogger("purge_user")

from gibh_agent.db.connection import SessionLocal  # noqa: E402
from gibh_agent.db.models import (  # noqa: E402
    Asset,
    Message,
    Session,
    Skill,
    User,
    UserSavedSkill,
    WorkflowTemplate,
)
from gibh_agent.plugin_system.registry import DynamicSkillPlugin  # noqa: E402


def _purge_sqlite_workflow_db(username: str, db_path: str) -> None:
    p = Path(db_path)
    if not p.is_file():
        logger.info("SQLite 跳过（不存在）: %s", p)
        return
    conn = sqlite3.connect(str(p))
    try:
        cur = conn.cursor()
        cur.execute("DELETE FROM saved_workflows WHERE user_id = ?", (username,))
        sw = cur.rowcount
        cur.execute("DELETE FROM job_history WHERE user_id = ?", (username,))
        jh = cur.rowcount
        conn.commit()
        logger.info("SQLite %s: saved_workflows 删 %s 行, job_history 删 %s 行", p, sw, jh)
    finally:
        conn.close()


def _rm_path(path: str | None, *, is_dir: bool = False) -> None:
    if not path or not str(path).strip():
        return
    p = Path(path)
    try:
        if not p.exists():
            return
        if is_dir:
            shutil.rmtree(p, ignore_errors=False)
        else:
            p.unlink(missing_ok=True)
        logger.info("已删磁盘: %s", p)
    except OSError as e:
        logger.warning("删磁盘失败（忽略继续）: %s — %s", p, e)


def purge(username: str, *, sqlite_paths: list[str]) -> None:
    db = SessionLocal()
    try:
        u = db.query(User).filter(User.username == username).one_or_none()
        if not u:
            logger.warning("MySQL users 中无此 username: %s（仍清理归属数据）", username)

        session_ids = [s.id for s in db.query(Session).filter(Session.owner_id == username).all()]
        if session_ids:
            n = db.query(Message).filter(Message.session_id.in_(session_ids)).delete(synchronize_session=False)
            logger.info("messages 删除 %s 条", n)
        n = db.query(Session).filter(Session.owner_id == username).delete(synchronize_session=False)
        logger.info("sessions 删除 %s 条", n)

        assets = db.query(Asset).filter(Asset.owner_id == username).all()
        for a in assets:
            _rm_path(a.file_path, is_dir=False)
        n = db.query(Asset).filter(Asset.owner_id == username).delete(synchronize_session=False)
        logger.info("assets 删除 %s 条", n)

        n = db.query(WorkflowTemplate).filter(WorkflowTemplate.owner_id == username).delete(
            synchronize_session=False
        )
        logger.info("workflow_templates 删除 %s 条", n)

        n = db.query(UserSavedSkill).filter(UserSavedSkill.owner_id == username).delete(
            synchronize_session=False
        )
        logger.info("user_saved_skills(owner) 删除 %s 条", n)

        skill_ids = [s.id for s in db.query(Skill).filter(Skill.author_id == username).all()]
        if skill_ids:
            n = db.query(UserSavedSkill).filter(UserSavedSkill.skill_id.in_(skill_ids)).delete(
                synchronize_session=False
            )
            logger.info("user_saved_skills(他人收藏该作者技能) 删除 %s 条", n)
            n = db.query(Skill).filter(Skill.author_id == username).delete(synchronize_session=False)
            logger.info("skills 删除 %s 条", n)

        plugins = db.query(DynamicSkillPlugin).filter(DynamicSkillPlugin.author_id == username).all()
        for pl in plugins:
            _rm_path(pl.script_path, is_dir=False)
            _rm_path(pl.extract_dir, is_dir=True)
        n = db.query(DynamicSkillPlugin).filter(DynamicSkillPlugin.author_id == username).delete(
            synchronize_session=False
        )
        logger.info("dynamic_skill_plugins 删除 %s 条", n)

        if u:
            db.delete(u)
            logger.info("users 行删除: id=%s username=%s", u.id, username)

        db.commit()
    except Exception:
        db.rollback()
        raise
    finally:
        db.close()

    for sp in sqlite_paths:
        _purge_sqlite_workflow_db(username, sp)


def main() -> int:
    ap = argparse.ArgumentParser(description="物理抹除指定 username 的全站关联数据")
    ap.add_argument("--username", required=True, help="注册用户名，如 mars")
    ap.add_argument("--yes", action="store_true", help="确认执行（必传）")
    ap.add_argument(
        "--sqlite",
        action="append",
        default=[],
        help="workflows.db 路径，可多次指定；默认 workflows.db 与 data/workflows.db",
    )
    args = ap.parse_args()
    if not args.yes:
        logger.error("拒绝执行：请附加 --yes 确认不可逆删除")
        return 2

    default_sqlite = ["workflows.db", "data/workflows.db"]
    sqlite_paths = list(dict.fromkeys(args.sqlite + default_sqlite))

    uname = args.username.strip()
    if not uname:
        logger.error("username 为空")
        return 2

    logger.warning("开始抹除用户 %r …", uname)
    purge(uname, sqlite_paths=sqlite_paths)
    logger.warning("完成: 用户 %r 相关数据已按脚本范围清理", uname)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
