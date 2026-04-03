#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""列出 MySQL users 表中最近注册的用户（需与 api-server 相同 MYSQL_* 环境变量）。"""
from __future__ import annotations

import argparse
import os
import sys

# 仓库根在 PYTHONPATH
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from gibh_agent.db.connection import SessionLocal, engine  # noqa: E402
from gibh_agent.db.models import User  # noqa: E402
from gibh_agent.db.user_approval_schema import ensure_users_approval_columns  # noqa: E402


def main() -> int:
    p = argparse.ArgumentParser(description="按 created_at 倒序列出注册用户")
    p.add_argument("-n", "--limit", type=int, default=30, help="最多条数，默认 30")
    args = p.parse_args()
    ensure_users_approval_columns(engine)
    db = SessionLocal()
    try:
        rows = (
            db.query(User)
            .order_by(User.created_at.desc())
            .limit(max(1, args.limit))
            .all()
        )
        if not rows:
            print("(users 表为空或不可读)")
            return 0
        for u in rows:
            ap = getattr(u, "approval_status", None) or "—"
            print(
                f"id={u.id}\tusername={u.username}\trole={u.role}\tapproval={ap}\tcreated_at={u.created_at}"
            )
    finally:
        db.close()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
