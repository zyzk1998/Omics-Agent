#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
将 Prompt 软技能写入存量 skills 表（幂等 UPSERT by name）。

用法（仓库根目录）:
  PYTHONPATH=. python scripts/patch_prompt_soft_skills.py --dry-run
  PYTHONPATH=. python scripts/patch_prompt_soft_skills.py
"""
from __future__ import annotations

import argparse
import sys
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from gibh_agent.db.prompt_soft_skill_templates import (  # noqa: E402
    PROMPT_SOFT_SKILL_SEEDS,
    PROMPT_SOFT_SKILL_TEMPLATES,
)
from gibh_agent.db.connection import SessionLocal
from gibh_agent.db.models import Skill as SkillModel


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()
    db = SessionLocal()
    try:
        now = datetime.now(timezone.utc)
        for row in PROMPT_SOFT_SKILL_SEEDS:
            name = row["name"]
            pt = PROMPT_SOFT_SKILL_TEMPLATES[name]
            existing = db.query(SkillModel).filter(SkillModel.name == name).first()
            if existing:
                print(f"UPDATE {name!r} main_category={row['main_category']!r}")
                if not args.dry_run:
                    existing.description = row["description"]
                    existing.main_category = row["main_category"]
                    existing.sub_category = row["sub_category"]
                    existing.prompt_template = pt
                    existing.status = "approved"
                    existing.created_at = now
            else:
                print(f"INSERT {name!r}")
                if not args.dry_run:
                    db.add(
                        SkillModel(
                            name=name,
                            description=row["description"],
                            main_category=row["main_category"],
                            sub_category=row["sub_category"],
                            prompt_template=pt,
                            author_id="system",
                            status="approved",
                            created_at=now,
                        )
                    )
        if args.dry_run:
            db.rollback()
            print("dry-run: 未写入")
        else:
            db.commit()
            print("已提交")
        return 0
    finally:
        db.close()


if __name__ == "__main__":
    raise SystemExit(main())
