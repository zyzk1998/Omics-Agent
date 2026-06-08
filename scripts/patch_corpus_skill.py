#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
将「科学语料数据加工」技能写入存量 skills 表（幂等 UPSERT by name）。

用法（仓库根目录）:
  PYTHONPATH=. python scripts/patch_corpus_skill.py --dry-run
  PYTHONPATH=. python scripts/patch_corpus_skill.py
"""
from __future__ import annotations

import argparse
import sys
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from gibh_agent.db.corpus_skill_prompt_templates import (  # noqa: E402
    CORPUS_SKILL_PROMPT_TEMPLATE,
    CORPUS_SKILL_SEED,
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
        name = CORPUS_SKILL_SEED["name"]
        existing = db.query(SkillModel).filter(SkillModel.name == name).first()
        if existing:
            print(f"UPDATE {name!r} main_category={CORPUS_SKILL_SEED['main_category']!r}")
            if not args.dry_run:
                existing.description = CORPUS_SKILL_SEED["description"]
                existing.main_category = CORPUS_SKILL_SEED["main_category"]
                existing.sub_category = CORPUS_SKILL_SEED["sub_category"]
                existing.prompt_template = CORPUS_SKILL_PROMPT_TEMPLATE
                existing.status = "approved"
                existing.created_at = now
        else:
            print(f"INSERT {name!r}")
            if not args.dry_run:
                db.add(
                    SkillModel(
                        name=name,
                        description=CORPUS_SKILL_SEED["description"],
                        main_category=CORPUS_SKILL_SEED["main_category"],
                        sub_category=CORPUS_SKILL_SEED["sub_category"],
                        prompt_template=CORPUS_SKILL_PROMPT_TEMPLATE,
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
