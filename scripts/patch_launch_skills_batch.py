#!/usr/bin/env python3
"""
一次性 PATCH：首发 10 项核心技能（gibh_agent/skills/ BaseSkill）写入 skills 表快车道模板。

用法:
  PYTHONPATH=. python scripts/patch_launch_skills_batch.py
  PYTHONPATH=. python scripts/patch_launch_skills_batch.py --dry-run
  docker exec gibh_v2_api python /app/scripts/patch_launch_skills_batch.py
"""
from __future__ import annotations

import argparse
import sys
from datetime import datetime
from pathlib import Path

_ROOT = Path(__file__).resolve().parents[1]
if str(_ROOT) not in sys.path:
    sys.path.insert(0, str(_ROOT))

from gibh_agent.db.connection import SessionLocal  # noqa: E402
from gibh_agent.db.launch_core_skill_prompt_templates import (  # noqa: E402
    LAUNCH_CORE_PROMPTS_BY_SKILL_NAME,
    LAUNCH_CORE_SKILL_SEEDS,
)
from gibh_agent.db.models import Skill  # noqa: E402

# 10 项全覆盖：8 项新增种子 + 2 项化学占位（3D / 格式转换）
_PATCH_META: dict[str, dict[str, str]] = {row["name"]: row for row in LAUNCH_CORE_SKILL_SEEDS}
_PATCH_META["3D分子结构渲染工具"] = {
    "name": "3D分子结构渲染工具",
    "main_category": "化学",
    "sub_category": "数据可视化",
    "description": "使用 RDKit 由 SMILES 生成 3D 构象（MOL/SDF 文本或 PNG）。",
}
_PATCH_META["分子格式转换工具"] = {
    "name": "分子格式转换工具",
    "main_category": "化学",
    "sub_category": "数据处理",
    "description": "使用 RDKit 在 SMILES、InChI、InChIKey、MOL 等格式之间互转。",
}

PATCH_ROWS: list[dict[str, str]] = []
for _name, _meta in _PATCH_META.items():
    PATCH_ROWS.append(
        {
            "name": _name,
            "main_category": _meta["main_category"],
            "sub_category": _meta["sub_category"],
            "description": _meta["description"],
            "prompt_template": LAUNCH_CORE_PROMPTS_BY_SKILL_NAME[_name],
        }
    )


def main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument("--dry-run", action="store_true")
    p.add_argument("--no-touch-created-at", action="store_true")
    args = p.parse_args()

    db = SessionLocal()
    try:
        n_up = 0
        for row in PATCH_ROWS:
            q = db.query(Skill).filter(Skill.name == row["name"]).one_or_none()
            if q is None:
                db.add(
                    Skill(
                        name=row["name"],
                        description=row["description"],
                        main_category=row["main_category"],
                        sub_category=row["sub_category"],
                        prompt_template=row["prompt_template"],
                        author_id="system",
                        status="approved",
                        created_at=datetime.utcnow(),
                    )
                )
                print("INSERT:", row["name"])
                n_up += 1
                continue
            print("UPDATE:", row["name"], (q.prompt_template or "")[:50], "...")
            q.main_category = row["main_category"]
            q.sub_category = row["sub_category"]
            q.description = row["description"]
            q.prompt_template = row["prompt_template"]
            if not args.no_touch_created_at:
                q.created_at = datetime.utcnow()
            n_up += 1

        if args.dry_run:
            db.rollback()
            print("Dry-run: rolled back.")
            return 0

        db.commit()
        print(f"OK: committed {n_up} launch skill row(s).")
        return 0
    except Exception as e:
        db.rollback()
        print(f"ERROR: {e}")
        raise
    finally:
        db.close()


if __name__ == "__main__":
    sys.exit(main())
