#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
一次性 PATCH：将 7 条组学旗舰管线的 detailed_spec 写回 skills 表（与代码注册表同源）。

用法:
  cd /home/ubuntu/GIBH-AGENT-V2 && PYTHONPATH=. python3 scripts/patch_omics_pipeline_detailed_specs.py
  docker exec gibh_v2_api python /app/scripts/patch_omics_pipeline_detailed_specs.py
"""
from __future__ import annotations

import sys
from pathlib import Path

_ROOT = Path(__file__).resolve().parents[1]
if str(_ROOT) not in sys.path:
    sys.path.insert(0, str(_ROOT))

from gibh_agent.db.connection import SessionLocal  # noqa: E402
from gibh_agent.db.models import Skill  # noqa: E402
from gibh_agent.db.skill_pipeline_specs import (  # noqa: E402
    OMICS_PIPELINE_SPECS_BY_TOOL_ID,
    SKILL_NAME_TO_PIPELINE_TOOL_ID,
)


def main() -> None:
    name_by_tid = {v: k for k, v in SKILL_NAME_TO_PIPELINE_TOOL_ID.items()}
    db = SessionLocal()
    try:
        n = 0
        for tid, spec in OMICS_PIPELINE_SPECS_BY_TOOL_ID.items():
            name = name_by_tid.get(tid)
            if not name:
                continue
            skill = db.query(Skill).filter(Skill.name == name).first()
            if not skill:
                print(f"skip missing skill: {name} ({tid})")
                continue
            skill.detailed_spec = dict(spec)
            n += 1
        db.commit()
        print(f"patch_omics_pipeline_detailed_specs: updated {n} rows")
    finally:
        db.close()


if __name__ == "__main__":
    main()
