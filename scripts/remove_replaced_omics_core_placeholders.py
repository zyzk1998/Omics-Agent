#!/usr/bin/env python3
"""
从 skills 表删除已被「三大组学全流程」快车道替换的旧版核心单项技能，避免多模态组学 tab 下出现 6+3=9/10 张重复卡片。

被替换名称（与历史 CORE_OMICS_SKILLS 一致）:
  - 基因组变异检测            → 由「基因组学全流程」承接
  - 表观遗传峰与 motif 分析   → 由「表观组学全流程」承接
  - 蛋白质互作网络分析      → 由「蛋白组学全流程」承接

用法:
  cd /home/ubuntu/GIBH-AGENT-V2 && PYTHONPATH=. python3 scripts/remove_replaced_omics_core_placeholders.py
"""
from __future__ import annotations

import sys
from pathlib import Path

_ROOT = Path(__file__).resolve().parents[1]
if str(_ROOT) not in sys.path:
    sys.path.insert(0, str(_ROOT))

REMOVED_NAMES = (
    "基因组变异检测",
    "表观遗传峰与 motif 分析",
    "蛋白质互作网络分析",
)


def main() -> None:
    from gibh_agent.db.connection import SessionLocal
    from gibh_agent.db.models import Skill

    db = SessionLocal()
    try:
        n = 0
        for name in REMOVED_NAMES:
            row = db.query(Skill).filter(Skill.name == name).first()
            if row:
                db.delete(row)
                n += 1
        db.commit()
        print(f"remove_replaced_omics_core_placeholders: deleted {n} row(s) (looked for {len(REMOVED_NAMES)} names).")
    finally:
        db.close()


if __name__ == "__main__":
    main()
