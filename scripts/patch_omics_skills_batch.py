#!/usr/bin/env python3
"""
一次性 PATCH：UPDATE skills 表中三大组学编排快车道技能，写入首行 [Omics_Route: ...]。
与 gibh_agent/db/omics_skill_prompt_templates.py 同源，禁止手写分叉模板。

用法:
  cd /home/ubuntu/GIBH-AGENT-V2 && PYTHONPATH=. python3 scripts/patch_omics_skills_batch.py
  docker exec gibh_v2_api python /app/scripts/patch_omics_skills_batch.py
"""
from __future__ import annotations

import sys
from datetime import datetime
from pathlib import Path

_ROOT = Path(__file__).resolve().parents[1]
if str(_ROOT) not in sys.path:
    sys.path.insert(0, str(_ROOT))

from gibh_agent.db.connection import SessionLocal  # noqa: E402
from gibh_agent.db.models import Skill  # noqa: E402
from gibh_agent.db.omics_skill_prompt_templates import OMICS_FAST_LANE_PROMPTS  # noqa: E402

PATCH_ROWS: list[dict[str, str]] = [
    {
        "name": "基因组学全流程",
        "main_category": "多模态组学",
        "sub_category": "基因组学",
        "description": "一键进入基因组胚系 DAG：含 [Omics_Route: genomics]，上传 FASTQ/BAM/VCF 注入首步 file_path。",
        "prompt_template": OMICS_FAST_LANE_PROMPTS["基因组学全流程"],
    },
    {
        "name": "蛋白组学全流程",
        "main_category": "多模态组学",
        "sub_category": "蛋白质组学",
        "description": "一键进入蛋白组搜库定量 DAG：含 [Omics_Route: proteomics]，上传 RAW/mzML。",
        "prompt_template": OMICS_FAST_LANE_PROMPTS["蛋白组学全流程"],
    },
    {
        "name": "表观组学全流程",
        "main_category": "多模态组学",
        "sub_category": "表观遗传组学",
        "description": "一键进入 ATAC/ChIP 主线：含 [Omics_Route: epigenomics]，上传 FASTQ/BAM。",
        "prompt_template": OMICS_FAST_LANE_PROMPTS["表观组学全流程"],
    },
]


def main() -> None:
    db = SessionLocal()
    try:
        now = datetime.utcnow()
        n = 0
        for row in PATCH_ROWS:
            name = row["name"]
            skill = db.query(Skill).filter(Skill.name == name).first()
            if skill:
                skill.description = row["description"]
                skill.main_category = row["main_category"]
                skill.sub_category = row["sub_category"]
                skill.prompt_template = row["prompt_template"]
                n += 1
            else:
                db.add(
                    Skill(
                        name=name,
                        description=row["description"],
                        main_category=row["main_category"],
                        sub_category=row["sub_category"],
                        prompt_template=row["prompt_template"],
                        author_id="system",
                        status="approved",
                        created_at=now,
                    )
                )
                n += 1
        db.commit()
        print(f"patch_omics_skills_batch: upserted {n} rows")
    finally:
        db.close()


if __name__ == "__main__":
    main()
