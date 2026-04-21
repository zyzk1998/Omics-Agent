#!/usr/bin/env python3
"""
一次性批量 PATCH：UPDATE（或 INSERT）skills 表中药物化学/RDKit 相关系统技能，
写入 [Skill_Route: ...]、分类「化学」、子类「药物筛选」/「分子计算」。

用法:
  python scripts/patch_chem_skills_batch.py
  python scripts/patch_chem_skills_batch.py --no-touch-created-at   # 不刷新排序时间戳
  docker exec gibh_v2_api python /app/scripts/patch_chem_skills_batch.py
"""
from __future__ import annotations

import argparse
import sys
from datetime import datetime
from pathlib import Path

_ROOT = Path(__file__).resolve().parents[1]
if str(_ROOT) not in sys.path:
    sys.path.insert(0, str(_ROOT))

from sqlalchemy import text  # noqa: E402

from gibh_agent.db.chem_skill_prompt_templates import CHEM_PROMPTS_BY_SKILL_NAME  # noqa: E402
from gibh_agent.db.connection import SessionLocal  # noqa: E402
from gibh_agent.db.models import Skill  # noqa: E402

# 与 seed_skills.py / chem_skill_prompt_templates.py 同源，禁止在此处手写分叉模板
PATCH_ROWS: list[dict[str, str]] = [
    {
        "name": "分子相似性评估工具",
        "main_category": "化学",
        "sub_category": "分子计算",
        "description": "双分子 Morgan 指纹 Tanimoto 相似度（2048-bit）；支持 SMILES 文本或文件。",
        "prompt_template": CHEM_PROMPTS_BY_SKILL_NAME["分子相似性评估工具"],
    },
    {
        "name": "分子假阳性片段检测工具",
        "main_category": "化学",
        "sub_category": "药物筛选",
        "description": "PAINS 泛测定干扰片段检测（RDKit FilterCatalog）。",
        "prompt_template": CHEM_PROMPTS_BY_SKILL_NAME["分子假阳性片段检测工具"],
    },
    {
        "name": "分子Morgan Fingerprint生成工具",
        "main_category": "化学",
        "sub_category": "分子计算",
        "description": "Morgan（ECFP 类）指纹摘要：开位数与哈希预览；附 2D 结构图。",
        "prompt_template": CHEM_PROMPTS_BY_SKILL_NAME["分子Morgan Fingerprint生成工具"],
    },
    {
        "name": "Brenk filter分子毒性检查工具",
        "main_category": "化学",
        "sub_category": "药物筛选",
        "description": "Brenk 不良片段筛查（RDKit BRENK 目录或 SMARTS 兜底）。",
        "prompt_template": CHEM_PROMPTS_BY_SKILL_NAME["Brenk filter分子毒性检查工具"],
    },
    {
        "name": "分子BBB评估工具",
        "main_category": "化学",
        "sub_category": "药物筛选",
        "description": "血脑屏障渗透启发式评估（理化描述符组合，非临床 PBPK）。",
        "prompt_template": CHEM_PROMPTS_BY_SKILL_NAME["分子BBB评估工具"],
    },
    {
        "name": "分子量计算工具",
        "main_category": "化学",
        "sub_category": "分子计算",
        "description": "平均分子量与精确分子量（RDKit Descriptors）。",
        "prompt_template": CHEM_PROMPTS_BY_SKILL_NAME["分子量计算工具"],
    },
]

INSERT_IF_MISSING = [
    {
        "name": "高通量药物相似性检索",
        "main_category": "化学",
        "sub_category": "药物筛选",
        "description": "对接 PubChem/ChEMBL 等的结构相似性检索与 HTML 报告（需 worker-pyskills + 第三方脚本；联网）。",
        "prompt_template": """[Skill_Route: drug_similarity_search]
您好。请提供查询分子的 **SMILES**，或上传含一行 SMILES 的文件；我将发起 **高通量相似性检索**。

**可选**：在后续参数中约定 **top_n**、**similarity_threshold**、目标数据库列表等（由 SkillAgent 从对话解析）。

（助手侧：**file_path** 必填若已上传附件；否则 **smiles_text**。）
""",
    },
]


def main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument("--dry-run", action="store_true")
    p.add_argument(
        "--no-touch-created-at",
        action="store_true",
        help="不更新 skills.created_at（默认会在 PATCH 时刷新，使本批技能在广场「已接入」分组内按最新时间置顶）",
    )
    args = p.parse_args()

    db = SessionLocal()
    try:
        n_up = 0
        for row in PATCH_ROWS:
            q = db.query(Skill).filter(Skill.name == row["name"]).one_or_none()
            if q is None:
                print(f"SKIP (no row): {row['name']}")
                continue
            print("Before:", row["name"], (q.prompt_template or "")[:60], "...")
            q.main_category = row["main_category"]
            q.sub_category = row["sub_category"]
            q.description = row["description"]
            q.prompt_template = row["prompt_template"]
            if not args.no_touch_created_at:
                q.created_at = datetime.utcnow()
            n_up += 1

        for ins in INSERT_IF_MISSING:
            ex = db.query(Skill).filter(Skill.name == ins["name"]).one_or_none()
            if ex:
                print("Exists, skip insert:", ins["name"])
                continue
            db.add(
                Skill(
                    name=ins["name"],
                    description=ins["description"],
                    main_category=ins["main_category"],
                    sub_category=ins["sub_category"],
                    prompt_template=ins["prompt_template"],
                    author_id="system",
                    status="approved",
                )
            )
            n_up += 1
            print("INSERT:", ins["name"])

        if args.dry_run:
            db.rollback()
            print("Dry-run: rolled back.")
            return 0

        db.commit()
        print(f"OK: committed updates/inserts (touched ~{n_up}).")
        return 0
    except Exception as e:
        db.rollback()
        print(f"ERROR: {e}")
        raise
    finally:
        db.close()


if __name__ == "__main__":
    sys.exit(main())
