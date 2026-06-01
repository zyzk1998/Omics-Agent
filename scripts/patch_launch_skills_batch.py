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
# 阵列二 Batch1：卡片名早已在 ADDITIONAL_BIOMED/CHEMISTRY 种子中，但须单独刷新 Prompt（否则仍为占位/旧文案）
for _extra in (
    {
        "name": "PubMed数据库查询",
        "main_category": "生物医药",
        "sub_category": "信息检索",
        "description": "PubMed E-utilities 文献题录检索，表格 + Markdown 摘要。",
    },
    {
        "name": "UniProt数据库查询",
        "main_category": "生物医药",
        "sub_category": "信息检索",
        "description": "UniProt REST 蛋白条目检索，表格 + 指标卡。",
    },
    {
        "name": "分子量计算工具",
        "main_category": "化学",
        "sub_category": "数据分析",
        "description": "RDKit 平均/精确分子量、分子式与 Lipinski 描述符（Metrics Board）。",
    },
    {
        "name": "ClinVar 数据库查询",
        "main_category": "生物医药",
        "sub_category": "信息检索",
        "description": "ClinVar E-utilities 变异临床意义检索，表格 + Markdown。",
    },
    {
        "name": "dbSNP 数据库查询",
        "main_category": "生物医药",
        "sub_category": "信息检索",
        "description": "dbSNP E-utilities SNP 位点与基因关联检索。",
    },
    {
        "name": "GWAS catalog 数据库查询",
        "main_category": "生物医药",
        "sub_category": "信息检索",
        "description": "GWAS Catalog REST 性状—SNP 关联与研究检索。",
    },
    {
        "name": "Reactome 数据库查询",
        "main_category": "生物医药",
        "sub_category": "信息检索",
        "description": "Reactome ContentService 人类通路检索。",
    },
    {
        "name": "InterPro 数据库查询",
        "main_category": "生物医药",
        "sub_category": "信息检索",
        "description": "InterPro REST 蛋白结构域与家族检索。",
    },
    {
        "name": "GEO 数据库查询",
        "main_category": "生物医药",
        "sub_category": "信息检索",
        "description": "NCBI GEO DataSets 元数据检索（不下载表达矩阵）。",
    },
    {
        "name": "基因蛋白信息查询器",
        "main_category": "生物医药",
        "sub_category": "信息检索",
        "description": "Ensembl REST 人类基因注释与转录本概览。",
    },
    {
        "name": "获取mRNA序列工具",
        "main_category": "生物医药",
        "sub_category": "信息检索",
        "description": "Ensembl REST 按基因符号获取代表性 cDNA 序列（FASTA）。",
    },
    {
        "name": "药物标识符交叉检索",
        "main_category": "化学",
        "sub_category": "信息检索",
        "description": "PubChem / ChEMBL / Registry 药物标识符交叉引用宽表。",
    },
    {
        "name": "代谢物名称检索",
        "main_category": "化学",
        "sub_category": "信息检索",
        "description": "RefMet 代谢物标准化名称与分类检索。",
    },
    {
        "name": "RNAcentral 按编号查询",
        "main_category": "生物医药",
        "sub_category": "信息检索",
        "description": "RNAcentral URS 条目与交叉 ID 检索。",
    },
    {
        "name": "Ensembl GO术语后代查询",
        "main_category": "生物医药",
        "sub_category": "信息检索",
        "description": "GO 术语子术语与 QuickGO 后代 ID 列表。",
    },
    {
        "name": "OpenTargets_按ChEMBL ID获取父分子和子分子",
        "main_category": "化学",
        "sub_category": "信息检索",
        "description": "Open Targets GraphQL 药物分子父子层级。",
    },
    {
        "name": "FDA药品标签字段检索",
        "main_category": "化学",
        "sub_category": "信息检索",
        "description": "openFDA 药品标签适应症/警告/用法字段检索。",
    },
):
    _nm = _extra["name"]
    if _nm in LAUNCH_CORE_PROMPTS_BY_SKILL_NAME:
        _PATCH_META[_nm] = _extra

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
