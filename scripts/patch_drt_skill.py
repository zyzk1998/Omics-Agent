#!/usr/bin/env python3
"""
一次性手术补丁：UPDATE skills 表中 DRTtools，写入 [Skill_Route: eis_drt_analysis] 等字段。

不插入新行，不放宽 total==0 的安全补种逻辑。执行后前端可推断 is_implemented 并走技能快车道。

用法（项目根目录，与 server 相同 MySQL 环境变量）:
  python scripts/patch_drt_skill.py
  docker exec gibh_v2_api python /app/scripts/patch_drt_skill.py
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

# 保证从任意 cwd 可 import
_ROOT = Path(__file__).resolve().parents[1]
if str(_ROOT) not in sys.path:
    sys.path.insert(0, str(_ROOT))

from sqlalchemy import text

from gibh_agent.db.connection import SessionLocal  # noqa: E402
from gibh_agent.db.models import Skill  # noqa: E402

SKILL_NAME = "DRTtools"

# 与 seed_skills.py ADDITIONAL_CHEMISTRY_SKILLS 中 DRTtools 保持一致（含暗号）
DRT_PROMPT_TEMPLATE = """[Skill_Route: eis_drt_analysis]
您好。我已上传 **电化学阻抗谱(EIS)** 数据，请使用 **弛豫时间分布（DRT）** 方法完成反演与可视化，用于分辨不同弛豫过程（如界面电荷转移、扩散等）。

**输入数据（请在此替换为您的文件）**
- 请在对话中 **上传 CSV**（或脚本支持的文本格式），包含 **frequency（Hz）**、**Z_real**、**Z_imag（Ω）** 等列；列名需与仪器导出一致。
- **占位说明**：若文件尚未上传，请先将「原始 EIS 数据文件」保存后拖入附件区；助手仅使用会话附件解析得到的 **`file_path`**，禁止臆造路径。

**可选参数**
- **regularization_lambda**：正则强度（默认 0.1，对应底层 `--lambda`）。
- **method**：反演/正则方法名（默认 `tikhonov`，与底层脚本一致）。

（助手侧：仅将附件列表中的绝对路径写入工具参数 `file_path`；未上传附件时请先提醒用户上传数据文件。）
"""

DRT_DESCRIPTION = (
    "EIS 阻抗谱 DRT 弛豫时间分析；上传含 frequency/Z_real/Z_imag 的 CSV，输出 JSON 与图。"
)


def main() -> int:
    p = argparse.ArgumentParser(description="Patch DRTtools row in skills table.")
    p.add_argument("--dry-run", action="store_true", help="Print SQL preview only, no commit.")
    args = p.parse_args()

    db = SessionLocal()
    try:
        row = db.query(Skill).filter(Skill.name == SKILL_NAME).one_or_none()
        if row is None:
            print(f"ERROR: no skill named {SKILL_NAME!r}; abort (use seed on empty DB first).")
            return 2

        before = {
            "main_category": row.main_category,
            "sub_category": row.sub_category,
            "description_len": len(row.description or ""),
            "prompt_has_route": "[Skill_Route: eis_drt_analysis]" in (row.prompt_template or ""),
        }
        print("Before:", before)

        row.main_category = "化学"
        row.sub_category = "电化学"
        row.description = DRT_DESCRIPTION
        row.prompt_template = DRT_PROMPT_TEMPLATE

        if args.dry_run:
            db.rollback()
            print("Dry-run: rolled back.")
            return 0

        db.commit()
        db.refresh(row)

        after = {
            "main_category": row.main_category,
            "sub_category": row.sub_category,
            "description_len": len(row.description or ""),
            "prompt_has_route": "[Skill_Route: eis_drt_analysis]" in (row.prompt_template or ""),
        }
        print("After:", after)

        # 二次确认受影响行（只读）
        chk = db.execute(
            text(
                "SELECT id, name, main_category, sub_category, "
                "LOCATE(:needle, COALESCE(prompt_template,'')) AS pos "
                "FROM skills WHERE name = :name LIMIT 1"
            ),
            {"needle": "[Skill_Route: eis_drt_analysis]", "name": SKILL_NAME},
        ).fetchone()
        print("Verify row:", dict(chk._mapping) if chk else None)

        print("OK: patch committed.")
        return 0
    except Exception as e:
        db.rollback()
        print(f"ERROR: {e}")
        raise
    finally:
        db.close()


if __name__ == "__main__":
    sys.exit(main())
