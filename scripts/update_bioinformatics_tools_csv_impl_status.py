#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""为 docs/生物信息分析工具.csv 追加/刷新最右列「实现情况」。

判定「已实现」：技能名与广场已落地条目一致（skill display_name、
PROMPT_SOFT_SKILL_TEMPLATES、LAUNCH_CORE_PROMPTS_BY_SKILL_NAME 等），
或下列 CSV 别名映射到已落地 display_name。
"""
from __future__ import annotations

import csv
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
CSV_PATH = ROOT / "docs" / "生物信息分析工具.csv"

ALIASES = {
    "核酸序列两两比对": "核酸序列比对",
    "蛋白质序列两两比对": "蛋白质序列比对",
}


def _collect_implemented_names() -> set[str]:
    names: set[str] = set()
    for p in (ROOT / "gibh_agent" / "skills").glob("skill_*.py"):
        text = p.read_text(encoding="utf-8", errors="ignore")
        m = re.search(r'display_name\s*=\s*["\']([^"\']+)["\']', text)
        if m:
            names.add(m.group(1).strip())

    from gibh_agent.db.prompt_soft_skill_templates import (
        PROMPT_SOFT_SKILL_SEEDS,
        PROMPT_SOFT_SKILL_TEMPLATES,
    )

    names.update(PROMPT_SOFT_SKILL_TEMPLATES.keys())
    for row in PROMPT_SOFT_SKILL_SEEDS:
        names.add(row["name"].strip())

    from gibh_agent.db.launch_core_skill_prompt_templates import (
        LAUNCH_CORE_PROMPTS_BY_SKILL_NAME,
        LAUNCH_CORE_SKILL_SEEDS,
    )

    names.update(LAUNCH_CORE_PROMPTS_BY_SKILL_NAME.keys())
    for row in LAUNCH_CORE_SKILL_SEEDS:
        names.add(row["name"].strip())

    try:
        from gibh_agent.db.chem_skill_prompt_templates import CHEM_PROMPT_TEMPLATES_BY_SKILL_NAME

        names.update(CHEM_PROMPT_TEMPLATES_BY_SKILL_NAME.keys())
    except ImportError:
        pass

    return names


def _is_implemented(csv_name: str, implemented: set[str]) -> bool:
    n = (csv_name or "").strip()
    if n in implemented:
        return True
    alias = ALIASES.get(n)
    return bool(alias and alias in implemented)


def main() -> int:
    if not CSV_PATH.is_file():
        print(f"缺少文件: {CSV_PATH}", file=sys.stderr)
        return 1

    sys.path.insert(0, str(ROOT))
    implemented = _collect_implemented_names()

    rows = list(csv.DictReader(CSV_PATH.open(encoding="utf-8-sig")))
    if not rows:
        print("CSV 为空", file=sys.stderr)
        return 1

    fieldnames = [k for k in rows[0].keys() if k != "实现情况"] + ["实现情况"]
    for row in rows:
        row["实现情况"] = "已实现" if _is_implemented(row.get("技能名称", ""), implemented) else "未实现"

    with CSV_PATH.open("w", encoding="utf-8-sig", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)

    done = sum(1 for r in rows if r["实现情况"] == "已实现")
    print(f"已写入 {CSV_PATH.name}：共 {len(rows)} 行，已实现 {done} 行，未实现 {len(rows) - done} 行")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
