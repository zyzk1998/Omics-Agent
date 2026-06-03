# -*- coding: utf-8 -*-
"""
化学大类 · 输出效果预览（demo_visualization）精修覆盖。

合并顺序：Phase1 → Phase2 → Bulk → 生物医药 → **本模块**。
"""
from __future__ import annotations

from typing import Any, Dict

from gibh_agent.db.skill_detail_demo_visualizations_chemistry_remainder import (
    CHEMISTRY_VIZ,
)

SKILL_DETAILED_SPECS_CHEMISTRY_VIZ: Dict[str, Dict[str, Any]] = {
    _tid: {"demo_visualization": _html} for _tid, _html in CHEMISTRY_VIZ.items()
}
