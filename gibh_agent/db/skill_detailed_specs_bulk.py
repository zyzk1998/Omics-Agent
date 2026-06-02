# -*- coding: utf-8 -*-
"""除多模态组学外、批量生成的技能 detailed_spec（由 seed 元数据自动构建）。"""
from __future__ import annotations

from gibh_agent.db.skill_detailed_spec_builder import build_bulk_detailed_specs
from gibh_agent.db.skill_detailed_specs_phase2 import SKILL_DETAILED_SPECS_PHASE2

# Phase 1 手工精修条目（勿被批量模板覆盖）
_PHASE1_TOOL_IDS = {
    "ppt_outline",
    "mindmap_gen",
    "diff_expr_interpreter",
    "weekly_report_writer",
    "deep_research",
    "blueprint_drafter",
    "go_kegg_narrative",
    "single_cell_checklist",
    "academic_abstract_refiner",
    "chembl_drug_search",
}

_SKIP_KEYS = set(SKILL_DETAILED_SPECS_PHASE2.keys()) | _PHASE1_TOOL_IDS

SKILL_DETAILED_SPECS_BULK: dict = build_bulk_detailed_specs(
    exclude_main_categories=("多模态组学",),
    skip_existing_keys=_SKIP_KEYS,
)
