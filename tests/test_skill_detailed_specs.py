# -*- coding: utf-8 -*-
"""技能广场 detailed_spec 解析与首批 10 项注册表。"""
from __future__ import annotations

from gibh_agent.core.skill_detailed_spec_utils import build_skill_plaza_payload
from gibh_agent.db.skill_detailed_specs import (
    SKILL_DETAILED_SPECS_BY_TOOL_ID,
    extract_tool_id_from_prompt,
    resolve_detailed_spec,
)


def test_phase1_has_ten_core_specs():
    assert len(SKILL_DETAILED_SPECS_BY_TOOL_ID) >= 20
    for tid in (
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
    ):
        spec = SKILL_DETAILED_SPECS_BY_TOOL_ID[tid]
        assert spec.get("description_long")
        assert spec.get("parameters_table")
        assert spec.get("query_examples")
        assert spec.get("demo_visualization"), f"missing demo_visualization for {tid}"


def test_bulk_covers_non_omics_routed_skills():
    from gibh_agent.db.seed_skills import get_all_system_skills_list
    from gibh_agent.db.skill_detailed_specs import SKILL_DETAILED_SPECS_BY_TOOL_ID, extract_tool_id_from_prompt

    missing = []
    for s in get_all_system_skills_list():
        mc = (s.get("main_category") or "").strip()
        if mc == "多模态组学":
            continue
        tid = extract_tool_id_from_prompt(s.get("prompt_template") or "")
        if not tid:
            continue
        if tid not in SKILL_DETAILED_SPECS_BY_TOOL_ID:
            missing.append(tid)
    assert not missing, f"missing detailed_spec: {missing[:10]}"


def test_omics_category_excluded_from_bulk_requirement():
    """多模态组学 seed 项无 [Skill_Route]，批量生成器亦排除该大类。"""
    from gibh_agent.db.seed_skills import get_all_system_skills_list
    from gibh_agent.db.skill_detailed_specs import extract_tool_id_from_prompt

    omics_routes = []
    for s in get_all_system_skills_list():
        if (s.get("main_category") or "").strip() != "多模态组学":
            continue
        tid = extract_tool_id_from_prompt(s.get("prompt_template") or "")
        if tid:
            omics_routes.append(tid)
    assert omics_routes == []


def test_extract_tool_id_from_prompt():
    pt = "[Skill_Route: ppt_outline]\nhello"
    assert extract_tool_id_from_prompt(pt) == "ppt_outline"


def test_build_skill_plaza_payload_attaches_spec():
    item = build_skill_plaza_payload(
        skill_id=1,
        name="科研汇报 PPT 大纲生成",
        description="短描述",
        main_category="其他技能",
        sub_category="文本处理",
        prompt_template="[Skill_Route: ppt_outline]\n```json\n{}\n```",
        is_implemented=True,
        author_id="system",
        created_at=None,
        saved=False,
    )
    assert item.get("tool_name") == "ppt_outline"
    assert item.get("detailed_spec", {}).get("tool_id") == "ppt_outline"


def test_db_override_wins():
    custom = {"tool_id": "ppt_outline", "description_long": "自定义长描述"}
    resolved = resolve_detailed_spec(
        "[Skill_Route: ppt_outline]",
        db_detailed_spec=custom,
    )
    assert resolved["description_long"] == "自定义长描述"
