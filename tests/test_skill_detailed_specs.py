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


def test_unrouted_biopharma_chemistry_have_rich_detailed_spec():
    """无 [Skill_Route] 的占位技能通过 name 映射获得完整说明书与预览。"""
    from gibh_agent.db.seed_skills import get_all_system_skills_list
    from gibh_agent.db.skill_detailed_specs import enrich_skill_plaza_item

    missing = []
    short_viz = []
    for s in get_all_system_skills_list():
        mc = (s.get("main_category") or "").strip()
        if mc not in ("生物医药", "化学"):
            continue
        name = s.get("name") or ""
        pt = s.get("prompt_template") or ""
        if extract_tool_id_from_prompt(pt):
            continue
        item = enrich_skill_plaza_item({"id": 1, "name": name, "prompt_template": pt})
        spec = item.get("detailed_spec")
        if not spec:
            missing.append(name)
            continue
        viz = spec.get("demo_visualization") or ""
        if len(viz) < 900:
            short_viz.append((name, len(viz)))
        assert spec.get("description_long")
        assert spec.get("query_examples")
    assert not missing, f"missing unrouted detailed_spec: {missing}"
    assert not short_viz, f"short demo_visualization: {short_viz[:5]}"


def test_omics_pipeline_specs_seven_flagship():
    from gibh_agent.db.skill_detailed_specs import SKILL_DETAILED_SPECS_BY_TOOL_ID
    from gibh_agent.db.skill_pipeline_specs import (
        OMICS_PIPELINE_SPECS_BY_TOOL_ID,
        SKILL_NAME_TO_PIPELINE_TOOL_ID,
        build_prompt_engineering_guide,
    )

    assert len(SKILL_NAME_TO_PIPELINE_TOOL_ID) == 7
    assert len(OMICS_PIPELINE_SPECS_BY_TOOL_ID) == 7
    for name, tid in SKILL_NAME_TO_PIPELINE_TOOL_ID.items():
        spec = SKILL_DETAILED_SPECS_BY_TOOL_ID[tid]
        assert spec.get("description_long"), name
        assert spec.get("demo_visualization"), name
        assert len(spec.get("demo_visualization") or "") >= 900, name
        qe = spec.get("query_examples") or []
        assert len(qe) == 3, name
        assert all(isinstance(x, dict) and x.get("prompt") for x in qe), name
        tiers = {x.get("tier") for x in qe}
        assert tiers == {"zero_shot", "dynamic_routing", "post_analysis"}, name
        labels = " ".join(x.get("label", "") for x in qe)
        assert "傻瓜" not in labels, name
        assert "一键化标准执行" in qe[0].get("label", ""), name

    from gibh_agent.db.skill_detail_demo_visualizations_omics import (
        OMICS_PIPELINE_REAL_RESULT_TOOL_IDS,
    )

    for tid in OMICS_PIPELINE_REAL_RESULT_TOOL_IDS:
        viz = SKILL_DETAILED_SPECS_BY_TOOL_ID[tid].get("demo_visualization") or ""
        assert "skill-viz-pipeline-gallery" in viz, tid

    scrna = SKILL_DETAILED_SPECS_BY_TOOL_ID["pipeline_transcriptomics"]
    scrna_viz = scrna["demo_visualization"]
    assert "/assets/images/demos/pipelines/transcriptomics/qc_violin.png" in scrna_viz
    assert "marker_dotplot.png" in scrna_viz
    spatial_viz = SKILL_DETAILED_SPECS_BY_TOOL_ID["pipeline_spatial"]["demo_visualization"]
    assert "spatial_multires.png" in spatial_viz
    radio_viz = SKILL_DETAILED_SPECS_BY_TOOL_ID["pipeline_radiomics"]["demo_visualization"]
    assert "original_firstorder_10Percentile" in radio_viz
    assert "Harmony" in scrna["query_examples"][1]["prompt"]
    assert "CD8+ T cells" in scrna["query_examples"][2]["prompt"]
    assert "```mermaid" in scrna["description_long"]
    assert 'Raw["Raw FASTQ"]' in scrna["description_long"]
    spatial = SKILL_DETAILED_SPECS_BY_TOOL_ID["pipeline_spatial"]
    assert "空间解卷积" in spatial["query_examples"][2]["prompt"]
    assert "```mermaid" in spatial["description_long"]
    radiomics = SKILL_DETAILED_SPECS_BY_TOOL_ID["pipeline_radiomics"]
    assert "LASSO" in radiomics["description_long"]
    assert "DimRed{" in radiomics["description_long"]
    for tid in OMICS_PIPELINE_SPECS_BY_TOOL_ID:
        dl = SKILL_DETAILED_SPECS_BY_TOOL_ID[tid]["description_long"]
        assert "```mermaid" in dl, tid
        assert "classDef core" in dl, tid

    guide = build_prompt_engineering_guide("a", "b", "c")
    assert [g["tier"] for g in guide] == ["zero_shot", "dynamic_routing", "post_analysis"]


def test_research_flow_specs_two_flagship():
    from gibh_agent.db.skill_detailed_specs import SKILL_DETAILED_SPECS_BY_TOOL_ID
    from gibh_agent.db.skill_research_flow_specs import (
        RESEARCH_FLOW_SPECS_BY_TOOL_ID,
        SKILL_NAME_TO_RESEARCH_FLOW_TOOL_ID,
    )

    assert len(SKILL_NAME_TO_RESEARCH_FLOW_TOOL_ID) == 2
    for name, tid in SKILL_NAME_TO_RESEARCH_FLOW_TOOL_ID.items():
        spec = SKILL_DETAILED_SPECS_BY_TOOL_ID[tid]
        assert spec.get("description_long"), name
        assert "```mermaid" in spec["description_long"], name
        assert "classDef core" in spec["description_long"], name
        assert len(spec.get("demo_visualization") or "") >= 900, name
        qe = spec.get("query_examples") or []
        assert len(qe) == 3, name
        assert qe[0].get("tier") == "zero_shot", name

    sted = SKILL_DETAILED_SPECS_BY_TOOL_ID["sted_ec_trajectory"]
    assert "Moscot" in sted["description_long"]
    full = SKILL_DETAILED_SPECS_BY_TOOL_ID["spatiotemporal_dynamics_sc"]
    assert "driver_gene" in full["description_long"] or "驱动基因" in full["description_long"]
    assert RESEARCH_FLOW_SPECS_BY_TOOL_ID["spatiotemporal_dynamics"] is full


def test_enrich_omics_pipeline_by_skill_name():
    from gibh_agent.db.skill_detailed_specs import enrich_skill_plaza_item

    item = enrich_skill_plaza_item(
        {
            "id": 99,
            "name": "转录组学标准全流程",
            "main_category": "多模态组学",
            "sub_category": "转录组学",
            "prompt_template": "（编排 Prompt，无 Skill_Route）",
        }
    )
    spec = item.get("detailed_spec") or {}
    assert item.get("tool_name") == "pipeline_transcriptomics"
    assert spec.get("tool_id") == "pipeline_transcriptomics"
    assert "Scrublet" in spec.get("description_long", "")


def test_name_to_tool_id_map_covers_drug_similarity():
    from gibh_agent.db.skill_detailed_specs import (
        SKILL_DETAILED_SPECS_BY_TOOL_ID,
        build_skill_name_to_tool_id_map,
    )

    name_map = build_skill_name_to_tool_id_map()
    tid = name_map.get("药物相似性评估工具")
    assert tid == "lipinski_druglikeness"
    spec = SKILL_DETAILED_SPECS_BY_TOOL_ID.get(tid) or {}
    assert spec.get("demo_visualization")
    assert len(spec.get("demo_visualization") or "") >= 900
