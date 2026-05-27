# -*- coding: utf-8 -*-
"""首发 10 项核心技能：ToolRegistry 注册与本地可执行冒烟（无网用例可跳过）。"""
from __future__ import annotations

import pytest

LAUNCH_TOOL_IDS = [
    "chembl_drug_search",
    "rdkit_substructure_search",
    "nucleotide_sequence_blast",
    "protein_sequence_blast",
    "pubchem_smiles_to_cid",
    "mhc_epitope_search",
    "chipatlas_experiment_search",
    "rdkit_3d_mol_render",
    "chembl_similar_molecules",
    "rdkit_mol_format_convert",
]

PROMPT_SKILL_NAMES = [
    "ChEMBL药物检索",
    "子结构搜索化合物",
    "核酸序列比对",
    "蛋白质序列比对",
    "通过SMILES获取CID",
    "MHC关联表位检索",
    "ChIPAtlas实验获取工具",
    "3D分子结构渲染工具",
    "检索相似小分子",
    "分子格式转换工具",
]


@pytest.fixture(scope="module", autouse=True)
def _load_skills():
    from gibh_agent.core.skill_registry import discover_and_register_skills

    discover_and_register_skills(force_reload=True)


@pytest.mark.parametrize("tool_id", LAUNCH_TOOL_IDS)
def test_tool_registry_has_launch_skill(tool_id: str):
    from gibh_agent.core.tool_registry import registry

    assert registry.get_tool(tool_id) is not None


def test_launch_prompts_have_skill_route():
    from gibh_agent.db.launch_core_skill_prompt_templates import LAUNCH_CORE_PROMPTS_BY_SKILL_NAME

    for name in PROMPT_SKILL_NAMES:
        tpl = LAUNCH_CORE_PROMPTS_BY_SKILL_NAME.get(name) or ""
        assert "[Skill_Route:" in tpl, name


def test_rdkit_substructure_local():
    pytest.importorskip("rdkit")
    from gibh_agent.skills.skill_rdkit_substructure import RdkitSubstructureSearchSkill

    r = RdkitSubstructureSearchSkill().execute(
        substructure_smiles="c1ccccc1",
        target_smiles="c1ccccc1CCO",
    )
    assert r["status"] == "success"
    assert r["data"]["hit_count"] >= 1


def test_launch_demo_defaults_fill_empty_params():
    from gibh_agent.skills.launch_skill_demos import apply_launch_demo_defaults

    chembl = apply_launch_demo_defaults("chembl_drug_search", {"query": "", "limit": 20})
    assert chembl["query"] == "aspirin"

    mhc = apply_launch_demo_defaults(
        "mhc_epitope_search", {"mhc_allele": "", "peptide_sequence": ""}
    )
    assert mhc["mhc_allele"] == "HLA-A*02:01"


def test_skill_agent_merge_launch_from_prompt_json():
    from gibh_agent.agents.skill_agent import _merge_launch_tool_args
    from gibh_agent.db.launch_core_skill_prompt_templates import LAUNCH_CORE_PROMPTS_BY_SKILL_NAME

    tpl = LAUNCH_CORE_PROMPTS_BY_SKILL_NAME["ChEMBL药物检索"]
    merged = _merge_launch_tool_args("chembl_drug_search", {}, tpl)
    assert merged.get("query") == "aspirin"


def test_mol_format_convert_local():
    pytest.importorskip("rdkit")
    from gibh_agent.skills.skill_mol_format_convert import RdkitMolFormatConvertSkill

    r = RdkitMolFormatConvertSkill().execute(
        input_text="CCO",
        input_format="smiles",
        output_format="inchi",
    )
    assert r["status"] == "success"
    assert r["data"]["result"]


def test_launch_skill_runtime_dependencies():
    """与 scripts/verify_launch_skill_deps.sh 同口径；宿主机无 BLAST 时跳过。"""
    import shutil

    pytest.importorskip("httpx")
    pytest.importorskip("chembl_webresource_client")
    pytest.importorskip("rdkit")
    for cli in ("blastn", "blastp"):
        if not shutil.which(cli):
            pytest.skip(
                f"未找到 {cli}；请在 api-server 容器执行 "
                "scripts/verify_launch_skill_deps.sh（见 docs/技能扩展规范文档.md §5.11.1）"
            )
