# -*- coding: utf-8 -*-
"""
首发核心技能 — 广场「一键体验」默认参数（化学/BLAST 10 项 + Prompt 软技能 10 项）。

- ``LAUNCH_ISOLATED_TOOL_IDS``：仅这些 tool_id 经 HTTP 委托 ``launch-skills``（BLAST/RDKit/ChEMBL）。
- ``LAUNCH_SKILL_DEMO_ARGS``：含 Prompt 软技能；演示默认值在 api-server 本地执行，**不**进 launch-skills 白名单。

SkillAgent 填参留空时回退；各 skill execute 在必填项为空时二次兜底。
"""
from __future__ import annotations

from typing import Any, Dict, FrozenSet

# 仅隔离栈（services/launch-skills/main.py _LAUNCH_SKILL_CLASSES）内执行的 tool_id
LAUNCH_ISOLATED_TOOL_IDS: FrozenSet[str] = frozenset(
    {
        "chembl_drug_search",
        "chembl_similar_molecules",
        "pubchem_smiles_to_cid",
        "mhc_epitope_search",
        "chipatlas_experiment_search",
        "rdkit_substructure_search",
        "rdkit_3d_mol_render",
        "rdkit_mol_format_convert",
        "nucleotide_sequence_blast",
        "protein_sequence_blast",
    }
)

# tool_id → 演示用 kwargs（须与 @registry / BaseSkill.execute 形参一致）
LAUNCH_SKILL_DEMO_ARGS: Dict[str, Dict[str, Any]] = {
    "chembl_drug_search": {
        "query": "aspirin",
        "search_type": "drug",
        "limit": 10,
    },
    "rdkit_substructure_search": {
        "substructure_smiles": "c1ccccc1",
        "target_smiles": "c1ccccc1CCO",
    },
    "nucleotide_sequence_blast": {
        "sequence_text": "ATGCGTACGTTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAG",
        "use_remote": True,
        "database": "core_nt",
        "blast_task": "blastn-short",
        "max_target_seqs": 5,
    },
    "protein_sequence_blast": {
        "sequence_text": "MKTIIALSYIFCLVFA",
        "use_remote": True,
        "database": "swissprot",
        "blast_task": "blastp-short",
        "max_target_seqs": 5,
    },
    "pubchem_smiles_to_cid": {
        "smiles": "CC(=O)Oc1ccccc1C(=O)O",
        "max_cids": 10,
    },
    "mhc_epitope_search": {
        "mhc_allele": "HLA-A*02:01",
        "limit": 10,
    },
    "chipatlas_experiment_search": {
        "expid": "SRX018625",
        "limit": 10,
    },
    "rdkit_3d_mol_render": {
        "smiles": "CC(=O)Oc1ccccc1C(=O)O",
        "output_format": "mol_block",
    },
    "chembl_similar_molecules": {
        "smiles": "CC(=O)Oc1ccccc1C(=O)O",
        "similarity_threshold": 70,
        "limit": 10,
    },
    "rdkit_mol_format_convert": {
        "input_text": "CC(=O)Oc1ccccc1C(=O)O",
        "input_format": "smiles",
        "output_format": "inchi",
    },
    "ppt_outline": {
        "user_request": (
            "单细胞时空动力学分析：研究背景、最优传输轨迹推断、"
            "驱动基因与通路富集、核心结论与展望"
        ),
        "context": "",
    },
    "mindmap_gen": {
        "user_request": (
            "单细胞时空动力学：数据校验→时序标准化→轨迹推断→"
            "驱动基因→通路富集→专家解读"
        ),
        "context": "",
    },
    "weekly_report_writer": {
        "user_request": (
            "周期 2026-05-19 至 2026-05-25；按 SKILL 工作流起草周报，"
            "含个人状态同步与团队对外同步"
        ),
        "context": "上一份报告：待用户补充路径；重点关注项目进展与 🟥/🟨 待办继承",
    },
    "tailored_resume": {
        "user_request": "Senior Data Analyst；SQL/Python/可视化/A/B 测试；医疗行业优先",
        "context": "5 年 RetailCo 数据分析；Tableau/Power BI；HealthPlus 实习；Business Analytics 学士",
    },
    "academic_poster_generator": {
        "user_request": "单细胞轨迹推断论文 → 会议海报故事板（四段式 + 至少 3 张配图方案）",
        "context": "摘要：最优传输轨迹推断；结果含 UMAP 与驱动基因（见用户附件）",
    },
    "academic_abstract_refiner": {
        "user_request": "将 context 中长文精炼为中英双语单段 SCI 风格摘要",
        "context": "请用户粘贴待精炼原文；禁止编造未出现的实验数据。",
    },
    "email_manager": {
        "user_request": "write email to client about project delay — professional, milestones + revised date",
        "context": "",
    },
    "pdf_extractor": {
        "user_request": "提取 PDF 章节大纲、表格摘要与关键结论",
        "context": "",
        "file_path": "",
    },
    "deep_research": {
        "user_request": "间歇性禁食对代谢健康的影响：益处、风险与适用人群",
        "context": "侧重近五年综述与临床指南；标注证据等级",
    },
    "blueprint_drafter": {
        "user_request": (
            "GIBH 技能快车道：用户 → 技能广场 → SkillAgent → ToolRegistry → "
            "ppt_outline/mindmap_gen 等工具 → 右栏 Markdown/Mermaid/HTML 可视化"
        ),
        "context": "扁平工程图；中文节点标签；含 launch-skills 委托支路",
    },
}


def _is_blank(val: Any) -> bool:
    if val is None:
        return True
    if isinstance(val, str):
        return not val.strip()
    return False


def apply_launch_demo_defaults(skill_id: str, params: Dict[str, Any]) -> Dict[str, Any]:
    """将演示默认值填入仍为空的字段（不覆盖用户已提供的非空值）。"""
    demo = LAUNCH_SKILL_DEMO_ARGS.get((skill_id or "").strip()) or {}
    out = dict(params)
    for key, default in demo.items():
        if _is_blank(out.get(key)):
            out[key] = default
    return out


def get_demo_arg(skill_id: str, param: str, default: Any = "") -> Any:
    demo = LAUNCH_SKILL_DEMO_ARGS.get((skill_id or "").strip()) or {}
    return demo.get(param, default)
