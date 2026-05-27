# -*- coding: utf-8 -*-
"""
首发 10 项核心技能 — 广场「一键体验」默认参数（与 launch_core_skill_prompt_templates 同源）。

SkillAgent 填参留空时回退；各 skill execute 在必填项为空时二次兜底。
"""
from __future__ import annotations

from typing import Any, Dict

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
