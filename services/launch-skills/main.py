# -*- coding: utf-8 -*-
"""首发 BaseSkill 隔离运行时：仅含 BLAST / ChEMBL / RDKit / httpx 等轻依赖。"""
from __future__ import annotations

import logging
import os
from typing import Any, Dict, Tuple  # noqa: F401 — Tuple used by _LAUNCH_SKILL_CLASSES

from fastapi import FastAPI
from pydantic import BaseModel, Field

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger("launch-skills")

# 在 worker 内直接执行 execute，禁止再 HTTP 回环
os.environ.setdefault("LAUNCH_SKILLS_IN_WORKER", "1")

app = FastAPI(title="GIBH Launch Skills", version="1.0.0")


class RunRequest(BaseModel):
    tool_id: str = Field(..., min_length=1)
    kwargs: Dict[str, Any] = Field(default_factory=dict)


@app.get("/health")
def health() -> Dict[str, str]:
    return {"status": "ok", "service": "launch-skills"}


_LAUNCH_SKILL_CLASSES: Dict[str, Tuple[str, str]] = {
    "chembl_drug_search": ("gibh_agent.skills.skill_chembl_search", "ChemblDrugSearchSkill"),
    "chembl_similar_molecules": ("gibh_agent.skills.skill_chembl_similarity", "ChemblSimilarMoleculesSkill"),
    "pubchem_smiles_to_cid": ("gibh_agent.skills.skill_pubchem_cid", "PubchemSmilesToCidSkill"),
    "mhc_epitope_search": ("gibh_agent.skills.skill_mhc_epitope", "MhcEpitopeSearchSkill"),
    "chipatlas_experiment_search": ("gibh_agent.skills.skill_chipatlas_experiment", "ChipatlasExperimentSearchSkill"),
    "rdkit_substructure_search": ("gibh_agent.skills.skill_rdkit_substructure", "RdkitSubstructureSearchSkill"),
    "rdkit_3d_mol_render": ("gibh_agent.skills.skill_rdkit_3d_render", "Rdkit3dMolRenderSkill"),
    "rdkit_mol_format_convert": ("gibh_agent.skills.skill_mol_format_convert", "RdkitMolFormatConvertSkill"),
    "nucleotide_sequence_blast": ("gibh_agent.skills.skill_nucleotide_alignment", "NucleotideSequenceAlignmentSkill"),
    "protein_sequence_blast": ("gibh_agent.skills.skill_protein_alignment", "ProteinSequenceAlignmentSkill"),
}


@app.post("/api/launch/run")
def run_launch_skill(req: RunRequest) -> Dict[str, Any]:
    import importlib

    tid = (req.tool_id or "").strip()
    spec = _LAUNCH_SKILL_CLASSES.get(tid)
    if not spec:
        return {"status": "error", "message": f"tool_id {tid!r} 不在首发技能白名单内"}

    mod_name, cls_name = spec
    try:
        mod = importlib.import_module(mod_name)
        cls = getattr(mod, cls_name)
        result = cls().execute(**dict(req.kwargs or {}))
    except Exception as exc:
        logger.exception("launch-skills execute failed: %s", tid)
        return {"status": "error", "message": str(exc)}

    if not isinstance(result, dict):
        return {"status": "success", "data": result, "message": str(result)}
    return result
