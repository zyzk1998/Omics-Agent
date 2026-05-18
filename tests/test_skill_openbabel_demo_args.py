# -*- coding: utf-8 -*-
"""Open Babel 技能广场一键演示：填参兜底须带 smiles_text + compute_properties。"""
from __future__ import annotations

from gibh_agent.agents.skill_agent import SkillAgent, _OPENBABEL_DEMO_SMILES
from gibh_agent.db.chem_skill_prompt_templates import CHEM_PROMPTS_BY_SKILL_NAME


def test_openbabel_prompt_template_has_demo_smiles() -> None:
    tpl = CHEM_PROMPTS_BY_SKILL_NAME["Open Babel"]
    assert "[Skill_Route: chem_openbabel]" in tpl
    assert _OPENBABEL_DEMO_SMILES in tpl
    assert "compute_properties" in tpl


def _agent() -> SkillAgent:
    return SkillAgent(llm_client=None, format_sse=lambda _t, _d: "")


def test_augment_openbabel_empty_llm_args_uses_demo() -> None:
    agent = _agent()
    prompt = CHEM_PROMPTS_BY_SKILL_NAME["Open Babel"]
    args = agent._augment_args_from_user_text("chem_openbabel", {}, prompt, [])
    assert args.get("smiles_text") == _OPENBABEL_DEMO_SMILES
    assert args.get("compute_properties") is True
    assert not (args.get("file_path") or "").strip()


def test_fallback_openbabel_no_upload() -> None:
    agent = _agent()
    prompt = CHEM_PROMPTS_BY_SKILL_NAME["Open Babel"]
    args = agent._fallback_args("chem_openbabel", prompt, [])
    assert args.get("smiles_text") == _OPENBABEL_DEMO_SMILES
    assert args.get("compute_properties") is True
