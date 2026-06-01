# -*- coding: utf-8 -*-
"""分子量计算工具 — RDKit Descriptors（Metrics Board）。"""
from __future__ import annotations

from typing import Any, Dict

from gibh_agent.skills._skill_payload import (
    error_result,
    metrics_cards_from_pairs,
    success_payload,
)
from gibh_agent.skills.base_skill import BaseSkill
from gibh_agent.skills.launch_skill_demos import apply_launch_demo_defaults


class CalcMolecularWeightSkill(BaseSkill):
    __abstractskill__ = False

    skill_id = "calc_molecular_weight"
    display_name = "分子量计算工具"
    description = "由 SMILES 计算平均分子量、精确分子量、分子式与 Lipinski 相关描述符（RDKit）。"
    category = "化学与分子信息学"
    sub_category = "数据分析"
    aliases = ["分子量", "MolWt", "精确质量", "mol_weight", "RDKit分子量"]
    required_parameters = ["smiles"]
    tool_chain_key = ""
    output_type = "mixed"
    __dependencies__ = ["pip:rdkit-pypi"]

    def execute(self, smiles: str = "", **kwargs: Any) -> Dict[str, Any]:
        filled = apply_launch_demo_defaults(
            self.skill_id,
            {"smiles": (smiles or kwargs.get("smiles_text") or "").strip()},
        )
        smi = str(filled.get("smiles") or "").strip()
        if not smi:
            return error_result("请提供 smiles")

        try:
            from rdkit import Chem
            from rdkit.Chem import Descriptors, Lipinski, rdMolDescriptors
        except ImportError:
            return error_result("未安装 RDKit（rdkit-pypi）")

        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            return error_result("无法解析 SMILES，请检查结构式")

        avg_mw = round(Descriptors.MolWt(mol), 4)
        exact_mw = round(rdMolDescriptors.CalcExactMolWt(mol), 4)
        formula = rdMolDescriptors.CalcMolFormula(mol)
        logp = round(Descriptors.MolLogP(mol), 3)
        hbd = Lipinski.NumHDonors(mol)
        hba = Lipinski.NumHAcceptors(mol)
        rot = Lipinski.NumRotatableBonds(mol)

        cards = metrics_cards_from_pairs(
            [
                ("平均分子量", avg_mw, "Da"),
                ("精确分子量", exact_mw, "Da"),
                ("分子式", formula, ""),
                ("LogP", logp, ""),
                ("氢键供体", hbd, ""),
                ("氢键受体", hba, ""),
                ("可旋转键", rot, ""),
            ]
        )
        md = (
            f"### 分子量与理化描述符\n\n"
            f"**SMILES**：`{smi}`\n\n"
            f"| 指标 | 数值 |\n|------|------|\n"
            f"| 平均分子量 (Da) | {avg_mw} |\n"
            f"| 精确分子量 (Da) | {exact_mw} |\n"
            f"| 分子式 | {formula} |\n"
            f"| LogP | {logp} |\n"
            f"| HBD / HBA | {hbd} / {hba} |\n"
            f"| 可旋转键 | {rot} |\n"
        )

        return success_payload(
            f"分子量计算完成（{formula}）",
            markdown=md,
            metrics_cards=cards,
            smiles=smi,
            average_molecular_weight=avg_mw,
            exact_molecular_weight=exact_mw,
            molecular_formula=formula,
            logp=logp,
        )
