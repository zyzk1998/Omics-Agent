# -*- coding: utf-8 -*-
"""RDKit 分子格式转换 — 首发核心技能。"""
from __future__ import annotations

from typing import Any, Dict

from gibh_agent.skills._skill_common import err, ok
from gibh_agent.skills.base_skill import BaseSkill
from gibh_agent.skills.launch_skill_demos import apply_launch_demo_defaults


class RdkitMolFormatConvertSkill(BaseSkill):
    __abstractskill__ = False

    """
    分子格式转换技能（RDKit）。

    在 SMILES、InChI、InChIKey、MOL 块等表示之间互转，并可选标准化分子。

    参数:
        input_text: 输入结构文本（SMILES、InChI 或 MOL 块）。
        input_format: ``smiles`` | ``inchi`` | ``mol``（默认 smiles）。
        output_format: ``smiles`` | ``inchi`` | ``inchikey`` | ``mol``（默认 inchi）。
        canonicalize: 是否输出规范 SMILES，默认 true（仅 output_format=smiles 时生效）。
    """

    skill_id = "rdkit_mol_format_convert"
    display_name = "分子格式转换工具"
    description = "使用 RDKit 在 SMILES、InChI、InChIKey、MOL 等化学结构格式之间转换。"
    category = "化学与分子信息学"
    sub_category = "数据处理"
    aliases = ["格式转换", "SMILES转InChI", "rdkitTransform", "分子格式"]
    required_parameters = ["input_text"]
    tool_chain_key = "rdkitTransform"
    __dependencies__ = ["pip:rdkit-pypi"]

    def execute(
        self,
        input_text: str = "",
        input_format: str = "smiles",
        output_format: str = "inchi",
        canonicalize: bool = True,
        **kwargs: Any,
    ) -> Dict[str, Any]:
        filled = apply_launch_demo_defaults(
            self.skill_id,
            {
                "input_text": (input_text or kwargs.get("smiles") or "").strip(),
                "input_format": (input_format or "smiles").strip(),
                "output_format": (output_format or "inchi").strip(),
            },
        )
        raw = str(filled.get("input_text") or "").strip()
        if not raw:
            return err("请提供 input_text（结构文本）")

        inp_fmt = str(filled.get("input_format") or "smiles").strip().lower()
        out_fmt = str(filled.get("output_format") or "inchi").strip().lower()

        try:
            from rdkit import Chem
        except ImportError:
            return err("未安装 RDKit")

        mol = None
        if inp_fmt == "smiles":
            mol = Chem.MolFromSmiles(raw)
        elif inp_fmt == "inchi":
            mol = Chem.MolFromInchi(raw)
        elif inp_fmt == "mol":
            mol = Chem.MolFromMolBlock(raw)
        else:
            return err(f"不支持的 input_format: {inp_fmt}")

        if mol is None:
            return err("无法解析输入结构")

        if out_fmt == "smiles":
            out = Chem.MolToSmiles(mol, canonical=bool(canonicalize))
        elif out_fmt == "inchi":
            out = Chem.MolToInchi(mol)
        elif out_fmt == "inchikey":
            out = Chem.MolToInchiKey(mol)
        elif out_fmt == "mol":
            out = Chem.MolToMolBlock(mol)
        else:
            return err(f"不支持的 output_format: {out_fmt}")

        return ok(
            f"格式转换完成：{inp_fmt} → {out_fmt}",
            input_format=inp_fmt,
            output_format=out_fmt,
            result=out,
        )
