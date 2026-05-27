# -*- coding: utf-8 -*-
"""RDKit 3D 分子结构生成 — 首发核心技能。"""
from __future__ import annotations

import os
from pathlib import Path
from typing import Any, Dict

from gibh_agent.skills._skill_common import err, ok
from gibh_agent.skills.base_skill import BaseSkill
from gibh_agent.skills.launch_skill_demos import apply_launch_demo_defaults


class Rdkit3dMolRenderSkill(BaseSkill):
    __abstractskill__ = False

    """
    3D 分子结构渲染/生成技能（RDKit）。

    由 SMILES 嵌入 3D 坐标并输出 MOL 块（SDF 文本）或可选 PNG 图像。
    适用于化学报告中的三维构象展示（轻量本地计算）。

    参数:
        smiles: 输入 SMILES（必填）。
        output_format: ``mol_block``（默认）或 ``sdf`` 或 ``png``。
        file_path: 当 output_format 为 png 时，可选输出 PNG 路径；默认写入临时目录。
        optimize: 是否执行 MMFF 优化，默认 true。
    """

    skill_id = "rdkit_3d_mol_render"
    display_name = "3D分子结构渲染工具"
    description = "使用 RDKit 由 SMILES 生成 3D 构象（MOL/SDF 文本或 PNG 图像）。"
    category = "化学与分子信息学"
    sub_category = "数据分析"
    aliases = ["3D分子", "分子3D", "RDKit3D", "molPreview"]
    required_parameters = ["smiles"]
    tool_chain_key = ""
    __dependencies__ = ["pip:rdkit-pypi"]

    def execute(
        self,
        smiles: str = "",
        output_format: str = "mol_block",
        file_path: str = "",
        optimize: bool = True,
        **kwargs: Any,
    ) -> Dict[str, Any]:
        filled = apply_launch_demo_defaults(
            self.skill_id,
            {"smiles": (smiles or "").strip(), "output_format": (output_format or "mol_block").strip()},
        )
        smi = str(filled.get("smiles") or "").strip()
        if not smi:
            return err("请提供 smiles")
        output_format = str(filled.get("output_format") or "mol_block")

        try:
            from rdkit import Chem
            from rdkit.Chem import AllChem
        except ImportError:
            return err("未安装 RDKit")

        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            return err("无法解析 SMILES")

        mol = Chem.AddHs(mol)
        params = AllChem.ETKDGv3()
        embed_code = AllChem.EmbedMolecule(mol, params)
        if embed_code != 0:
            return err("3D 坐标嵌入失败（构象生成未收敛）")

        if optimize:
            try:
                AllChem.MMFFOptimizeMolecule(mol)
            except Exception:
                pass

        fmt = (output_format or "mol_block").strip().lower()
        if fmt == "png":
            try:
                from rdkit.Chem import Draw
            except ImportError:
                return err("RDKit Draw 模块不可用")
            out = (file_path or "").strip()
            if not out:
                import tempfile

                fd, out = tempfile.mkstemp(suffix=".png", prefix="mol3d_")
                os.close(fd)
            Draw.MolToFile(mol, out, size=(400, 400))
            return ok("3D 分子 PNG 已生成", smiles=smi, image_path=out, format="png")

        if fmt == "sdf":
            block = Chem.MolToMolBlock(mol)
            sdf_text = block + "\n$$$$\n"
            return ok("3D SDF 文本已生成", smiles=smi, sdf_text=sdf_text, format="sdf")

        block = Chem.MolToMolBlock(mol)
        return ok("3D MOL 块已生成", smiles=smi, mol_block=block, format="mol_block")
