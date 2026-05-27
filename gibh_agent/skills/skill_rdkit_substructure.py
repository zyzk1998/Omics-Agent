# -*- coding: utf-8 -*-
"""RDKit 子结构匹配 — 首发核心技能。"""
from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, List

from gibh_agent.skills._skill_common import err, ok
from gibh_agent.skills.base_skill import BaseSkill
from gibh_agent.skills.launch_skill_demos import apply_launch_demo_defaults


class RdkitSubstructureSearchSkill(BaseSkill):
    __abstractskill__ = False

    """
    RDKit 子结构搜索技能。

    在给定目标分子（或分子库）中检测是否包含指定子结构（SMARTS/SMILES）。
    适用于药效团筛选、骨架匹配与本地化合物库过滤。

    参数:
        substructure_smiles: 子结构 SMILES 或 SMARTS。
        target_smiles: 单个待检分子 SMILES（与 file_path / smiles_list 二选一）。
        smiles_list: 逗号分隔的多个 SMILES，批量筛选命中项。
        file_path: 每行一条 SMILES 的文本文件路径。
        use_smarts: 为 true 时将 substructure_smiles 按 SMARTS 解析（默认 false，按 SMILES）。
    """

    skill_id = "rdkit_substructure_search"
    display_name = "子结构搜索化合物"
    description = (
        "使用 RDKit 在目标分子或 SMILES 列表中执行子结构匹配；"
        "适用于药效团与骨架筛选。"
    )
    category = "化学与分子信息学"
    sub_category = "数据分析"
    aliases = ["子结构搜索", "子结构匹配", "SMARTS", "RDKit子结构"]
    required_parameters = ["substructure_smiles"]
    tool_chain_key = ""
    __dependencies__ = ["pip:rdkit-pypi"]

    def execute(
        self,
        substructure_smiles: str = "",
        target_smiles: str = "",
        smiles_list: str = "",
        file_path: str = "",
        use_smarts: bool = False,
        **kwargs: Any,
    ) -> Dict[str, Any]:
        demo = apply_launch_demo_defaults(
            self.skill_id,
            {
                "substructure_smiles": (substructure_smiles or kwargs.get("query_smiles") or "").strip(),
                "target_smiles": (target_smiles or "").strip(),
                "smiles_list": (smiles_list or "").strip(),
                "file_path": (file_path or "").strip(),
            },
        )
        pattern = str(demo.get("substructure_smiles") or "").strip()
        if not pattern:
            return err("请提供 substructure_smiles（子结构 SMILES/SMARTS）")

        targets: List[str] = []
        target_smiles = str(demo.get("target_smiles") or "")
        smiles_list = str(demo.get("smiles_list") or "")
        file_path = str(demo.get("file_path") or "")
        if target_smiles.strip():
            targets.append(target_smiles.strip())
        if smiles_list.strip():
            targets.extend(s.strip() for s in smiles_list.split(",") if s.strip())
        fp = (file_path or "").strip()
        if fp:
            p = Path(fp)
            if not p.is_file():
                return err(f"SMILES 列表文件不存在: {fp}")
            for line in p.read_text(encoding="utf-8", errors="replace").splitlines():
                line = line.strip()
                if line and not line.startswith("#"):
                    targets.append(line.split()[0])

        if not targets:
            return err("请提供 target_smiles、smiles_list 或 file_path 之一作为待检分子库")

        try:
            from rdkit import Chem
        except ImportError:
            return err("未安装 RDKit（pip install rdkit-pypi）")

        if use_smarts:
            query_mol = Chem.MolFromSmarts(pattern)
        else:
            query_mol = Chem.MolFromSmiles(pattern)
            if query_mol is None:
                query_mol = Chem.MolFromSmarts(pattern)
        if query_mol is None:
            return err("无法解析子结构 SMILES/SMARTS")

        hits: List[Dict[str, Any]] = []
        for idx, smi in enumerate(targets):
            mol = Chem.MolFromSmiles(smi)
            if mol is None:
                continue
            if mol.HasSubstructMatch(query_mol):
                hits.append({"index": idx, "smiles": smi})

        return ok(
            f"子结构筛选完成：{len(hits)}/{len(targets)} 条命中",
            substructure=pattern,
            total=len(targets),
            hit_count=len(hits),
            hits=hits,
        )
