# -*- coding: utf-8 -*-
"""ChEMBL 药物/分子检索 — 首发核心技能。"""
from __future__ import annotations

from typing import Any, Dict, List

from gibh_agent.skills._skill_common import err, ok
from gibh_agent.skills.base_skill import BaseSkill
from gibh_agent.skills.launch_skill_demos import apply_launch_demo_defaults


class ChemblDrugSearchSkill(BaseSkill):
    __abstractskill__ = False

    """
    ChEMBL 药物检索技能。

    根据药物通用名、同义词或 ChEMBL ID 在 ChEMBL 数据库中检索药物与关联分子记录。
    适用于药物发现、靶点调研、化合物标识符解析等场景。

    参数:
        query: 检索关键词（药物名称片段或 ChEMBL ID，如 CHEMBL25）。
        search_type: ``drug``（默认，查 drug 表）或 ``molecule``（查 molecule 表）。
        limit: 返回条数上限，默认 20，最大 50。
    """

    skill_id = "chembl_drug_search"
    display_name = "ChEMBL药物检索"
    description = (
        "根据药物名称、同义词或 ChEMBL ID 检索 ChEMBL 药物与分子信息；"
        "适用于药物发现与化合物标识符解析。"
    )
    category = "化学与分子信息学"
    sub_category = "数据分析"
    aliases = ["ChEMBL", "药物检索", "chembl", "ChEMBL药物检索"]
    required_parameters = ["query"]
    tool_chain_key = "chembl"
    __dependencies__ = ["pip:chembl_webresource_client"]

    def execute(
        self,
        query: str = "",
        search_type: str = "drug",
        limit: int = 20,
        **kwargs: Any,
    ) -> Dict[str, Any]:
        filled = apply_launch_demo_defaults(
            self.skill_id,
            {
                "query": (query or kwargs.get("drug_name") or "").strip(),
                "search_type": (search_type or "drug").strip(),
                "limit": limit,
            },
        )
        q = str(filled.get("query") or "").strip()
        if not q:
            return err("请提供 query（药物名称或 ChEMBL ID）")

        lim = max(1, min(int(filled.get("limit") or 20), 50))
        st = str(filled.get("search_type") or "drug").strip().lower()

        try:
            from chembl_webresource_client.new_client import new_client
        except ImportError:
            return err("未安装 chembl_webresource_client，请在镜像中 pip install chembl_webresource_client")

        try:
            if st == "molecule" or q.upper().startswith("CHEMBL"):
                resource = new_client.molecule
                if q.upper().startswith("CHEMBL"):
                    rows = list(
                        resource.filter(molecule_chembl_id__iexact=q.upper())
                        .only(
                            "molecule_chembl_id",
                            "pref_name",
                            "max_phase",
                            "molecule_type",
                        )[:lim]
                    )
                else:
                    rows = list(
                        resource.filter(pref_name__icontains=q)
                        .only(
                            "molecule_chembl_id",
                            "pref_name",
                            "max_phase",
                            "molecule_type",
                        )[:lim]
                    )
            else:
                resource = new_client.drug
                if q.upper().startswith("CHEMBL"):
                    rows = list(
                        resource.filter(molecule_chembl_id__iexact=q.upper())
                        .only(
                            "molecule_chembl_id",
                            "pref_name",
                            "drug_type",
                            "max_phase",
                        )[:lim]
                    )
                else:
                    rows = list(
                        resource.filter(pref_name__icontains=q)
                        .only(
                            "molecule_chembl_id",
                            "pref_name",
                            "drug_type",
                            "max_phase",
                        )[:lim]
                    )
        except Exception as exc:
            return err(f"ChEMBL API 调用失败: {exc}")

        hits: List[Dict[str, Any]] = [dict(r) for r in rows]
        if not hits:
            return ok(f"未找到与「{q}」匹配的 ChEMBL 记录", query=q, hits=[], count=0)

        return ok(
            f"ChEMBL 检索完成，共 {len(hits)} 条",
            query=q,
            search_type=st,
            hits=hits,
            count=len(hits),
        )
