# -*- coding: utf-8 -*-
"""ChEMBL 相似分子检索 — 首发核心技能。"""
from __future__ import annotations

from typing import Any, Dict, List
from urllib.parse import quote

import httpx

from gibh_agent.skills._skill_common import DEFAULT_HTTP_TIMEOUT, err, ok
from gibh_agent.skills.base_skill import BaseSkill
from gibh_agent.skills.launch_skill_demos import apply_launch_demo_defaults

CHEMBL_SIMILARITY_API = "https://www.ebi.ac.uk/chembl/api/data/similarity/{smiles}/{threshold}.json"


class ChemblSimilarMoleculesSkill(BaseSkill):
    __abstractskill__ = False

    """
    ChEMBL 相似小分子检索。

    根据 SMILES 在 ChEMBL 网络服务中检索 Tanimoto 相似化合物（仅适用于小分子）。
    也可通过 ChEMBL ID 先解析 SMILES 再检索（需提供 chembl_id）。

    参数:
        smiles: 查询 SMILES（与 chembl_id 二选一）。
        chembl_id: ChEMBL 分子 ID（可选，用于拉取结构后再相似性搜索）。
        similarity_threshold: Tanimoto 阈值 0–100，默认 70。
        limit: 返回条数上限，默认 20。
    """

    skill_id = "chembl_similar_molecules"
    display_name = "检索相似小分子"
    description = (
        "通过 ChEMBL 相似性 API 根据 SMILES 或 ChEMBL ID 检索结构类似的小分子化合物。"
    )
    category = "化学与分子信息学"
    sub_category = "数据分析"
    aliases = ["相似分子", "Tanimoto", "ChEMBL相似", "结构相似"]
    required_parameters = []
    tool_chain_key = "chembl"
    __dependencies__ = ["pip:httpx", "pip:chembl_webresource_client"]

    def _smiles_from_chembl_id(self, chembl_id: str) -> str:
        from chembl_webresource_client.new_client import new_client

        rows = list(
            new_client.molecule.filter(molecule_chembl_id__iexact=chembl_id.upper()).only(
                "molecule_structures"
            )[:1]
        )
        if not rows:
            raise ValueError(f"未找到 ChEMBL ID: {chembl_id}")
        structs = rows[0].get("molecule_structures") or {}
        smi = structs.get("canonical_smiles") if isinstance(structs, dict) else None
        if not smi:
            raise ValueError(f"ChEMBL {chembl_id} 无可用 SMILES")
        return smi

    def execute(
        self,
        smiles: str = "",
        chembl_id: str = "",
        similarity_threshold: int = 70,
        limit: int = 20,
        **kwargs: Any,
    ) -> Dict[str, Any]:
        filled = apply_launch_demo_defaults(
            self.skill_id,
            {
                "smiles": (smiles or "").strip(),
                "chembl_id": (chembl_id or kwargs.get("molecule_chembl_id") or "").strip(),
                "similarity_threshold": similarity_threshold,
                "limit": limit,
            },
        )
        smi = str(filled.get("smiles") or "").strip()
        cid = str(filled.get("chembl_id") or "").strip()

        if not smi and cid:
            try:
                smi = self._smiles_from_chembl_id(cid)
            except ImportError:
                return err("需要 chembl_webresource_client 以根据 ChEMBL ID 解析 SMILES")
            except ValueError as exc:
                return err(str(exc))
        if not smi:
            return err("请提供 smiles 或 chembl_id")

        thresh = max(40, min(int(filled.get("similarity_threshold") or 70), 100))
        lim = max(1, min(int(filled.get("limit") or 20), 50))
        url = CHEMBL_SIMILARITY_API.format(
            smiles=quote(smi, safe=""),
            threshold=thresh,
        )

        try:
            with httpx.Client(timeout=DEFAULT_HTTP_TIMEOUT, follow_redirects=True) as client:
                resp = client.get(url, params={"limit": lim})
                resp.raise_for_status()
                payload = resp.json()
        except httpx.HTTPError as exc:
            return err(f"ChEMBL 相似性 API 失败: {exc}")

        mols = payload.get("molecules") or []
        hits: List[Dict[str, Any]] = []
        for item in mols[:lim]:
            if isinstance(item, dict):
                hits.append(
                    {
                        "molecule_chembl_id": item.get("molecule_chembl_id"),
                        "pref_name": item.get("pref_name"),
                        "similarity": item.get("similarity"),
                    }
                )

        return ok(
            f"ChEMBL 相似性检索完成，返回 {len(hits)} 条（阈值 {thresh}%）",
            query_smiles=smi,
            threshold=thresh,
            hits=hits,
            count=len(hits),
        )
