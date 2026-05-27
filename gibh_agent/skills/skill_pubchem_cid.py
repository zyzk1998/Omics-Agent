# -*- coding: utf-8 -*-
"""PubChem SMILES → CID — 首发核心技能。"""
from __future__ import annotations

from typing import Any, Dict, List

import httpx

from gibh_agent.skills._skill_common import DEFAULT_HTTP_TIMEOUT, err, ok, pubchem_smiles_url
from gibh_agent.skills.base_skill import BaseSkill
from gibh_agent.skills.launch_skill_demos import apply_launch_demo_defaults


class PubchemSmilesToCidSkill(BaseSkill):
    __abstractskill__ = False

    """
    通过 SMILES 获取 PubChem CID。

    调用 PubChem PUG REST API，将 SMILES 解析为化合物 CID 列表。
    适用于化合物标识符桥接、后续 PubChem/ChEMBL 联合检索。

    参数:
        smiles: 输入 SMILES 字符串（必填）。
        max_cids: 最多返回 CID 数量，默认 20。
    """

    skill_id = "pubchem_smiles_to_cid"
    display_name = "通过SMILES获取CID"
    description = "通过 SMILES 字符串查询 PubChem 化合物 CID 列表（PUG REST API）。"
    category = "化学与分子信息学"
    sub_category = "数据分析"
    aliases = ["PubChem CID", "SMILES转CID", "smilesCid", "PubChem"]
    required_parameters = ["smiles"]
    tool_chain_key = "smilesCid"
    __dependencies__ = ["pip:httpx"]

    def execute(self, smiles: str = "", max_cids: int = 20, **kwargs: Any) -> Dict[str, Any]:
        filled = apply_launch_demo_defaults(
            self.skill_id,
            {"smiles": (smiles or kwargs.get("smiles_text") or "").strip(), "max_cids": max_cids},
        )
        smi = str(filled.get("smiles") or "").strip()
        if not smi:
            return err("请提供 smiles")

        url = pubchem_smiles_url(smi, "cids/JSON")
        lim = max(1, min(int(filled.get("max_cids") or 20), 100))

        try:
            with httpx.Client(timeout=DEFAULT_HTTP_TIMEOUT, follow_redirects=True) as client:
                resp = client.get(url)
                if resp.status_code == 404:
                    return ok("PubChem 未找到匹配化合物", smiles=smi, cids=[], count=0)
                resp.raise_for_status()
                payload = resp.json()
        except httpx.HTTPError as exc:
            return err(f"PubChem API 请求失败: {exc}")

        cids: List[int] = []
        block = payload.get("IdentifierList") or {}
        for cid in block.get("CID") or []:
            try:
                cids.append(int(cid))
            except (TypeError, ValueError):
                continue
        cids = cids[:lim]

        return ok(
            f"已解析 {len(cids)} 个 PubChem CID",
            smiles=smi,
            cids=cids,
            count=len(cids),
        )
