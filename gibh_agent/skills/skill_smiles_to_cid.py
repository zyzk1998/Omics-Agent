# -*- coding: utf-8 -*-
"""通过 SMILES 获取 PubChem CID — PUG REST（tool_id: smiles_to_cid）。"""
from __future__ import annotations

from typing import Any, Dict, List

import httpx

from gibh_agent.skills._skill_common import DEFAULT_HTTP_TIMEOUT, pubchem_smiles_url
from gibh_agent.skills._skill_payload import (
    empty_result,
    error_result,
    metrics_cards_from_pairs,
    success_payload,
    table_data_from_rows,
    timeout_result,
)
from gibh_agent.skills.base_skill import BaseSkill
from gibh_agent.skills.launch_skill_demos import apply_launch_demo_defaults


class SmilesToCidSkill(BaseSkill):
    __abstractskill__ = False

    skill_id = "smiles_to_cid"
    display_name = "通过SMILES获取CID"
    description = "SMILES → PubChem CID 列表（PUG REST）；指标卡 + 表格。"
    category = "化学与分子信息学"
    sub_category = "数据分析"
    aliases = ["SMILES转CID", "PubChem CID", "smilesCid", "smiles_to_cid"]
    required_parameters = ["smiles"]
    tool_chain_key = "smilesCid"
    output_type = "mixed"
    __dependencies__ = ["pip:httpx"]

    def execute(self, smiles: str = "", max_cids: int = 20, **kwargs: Any) -> Dict[str, Any]:
        filled = apply_launch_demo_defaults(
            self.skill_id,
            {"smiles": (smiles or kwargs.get("smiles_text") or "").strip(), "max_cids": max_cids},
        )
        smi = str(filled.get("smiles") or "").strip()
        if not smi:
            return error_result("请提供 smiles")

        lim = max(1, min(int(filled.get("max_cids") or 20), 100))
        url = pubchem_smiles_url(smi, "cids/JSON")
        timeout = httpx.Timeout(DEFAULT_HTTP_TIMEOUT)

        try:
            with httpx.Client(timeout=timeout, follow_redirects=True) as client:
                resp = client.get(url)
                if resp.status_code == 404:
                    return empty_result(smiles=smi)
                resp.raise_for_status()
                payload = resp.json()
        except httpx.TimeoutException:
            return timeout_result(smiles=smi)
        except httpx.HTTPError as exc:
            return error_result(f"PubChem API 请求失败：{exc}")

        cids: List[int] = []
        block = payload.get("IdentifierList") or {}
        for cid in block.get("CID") or []:
            try:
                cids.append(int(cid))
            except (TypeError, ValueError):
                continue
        cids = cids[:lim]

        if not cids:
            return empty_result(smiles=smi)

        hits = [{"CID": cid, "SMILES": smi} for cid in cids]
        cards = metrics_cards_from_pairs(
            [
                ("SMILES", smi[:80] + ("…" if len(smi) > 80 else ""), ""),
                ("命中 CID 数", len(cids), "个"),
                ("首选 CID", cids[0], ""),
            ]
        )
        md = (
            f"### PubChem CID 解析\n\n"
            f"**SMILES**：`{smi}`\n\n"
            f"**CID 列表**：{', '.join(str(c) for c in cids[:15])}"
            + (f" … 共 {len(cids)} 个" if len(cids) > 15 else "")
            + "\n"
        )

        return success_payload(
            f"已解析 {len(cids)} 个 PubChem CID",
            markdown=md,
            hits=hits,
            table_data=table_data_from_rows(["CID", "SMILES"], hits),
            metrics_cards=cards,
            smiles=smi,
            cids=cids,
            count=len(cids),
        )
