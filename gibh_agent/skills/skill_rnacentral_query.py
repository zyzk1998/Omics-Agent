# -*- coding: utf-8
"""RNAcentral 按编号查询 — REST + EBI Search。"""
from __future__ import annotations

import re
from typing import Any, Dict, List

import httpx

from gibh_agent.skills._skill_common import DEFAULT_HTTP_TIMEOUT
from gibh_agent.skills._skill_payload import (
    empty_result,
    error_result,
    success_payload,
    table_data_from_rows,
    timeout_result,
)
from gibh_agent.skills.base_skill import BaseSkill
from gibh_agent.skills.launch_skill_demos import apply_launch_demo_defaults

_RNACENTRAL_RNA = "https://rnacentral.org/api/v1/rna/{rna_id}"
_EBI_RNACENTRAL = "https://www.ebi.ac.uk/ebisearch/ws/rest/rnacentral"


def _parse_rnacentral_entry(data: Dict[str, Any]) -> Dict[str, str]:
    xrefs = data.get("xrefs") or []
    xref_str = "—"
    if isinstance(xrefs, list) and xrefs:
        parts: List[str] = []
        for x in xrefs[:3]:
            if isinstance(x, dict):
                parts.append(f"{x.get('database', '')}:{x.get('accession', '')}")
        xref_str = "; ".join(p for p in parts if p.strip(":")) or "—"
    return {
        "RNAcentralID": str(data.get("rnacentral_id") or data.get("id") or "—"),
        "Length": str(data.get("length") or "—"),
        "MD5": str(data.get("md5") or "—")[:12],
        "CrossRefs": xref_str[:120],
        "URL": str(data.get("url") or "—"),
    }


class RnacentralQuerySkill(BaseSkill):
    __abstractskill__ = False

    skill_id = "rnacentral_query"
    display_name = "RNAcentral 按编号查询"
    description = "按 RNAcentral URS 编号或关键词检索 RNA 家族条目与交叉 ID。"
    category = "生物医药"
    sub_category = "信息检索"
    aliases = ["RNAcentral", "RNA", "rnacentral", "非编码RNA"]
    required_parameters = ["query"]
    tool_chain_key = "rnacentral"
    output_type = "mixed"
    __dependencies__ = ["pip:httpx"]

    def execute(self, query: str = "", limit: int = 10, **kwargs: Any) -> Dict[str, Any]:
        filled = apply_launch_demo_defaults(
            self.skill_id,
            {"query": (query or kwargs.get("rna_id") or "").strip(), "limit": limit},
        )
        q = str(filled.get("query") or "").strip()
        if not q:
            return error_result("请提供 query（URS 编号或 RNA 关键词）")

        lim = max(1, min(int(filled.get("limit") or 10), 15))
        is_urs = bool(re.match(r"^URS[0-9A-F]{10}$", q, re.I))

        try:
            with httpx.Client(timeout=DEFAULT_HTTP_TIMEOUT, follow_redirects=True) as client:
                if is_urs:
                    resp = client.get(_RNACENTRAL_RNA.format(rna_id=q.upper()))
                    if resp.status_code == 404:
                        return empty_result(query=q)
                    resp.raise_for_status()
                    hits = [_parse_rnacentral_entry(resp.json())]
                else:
                    resp = client.get(
                        _EBI_RNACENTRAL,
                        params={"query": q, "size": lim, "format": "json"},
                    )
                    resp.raise_for_status()
                    entries = (resp.json().get("entries") or [])
                    hits = []
                    for ent in entries[:lim]:
                        if not isinstance(ent, dict):
                            continue
                        rid = str(ent.get("id") or "").split("_")[0]
                        if not rid.startswith("URS"):
                            continue
                        r2 = client.get(_RNACENTRAL_RNA.format(rna_id=rid))
                        if r2.status_code == 200:
                            hits.append(_parse_rnacentral_entry(r2.json()))
                        if len(hits) >= lim:
                            break
        except httpx.TimeoutException:
            return timeout_result(query=q)
        except httpx.HTTPError as exc:
            return error_result(f"RNAcentral API 请求失败：{exc}")

        if not hits:
            return empty_result(query=q)

        cols = ["RNAcentralID", "Length", "MD5", "CrossRefs", "URL"]
        md = (
            f"### RNAcentral 检索（{len(hits)} 条）\n\n"
            f"**查询**：`{q}`\n\n"
            f"**首条**：{hits[0]['RNAcentralID']}（{hits[0]['Length']} nt）\n"
        )

        return success_payload(
            f"RNAcentral 命中 {len(hits)} 条",
            markdown=md,
            hits=hits,
            table_data=table_data_from_rows(cols, hits),
            query=q,
            result_count=len(hits),
        )
