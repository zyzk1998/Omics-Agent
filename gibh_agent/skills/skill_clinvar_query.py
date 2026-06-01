# -*- coding: utf-8 -*-
"""ClinVar 数据库查询 — NCBI E-utilities（JSON esearch + esummary）。"""
from __future__ import annotations

from typing import Any, Dict, List

import httpx

from gibh_agent.skills._eutils_json import eutils_esearch, eutils_esummary
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


def _parse_clinvar_row(rec: Dict[str, Any]) -> Dict[str, str]:
    germ = rec.get("germline_classification") or {}
    if not isinstance(germ, dict):
        germ = {}
    clin_sig = str(germ.get("description") or "").strip()
    review = str(germ.get("review_status") or "").strip()
    traits: List[str] = []
    for ts in germ.get("trait_set") or []:
        if isinstance(ts, dict) and ts.get("trait_name"):
            traits.append(str(ts["trait_name"]))
    loc = ""
    vs = rec.get("variation_set") or []
    if isinstance(vs, list) and vs and isinstance(vs[0], dict):
        vl = vs[0].get("variation_loc") or []
        if isinstance(vl, list) and vl and isinstance(vl[0], dict):
            loc = f"chr{vl[0].get('chr', '')}:{vl[0].get('start', '')}"
    return {
        "Accession": str(rec.get("accession") or rec.get("uid") or "—"),
        "Title": str(rec.get("title") or "—")[:120],
        "Gene": str(rec.get("gene_sort") or "—"),
        "ClinicalSignificance": clin_sig or "—",
        "ReviewStatus": (review[:60] + "…") if len(review) > 60 else (review or "—"),
        "Location": loc or "—",
        "Traits": "; ".join(traits[:2])[:100] if traits else "—",
    }


class ClinvarQuerySkill(BaseSkill):
    __abstractskill__ = False

    skill_id = "clinvar_query"
    display_name = "ClinVar 数据库查询"
    description = "检索 ClinVar 变异临床意义、基因与提交记录（NCBI E-utilities）。"
    category = "生物医药"
    sub_category = "信息检索"
    aliases = ["ClinVar", "临床变异", "clinvar", "变异临床意义"]
    required_parameters = ["query"]
    tool_chain_key = "clinvar"
    output_type = "mixed"
    __dependencies__ = ["pip:httpx"]

    def execute(self, query: str = "", limit: int = 10, **kwargs: Any) -> Dict[str, Any]:
        filled = apply_launch_demo_defaults(
            self.skill_id,
            {"query": (query or kwargs.get("search_query") or "").strip(), "limit": limit},
        )
        q = str(filled.get("query") or "").strip()
        if not q:
            return error_result("请提供 query（基因名、rsID 或 HGVS 描述）")

        lim = max(1, min(int(filled.get("limit") or 10), 25))
        try:
            ids = eutils_esearch("clinvar", q, retmax=lim, timeout=DEFAULT_HTTP_TIMEOUT)
            if not ids:
                return empty_result(query=q)
            summaries = eutils_esummary("clinvar", ids[:lim], timeout=DEFAULT_HTTP_TIMEOUT)
            hits = [_parse_clinvar_row(summaries[uid]) for uid in ids[:lim] if uid in summaries]
        except httpx.TimeoutException:
            return timeout_result(query=q)
        except httpx.HTTPError as exc:
            return error_result(f"ClinVar API 请求失败：{exc}")

        if not hits:
            return empty_result(query=q)

        cols = ["Accession", "Title", "Gene", "ClinicalSignificance", "ReviewStatus", "Location", "Traits"]
        md_lines = [
            f"### ClinVar 检索结果（{len(hits)} 条）\n",
            f"**检索式**：`{q}`\n",
        ]
        for i, row in enumerate(hits[:6], 1):
            sig = row.get("ClinicalSignificance", "—")
            md_lines.append(
                f"{i}. **{row.get('Gene')}** — {row.get('Title')} "
                f"（`{sig}`）[Acc:{row.get('Accession')}]\n"
            )

        return success_payload(
            f"ClinVar 命中 {len(hits)} 条变异",
            markdown="\n".join(md_lines),
            hits=hits,
            table_data=table_data_from_rows(cols, hits),
            query=q,
            result_count=len(hits),
        )
