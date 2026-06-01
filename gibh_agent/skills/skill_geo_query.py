# -*- coding: utf-8 -*-
"""GEO 数据库查询 — NCBI GEO DataSets 元数据（E-utilities JSON，不下载表达矩阵）。"""
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


def _parse_geo_row(rec: Dict[str, Any]) -> Dict[str, str]:
    acc = str(rec.get("accession") or f"GSE{rec.get('gse')}" or rec.get("uid") or "—")
    gse_num = rec.get("gse")
    if acc == "—" and gse_num:
        acc = f"GSE{gse_num}"
    return {
        "Accession": acc,
        "Title": str(rec.get("title") or "—")[:120],
        "Type": str(rec.get("gdstype") or rec.get("entrytype") or "—")[:80],
        "Organism": str(rec.get("taxon") or rec.get("platformtaxa") or "—"),
        "Samples": str(rec.get("n_samples") or "—"),
        "Platform": str(rec.get("gpl") or "—"),
        "PubDate": str(rec.get("pdat") or "—"),
        "GEO_URL": f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={acc}" if acc != "—" else "—",
    }


class GeoQuerySkill(BaseSkill):
    __abstractskill__ = False

    skill_id = "geo_query"
    display_name = "GEO 数据库查询"
    description = "检索 NCBI GEO 表达数据集 accession 与实验设计摘要（仅元数据，不下载矩阵）。"
    category = "生物医药"
    sub_category = "信息检索"
    aliases = ["GEO", "基因表达", "geo", "GSE"]
    required_parameters = ["query"]
    tool_chain_key = "geo"
    output_type = "mixed"
    __dependencies__ = ["pip:httpx"]

    def execute(self, query: str = "", limit: int = 10, **kwargs: Any) -> Dict[str, Any]:
        filled = apply_launch_demo_defaults(
            self.skill_id,
            {"query": (query or kwargs.get("accession") or "").strip(), "limit": limit},
        )
        q = str(filled.get("query") or "").strip()
        if not q:
            return error_result("请提供 query（关键词、GSE accession 或 MeSH 主题）")

        lim = max(1, min(int(filled.get("limit") or 10), 20))
        term = q
        if q.upper().startswith("GSE") and q[3:].isdigit():
            term = f"{q}[Accession]"
        try:
            ids = eutils_esearch("gds", term, retmax=lim, timeout=DEFAULT_HTTP_TIMEOUT)
            if not ids:
                return empty_result(query=q)
            summaries = eutils_esummary("gds", ids[:lim], timeout=DEFAULT_HTTP_TIMEOUT)
            hits = [_parse_geo_row(summaries[uid]) for uid in ids[:lim] if uid in summaries]
        except httpx.TimeoutException:
            return timeout_result(query=q)
        except httpx.HTTPError as exc:
            return error_result(f"GEO API 请求失败：{exc}")

        if not hits:
            return empty_result(query=q)

        cols = ["Accession", "Title", "Type", "Organism", "Samples", "Platform", "PubDate", "GEO_URL"]
        md_lines = [
            f"### GEO 数据集检索（{len(hits)} 条 · 仅元数据）\n",
            f"**检索式**：`{q}`\n",
        ]
        for i, row in enumerate(hits[:6], 1):
            md_lines.append(
                f"{i}. **{row.get('Accession')}** — {row.get('Title')} "
                f"（{row.get('Samples')} samples）[GEO]({row.get('GEO_URL')})\n"
            )

        return success_payload(
            f"GEO 命中 {len(hits)} 个数据集",
            markdown="\n".join(md_lines),
            hits=hits,
            table_data=table_data_from_rows(cols, hits),
            query=q,
            result_count=len(hits),
        )
