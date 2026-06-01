# -*- coding: utf-8 -*-
"""PubMed 数据库查询 — NCBI E-utilities（Biopython Entrez）。"""
from __future__ import annotations

import os
from typing import Any, Dict, List
from urllib.error import URLError

from gibh_agent.skills._entrez_util import entrez_as_list, entrez_get, iter_pubmed_esummary_docs, parse_pubmed_esummary_doc
from gibh_agent.skills._skill_payload import (
    empty_result,
    error_result,
    success_payload,
    table_data_from_rows,
    timeout_result,
)
from gibh_agent.skills.base_skill import BaseSkill
from gibh_agent.skills.launch_skill_demos import apply_launch_demo_defaults

_ENTREZ_EMAIL = (os.getenv("ENTREZ_EMAIL") or "omics-agent@example.com").strip()
_ENTREZ_TOOL = (os.getenv("ENTREZ_TOOL") or "GIBH_OmicsAgent").strip()
_ENTREZ_TIMEOUT = float(os.getenv("GIBH_ENTREZ_TIMEOUT", "45"))


def _fetch_pubmed_hits(query: str, limit: int) -> List[Dict[str, Any]]:
    from Bio import Entrez

    Entrez.email = _ENTREZ_EMAIL
    Entrez.tool = _ENTREZ_TOOL
    Entrez.timeout = _ENTREZ_TIMEOUT

    with Entrez.esearch(db="pubmed", term=query, retmax=limit, sort="relevance") as h_search:
        search_rec = Entrez.read(h_search)
    id_list = entrez_as_list(entrez_get(search_rec, "IdList"))
    if not id_list:
        return []

    id_strs = [str(i) for i in id_list if i is not None]
    with Entrez.esummary(db="pubmed", id=",".join(id_strs)) as h_sum:
        summaries = Entrez.read(h_sum)

    hits: List[Dict[str, Any]] = []
    for doc in iter_pubmed_esummary_docs(summaries):
        if doc is None:
            continue
        row = parse_pubmed_esummary_doc(doc)
        if row.get("PMID"):
            hits.append(row)
    return hits


class PubmedQuerySkill(BaseSkill):
    __abstractskill__ = False

    skill_id = "pubmed_query"
    display_name = "PubMed数据库查询"
    description = "按关键词检索 PubMed 文献题录（E-utilities），返回表格与 Markdown 摘要列表。"
    category = "生物医药"
    sub_category = "信息检索"
    aliases = ["PubMed", "文献检索", "pubmed", "NCBI PubMed"]
    required_parameters = ["query"]
    tool_chain_key = "pubmed"
    output_type = "mixed"
    __dependencies__ = ["pip:biopython", "pip:httpx"]

    def execute(self, query: str = "", limit: int = 10, **kwargs: Any) -> Dict[str, Any]:
        filled = apply_launch_demo_defaults(
            self.skill_id,
            {
                "query": (query or kwargs.get("search_query") or "").strip(),
                "limit": limit,
            },
        )
        q = str(filled.get("query") or "").strip()
        if not q:
            return error_result("请提供 query（检索关键词或 PMID 列表）")

        lim = max(1, min(int(filled.get("limit") or 10), 30))

        try:
            hits = _fetch_pubmed_hits(q, lim)
        except TimeoutError:
            return timeout_result(query=q)
        except URLError as exc:
            if "timed out" in str(exc).lower():
                return timeout_result(query=q)
            return error_result(f"PubMed 网络请求失败：{exc}")
        except Exception as exc:
            err_s = str(exc).lower()
            if "timed out" in err_s or "timeout" in err_s:
                return timeout_result(query=q)
            return error_result(f"PubMed 检索失败：{exc}")

        if not hits:
            return empty_result(query=q)

        cols = ["PMID", "Title", "Journal", "PubDate", "Authors"]
        md_lines = [
            f"### PubMed 检索结果（{len(hits)} 条）\n",
            f"**检索式**：`{q}`\n",
        ]
        for i, row in enumerate(hits[:8], 1):
            md_lines.append(
                f"{i}. **{row.get('Title', '—')}** — *{row.get('Journal', '')}* "
                f"({row.get('PubDate', '')}) [PMID:{row.get('PMID', '')}]\n"
            )
        if len(hits) > 8:
            md_lines.append(f"\n*另有 {len(hits) - 8} 条见下表。*\n")

        return success_payload(
            f"PubMed 命中 {len(hits)} 条文献",
            markdown="\n".join(md_lines),
            hits=hits,
            table_data=table_data_from_rows(cols, hits),
            query=q,
            result_count=len(hits),
        )
