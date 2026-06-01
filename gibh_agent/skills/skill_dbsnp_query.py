# -*- coding: utf-8 -*-
"""dbSNP 数据库查询 — NCBI E-utilities（JSON esearch + esummary）。"""
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


def _format_genes(genes: Any) -> str:
    if not isinstance(genes, list):
        return "—"
    names: List[str] = []
    for g in genes:
        if isinstance(g, dict) and g.get("name"):
            names.append(str(g["name"]))
    return ", ".join(names[:3]) if names else "—"


def _parse_dbsnp_row(rec: Dict[str, Any]) -> Dict[str, str]:
    snp_num = rec.get("snp_id")
    rs = f"rs{snp_num}" if snp_num else str(rec.get("uid") or "—")
    return {
        "RsID": rs,
        "Chr": str(rec.get("chr") or "—"),
        "ClinicalSignificance": str(rec.get("clinical_significance") or "—")[:80],
        "Genes": _format_genes(rec.get("genes")),
        "FunctionClass": str(rec.get("fxn_class") or "—")[:80],
        "GlobalMAF": str(rec.get("global_maf") or "—"),
    }


class DbsnpQuerySkill(BaseSkill):
    __abstractskill__ = False

    skill_id = "dbsnp_query"
    display_name = "dbSNP 数据库查询"
    description = "检索 dbSNP 位点、基因关联与临床意义概览（NCBI E-utilities）。"
    category = "生物医药"
    sub_category = "信息检索"
    aliases = ["dbSNP", "SNP", "dbsnp", "rsID"]
    required_parameters = ["query"]
    tool_chain_key = "dbsnp"
    output_type = "mixed"
    __dependencies__ = ["pip:httpx"]

    def execute(self, query: str = "", limit: int = 10, **kwargs: Any) -> Dict[str, Any]:
        filled = apply_launch_demo_defaults(
            self.skill_id,
            {"query": (query or kwargs.get("rsid") or "").strip(), "limit": limit},
        )
        q = str(filled.get("query") or "").strip()
        if not q:
            return error_result("请提供 query（rsID、基因名或染色体位置）")

        lim = max(1, min(int(filled.get("limit") or 10), 25))
        term = q if q.lower().startswith("rs") else q
        try:
            ids = eutils_esearch("snp", term, retmax=lim, timeout=DEFAULT_HTTP_TIMEOUT)
            if not ids:
                return empty_result(query=q)
            summaries = eutils_esummary("snp", ids[:lim], timeout=DEFAULT_HTTP_TIMEOUT)
            hits = [_parse_dbsnp_row(summaries[uid]) for uid in ids[:lim] if uid in summaries]
        except httpx.TimeoutException:
            return timeout_result(query=q)
        except httpx.HTTPError as exc:
            return error_result(f"dbSNP API 请求失败：{exc}")

        if not hits:
            return empty_result(query=q)

        cols = ["RsID", "Chr", "ClinicalSignificance", "Genes", "FunctionClass", "GlobalMAF"]
        md = (
            f"### dbSNP 检索结果（{len(hits)} 条）\n\n"
            f"**查询**：`{q}`\n\n"
            f"**首条**：{hits[0].get('RsID')} — {hits[0].get('Genes')} "
            f"（chr{hits[0].get('Chr')}）\n"
        )

        return success_payload(
            f"dbSNP 命中 {len(hits)} 条",
            markdown=md,
            hits=hits,
            table_data=table_data_from_rows(cols, hits),
            query=q,
            result_count=len(hits),
        )
