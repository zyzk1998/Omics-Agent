# -*- coding: utf-8 -*-
"""GWAS Catalog 数据库查询 — EBI REST API。"""
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

_GWAS_BASE = "https://www.ebi.ac.uk/gwas/rest/api"


def _normalize_rsid(q: str) -> str:
    s = q.strip()
    if s.lower().startswith("rs"):
        return s
    if re.fullmatch(r"\d+", s):
        return f"rs{s}"
    return s


def _fetch_by_rsid(rsid: str, limit: int, timeout: float) -> List[Dict[str, Any]]:
    with httpx.Client(timeout=timeout, follow_redirects=True) as client:
        resp = client.get(
            f"{_GWAS_BASE}/associations/search/findByRsId",
            params={"rsId": rsid, "page": 0, "size": limit},
        )
        resp.raise_for_status()
        data = resp.json()
    assocs = (data.get("_embedded") or {}).get("associations") or []
    hits: List[Dict[str, Any]] = []
    for a in assocs:
        if not isinstance(a, dict):
            continue
        locus = a.get("loci") or []
        rs = rsid
        trait = ""
        if isinstance(locus, list) and locus and isinstance(locus[0], dict):
            rs = str(locus[0].get("strongestRiskAlleles") or rsid)
            trait = str(locus[0].get("authorReportedGene") or "")
        pval = a.get("pvalue")
        pval_s = f"{pval:.2e}" if isinstance(pval, (int, float)) else str(pval or "—")
        hits.append(
            {
                "RsID": rs,
                "PValue": pval_s,
                "RiskFrequency": str(a.get("riskFrequency") or "—"),
                "OR": str(a.get("orPerCopyNum") or "—"),
                "ReportedGene": trait or "—",
            }
        )
    return hits


def _fetch_by_trait(trait: str, limit: int, timeout: float) -> List[Dict[str, Any]]:
    with httpx.Client(timeout=timeout, follow_redirects=True) as client:
        resp = client.get(
            f"{_GWAS_BASE}/studies/search/findByDiseaseTrait",
            params={"diseaseTrait": trait, "page": 0, "size": limit},
        )
        resp.raise_for_status()
        data = resp.json()
    studies = (data.get("_embedded") or {}).get("studies") or []
    hits: List[Dict[str, Any]] = []
    for st in studies:
        if not isinstance(st, dict):
            continue
        pubmed = ""
        pubs = st.get("publicationInfo") or {}
        if isinstance(pubs, dict):
            pubmed = str(pubs.get("pubmedId") or "")
        hits.append(
            {
                "StudyID": str(st.get("accessionId") or st.get("gwasId") or "—"),
                "Trait": str(st.get("diseaseTrait") or trait)[:80],
                "InitialSampleSize": str(st.get("initialSampleSize") or "—")[:100],
                "PubmedID": pubmed or "—",
                "SnpCount": str(st.get("snpCount") or "—"),
            }
        )
    return hits


class GwasCatalogQuerySkill(BaseSkill):
    __abstractskill__ = False

    skill_id = "gwas_catalog_query"
    display_name = "GWAS catalog 数据库查询"
    description = "检索 GWAS Catalog 性状—SNP 关联或研究元数据（EBI REST）。"
    category = "生物医药"
    sub_category = "信息检索"
    aliases = ["GWAS", "GWAS Catalog", "gwas", "全基因组关联"]
    required_parameters = ["query"]
    tool_chain_key = "gwasCatalog"
    output_type = "mixed"
    __dependencies__ = ["pip:httpx"]

    def execute(self, query: str = "", limit: int = 10, **kwargs: Any) -> Dict[str, Any]:
        filled = apply_launch_demo_defaults(
            self.skill_id,
            {"query": (query or kwargs.get("trait") or "").strip(), "limit": limit},
        )
        q = str(filled.get("query") or "").strip()
        if not q:
            return error_result("请提供 query（rsID 或疾病/性状名称）")

        lim = max(1, min(int(filled.get("limit") or 10), 20))
        is_rsid = bool(re.match(r"^rs?\d+$", q, re.I))
        try:
            if is_rsid:
                hits = _fetch_by_rsid(_normalize_rsid(q), lim, DEFAULT_HTTP_TIMEOUT)
                cols = ["RsID", "PValue", "RiskFrequency", "OR", "ReportedGene"]
                mode = "association"
            else:
                hits = _fetch_by_trait(q, lim, DEFAULT_HTTP_TIMEOUT)
                cols = ["StudyID", "Trait", "InitialSampleSize", "PubmedID", "SnpCount"]
                mode = "study"
        except httpx.TimeoutException:
            return timeout_result(query=q)
        except httpx.HTTPError as exc:
            return error_result(f"GWAS Catalog API 请求失败：{exc}")

        if not hits:
            return empty_result(query=q)

        md = (
            f"### GWAS Catalog 检索结果（{len(hits)} 条 · {mode}）\n\n"
            f"**查询**：`{q}`\n"
        )

        return success_payload(
            f"GWAS Catalog 命中 {len(hits)} 条",
            markdown=md,
            hits=hits,
            table_data=table_data_from_rows(cols, hits),
            query=q,
            result_count=len(hits),
            search_mode=mode,
        )
