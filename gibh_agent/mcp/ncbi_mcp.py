# -*- coding: utf-8 -*-
"""
NCBI 文献与基因检索 MCP（PubMed / Gene，E-utilities + 可选 Datasets v2）。

- E-utilities：在 URL 查询参数中附加 `api_key`（NCBI 官方方式，提升速率限额）。
- NCBI Datasets API v2：在请求头中附加 `api-key`（基因符号精准解析）。
"""
from __future__ import annotations

import logging
import os
from typing import Any, Dict, List, Optional

import requests

from gibh_agent.core.tool_registry import registry

logger = logging.getLogger(__name__)

EUTIL_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
DATASETS_GENE_SYMBOL = "https://api.ncbi.nlm.nih.gov/datasets/v2/gene/symbol"


def _eutils_params(extra: Dict[str, Any]) -> Dict[str, Any]:
    p: Dict[str, Any] = dict(extra)
    p.setdefault("retmode", "json")
    p["tool"] = "gibh_agent"
    p["email"] = os.getenv("NCBI_CONTACT_EMAIL", "gibh-agent@localhost")
    key = os.getenv("NCBI_API_KEY", "").strip()
    if key:
        p["api_key"] = key
    return p


def _pubmed_search(query: str, max_results: int) -> List[Dict[str, Any]]:
    q = (query or "").strip()
    if not q:
        return []
    max_results = max(1, min(int(max_results or 10), 50))
    r = requests.get(
        f"{EUTIL_BASE}/esearch.fcgi",
        params=_eutils_params({"db": "pubmed", "term": q, "retmax": max_results}),
        timeout=45,
    )
    r.raise_for_status()
    data = r.json()
    idlist = (data.get("esearchresult") or {}).get("idlist") or []
    if not idlist:
        return []
    r2 = requests.get(
        f"{EUTIL_BASE}/esummary.fcgi",
        params=_eutils_params({"db": "pubmed", "id": ",".join(idlist)}),
        timeout=45,
    )
    r2.raise_for_status()
    summ = r2.json()
    result = (summ.get("result") or {})
    uids = result.get("uids") or []
    out: List[Dict[str, Any]] = []
    for uid in uids:
        rec = result.get(uid) or {}
        title = rec.get("title") or rec.get("sorttitle") or ""
        authors = rec.get("authors") or []
        auth_str = ""
        if isinstance(authors, list) and authors:
            first = authors[0]
            if isinstance(first, dict):
                auth_str = first.get("name", "")
        src = rec.get("source") or rec.get("fulljournalname") or ""
        pubdate = rec.get("pubdate") or rec.get("epubdate") or ""
        body = " ".join(x for x in (src, pubdate, auth_str) if x).strip()
        out.append(
            {
                "pmid": str(uid),
                "title": (title or "")[:800],
                "href": f"https://pubmed.ncbi.nlm.nih.gov/{uid}/",
                "body": (body or "")[:1500],
            }
        )
    return out


def _gene_search_eutils(query: str, max_results: int) -> List[Dict[str, Any]]:
    q = (query or "").strip()
    if not q:
        return []
    max_results = max(1, min(int(max_results or 10), 50))
    r = requests.get(
        f"{EUTIL_BASE}/esearch.fcgi",
        params=_eutils_params({"db": "gene", "term": q, "retmax": max_results}),
        timeout=45,
    )
    r.raise_for_status()
    data = r.json()
    idlist = (data.get("esearchresult") or {}).get("idlist") or []
    if not idlist:
        return []
    r2 = requests.get(
        f"{EUTIL_BASE}/esummary.fcgi",
        params=_eutils_params({"db": "gene", "id": ",".join(idlist)}),
        timeout=45,
    )
    r2.raise_for_status()
    summ = r2.json()
    result = (summ.get("result") or {})
    uids = result.get("uids") or []
    out: List[Dict[str, Any]] = []
    for uid in uids:
        rec = result.get(uid) or {}
        name = rec.get("name") or rec.get("description") or ""
        desc = rec.get("description") or rec.get("summary") or ""
        org = rec.get("organism")
        organism = org.get("scientificname") if isinstance(org, dict) else (org if isinstance(org, str) else "")
        body = " ".join(x for x in (desc, organism) if x)[:1500]
        out.append(
            {
                "gene_id": str(uid),
                "title": (name or f"Gene {uid}")[:500],
                "href": f"https://www.ncbi.nlm.nih.gov/gene/{uid}",
                "body": body,
            }
        )
    return out


def _gene_datasets_by_symbol(symbol: str) -> Optional[Dict[str, Any]]:
    """Datasets v2：请求头 `api-key`（需 NCBI_API_KEY）。"""
    key = os.getenv("NCBI_API_KEY", "").strip()
    if not key:
        return None
    sym = symbol.strip()
    if not sym or len(sym) > 32 or any(c in sym for c in " \t\n\r"):
        return None
    url = f"{DATASETS_GENE_SYMBOL}/{sym}"
    resp = requests.get(url, headers={"api-key": key}, timeout=45)
    if resp.status_code != 200:
        logger.debug("[mcp_ncbi_search] Datasets symbol %s HTTP %s", sym, resp.status_code)
        return None
    return resp.json()


@registry.register(
    name="mcp_ncbi_search",
    description=(
        "在 NCBI 检索 PubMed 文献摘要或 Gene 基因条目。"
        "适合精准文献、PMID/关键词、基因符号与基因 ID 查询。"
        "参数 database 取 pubmed、gene 或 auto（短词偏基因、否则偏文献）。"
    ),
    category="MCP",
    output_type="json",
)
def mcp_ncbi_search(
    query: str,
    database: str = "pubmed",
    max_results: int = 10,
) -> Dict[str, Any]:
    """
    NCBI 检索（E-utilities；配置 NCBI_API_KEY 时自动附加 api_key / api-key）。

    Args:
        query: 检索式（PubMed 支持 MeSH/字段；Gene 支持符号或描述）
        database: pubmed | gene | auto
        max_results: 返回条数上限（1–50）
    """
    q = (query or "").strip()
    if not q:
        return {"status": "error", "error": "query 不能为空", "results": [], "database": database}

    db = (database or "pubmed").strip().lower()
    if db not in ("pubmed", "gene", "auto"):
        db = "pubmed"

    max_results = max(1, min(int(max_results or 10), 50))
    has_key = bool(os.getenv("NCBI_API_KEY", "").strip())

    resolved = db
    if db == "auto":
        resolved = "gene" if len(q) <= 12 and q.isalnum() else "pubmed"

    results: List[Dict[str, Any]] = []
    datasets_note: Optional[str] = None

    try:
        if resolved == "gene":
            if has_key:
                ds = _gene_datasets_by_symbol(q)
                if ds:
                    datasets_note = "datasets_v2_symbol"
                    genes = ds.get("genes") or ds.get("gene_reports") or ds.get("reports") or []
                    if isinstance(genes, list):
                        for g in genes[:max_results]:
                            if not isinstance(g, dict):
                                continue
                            node = g.get("gene") if isinstance(g.get("gene"), dict) else g
                            if not isinstance(node, dict):
                                continue
                            gid = node.get("gene_id") or node.get("geneId") or ""
                            sym = node.get("symbol") or node.get("gene_symbol") or q
                            desc = node.get("description") or node.get("name") or ""
                            tax = node.get("taxname") or node.get("tax_name") or ""
                            results.append(
                                {
                                    "gene_id": str(gid),
                                    "title": str(sym)[:200],
                                    "href": f"https://www.ncbi.nlm.nih.gov/gene/{gid}" if gid else "",
                                    "body": " ".join(x for x in (desc, tax) if x)[:1500],
                                }
                            )
            if not results:
                results = _gene_search_eutils(q, max_results)
        else:
            results = _pubmed_search(q, max_results)

        return {
            "status": "success",
            "database": resolved,
            "query": q,
            "api_key_configured": has_key,
            "datasets_enrichment": datasets_note,
            "results": results,
        }
    except Exception as e:
        logger.exception("❌ [mcp_ncbi_search] 失败: %s", e)
        return {
            "status": "error",
            "error": str(e),
            "database": resolved,
            "query": q,
            "api_key_configured": has_key,
            "results": [],
        }
