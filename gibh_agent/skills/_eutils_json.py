# -*- coding: utf-8 -*-
"""NCBI E-utilities JSON API 辅助（ClinVar / dbSNP / GEO 等，避免 Entrez XML DTD 校验失败）。"""
from __future__ import annotations

import os
from typing import Any, Dict, List, Optional

import httpx

_EUTIL_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
_EUTIL_EMAIL = (os.getenv("ENTREZ_EMAIL") or "omics-agent@example.com").strip()
_EUTIL_TOOL = (os.getenv("ENTREZ_TOOL") or "GIBH_OmicsAgent").strip()
_EUTIL_TIMEOUT = float(os.getenv("GIBH_ENTREZ_TIMEOUT", "45"))


def _base_params(**extra: Any) -> Dict[str, Any]:
    p: Dict[str, Any] = {"retmode": "json", "email": _EUTIL_EMAIL, "tool": _EUTIL_TOOL}
    p.update(extra)
    key = (os.getenv("NCBI_API_KEY") or "").strip()
    if key:
        p["api_key"] = key
    return p


def eutils_esearch(
    db: str,
    term: str,
    *,
    retmax: int = 10,
    timeout: Optional[float] = None,
) -> List[str]:
    """esearch.fcgi → uid 列表。"""
    with httpx.Client(timeout=timeout or _EUTIL_TIMEOUT, follow_redirects=True) as client:
        resp = client.get(
            f"{_EUTIL_BASE}/esearch.fcgi",
            params=_base_params(db=db, term=term, retmax=max(1, retmax)),
        )
        resp.raise_for_status()
        data = resp.json()
    idlist = (data.get("esearchresult") or {}).get("idlist") or []
    return [str(i) for i in idlist if i]


def eutils_esummary(
    db: str,
    ids: List[str],
    *,
    timeout: Optional[float] = None,
) -> Dict[str, Dict[str, Any]]:
    """esummary.fcgi → {uid: record_dict}。"""
    if not ids:
        return {}
    with httpx.Client(timeout=timeout or _EUTIL_TIMEOUT, follow_redirects=True) as client:
        resp = client.get(
            f"{_EUTIL_BASE}/esummary.fcgi",
            params=_base_params(db=db, id=",".join(ids)),
        )
        resp.raise_for_status()
        data = resp.json()
    result = (data.get("result") or {})
    out: Dict[str, Dict[str, Any]] = {}
    for uid in result.get("uids") or []:
        rec = result.get(uid)
        if isinstance(rec, dict):
            out[str(uid)] = rec
    return out
