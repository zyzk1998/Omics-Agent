# -*- coding: utf-8 -*-
"""UniProt 数据库查询 — REST API。"""
from __future__ import annotations

from typing import Any, Dict, List

import httpx

from gibh_agent.skills._skill_common import DEFAULT_HTTP_TIMEOUT
from gibh_agent.skills._skill_payload import (
    empty_result,
    error_result,
    success_payload,
    table_data_from_rows,
    timeout_result,
    metrics_cards_from_pairs,
)
from gibh_agent.skills.base_skill import BaseSkill
from gibh_agent.skills.launch_skill_demos import apply_launch_demo_defaults

_UNIPROT_SEARCH = "https://rest.uniprot.org/uniprotkb/search"
_UNIPROT_ENTRY = "https://rest.uniprot.org/uniprotkb/{acc}.json"


def _nested_value(obj: Any) -> str:
    """UniProt JSON 中 value 字段或扁平字符串。"""
    if obj is None:
        return ""
    if isinstance(obj, str):
        return obj.strip()
    if isinstance(obj, dict):
        return str(obj.get("value") or "").strip()
    return str(obj).strip()


def _parse_uniprot_entry(entry: Dict[str, Any]) -> Dict[str, Any]:
    if not isinstance(entry, dict):
        return {
            "Accession": "—",
            "Protein": "—",
            "Organism": "—",
            "Length": "—",
            "Gene": "—",
        }

    acc = (entry.get("primaryAccession") or entry.get("uniProtkbId") or "").strip()
    rec_name = ""
    prot = entry.get("proteinDescription") or {}
    if isinstance(prot, dict):
        rn = prot.get("recommendedName") or {}
        if isinstance(rn, dict):
            rec_name = _nested_value(rn.get("fullName"))
        if not rec_name:
            sub = prot.get("submissionNames") or []
            if isinstance(sub, list) and sub:
                first = sub[0]
                rec_name = _nested_value(first.get("fullName") if isinstance(first, dict) else first)
        if not rec_name:
            rec_name = _nested_value(prot.get("proteinName"))

    org = ""
    organism = entry.get("organism") or {}
    if isinstance(organism, dict):
        sci = organism.get("scientificName")
        if isinstance(sci, str):
            org = sci.strip()
        else:
            org = _nested_value(sci)

    length: Any = ""
    seq = entry.get("sequence") or {}
    if isinstance(seq, dict):
        length = seq.get("length") or ""
    if length == "" and entry.get("sequenceLength") is not None:
        length = entry.get("sequenceLength")

    genes: List[str] = []
    raw_genes = entry.get("genes") or []
    if isinstance(raw_genes, list):
        for g in raw_genes:
            if not isinstance(g, dict):
                continue
            gn = g.get("geneName") or {}
            val = _nested_value(gn)
            if val:
                genes.append(val)

    return {
        "Accession": acc or "—",
        "Protein": rec_name or "—",
        "Organism": org or "—",
        "Length": length if length != "" else "—",
        "Gene": ", ".join(genes[:3]) if genes else "—",
    }


def _looks_like_uniprot_accession(query: str) -> bool:
    """UniProtKB accession 形如 P04637、Q12888、A0A0B4J2F0。"""
    import re

    q = (query or "").strip().upper()
    return bool(re.fullmatch(r"[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}", q))


class UniprotQuerySkill(BaseSkill):
    __abstractskill__ = False

    skill_id = "uniprot_query"
    display_name = "UniProt数据库查询"
    description = "检索 UniProt 蛋白条目：功能名、物种、序列长度与基因名；表格 + 指标卡。"
    category = "生物医药"
    sub_category = "信息检索"
    aliases = ["UniProt", "蛋白数据库", "uniprot", "蛋白质注释"]
    required_parameters = ["query"]
    tool_chain_key = "uniprot"
    output_type = "mixed"
    __dependencies__ = ["pip:httpx"]

    def execute(self, query: str = "", limit: int = 10, **kwargs: Any) -> Dict[str, Any]:
        filled = apply_launch_demo_defaults(
            self.skill_id,
            {
                "query": (query or kwargs.get("accession") or "").strip(),
                "limit": limit,
            },
        )
        q = str(filled.get("query") or "").strip()
        if not q:
            return error_result("请提供 query（蛋白名、基因名或 UniProt accession）")

        lim = max(1, min(int(filled.get("limit") or 10), 25))
        timeout = httpx.Timeout(DEFAULT_HTTP_TIMEOUT)

        try:
            with httpx.Client(timeout=timeout, follow_redirects=True) as client:
                if _looks_like_uniprot_accession(q):
                    resp = client.get(_UNIPROT_ENTRY.format(acc=q.upper()))
                    if resp.status_code == 404:
                        hits: List[Dict[str, Any]] = []
                    else:
                        resp.raise_for_status()
                        hits = [_parse_uniprot_entry(resp.json())]
                else:
                    params = {
                        "query": q,
                        "size": lim,
                        "fields": "accession,protein_name,organism_name,length,gene_names",
                    }
                    resp = client.get(_UNIPROT_SEARCH, params=params)
                    resp.raise_for_status()
                    payload = resp.json()
                    results = payload.get("results") or []
                    hits = [
                        _parse_uniprot_entry(r)
                        for r in results
                        if isinstance(r, dict) and (r.get("primaryAccession") or r.get("uniProtkbId"))
                    ]
        except httpx.TimeoutException:
            return timeout_result(query=q)
        except httpx.HTTPError as exc:
            return error_result(f"UniProt API 请求失败：{exc}")

        if not hits:
            return empty_result(query=q)

        cols = ["Accession", "Protein", "Organism", "Length", "Gene"]
        top = hits[0]
        cards = metrics_cards_from_pairs(
            [
                ("UniProt Accession", top.get("Accession"), ""),
                ("蛋白名称", top.get("Protein"), ""),
                ("序列长度", top.get("Length"), "aa"),
                ("物种", top.get("Organism"), ""),
            ]
        )
        md = (
            f"### UniProt 检索结果（{len(hits)} 条）\n\n"
            f"**查询**：`{q}`\n\n"
            f"**首条**：{top.get('Protein')}（`{top.get('Accession')}`）— "
            f"{top.get('Organism')}，{top.get('Length')} aa\n"
        )

        return success_payload(
            f"UniProt 命中 {len(hits)} 条",
            markdown=md,
            hits=hits,
            table_data=table_data_from_rows(cols, hits),
            metrics_cards=cards,
            query=q,
            result_count=len(hits),
        )
