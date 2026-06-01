# -*- coding: utf-8 -*-
"""基因蛋白信息查询器 — Ensembl REST 轻聚合。"""
from __future__ import annotations

from typing import Any, Dict, List

import httpx

from gibh_agent.skills._skill_common import DEFAULT_HTTP_TIMEOUT
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

_ENSEMBL_LOOKUP = "https://rest.ensembl.org/lookup/symbol/homo_sapiens/{symbol}"
_ENSEMBL_HEADERS = {"Content-Type": "application/json", "Accept": "application/json"}


def _parse_gene_info(data: Dict[str, Any], symbol: str) -> Dict[str, Any]:
    desc = data.get("description") or "—"
    if isinstance(desc, str) and len(desc) > 150:
        desc = desc[:150] + "…"
    transcripts = data.get("Transcript") or []
    tx_count = len(transcripts) if isinstance(transcripts, list) else 0
    canonical = str(data.get("canonical_transcript") or "—")
    biotype = str(data.get("biotype") or "—")
    chrom = str(data.get("seq_region_name") or "—")
    start = data.get("start")
    end = data.get("end")
    loc = f"chr{chrom}:{start}-{end}" if chrom != "—" and start and end else "—"
    return {
        "GeneSymbol": symbol,
        "EnsemblID": str(data.get("id") or "—"),
        "Biotype": biotype,
        "Location": loc,
        "Transcripts": tx_count,
        "CanonicalTranscript": canonical,
        "Description": str(desc),
        "EnsemblURL": (
            f"https://www.ensembl.org/Homo_sapiens/Gene/Summary?g={data.get('id')}"
            if data.get("id")
            else "—"
        ),
    }


class GeneProteinInfoQuerySkill(BaseSkill):
    __abstractskill__ = False

    skill_id = "gene_protein_info_query"
    display_name = "基因蛋白信息查询器"
    description = "按人类基因符号查询 Ensembl 基因注释、定位与转录本概览（REST 轻聚合）。"
    category = "生物医药"
    sub_category = "信息检索"
    aliases = ["基因查询", "Ensembl", "gene info", "基因注释"]
    required_parameters = ["query"]
    tool_chain_key = "gene"
    output_type = "mixed"
    __dependencies__ = ["pip:httpx"]

    def execute(self, query: str = "", **kwargs: Any) -> Dict[str, Any]:
        filled = apply_launch_demo_defaults(
            self.skill_id,
            {"query": (query or kwargs.get("gene_symbol") or "").strip()},
        )
        q = str(filled.get("query") or "").strip()
        if not q:
            return error_result("请提供 query（人类基因符号，如 TP53、BRCA1）")

        url = _ENSEMBL_LOOKUP.format(symbol=q.upper())
        try:
            with httpx.Client(timeout=DEFAULT_HTTP_TIMEOUT, follow_redirects=True) as client:
                resp = client.get(url, params={"expand": 1}, headers=_ENSEMBL_HEADERS)
                if resp.status_code == 404:
                    return empty_result(query=q)
                resp.raise_for_status()
                data = resp.json()
        except httpx.TimeoutException:
            return timeout_result(query=q)
        except httpx.HTTPError as exc:
            return error_result(f"Ensembl API 请求失败：{exc}")

        if not isinstance(data, dict) or not data.get("id"):
            return empty_result(query=q)

        info = _parse_gene_info(data, q.upper())
        hits: List[Dict[str, Any]] = [info]
        cols = [
            "GeneSymbol",
            "EnsemblID",
            "Biotype",
            "Location",
            "Transcripts",
            "CanonicalTranscript",
            "Description",
            "EnsemblURL",
        ]
        cards = metrics_cards_from_pairs(
            [
                ("基因符号", info["GeneSymbol"], ""),
                ("Ensembl ID", info["EnsemblID"], ""),
                ("定位", info["Location"], ""),
                ("转录本数", info["Transcripts"], ""),
            ]
        )
        md = (
            f"### 基因蛋白信息（{info['GeneSymbol']}）\n\n"
            f"**Ensembl**：`{info['EnsemblID']}` · {info['Biotype']}\n\n"
            f"**定位**：{info['Location']}\n\n"
            f"**描述**：{info['Description']}\n\n"
            f"[Ensembl 详情]({info['EnsemblURL']})\n"
        )

        return success_payload(
            f"Ensembl 基因 {info['GeneSymbol']} 注释已获取",
            markdown=md,
            hits=hits,
            table_data=table_data_from_rows(cols, hits),
            metrics_cards=cards,
            query=q,
            result_count=1,
        )
