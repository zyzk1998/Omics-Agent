# -*- coding: utf-8 -*-
"""获取 mRNA 序列工具 — Ensembl REST（cDNA / 转录本序列，FASTA 文本）。"""
from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, Optional, Tuple

import httpx

from gibh_agent.skills._skill_common import DEFAULT_HTTP_TIMEOUT
from gibh_agent.skills._skill_payload import (
    empty_result,
    error_result,
    metrics_cards_from_pairs,
    success_payload,
    timeout_result,
)
from gibh_agent.skills.base_skill import BaseSkill
from gibh_agent.skills.launch_skill_demos import apply_launch_demo_defaults

_ENSEMBL_LOOKUP = "https://rest.ensembl.org/lookup/symbol/homo_sapiens/{symbol}"
_ENSEMBL_SEQ = "https://rest.ensembl.org/sequence/id/{transcript_id}"
_ENSEMBL_HEADERS = {"Content-Type": "application/json", "Accept": "application/json"}


def _resolve_gene_symbol(
    sequence_or_path: str = "",
    *,
    query: str = "",
    gene_symbol: str = "",
) -> Tuple[Optional[str], Optional[str]]:
    """从基因符号文本或 sequence_or_path（非文件路径时）解析符号。"""
    for raw in (gene_symbol, query, sequence_or_path):
        s = (raw or "").strip()
        if not s:
            continue
        p = Path(s)
        if p.is_file():
            try:
                text = p.read_text(encoding="utf-8", errors="replace").strip()
            except OSError as exc:
                return None, f"无法读取文件: {exc}"
            first = text.splitlines()[0].strip() if text else ""
            if first.startswith(">"):
                hdr = first[1:].split()[0]
                return hdr, None
            return None, "FASTA 文件未含有效基因符号头"
        if len(s) <= 30 and " " not in s and "/" not in s:
            return s.upper(), None
    return None, "请提供 sequence_or_path（基因符号，如 TP53）或 query"


def _pick_transcript_id(gene_data: Dict[str, Any]) -> str:
    canonical = (gene_data.get("canonical_transcript") or "").strip()
    if canonical:
        return canonical.split(".")[0]
    transcripts = gene_data.get("Transcript") or []
    if isinstance(transcripts, list):
        for tx in transcripts:
            if isinstance(tx, dict) and tx.get("biotype") == "protein_coding" and tx.get("id"):
                return str(tx["id"]).split(".")[0]
        if transcripts and isinstance(transcripts[0], dict) and transcripts[0].get("id"):
            return str(transcripts[0]["id"]).split(".")[0]
    return ""


class MrnaSequenceFetchSkill(BaseSkill):
    __abstractskill__ = False

    skill_id = "mrna_sequence_fetch"
    display_name = "获取mRNA序列工具"
    description = "按人类基因符号获取 Ensembl 代表性 mRNA（cDNA）序列，返回 FASTA 文本。"
    category = "生物医药"
    sub_category = "信息检索"
    aliases = ["mRNA", "cDNA", "mrna", "转录本序列"]
    required_parameters = ["sequence_or_path"]
    tool_chain_key = "mrna"
    output_type = "mixed"
    __dependencies__ = ["pip:httpx"]

    def execute(
        self,
        sequence_or_path: str = "",
        query: str = "",
        gene_symbol: str = "",
        **kwargs: Any,
    ) -> Dict[str, Any]:
        filled = apply_launch_demo_defaults(
            self.skill_id,
            {
                "sequence_or_path": (sequence_or_path or "").strip(),
                "query": (query or gene_symbol or "").strip(),
            },
        )
        symbol, err = _resolve_gene_symbol(
            str(filled.get("sequence_or_path") or ""),
            query=str(filled.get("query") or ""),
            gene_symbol=gene_symbol,
        )
        if err:
            return error_result(err)
        if not symbol:
            return error_result("请提供 sequence_or_path（基因符号，如 TP53）")

        try:
            with httpx.Client(timeout=DEFAULT_HTTP_TIMEOUT, follow_redirects=True) as client:
                lookup_resp = client.get(
                    _ENSEMBL_LOOKUP.format(symbol=symbol),
                    params={"expand": 1},
                    headers=_ENSEMBL_HEADERS,
                )
                if lookup_resp.status_code == 404:
                    return empty_result(gene_symbol=symbol)
                lookup_resp.raise_for_status()
                gene_data = lookup_resp.json()

                tx_id = _pick_transcript_id(gene_data if isinstance(gene_data, dict) else {})
                if not tx_id:
                    return empty_result(gene_symbol=symbol)

                seq_resp = client.get(
                    _ENSEMBL_SEQ.format(transcript_id=tx_id),
                    params={"type": "cdna"},
                    headers=_ENSEMBL_HEADERS,
                )
                seq_resp.raise_for_status()
                seq_payload = seq_resp.json()
        except httpx.TimeoutException:
            return timeout_result(gene_symbol=symbol)
        except httpx.HTTPError as exc:
            return error_result(f"Ensembl 序列请求失败：{exc}")

        seq = ""
        if isinstance(seq_payload, dict):
            seq = str(seq_payload.get("seq") or "")
        elif isinstance(seq_payload, str):
            seq = seq_payload
        if not seq:
            return empty_result(gene_symbol=symbol)

        fasta_id = f"{symbol}|{tx_id}|cdna"
        fasta = f">{fasta_id}\n{seq}\n"
        preview = seq[:80] + ("…" if len(seq) > 80 else "")
        cards = metrics_cards_from_pairs(
            [
                ("基因符号", symbol, ""),
                ("转录本 ID", tx_id, ""),
                ("序列长度", len(seq), "bp"),
            ]
        )
        md = (
            f"### mRNA（cDNA）序列 — {symbol}\n\n"
            f"**转录本**：`{tx_id}`\n\n"
            f"**长度**：{len(seq)} bp\n\n"
            f"```fasta\n>{fasta_id}\n{preview}\n```\n\n"
            f"*完整序列见下方 data.fasta 字段。*\n"
        )

        return success_payload(
            f"{symbol} mRNA 序列已获取（{len(seq)} bp）",
            markdown=md,
            metrics_cards=cards,
            gene_symbol=symbol,
            transcript_id=tx_id,
            sequence_length=len(seq),
            fasta=fasta,
            sequence_preview=preview,
        )
