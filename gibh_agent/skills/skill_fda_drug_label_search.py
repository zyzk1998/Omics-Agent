# -*- coding: utf-8
"""FDA 药品标签字段检索 — openFDA drug/label API。"""
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
)
from gibh_agent.skills.base_skill import BaseSkill
from gibh_agent.skills.launch_skill_demos import apply_launch_demo_defaults

_OPENFDA_LABEL = "https://api.fda.gov/drug/label.json"


def _first_text(val: Any, max_len: int = 300) -> str:
    if val is None:
        return "—"
    if isinstance(val, list):
        val = val[0] if val else ""
    s = str(val).strip()
    if len(s) > max_len:
        return s[:max_len] + "…"
    return s or "—"


def _parse_fda_label(rec: Dict[str, Any]) -> Dict[str, str]:
    openfda = rec.get("openfda") or {}
    brand = "—"
    generic = "—"
    if isinstance(openfda, dict):
        brands = openfda.get("brand_name") or []
        generics = openfda.get("generic_name") or []
        if brands:
            brand = str(brands[0])
        if generics:
            generic = str(generics[0])
    return {
        "BrandName": brand,
        "GenericName": generic,
        "Indications": _first_text(rec.get("indications_and_usage")),
        "Warnings": _first_text(rec.get("warnings")),
        "Dosage": _first_text(rec.get("dosage_and_administration")),
        "AdverseReactions": _first_text(rec.get("adverse_reactions")),
        "SetID": str(rec.get("set_id") or "—"),
    }


class FdaDrugLabelSearchSkill(BaseSkill):
    __abstractskill__ = False

    skill_id = "fda_drug_label_search"
    display_name = "FDA药品标签字段检索"
    description = "openFDA 按品牌名/通用名检索药品标签关键字段（适应症、警告、用法等）。"
    category = "化学"
    sub_category = "信息检索"
    aliases = ["FDA", "药品标签", "openFDA", "说明书"]
    required_parameters = ["query"]
    tool_chain_key = "fda"
    output_type = "mixed"
    __dependencies__ = ["pip:httpx"]

    def execute(self, query: str = "", limit: int = 5, **kwargs: Any) -> Dict[str, Any]:
        filled = apply_launch_demo_defaults(
            self.skill_id,
            {
                "query": (query or kwargs.get("drug_name") or "").strip(),
                "limit": limit,
            },
        )
        q = str(filled.get("query") or "").strip()
        if not q:
            return error_result("请提供 query（药品品牌名或通用名）")

        lim = max(1, min(int(filled.get("limit") or 5), 10))
        search = f'openfda.brand_name:"{q}" OR openfda.generic_name:"{q}"'
        try:
            with httpx.Client(timeout=DEFAULT_HTTP_TIMEOUT, follow_redirects=True) as client:
                resp = client.get(
                    _OPENFDA_LABEL,
                    params={"search": search, "limit": lim},
                )
                if resp.status_code == 404:
                    return empty_result(query=q)
                resp.raise_for_status()
                data = resp.json()
        except httpx.TimeoutException:
            return timeout_result(query=q)
        except httpx.HTTPError as exc:
            return error_result(f"openFDA API 请求失败：{exc}")

        results = data.get("results") or []
        hits = [_parse_fda_label(r) for r in results if isinstance(r, dict)]
        if not hits:
            return empty_result(query=q)

        cols = [
            "BrandName",
            "GenericName",
            "Indications",
            "Warnings",
            "Dosage",
            "AdverseReactions",
            "SetID",
        ]
        top = hits[0]
        md_lines = [
            f"### FDA 药品标签检索（{len(hits)} 条）\n",
            f"**查询**：`{q}`\n",
            f"**首条**：{top.get('BrandName')} / {top.get('GenericName')}\n",
            f"\n**适应症节选**：\n> {top.get('Indications')}\n",
        ]

        return success_payload(
            f"openFDA 命中 {len(hits)} 条标签",
            markdown="\n".join(md_lines),
            hits=hits,
            table_data=table_data_from_rows(cols, hits),
            query=q,
            result_count=len(hits),
        )
