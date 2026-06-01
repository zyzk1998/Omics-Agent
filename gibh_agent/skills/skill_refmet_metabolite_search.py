# -*- coding: utf-8
"""代谢物名称检索 — RefMet API（Metabolomics Workbench）。"""
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

_REFMET_MATCH = "https://www.metabolomicsworkbench.org/rest/refmet/match/{name}/{threshold}/json"
_REFMET_NAME = "https://www.metabolomicsworkbench.org/rest/refmet/name/{name}/all/json"


def _parse_refmet_item(item: Dict[str, Any]) -> Dict[str, str]:
    return {
        "RefMetName": str(item.get("refmet_name") or item.get("name") or "—"),
        "Formula": str(item.get("formula") or "—"),
        "ExactMass": str(item.get("exactmass") or "—"),
        "SuperClass": str(item.get("super_class") or "—"),
        "MainClass": str(item.get("main_class") or "—"),
        "SubClass": str(item.get("sub_class") or "—"),
        "PubChemCID": str(item.get("pubchem_cid") or "—"),
        "InChIKey": str(item.get("inchi_key") or "—"),
    }


class RefmetMetaboliteSearchSkill(BaseSkill):
    __abstractskill__ = False

    skill_id = "refmet_metabolite_search"
    display_name = "代谢物名称检索"
    description = "RefMet 标准化代谢物名称、分子式与分类（Metabolomics Workbench REST）。"
    category = "化学"
    sub_category = "信息检索"
    aliases = ["RefMet", "代谢物", "refmet", "代谢组"]
    required_parameters = ["query"]
    tool_chain_key = "refmet"
    output_type = "mixed"
    __dependencies__ = ["pip:httpx"]

    def execute(self, query: str = "", match_threshold: float = 0.8, **kwargs: Any) -> Dict[str, Any]:
        filled = apply_launch_demo_defaults(
            self.skill_id,
            {
                "query": (query or kwargs.get("metabolite_name") or "").strip(),
                "match_threshold": match_threshold,
            },
        )
        q = str(filled.get("query") or "").strip()
        if not q:
            return error_result("请提供 query（代谢物名称）")

        threshold = max(0.5, min(float(filled.get("match_threshold") or 0.8), 1.0))
        try:
            with httpx.Client(timeout=DEFAULT_HTTP_TIMEOUT, follow_redirects=True) as client:
                resp = client.get(_REFMET_MATCH.format(name=q, threshold=threshold))
                if resp.status_code == 404:
                    resp = client.get(_REFMET_NAME.format(name=q))
                if resp.status_code == 404:
                    return empty_result(query=q)
                resp.raise_for_status()
                data = resp.json()
        except httpx.TimeoutException:
            return timeout_result(query=q)
        except httpx.HTTPError as exc:
            return error_result(f"RefMet API 请求失败：{exc}")

        items: List[Dict[str, Any]] = []
        if isinstance(data, list):
            items = [x for x in data if isinstance(x, dict)]
        elif isinstance(data, dict):
            items = [data]

        hits = [_parse_refmet_item(it) for it in items[:15]]
        if not hits:
            return empty_result(query=q)

        cols = [
            "RefMetName",
            "Formula",
            "ExactMass",
            "SuperClass",
            "MainClass",
            "SubClass",
            "PubChemCID",
            "InChIKey",
        ]
        md = (
            f"### RefMet 代谢物检索（{len(hits)} 条）\n\n"
            f"**查询**：`{q}`\n\n"
            f"**首条**：{hits[0]['RefMetName']} · {hits[0]['Formula']}\n"
        )

        return success_payload(
            f"RefMet 命中 {len(hits)} 条",
            markdown=md,
            hits=hits,
            table_data=table_data_from_rows(cols, hits),
            query=q,
            result_count=len(hits),
        )
