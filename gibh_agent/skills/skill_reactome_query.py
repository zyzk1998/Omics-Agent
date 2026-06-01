# -*- coding: utf-8 -*-
"""Reactome 数据库查询 — ContentService REST API。"""
from __future__ import annotations

import re
from html import unescape
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

_REACTOME_SEARCH = "https://reactome.org/ContentService/search/query"


def _strip_html(text: str) -> str:
    s = unescape(re.sub(r"<[^>]+>", "", text or ""))
    return re.sub(r"\s+", " ", s).strip()


def _parse_reactome_entry(entry: Dict[str, Any]) -> Dict[str, str]:
    name = _strip_html(str(entry.get("name") or ""))
    summ = _strip_html(str(entry.get("summation") or ""))[:200]
    st_id = str(entry.get("stId") or entry.get("id") or "—")
    species = entry.get("species") or []
    sp = species[0] if isinstance(species, list) and species else "Homo sapiens"
    return {
        "StableID": st_id,
        "Name": name or "—",
        "Type": str(entry.get("type") or entry.get("exactType") or "Pathway"),
        "Species": str(sp),
        "Summary": summ or "—",
        "URL": f"https://reactome.org/content/detail/{st_id}" if st_id != "—" else "—",
    }


class ReactomeQuerySkill(BaseSkill):
    __abstractskill__ = False

    skill_id = "reactome_query"
    display_name = "Reactome 数据库查询"
    description = "检索 Reactome 人类通路层级、稳定 ID 与摘要（ContentService REST）。"
    category = "生物医药"
    sub_category = "信息检索"
    aliases = ["Reactome", "通路", "reactome", "信号通路"]
    required_parameters = ["query"]
    tool_chain_key = "reactome"
    output_type = "mixed"
    __dependencies__ = ["pip:httpx"]

    def execute(self, query: str = "", limit: int = 10, **kwargs: Any) -> Dict[str, Any]:
        filled = apply_launch_demo_defaults(
            self.skill_id,
            {"query": (query or kwargs.get("pathway") or "").strip(), "limit": limit},
        )
        q = str(filled.get("query") or "").strip()
        if not q:
            return error_result("请提供 query（通路名或关键词，如 apoptosis）")

        lim = max(1, min(int(filled.get("limit") or 10), 20))
        params = {
            "query": q,
            "species": "Homo sapiens",
            "types": "Pathway",
            "cluster": "true",
        }
        try:
            with httpx.Client(timeout=DEFAULT_HTTP_TIMEOUT, follow_redirects=True) as client:
                resp = client.get(_REACTOME_SEARCH, params=params)
                resp.raise_for_status()
                data = resp.json()
        except httpx.TimeoutException:
            return timeout_result(query=q)
        except httpx.HTTPError as exc:
            return error_result(f"Reactome API 请求失败：{exc}")

        results = data.get("results") or []
        hits: List[Dict[str, Any]] = []
        for group in results:
            if not isinstance(group, dict):
                continue
            for entry in group.get("entries") or []:
                if isinstance(entry, dict):
                    row = _parse_reactome_entry(entry)
                    if row.get("StableID") != "—":
                        hits.append(row)
                if len(hits) >= lim:
                    break
            if len(hits) >= lim:
                break

        if not hits:
            return empty_result(query=q)

        cols = ["StableID", "Name", "Type", "Species", "Summary", "URL"]
        md_lines = [
            f"### Reactome 通路检索（{len(hits)} 条）\n",
            f"**关键词**：`{q}`\n",
        ]
        for i, row in enumerate(hits[:8], 1):
            md_lines.append(
                f"{i}. **{row.get('Name')}** — `{row.get('StableID')}` "
                f"([详情]({row.get('URL')})\n"
            )
            if row.get("Summary") and row["Summary"] != "—":
                md_lines.append(f"   {row['Summary'][:160]}…\n")

        return success_payload(
            f"Reactome 命中 {len(hits)} 条通路",
            markdown="\n".join(md_lines),
            hits=hits,
            table_data=table_data_from_rows(cols, hits),
            query=q,
            result_count=len(hits),
        )
