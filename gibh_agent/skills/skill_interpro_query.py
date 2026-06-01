# -*- coding: utf-8 -*-
"""InterPro 数据库查询 — EBI InterPro REST API。"""
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

_INTERPRO_SEARCH = "https://www.ebi.ac.uk/interpro/api/entry/InterPro"


def _parse_interpro_entry(item: Dict[str, Any]) -> Dict[str, str]:
    meta = item.get("metadata") or {}
    if not isinstance(meta, dict):
        meta = {}
    acc = str(meta.get("accession") or "—")
    name_obj = meta.get("name") or {}
    name = name_obj.get("name") if isinstance(name_obj, dict) else str(name_obj or "—")
    entry_type = str(meta.get("type") or meta.get("entry_type") or "—")
    desc = meta.get("description") or []
    desc_text = "—"
    if isinstance(desc, list) and desc:
        first = desc[0]
        if isinstance(first, dict):
            desc_text = str(first.get("text") or "—")[:120]
        else:
            desc_text = str(first)[:120]
    return {
        "Accession": acc,
        "Name": str(name or "—"),
        "Type": entry_type,
        "Description": desc_text,
        "URL": f"https://www.ebi.ac.uk/interpro/entry/InterPro/{acc}/" if acc != "—" else "—",
    }


class InterproQuerySkill(BaseSkill):
    __abstractskill__ = False

    skill_id = "interpro_query"
    display_name = "InterPro 数据库查询"
    description = "检索 InterPro 蛋白结构域、家族与功能分类（EBI REST）。"
    category = "生物医药"
    sub_category = "信息检索"
    aliases = ["InterPro", "结构域", "interpro", "蛋白结构域"]
    required_parameters = ["query"]
    tool_chain_key = "interpro"
    output_type = "mixed"
    __dependencies__ = ["pip:httpx"]

    def execute(self, query: str = "", limit: int = 10, **kwargs: Any) -> Dict[str, Any]:
        filled = apply_launch_demo_defaults(
            self.skill_id,
            {"query": (query or kwargs.get("domain") or "").strip(), "limit": limit},
        )
        q = str(filled.get("query") or "").strip()
        if not q:
            return error_result("请提供 query（结构域名、IPR 编号或关键词）")

        lim = max(1, min(int(filled.get("limit") or 10), 20))
        q_upper = q.upper()
        try:
            with httpx.Client(timeout=DEFAULT_HTTP_TIMEOUT, follow_redirects=True) as client:
                if q_upper.startswith("IPR") and q_upper[3:].isdigit():
                    resp = client.get(
                        f"{_INTERPRO_SEARCH}/{q_upper}",
                        headers={"Accept": "application/json"},
                    )
                    if resp.status_code == 404:
                        return empty_result(query=q)
                    resp.raise_for_status()
                    data = {"results": [{"metadata": resp.json().get("metadata") or resp.json()}]}
                else:
                    resp = client.get(
                        _INTERPRO_SEARCH,
                        params={"search": q, "page": 1, "page_size": lim},
                        headers={"Accept": "application/json"},
                    )
                    resp.raise_for_status()
                    data = resp.json()
        except httpx.TimeoutException:
            return timeout_result(query=q)
        except httpx.HTTPError as exc:
            return error_result(f"InterPro API 请求失败：{exc}")

        hits = [
            _parse_interpro_entry(item)
            for item in (data.get("results") or [])
            if isinstance(item, dict)
        ]

        if not hits:
            return empty_result(query=q)

        cols = ["Accession", "Name", "Type", "Description", "URL"]
        md = (
            f"### InterPro 检索结果（{len(hits)} 条）\n\n"
            f"**查询**：`{q}`\n\n"
            f"**首条**：{hits[0].get('Name')}（`{hits[0].get('Accession')}`）\n"
        )

        return success_payload(
            f"InterPro 命中 {len(hits)} 条",
            markdown=md,
            hits=hits,
            table_data=table_data_from_rows(cols, hits),
            query=q,
            result_count=len(hits),
        )
