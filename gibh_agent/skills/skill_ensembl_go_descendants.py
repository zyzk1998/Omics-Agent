# -*- coding: utf-8
"""Ensembl GO 术语后代查询 — QuickGO + Ensembl Ontology REST。"""
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

_QUICKGO_DESC = "https://www.ebi.ac.uk/QuickGO/services/ontology/go/terms/{go_id}/descendants"
_ENSEMBL_ONT = "https://rest.ensembl.org/ontology/id/{go_id}"
_HEADERS = {"Content-Type": "application/json", "Accept": "application/json"}


def _normalize_go_id(q: str) -> str:
    s = q.strip().upper()
    if s.startswith("GO:"):
        return s
    if re.fullmatch(r"\d{7}", s):
        return f"GO:{s}"
    return s


class EnsemblGoDescendantsSkill(BaseSkill):
    __abstractskill__ = False

    skill_id = "ensembl_go_descendants"
    display_name = "Ensembl GO术语后代查询"
    description = "按 GO ID 查询术语定义与后代/子术语列表（QuickGO + Ensembl REST）。"
    category = "生物医药"
    sub_category = "信息检索"
    aliases = ["GO", "Gene Ontology", "GO术语", "ensemblGo"]
    required_parameters = ["query"]
    tool_chain_key = "ensemblGo"
    output_type = "mixed"
    __dependencies__ = ["pip:httpx"]

    def execute(self, query: str = "", limit: int = 20, **kwargs: Any) -> Dict[str, Any]:
        filled = apply_launch_demo_defaults(
            self.skill_id,
            {"query": (query or kwargs.get("go_id") or "").strip(), "limit": limit},
        )
        q = str(filled.get("query") or "").strip()
        if not q:
            return error_result("请提供 query（GO ID，如 GO:0006915）")

        go_id = _normalize_go_id(q)
        if not go_id.startswith("GO:"):
            return error_result("GO ID 格式无效，示例：GO:0006915")

        lim = max(1, min(int(filled.get("limit") or 20), 50))
        try:
            with httpx.Client(timeout=DEFAULT_HTTP_TIMEOUT, follow_redirects=True) as client:
                term_resp = client.get(_ENSEMBL_ONT.format(go_id=go_id), headers=_HEADERS)
                if term_resp.status_code == 404:
                    return empty_result(query=go_id)
                term_resp.raise_for_status()
                term = term_resp.json()

                desc_resp = client.get(_QUICKGO_DESC.format(go_id=go_id))
                desc_resp.raise_for_status()
                desc_data = desc_resp.json()
        except httpx.TimeoutException:
            return timeout_result(query=go_id)
        except httpx.HTTPError as exc:
            return error_result(f"GO 术语 API 请求失败：{exc}")

        term_name = str(term.get("name") or "—")
        namespace = str(term.get("namespace") or "—")
        desc_ids: List[str] = []
        results = (desc_data.get("results") or [])
        if results and isinstance(results[0], dict):
            raw = results[0].get("descendants") or []
            desc_ids = [str(x) for x in raw if x][:lim]

        direct_children = term.get("children") or []
        child_rows: List[Dict[str, str]] = []
        for ch in direct_children[:lim]:
            if isinstance(ch, dict):
                child_rows.append(
                    {
                        "GOID": str(ch.get("accession") or "—"),
                        "Name": str(ch.get("name") or "—")[:100],
                        "Namespace": str(ch.get("namespace") or namespace),
                        "Relation": "direct_child",
                    }
                )

        for did in desc_ids:
            if any(r["GOID"] == did for r in child_rows):
                continue
            child_rows.append(
                {"GOID": did, "Name": "—", "Namespace": namespace, "Relation": "descendant"}
            )
            if len(child_rows) >= lim:
                break

        if not child_rows:
            return empty_result(query=go_id)

        cols = ["GOID", "Name", "Namespace", "Relation"]
        md_lines = [
            f"### GO 术语后代查询\n",
            f"**术语**：{term_name}（`{go_id}`）· {namespace}\n",
            f"**后代/子术语**：展示 {len(child_rows)} 条（含直接子术语与 QuickGO 后代 ID）\n",
        ]
        for i, row in enumerate(child_rows[:10], 1):
            md_lines.append(f"{i}. `{row['GOID']}` — {row['Name']} ({row['Relation']})\n")

        return success_payload(
            f"GO {go_id} 后代/子术语 {len(child_rows)} 条",
            markdown="\n".join(md_lines),
            hits=child_rows,
            table_data=table_data_from_rows(cols, child_rows),
            query=go_id,
            term_name=term_name,
            descendant_count=len(desc_ids),
            result_count=len(child_rows),
        )
