# -*- coding: utf-8
"""OpenTargets 按 ChEMBL ID 获取父分子和子分子 — GraphQL API。"""
from __future__ import annotations

import re
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

_OT_GRAPHQL = "https://api.platform.opentargets.org/api/v4/graphql"
_OT_QUERY = """
query DrugHierarchy($chemblId: String!) {
  drug(chemblId: $chemblId) {
    id
    name
    parentMolecule { id name }
    childMolecules { id name }
  }
}
"""


def _normalize_chembl_id(q: str) -> str:
    s = q.strip().upper()
    if s.startswith("CHEMBL"):
        return s
    if re.fullmatch(r"\d+", s):
        return f"CHEMBL{s}"
    return s


class OpentargetsChemblHierarchySkill(BaseSkill):
    __abstractskill__ = False

    skill_id = "opentargets_chembl_hierarchy"
    display_name = "OpenTargets_按ChEMBL ID获取父分子和子分子"
    description = "Open Targets GraphQL：按 ChEMBL ID 查询药物分子父子层级关系。"
    category = "化学"
    sub_category = "信息检索"
    aliases = ["OpenTargets", "ChEMBL层级", "opentargets", "药物层级"]
    required_parameters = ["query"]
    tool_chain_key = "opentargets"
    output_type = "mixed"
    __dependencies__ = ["pip:httpx"]

    def execute(self, query: str = "", chembl_id: str = "", **kwargs: Any) -> Dict[str, Any]:
        filled = apply_launch_demo_defaults(
            self.skill_id,
            {"query": (query or chembl_id or "").strip()},
        )
        q = _normalize_chembl_id(str(filled.get("query") or ""))
        if not q or not q.startswith("CHEMBL"):
            return error_result("请提供 query（ChEMBL ID，如 CHEMBL25）")

        try:
            with httpx.Client(timeout=DEFAULT_HTTP_TIMEOUT, follow_redirects=True) as client:
                resp = client.post(
                    _OT_GRAPHQL,
                    json={"query": _OT_QUERY, "variables": {"chemblId": q}},
                )
                resp.raise_for_status()
                payload = resp.json()
        except httpx.TimeoutException:
            return timeout_result(chembl_id=q)
        except httpx.HTTPError as exc:
            return error_result(f"Open Targets API 请求失败：{exc}")

        if payload.get("errors"):
            return error_result(str(payload["errors"][0].get("message") or payload["errors"]))

        drug = (payload.get("data") or {}).get("drug")
        if not drug:
            return empty_result(chembl_id=q)

        hits: List[Dict[str, str]] = []
        parent = drug.get("parentMolecule")
        if isinstance(parent, dict) and parent.get("id"):
            hits.append(
                {
                    "Relation": "parent",
                    "ChEMBLID": str(parent.get("id")),
                    "Name": str(parent.get("name") or "—"),
                }
            )
        for child in drug.get("childMolecules") or []:
            if isinstance(child, dict) and child.get("id"):
                hits.append(
                    {
                        "Relation": "child",
                        "ChEMBLID": str(child.get("id")),
                        "Name": str(child.get("name") or "—"),
                    }
                )

        drug_name = str(drug.get("name") or q)
        if not hits:
            hits.append({"Relation": "self", "ChEMBLID": q, "Name": drug_name})

        cols = ["Relation", "ChEMBLID", "Name"]
        cards = metrics_cards_from_pairs(
            [
                ("查询 ChEMBL", q, ""),
                ("药物名称", drug_name, ""),
                ("子分子数", sum(1 for h in hits if h["Relation"] == "child"), ""),
            ]
        )
        md_lines = [
            f"### Open Targets 分子层级 — {drug_name}\n",
            f"**ChEMBL**：`{q}`\n",
        ]
        if parent and isinstance(parent, dict) and parent.get("id"):
            md_lines.append(f"- **父分子**：{parent.get('name')} (`{parent.get('id')}`)\n")
        children = [h for h in hits if h["Relation"] == "child"]
        if children:
            md_lines.append(f"- **子分子**（{len(children)}）：\n")
            for c in children[:8]:
                md_lines.append(f"  - {c['Name']} (`{c['ChEMBLID']}`)\n")

        return success_payload(
            f"Open Targets {q} 层级已获取",
            markdown="\n".join(md_lines),
            hits=hits,
            table_data=table_data_from_rows(cols, hits),
            metrics_cards=cards,
            chembl_id=q,
            drug_name=drug_name,
            result_count=len(hits),
        )
