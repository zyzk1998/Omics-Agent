# -*- coding: utf-8
"""药物标识符交叉检索 — PubChem / ChEMBL / Registry 多源 HTTP 宽表（超时友好）。"""
from __future__ import annotations

import re
from typing import Any, Dict, List, Optional

import httpx

from gibh_agent.skills._skill_common import DEFAULT_HTTP_TIMEOUT, pubchem_smiles_url
from gibh_agent.skills._skill_payload import (
    empty_result,
    error_result,
    success_payload,
    table_data_from_rows,
    timeout_result,
)
from gibh_agent.skills.base_skill import BaseSkill
from gibh_agent.skills.launch_skill_demos import apply_launch_demo_defaults

_PUBCHEM_XREFS = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/xrefs/RegistryID,RN/JSON"
_CHEMBL_SEARCH = "https://www.ebi.ac.uk/chembl/api/data/molecule/search.json"


def _resolve_cid(client: httpx.Client, query: str, smiles: str) -> Optional[int]:
    from urllib.parse import quote

    if smiles:
        resp = client.get(pubchem_smiles_url(smiles, "cids/JSON"))
        if resp.status_code == 404:
            return None
        resp.raise_for_status()
        cids = (resp.json().get("IdentifierList") or {}).get("CID") or []
        return int(cids[0]) if cids else None
    q = query.strip()
    if re.fullmatch(r"\d+", q):
        return int(q)
    resp = client.get(
        f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{quote(q, safe='')}/cids/JSON"
    )
    if resp.status_code == 404:
        return None
    resp.raise_for_status()
    cids = (resp.json().get("IdentifierList") or {}).get("CID") or []
    return int(cids[0]) if cids else None


def _fetch_chembl_id(client: httpx.Client, query: str) -> str:
    resp = client.get(_CHEMBL_SEARCH, params={"q": query, "limit": 1})
    if resp.status_code != 200:
        return "—"
    mols = (resp.json().get("molecules") or [])
    if not mols:
        return "—"
    return str(mols[0].get("molecule_chembl_id") or "—")


class DrugIdCrossrefSkill(BaseSkill):
    __abstractskill__ = False

    skill_id = "drug_id_crossref"
    display_name = "药物标识符交叉检索"
    description = "按药名或 SMILES 聚合 PubChem CID、ChEMBL ID 与 Registry 交叉引用（多源 HTTP 宽表）。"
    category = "化学"
    sub_category = "信息检索"
    aliases = ["药物ID", "交叉引用", "DrugBank", "标识符互证"]
    required_parameters = ["query"]
    tool_chain_key = "drug"
    output_type = "mixed"
    __dependencies__ = ["pip:httpx"]

    def execute(
        self,
        query: str = "",
        smiles: str = "",
        **kwargs: Any,
    ) -> Dict[str, Any]:
        filled = apply_launch_demo_defaults(
            self.skill_id,
            {
                "query": (query or kwargs.get("drug_name") or "").strip(),
                "smiles": (smiles or "").strip(),
            },
        )
        q = str(filled.get("query") or "").strip()
        smi = str(filled.get("smiles") or "").strip()
        if not q and not smi:
            return error_result("请提供 query（药名）或 smiles")

        label = smi or q
        timeout = httpx.Timeout(min(DEFAULT_HTTP_TIMEOUT, 45))
        try:
            with httpx.Client(timeout=timeout, follow_redirects=True) as client:
                cid = _resolve_cid(client, q, smi)
                chembl_id = _fetch_chembl_id(client, q or "aspirin")
                registry_ids: List[str] = []
                if cid:
                    xr = client.get(_PUBCHEM_XREFS.format(cid=cid))
                    if xr.status_code == 200:
                        info = (xr.json().get("InformationList") or {}).get("Information") or []
                        if info and isinstance(info[0], dict):
                            registry_ids = [str(x) for x in (info[0].get("RegistryID") or [])[:5]]
        except httpx.TimeoutException:
            return timeout_result(query=label)
        except httpx.HTTPError as exc:
            return error_result(f"药物标识符交叉检索失败：{exc}")

        if not cid and chembl_id == "—":
            return empty_result(query=label)

        row = {
            "Query": label,
            "PubChemCID": str(cid) if cid else "—",
            "ChEMBLID": chembl_id,
            "RegistryIDs": ", ".join(registry_ids[:3]) if registry_ids else "—",
            "PubChemURL": f"https://pubchem.ncbi.nlm.nih.gov/compound/{cid}" if cid else "—",
            "ChEMBLURL": (
                f"https://www.ebi.ac.uk/chembl/compound_report_card/{chembl_id}/"
                if chembl_id != "—"
                else "—"
            ),
        }
        hits = [row]
        cols = ["Query", "PubChemCID", "ChEMBLID", "RegistryIDs", "PubChemURL", "ChEMBLURL"]
        md = (
            f"### 药物标识符交叉检索\n\n"
            f"**查询**：`{label}`\n\n"
            f"| 来源 | ID |\n|------|----|\n"
            f"| PubChem | {row['PubChemCID']} |\n"
            f"| ChEMBL | {row['ChEMBLID']} |\n"
        )

        return success_payload(
            f"药物标识符交叉引用已聚合",
            markdown=md,
            hits=hits,
            table_data=table_data_from_rows(cols, hits),
            query=label,
            pubchem_cid=cid,
            chembl_id=chembl_id,
        )
