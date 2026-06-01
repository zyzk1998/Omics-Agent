# -*- coding: utf-8 -*-
"""阵列二 Batch1：注册、委托与异常载荷单元测试（无真实外网）。"""
from __future__ import annotations

import importlib
from unittest.mock import MagicMock, patch

import httpx
import pytest

from gibh_agent.core.skill_registry import discover_and_register_skills
from gibh_agent.core.tool_registry import registry
from gibh_agent.skills.launch_skill_demos import LAUNCH_ISOLATED_TOOL_IDS, apply_launch_demo_defaults
from gibh_agent.skills.launch_skill_isolated import should_delegate_to_launch_worker


@pytest.fixture(scope="module", autouse=True)
def _register():
    discover_and_register_skills(force_reload=True)


@pytest.mark.parametrize(
    "tid",
    [
        "pubmed_query",
        "uniprot_query",
        "smiles_to_cid",
        "seq_format_converter",
        "calc_molecular_weight",
    ],
)
def test_array2_tools_registered_and_isolated(tid: str, monkeypatch):
    assert registry.get_tool(tid) is not None
    assert tid in LAUNCH_ISOLATED_TOOL_IDS
    monkeypatch.setenv("LAUNCH_SKILLS_BASE_URL", "http://127.0.0.1:8010")
    assert should_delegate_to_launch_worker(tid) is True


def test_demo_defaults_fill_empty():
    for tid in (
        "pubmed_query",
        "uniprot_query",
        "smiles_to_cid",
        "seq_format_converter",
        "calc_molecular_weight",
    ):
        out = apply_launch_demo_defaults(tid, {})
        assert out


def _smiles_mod():
    return importlib.import_module("gibh_agent.skills.skill_smiles_to_cid")


def test_smiles_to_cid_timeout(monkeypatch):
    mod = _smiles_mod()

    class _FakeClient:
        def __enter__(self):
            return self

        def __exit__(self, *args):
            return False

        def get(self, *args, **kwargs):
            raise httpx.TimeoutException("timeout")

    monkeypatch.setattr(mod.httpx, "Client", lambda **kw: _FakeClient())
    res = mod.SmilesToCidSkill().execute(smiles="CCO")
    assert res["status"] == "error"
    assert "超时" in res["message"]


def test_smiles_to_cid_empty(monkeypatch):
    mod = _smiles_mod()

    class _FakeClient:
        def __enter__(self):
            return self

        def __exit__(self, *args):
            return False

        def get(self, *args, **kwargs):
            resp = MagicMock()
            resp.status_code = 404
            return resp

    monkeypatch.setattr(mod.httpx, "Client", lambda **kw: _FakeClient())
    res = mod.SmilesToCidSkill().execute(smiles="invalid_smiles_xyz")
    assert res["status"] == "empty"
    assert "未能根据" in res["message"]


def test_smiles_to_cid_success_payload(monkeypatch):
    mod = _smiles_mod()

    class _FakeClient:
        def __enter__(self):
            return self

        def __exit__(self, *args):
            return False

        def get(self, *args, **kwargs):
            resp = MagicMock()
            resp.status_code = 200
            resp.json.return_value = {"IdentifierList": {"CID": [2244, 338]}}
            return resp

    monkeypatch.setattr(mod.httpx, "Client", lambda **kw: _FakeClient())
    res = mod.SmilesToCidSkill().execute(smiles="CC(=O)Oc1ccccc1C(=O)O")
    assert res["status"] == "success"
    assert res["hits"]
    assert res["table_data"]["columns"]
    assert res["metrics_cards"]


def test_pubmed_esummary_list_root_parsing():
    """NCBI 当前 esummary 根节点为 ListElement，非 DocumentSummarySet 包裹。"""
    from gibh_agent.skills._entrez_util import iter_pubmed_esummary_docs, parse_pubmed_esummary_doc
    from gibh_agent.skills.skill_pubmed_query import _fetch_pubmed_hits

    mock_summaries = [
        {
            "Id": "12345",
            "Title": "Demo paper",
            "FullJournalName": "Demo Journal",
            "PubDate": "2024 Jan",
            "AuthorList": ["Zhang A", "Li B"],
        }
    ]
    docs = list(iter_pubmed_esummary_docs(mock_summaries))
    assert len(docs) == 1
    row = parse_pubmed_esummary_doc(docs[0])
    assert row["PMID"] == "12345"
    assert "Zhang A" in row["Authors"]

    from Bio import Entrez

    with patch.object(Entrez, "read") as mock_read, patch.object(Entrez, "esearch") as mock_esearch, patch.object(
        Entrez, "esummary"
    ) as mock_esummary:
        mock_read.side_effect = lambda h: (
            {"IdList": ["12345"]} if mock_read.call_count == 1 else mock_summaries
        )
        mock_esearch.return_value.__enter__ = lambda s: s
        mock_esearch.return_value.__exit__ = lambda *a: None
        mock_esummary.return_value.__enter__ = lambda s: s
        mock_esummary.return_value.__exit__ = lambda *a: None

        hits = _fetch_pubmed_hits("demo", 5)
    assert len(hits) == 1
    assert hits[0]["Authors"] == "Zhang A, Li B"


def test_uniprot_parse_search_shape():
    from gibh_agent.skills.skill_uniprot_query import _parse_uniprot_entry

    row = _parse_uniprot_entry(
        {
            "primaryAccession": "P04637",
            "organism": {"scientificName": "Homo sapiens"},
            "proteinDescription": {
                "recommendedName": {"fullName": {"value": "Cellular tumor antigen p53"}}
            },
            "sequence": {"length": 393},
            "genes": [{"geneName": {"value": "TP53"}}],
        }
    )
    assert row["Accession"] == "P04637"
    assert row["Organism"] == "Homo sapiens"
    assert row["Length"] == 393
    assert "TP53" in row["Gene"]


def test_calc_molecular_weight_metrics():
    from gibh_agent.skills.skill_calc_molecular_weight import CalcMolecularWeightSkill

    res = CalcMolecularWeightSkill().execute(smiles="CCO")
    assert res["status"] == "success"
    assert res["metrics_cards"]
    assert res["data"]["average_molecular_weight"]
