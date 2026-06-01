# -*- coding: utf-8 -*-
"""阵列二 Batch2a + Batch2b：14 项 DB 技能 — 注册、委托与 mock 解析单测（无真实外网）。"""
from __future__ import annotations

import importlib
from unittest.mock import MagicMock

import httpx
import pytest

from gibh_agent.core.skill_registry import discover_and_register_skills
from gibh_agent.core.tool_registry import registry
from gibh_agent.skills.launch_skill_demos import LAUNCH_ISOLATED_TOOL_IDS, apply_launch_demo_defaults
from gibh_agent.skills.launch_skill_isolated import should_delegate_to_launch_worker

_BATCH2A_IDS = (
    "clinvar_query",
    "dbsnp_query",
    "gwas_catalog_query",
    "reactome_query",
    "interpro_query",
    "geo_query",
    "gene_protein_info_query",
    "mrna_sequence_fetch",
)

_BATCH2B_IDS = (
    "drug_id_crossref",
    "refmet_metabolite_search",
    "rnacentral_query",
    "ensembl_go_descendants",
    "opentargets_chembl_hierarchy",
    "fda_drug_label_search",
)

_ALL_BATCH2 = _BATCH2A_IDS + _BATCH2B_IDS


@pytest.fixture(scope="module", autouse=True)
def _register():
    discover_and_register_skills(force_reload=True)


@pytest.mark.parametrize("tid", _BATCH2A_IDS)
def test_batch2a_tools_registered_and_isolated(tid: str, monkeypatch):
    assert registry.get_tool(tid) is not None
    assert tid in LAUNCH_ISOLATED_TOOL_IDS
    monkeypatch.setenv("LAUNCH_SKILLS_BASE_URL", "http://127.0.0.1:8010")
    assert should_delegate_to_launch_worker(tid) is True


def test_batch2a_demo_defaults_fill_empty():
    for tid in _BATCH2A_IDS:
        out = apply_launch_demo_defaults(tid, {})
        assert out


@pytest.mark.parametrize("tid", _BATCH2B_IDS)
def test_batch2b_tools_registered_and_isolated(tid: str, monkeypatch):
    assert registry.get_tool(tid) is not None
    assert tid in LAUNCH_ISOLATED_TOOL_IDS
    monkeypatch.setenv("LAUNCH_SKILLS_BASE_URL", "http://127.0.0.1:8010")
    assert should_delegate_to_launch_worker(tid) is True


def test_batch2b_demo_defaults_fill_empty():
    for tid in _BATCH2B_IDS:
        out = apply_launch_demo_defaults(tid, {})
        assert out


def test_clinvar_parse_row():
    from gibh_agent.skills.skill_clinvar_query import _parse_clinvar_row

    row = _parse_clinvar_row(
        {
            "accession": "VCV000012345",
            "title": "NM_007294.4(BRCA1):c.3547A>T",
            "gene_sort": "BRCA1",
            "germline_classification": {
                "description": "Pathogenic",
                "review_status": "criteria provided",
                "trait_set": [{"trait_name": "Breast cancer"}],
            },
            "variation_set": [{"variation_loc": [{"chr": "17", "start": "43091984"}]}],
        }
    )
    assert row["Gene"] == "BRCA1"
    assert row["ClinicalSignificance"] == "Pathogenic"
    assert "17" in row["Location"]


def test_clinvar_success_mock(monkeypatch):
    mod = importlib.import_module("gibh_agent.skills.skill_clinvar_query")

    monkeypatch.setattr(
        mod,
        "eutils_esearch",
        lambda db, term, **kw: ["4851990"],
    )
    monkeypatch.setattr(
        mod,
        "eutils_esummary",
        lambda db, ids, **kw: {
            "4851990": {
                "accession": "VCV004851990",
                "title": "Demo variant",
                "gene_sort": "BRCA1",
                "germline_classification": {"description": "Pathogenic", "review_status": "single"},
                "variation_set": [],
            }
        },
    )
    res = mod.ClinvarQuerySkill().execute(query="BRCA1")
    assert res["status"] == "success"
    assert res["hits"]
    assert res["table_data"]["columns"]


def test_dbsnp_success_mock(monkeypatch):
    mod = importlib.import_module("gibh_agent.skills.skill_dbsnp_query")

    monkeypatch.setattr(mod, "eutils_esearch", lambda db, term, **kw: ["386606420"])
    monkeypatch.setattr(
        mod,
        "eutils_esummary",
        lambda db, ids, **kw: {
            "386606420": {
                "snp_id": 699,
                "chr": "1",
                "clinical_significance": "benign",
                "genes": [{"name": "AGT"}],
                "fxn_class": "missense_variant",
            }
        },
    )
    res = mod.DbsnpQuerySkill().execute(query="rs699")
    assert res["status"] == "success"
    assert res["hits"][0]["RsID"] == "rs699"


def test_gwas_trait_mode_mock(monkeypatch):
    mod = importlib.import_module("gibh_agent.skills.skill_gwas_catalog_query")

    def _fake_trait(trait, limit, timeout):
        return [
            {
                "StudyID": "GCST000001",
                "Trait": "type 2 diabetes",
                "InitialSampleSize": "1000 cases",
                "PubmedID": "12345",
                "SnpCount": "500000",
            }
        ]

    monkeypatch.setattr(mod, "_fetch_by_trait", _fake_trait)
    res = mod.GwasCatalogQuerySkill().execute(query="type 2 diabetes")
    assert res["status"] == "success"
    assert res["data"]["search_mode"] == "study"


def test_reactome_parse_entry():
    from gibh_agent.skills.skill_reactome_query import _parse_reactome_entry

    row = _parse_reactome_entry(
        {
            "stId": "R-HSA-109581",
            "name": "<span>Apoptosis</span>",
            "type": "Pathway",
            "species": ["Homo sapiens"],
            "summation": "Apoptosis pathway summary.",
        }
    )
    assert row["StableID"] == "R-HSA-109581"
    assert "Apoptosis" in row["Name"]


def test_gene_protein_info_mock(monkeypatch):
    mod = importlib.import_module("gibh_agent.skills.skill_gene_protein_info_query")

    class _FakeClient:
        def __enter__(self):
            return self

        def __exit__(self, *args):
            return False

        def get(self, url, **kwargs):
            resp = MagicMock()
            resp.status_code = 200
            resp.json.return_value = {
                "id": "ENSG00000141510",
                "display_name": "TP53",
                "biotype": "protein_coding",
                "seq_region_name": "17",
                "start": 7668402,
                "end": 7687550,
                "description": "tumor protein p53",
                "canonical_transcript": "ENST00000269305",
                "Transcript": [{"id": "ENST00000269305", "biotype": "protein_coding"}],
            }
            return resp

    monkeypatch.setattr(mod.httpx, "Client", lambda **kw: _FakeClient())
    res = mod.GeneProteinInfoQuerySkill().execute(query="TP53")
    assert res["status"] == "success"
    assert res["metrics_cards"]


def test_mrna_sequence_fetch_mock(monkeypatch):
    mod = importlib.import_module("gibh_agent.skills.skill_mrna_sequence_fetch")

    class _FakeClient:
        def __enter__(self):
            return self

        def __exit__(self, *args):
            return False

        def get(self, url, **kwargs):
            resp = MagicMock()
            resp.status_code = 200
            if "lookup" in url:
                resp.json.return_value = {
                    "id": "ENSG00000141510",
                    "canonical_transcript": "ENST00000269305",
                    "Transcript": [{"id": "ENST00000269305", "biotype": "protein_coding"}],
                }
            else:
                resp.json.return_value = {"seq": "ATGC" * 20}
            return resp

    monkeypatch.setattr(mod.httpx, "Client", lambda **kw: _FakeClient())
    res = mod.MrnaSequenceFetchSkill().execute(sequence_or_path="TP53")
    assert res["status"] == "success"
    assert res["data"]["fasta"].startswith(">TP53|")
    assert res["data"]["sequence_length"] == 80


def test_clinvar_timeout(monkeypatch):
    mod = importlib.import_module("gibh_agent.skills.skill_clinvar_query")

    def _raise(*a, **k):
        raise httpx.TimeoutException("timeout")

    monkeypatch.setattr(mod, "eutils_esearch", _raise)
    res = mod.ClinvarQuerySkill().execute(query="BRCA1")
    assert res["status"] == "error"
    assert "超时" in res["message"]


def test_opentargets_hierarchy_mock(monkeypatch):
    mod = importlib.import_module("gibh_agent.skills.skill_opentargets_chembl_hierarchy")

    class _FakeClient:
        def __enter__(self):
            return self

        def __exit__(self, *args):
            return False

        def post(self, url, **kwargs):
            resp = MagicMock()
            resp.status_code = 200
            resp.json.return_value = {
                "data": {
                    "drug": {
                        "id": "CHEMBL25",
                        "name": "ASPIRIN",
                        "parentMolecule": None,
                        "childMolecules": [{"id": "CHEMBL1697753", "name": "ASPIRIN DL-LYSINE"}],
                    }
                }
            }
            return resp

    monkeypatch.setattr(mod.httpx, "Client", lambda **kw: _FakeClient())
    res = mod.OpentargetsChemblHierarchySkill().execute(query="CHEMBL25")
    assert res["status"] == "success"
    assert any(h["Relation"] == "child" for h in res["hits"])


def test_ensembl_go_normalize():
    from gibh_agent.skills.skill_ensembl_go_descendants import _normalize_go_id

    assert _normalize_go_id("0006915") == "GO:0006915"
    assert _normalize_go_id("GO:0006915") == "GO:0006915"
