"""Planner fallback / metadata overrides for genomics / proteomics / epigenomics."""

from unittest.mock import MagicMock

import pytest

from gibh_agent.core.planner import SOPPlanner


@pytest.fixture
def planner() -> SOPPlanner:
    return SOPPlanner(MagicMock(), MagicMock())


def test_fallback_mzml_proteomics(planner: SOPPlanner) -> None:
    r = planner._fallback_intent_classification(
        "蛋白质鉴定",
        True,
        {"file_path": "/data/sample.mzML", "file_type": "proteomics_ms"},
    )
    assert r["domain_name"] == "proteomics"


def test_fallback_fastq_epigenomics(planner: SOPPlanner) -> None:
    r = planner._fallback_intent_classification(
        "ATAC-seq peak calling macs2",
        True,
        {"file_path": "/data/a.fastq.gz", "file_type": "fastq"},
    )
    assert r["domain_name"] == "epigenomics"


def test_fallback_fastq_genomics(planner: SOPPlanner) -> None:
    r = planner._fallback_intent_classification(
        "WGS variant calling with GATK",
        True,
        {"file_path": "/reads.fastq.gz", "file_type": "fastq"},
    )
    assert r["domain_name"] == "genomics"


def test_fallback_fastq_default_rna(planner: SOPPlanner) -> None:
    r = planner._fallback_intent_classification(
        "帮我分析一下",
        True,
        {"file_path": "/reads.fastq.gz", "file_type": "fastq"},
    )
    assert r["domain_name"] == "RNA"
