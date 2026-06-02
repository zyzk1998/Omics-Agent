"""omics_io_registry 与全局 I/O 自愈单元测试。"""
from __future__ import annotations

import tempfile
from pathlib import Path

import pytest

from gibh_agent.core.omics_io_registry import (
    fuzzy_heal_path_by_basename,
    get_extension_warning,
    is_known_omics_extension,
)
from gibh_agent.core.tool_input_validator import validate_and_normalize_inputs
from gibh_agent.utils.path_resolver import validate_inspector_file_format as gate_format


def _fake_tool(file_path: str) -> dict:
    return {"status": "success", "file_path": file_path}


def test_fastq_extension_is_known():
    assert is_known_omics_extension("SRR622461_1.fastq.gz")
    assert get_extension_warning("SRR622461_1.fastq.gz") is None


def test_unknown_extension_returns_warning_not_error():
    warn = get_extension_warning("sample.weirdomics")
    assert warn is not None
    assert gate_format("sample.weirdomics") is None


def test_fuzzy_heal_wrong_absolute_path_by_basename():
    with tempfile.TemporaryDirectory() as td:
        real_dir = Path(td) / "user_a" / "20250602_batch"
        real_dir.mkdir(parents=True)
        fq = real_dir / "SRR622461_1.fastq.gz"
        fq.write_bytes(b"mock")
        wrong_abs = f"/app/uploads/user_a/old_timestamp/SRR622461_1.fastq.gz"
        healed, _ = fuzzy_heal_path_by_basename(
            wrong_abs,
            upload_dir=Path(td),
            owner_id="user_a",
            context_paths=[str(fq)],
            uploaded_entries=[{"name": "SRR622461_1.fastq.gz", "path": str(fq)}],
        )
        assert healed == str(fq.resolve())


def test_validator_fuzzy_heals_llm_wrong_path():
    with tempfile.TemporaryDirectory() as td:
        real = Path(td) / "batch" / "SRR622461_2.fastq.gz"
        real.parent.mkdir(parents=True)
        real.write_bytes(b"x")
        wrong = f"/app/uploads/stale_dir/SRR622461_2.fastq.gz"
        params = {"file_path": wrong}
        ctx = {
            "upload_dir": td,
            "owner_id": "batch",
            "uploaded_files": [{"path": str(real), "name": real.name}],
        }
        out, err = validate_and_normalize_inputs(
            "genomics_fastp_qc",
            params,
            tool_func=_fake_tool,
            context=ctx,
        )
        assert err is None
        assert out["file_path"] == str(real.resolve())
