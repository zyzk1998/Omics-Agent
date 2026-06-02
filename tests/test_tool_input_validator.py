"""tool_input_validator 单元测试。"""
from __future__ import annotations

import os
import tempfile
from pathlib import Path

import pytest

from gibh_agent.core.tool_input_validator import validate_and_normalize_inputs


def _fake_rna_validation(adata_path: str) -> dict:
    return {"status": "success", "adata_path": adata_path}


_fake_rna_validation.__tool_name__ = "rna_data_validation"


def test_fallback_fills_empty_adata_path_from_uploaded_files():
    with tempfile.TemporaryDirectory() as td:
        tenx_dir = Path(td) / "10x_data_test"
        tenx_dir.mkdir()
        (tenx_dir / "matrix.mtx").write_text("mock")
        params = {"adata_path": ""}
        ctx = {
            "upload_dir": td,
            "uploaded_files": [
                {"is_10x": True, "group_dir": str(tenx_dir), "path": str(tenx_dir)}
            ],
        }
        out, err = validate_and_normalize_inputs(
            "rna_data_validation",
            params,
            tool_func=_fake_rna_validation,
            context=ctx,
        )
        assert err is None
        assert out["adata_path"] == str(tenx_dir.resolve())


def test_file_path_mapped_to_adata_path():
    with tempfile.NamedTemporaryFile(suffix=".h5ad", delete=False) as f:
        f.write(b"x")
        h5 = f.name
    try:
        params = {"file_path": h5}
        out, err = validate_and_normalize_inputs(
            "rna_data_validation",
            params,
            tool_func=_fake_rna_validation,
            context={"upload_dir": os.path.dirname(h5)},
        )
        assert err is None
        assert out["adata_path"] == h5
    finally:
        os.unlink(h5)


def test_missing_path_returns_structured_error():
    params = {"adata_path": "/app/uploads/nonexistent_dir_abc123/matrix.mtx"}
    out, err = validate_and_normalize_inputs(
        "rna_data_validation",
        params,
        tool_func=_fake_rna_validation,
        context={"upload_dir": "/app/uploads"},
    )
    assert err is not None
    assert err["status"] == "error"
    assert err["error_category"] == "input_path_not_found"
    assert "adata_path" in err["user_message"]


def test_chain_h5ad_replaces_upload_dir_for_downstream_tool():
    with tempfile.TemporaryDirectory() as td:
        h5 = Path(td) / "filtered.h5ad"
        h5.write_bytes(b"mock")
        tenx_dir = Path(td) / "10x_data_test"
        tenx_dir.mkdir()
        params = {"adata_path": str(tenx_dir)}
        ctx = {
            "upload_dir": td,
            "current_file_path": str(h5),
            "uploaded_files": [{"is_10x": True, "group_dir": str(tenx_dir), "path": str(tenx_dir)}],
        }
        out, err = validate_and_normalize_inputs(
            "rna_hvg",
            params,
            tool_func=_fake_rna_validation,
            context=ctx,
        )
        assert err is None
        assert out["adata_path"] == str(h5.resolve())
