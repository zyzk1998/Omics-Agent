"""asset_locator 单元测试。"""
from __future__ import annotations

import os
import tempfile
from pathlib import Path

import pytest

from gibh_agent.core.asset_locator import (
    AssetLocatorError,
    normalize_uploaded_file_entry,
    resolve_asset_path,
)


def test_resolve_relative_path_under_upload_dir():
    with tempfile.TemporaryDirectory() as td:
        f = Path(td) / "owner" / "batch" / "brain_mri.nii.gz"
        f.parent.mkdir(parents=True)
        f.write_bytes(b"nii")
        rel = "owner/batch/brain_mri.nii.gz"
        got = resolve_asset_path(rel, upload_dir=td, must_exist=True)
        assert got == str(f.resolve())


def test_normalize_entry_preserves_asset_id():
    with tempfile.TemporaryDirectory() as td:
        f = Path(td) / "sample.csv"
        f.write_text("a,b\n1,2\n")
        out = normalize_uploaded_file_entry(
            {"path": str(f), "name": "sample.csv", "asset_id": 42, "source": "data_asset"},
            upload_dir=td,
            strict=True,
        )
        assert out["asset_id"] == 42
        assert out["path"] == str(f.resolve())


def test_missing_file_strict_raises():
    with pytest.raises(AssetLocatorError):
        resolve_asset_path(
            "/app/uploads/no_such_dir_abc/file.nii.gz",
            upload_dir=tempfile.gettempdir(),
            must_exist=True,
        )


def test_basename_rglob_scoped_to_owner_only():
    """同名文件存在于不同 owner 目录时，不得跨用户解析。"""
    with tempfile.TemporaryDirectory() as td:
        ud = Path(td)
        user_a = ud / "user_a" / "batch"
        user_b = ud / "user_b" / "batch"
        user_a.mkdir(parents=True)
        user_b.mkdir(parents=True)
        file_a = user_a / "brain.nii.gz"
        file_b = user_b / "brain.nii.gz"
        file_a.write_bytes(b"A")
        file_b.write_bytes(b"B")

        got_a = resolve_asset_path(
            "brain.nii.gz",
            upload_dir=td,
            owner_id="user_a",
            must_exist=True,
        )
        got_b = resolve_asset_path(
            "brain.nii.gz",
            upload_dir=td,
            owner_id="user_b",
            must_exist=True,
        )
        assert got_a == str(file_a.resolve())
        assert got_b == str(file_b.resolve())
        assert got_a != got_b

        with pytest.raises(AssetLocatorError):
            resolve_asset_path(
                "brain.nii.gz",
                upload_dir=td,
                owner_id=None,
                must_exist=True,
            )
