"""数据资产 asset_id 查库桥接测试。"""
from __future__ import annotations

import tempfile
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import MagicMock

from gibh_agent.core.asset_locator import (
    bridge_data_asset_to_physical_path,
    bridge_uploaded_files_at_entry,
)


def test_asset_id_db_path_replaces_zombie_frontend_path():
    with tempfile.TemporaryDirectory() as td:
        real = Path(td) / "user1" / "current_batch" / "SRR622461_1.fastq.gz"
        real.parent.mkdir(parents=True)
        real.write_bytes(b"fq")
        zombie = f"/app/uploads/user1/old_ts_batch/SRR622461_1.fastq.gz"

        mock_row = SimpleNamespace(
            id=42,
            owner_id="user1",
            file_name="SRR622461_1.fastq.gz",
            file_path=str(real),
        )
        mock_db = MagicMock()
        mock_db.query.return_value.filter.return_value.first.return_value = mock_row

        resolved, warnings = bridge_data_asset_to_physical_path(
            {
                "path": zombie,
                "file_path": zombie,
                "name": "SRR622461_1.fastq.gz",
                "source": "data_asset",
                "asset_id": 42,
            },
            upload_dir=td,
            owner_id="user1",
            db=mock_db,
        )
        assert resolved == str(real.resolve())
        assert any("asset_id" in w or "映射" in w for w in warnings)


def test_bridge_uploaded_files_list_at_entry():
    with tempfile.TemporaryDirectory() as td:
        real = Path(td) / "u" / "data.h5ad"
        real.parent.mkdir(parents=True)
        real.write_bytes(b"h5")
        mock_row = SimpleNamespace(
            id=7,
            owner_id="u",
            file_name="data.h5ad",
            file_path=str(real),
        )
        mock_db = MagicMock()
        mock_db.query.return_value.filter.return_value.first.return_value = mock_row

        out, warns = bridge_uploaded_files_at_entry(
            [
                {
                    "path": "/app/uploads/u/stale/data.h5ad",
                    "name": "data.h5ad",
                    "source": "data_asset",
                    "asset_id": 7,
                }
            ],
            upload_dir=td,
            owner_id="u",
            db=mock_db,
            strict=False,
        )
        assert len(out) == 1
        assert out[0]["path"] == str(real.resolve())
        assert out[0]["asset_id"] == 7
