# -*- coding: utf-8 -*-
"""Phase 3 最终闭环：入库物理路由单元测试。"""
from __future__ import annotations

import json
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

from gibh_agent.core.ingestion_router import (
    IngestionRouter,
    IngestionRouterError,
    deliver_archive,
)


def test_local_volume_deliver_success(tmp_path):
    mount = tmp_path / "db_mount"
    mount.mkdir()
    archive = tmp_path / "bundle.tar.gz"
    archive.write_bytes(b"fake-archive")

    settings = {
        "mount_type": "local_volume",
        "local_volume": {"mount_path": str(mount)},
    }
    result = deliver_archive(str(archive), settings, session_id="sess-1", owner_id="owner-x")

    assert result["strategy"] == "local_volume"
    dest = Path(result["destination"])
    assert dest.is_file()
    assert dest.read_bytes() == b"fake-archive"
    assert "gibh_ingest/sess-1" in result["destination_dir"]


def test_local_volume_permission_error(tmp_path):
    mount = tmp_path / "readonly"
    mount.mkdir()
    archive = tmp_path / "bundle.tar.gz"
    archive.write_bytes(b"x")

    settings = {
        "mount_type": "local_volume",
        "local_volume": {"mount_path": str(mount)},
    }

    with patch.object(Path, "mkdir", side_effect=PermissionError("denied")):
        with pytest.raises(IngestionRouterError) as exc:
            deliver_archive(str(archive), settings, session_id="s1")
    assert exc.value.mount_type == "local_volume"
    assert "权限" in exc.value.message


def test_hpc_sftp_upload_ioerror(tmp_path):
    archive = tmp_path / "bundle.tar.gz"
    archive.write_bytes(b"x")

    settings = {
        "mount_type": "hpc_slurm",
        "hpc_slurm": {
            "host": "127.0.0.1",
            "port": 22,
            "username": "ubuntu",
            "base_path": "/data/inbox",
        },
    }

    mock_transport = MagicMock()
    mock_sftp = MagicMock()
    mock_sftp.put.side_effect = IOError("Permission denied")

    with patch("paramiko.Transport", return_value=mock_transport):
        with patch("paramiko.SFTPClient.from_transport", return_value=mock_sftp):
            with pytest.raises(IngestionRouterError) as exc:
                deliver_archive(str(archive), settings, session_id="s2")

    assert exc.value.mount_type == "hpc_slurm"
    assert "SFTP 上传失败" in exc.value.message
    mock_sftp.close.assert_called_once()
    mock_transport.close.assert_called_once()


def test_hpc_ssh_auth_failure(tmp_path):
    import paramiko

    archive = tmp_path / "bundle.tar.gz"
    archive.write_bytes(b"x")
    settings = {
        "mount_type": "hpc_slurm",
        "hpc_slurm": {"host": "127.0.0.1", "port": 22, "username": "u", "base_path": ""},
    }

    mock_transport = MagicMock()
    mock_transport.connect.side_effect = paramiko.AuthenticationException("auth fail")

    with patch("paramiko.Transport", return_value=mock_transport):
        with pytest.raises(IngestionRouterError) as exc:
            deliver_archive(str(archive), settings)

    assert "SSH 认证失败" in exc.value.message


def test_api_url_upload_success(tmp_path):
    archive = tmp_path / "bundle.tar.gz"
    archive.write_bytes(b"payload")
    endpoint = "https://example.com/ingest"

    mock_resp = MagicMock()
    mock_resp.status_code = 200
    mock_resp.content = b'{"ok": true}'
    mock_resp.json.return_value = {"ok": True}
    mock_resp.text = '{"ok": true}'

    settings = {
        "mount_type": "api_url",
        "api_url": {"endpoint": endpoint, "token": "tok"},
    }

    with patch("requests.post", return_value=mock_resp) as mock_post:
        result = deliver_archive(str(archive), settings, session_id="s3")

    assert result["http_status"] == 200
    assert result["strategy"] == "api_url"
    mock_post.assert_called_once()
    call_kwargs = mock_post.call_args.kwargs
    assert call_kwargs["headers"]["Authorization"] == "Bearer tok"


def test_missing_archive_raises():
    with pytest.raises(IngestionRouterError, match="归档文件不存在"):
        IngestionRouter().deliver("/nonexistent/archive.tar.gz", {"mount_type": "local_volume"})
