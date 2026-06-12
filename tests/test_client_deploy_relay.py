# -*- coding: utf-8 -*-
"""浏览器中继落盘 staging / inline 分流测试。"""
from __future__ import annotations

from pathlib import Path

import pytest

from gibh_agent.core.client_deploy_relay import (
    INLINE_MAX_BYTES,
    build_client_relay_deploy_package,
    collect_ingestion_deploy_entries,
    register_staging_record,
    resolve_staged_file,
)


def test_inline_small_archive_only(tmp_path: Path):
    archive = tmp_path / "bundle_small.tar.gz"
    archive.write_bytes(b"x" * 512)
    entries = collect_ingestion_deploy_entries(archive_path=str(archive))
    pkg = build_client_relay_deploy_package(
        owner_id="owner1",
        session_id="sess1",
        session_title="T",
        host_mount_path=r"C:\proj\rawdata",
        entries=entries,
        folder_timestamp="202606111628",
    )
    assert pkg["transport"] == "inline"
    assert len(pkg["files"]) == 1
    assert not pkg["download_items"]
    assert pkg["files"][0]["dest_name"] == archive.name
    assert pkg["session_folder"] == "T-202606111628"
    assert pkg["folder_timestamp"] == "202606111628"


def test_url_mode_for_large_archive(tmp_path: Path, monkeypatch):
    monkeypatch.setenv("CLIENT_DEPLOY_STAGING_DIR", str(tmp_path / "staging"))
    archive = tmp_path / "bundle_large.tar.gz"
    archive.write_bytes(b"L" * (INLINE_MAX_BYTES + 1))
    entries = collect_ingestion_deploy_entries(archive_path=str(archive))
    pkg = build_client_relay_deploy_package(
        owner_id="owner1",
        session_id="sess1",
        session_title="T",
        host_mount_path=r"C:\proj\rawdata",
        entries=entries,
        folder_timestamp="202606111628",
    )
    assert pkg["transport"] == "url"
    assert not pkg["files"]
    assert len(pkg["download_items"]) == 1
    assert pkg["download_items"][0]["download_url"].startswith("/api/ingestion/deploy-staging/")
    assert pkg["staging_token"]


def test_resolve_staged_file_owner_bound(tmp_path: Path, monkeypatch):
    monkeypatch.setenv("CLIENT_DEPLOY_STAGING_DIR", str(tmp_path / "staging"))
    src = tmp_path / "big.bin"
    src.write_bytes(b"Z" * 10)
    token = "tok123"
    root = tmp_path / "staging" / "owner1" / "sess1" / token
    root.mkdir(parents=True)
    staged = root / "bundle_large.tar.gz"
    staged.write_bytes(src.read_bytes())
    register_staging_record(
        token=token,
        owner_id="owner1",
        session_id="sess1",
        root=root,
        file_map={"bundle_large.tar.gz": staged},
    )
    resolved = resolve_staged_file(token, "bundle_large.tar.gz", owner_id="owner1")
    assert resolved.read_bytes() == b"Z" * 10
    with pytest.raises(PermissionError):
        resolve_staged_file(token, "bundle_large.tar.gz", owner_id="other")
