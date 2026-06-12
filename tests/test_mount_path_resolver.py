# -*- coding: utf-8 -*-
from gibh_agent.core.storage.mount_path_resolver import (
    normalize_workspace_mount_root,
    resolve_container_writable_mount,
    resolve_ingestion_mount_target,
)


def test_windows_mount_does_not_fallback_to_container_default(tmp_path, monkeypatch):
    data_dir = tmp_path / "data"
    data_dir.mkdir()
    monkeypatch.setenv("DEFAULT_INGESTION_DIR", str(data_dir))
    win_path = r"C:\Users\test\Desktop\rawdata\result"
    assert resolve_container_writable_mount(win_path) is None
    mt, path, err = resolve_ingestion_mount_target(win_path)
    assert err is None
    assert mt == "local_sidecar"
    assert path == r"C:\Users\test\Desktop\rawdata"


def test_normalize_workspace_mount_root_strips_result_suffix():
    assert normalize_workspace_mount_root(r"C:\proj\rawdata\result") == r"C:\proj\rawdata"
    assert normalize_workspace_mount_root("/app/data/result") == "/app/data"
