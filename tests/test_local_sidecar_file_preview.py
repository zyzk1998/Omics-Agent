"""本地 Sidecar 文件预览：POST /api/files/download 须正确传递 Windows 风格路径。"""
from __future__ import annotations

import importlib.util
import os
import sys
from pathlib import Path

import pytest
from fastapi.testclient import TestClient


def _load_local_server():
    repo = Path(__file__).resolve().parents[1]
    mod_path = repo / "gibh-desktop-app" / "local_sidecar" / "local_server.py"
    spec = importlib.util.spec_from_file_location("omics_local_sidecar_test", mod_path)
    assert spec and spec.loader
    mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mod
    spec.loader.exec_module(mod)
    return mod


@pytest.fixture()
def sidecar_client(tmp_path: Path):
    mod = _load_local_server()
    mod._workspace_context.clear()
    mod._set_workspace_context(str(tmp_path))
    mod._assert_deploy_client = lambda _request: None
    img = tmp_path / "volcano_plot.png"
    img.write_bytes(b"\x89PNG\r\n\x1a\n" + b"\x00" * 32)
    with TestClient(mod.app) as client:
        yield client, img, tmp_path, mod


@pytest.mark.skipif(os.name != "nt", reason="Windows 反斜杠路径仅在 nt 平台可 resolve")
def test_download_post_windows_backslash_path(sidecar_client):
    client, img, root, _mod = sidecar_client
    win_style = str(root).replace("/", "\\") + "\\volcano_plot.png"
    res = client.post("/api/files/download", json={"path": win_style})
    assert res.status_code == 200, res.text
    assert res.headers.get("content-type", "").startswith("image/png")
    assert res.content[:8] == b"\x89PNG\r\n\x1a\n"


def test_show_in_explorer_accepts_workspace_directory(sidecar_client):
    client, _img, root, mod = sidecar_client
    mod._assert_loopback_client = lambda _request: None
    res = client.post("/api/fs/show_in_explorer", json={"path": str(root)})
    assert res.status_code == 200, res.text
    data = res.json()
    assert data.get("status") == "success"
    assert data.get("action") == "show_in_explorer"


def test_health_lists_files_download_capability(sidecar_client):
    client, _img, _root, _mod = sidecar_client
    res = client.get("/health")
    assert res.status_code == 200
    data = res.json()
    caps = data.get("capabilities") or []
    assert "files_download" in caps
    assert "silent_deploy" in caps
    assert "silent_deploy_from_url" in caps
    assert "write_bytes" in caps


def test_silent_deploy_writes_to_session_result(sidecar_client):
    import base64

    client, _img, root, _mod = sidecar_client
    payload = {
        "session_id": "sess-a",
        "session_title": "My Session",
        "host_mount_path": str(root),
        "files": [
            {
                "dest_name": "bundle_test.tar.gz",
                "content_b64": base64.b64encode(b"fake-tar").decode("ascii"),
            }
        ],
    }
    res = client.post("/api/workspace/silent_deploy", json=payload)
    assert res.status_code == 200, res.text
    data = res.json()
    assert data.get("status") == "success"
    assert data.get("changed_paths")


def test_write_bytes_via_write_file_b64(sidecar_client, tmp_path):
    import base64

    client, _img, root, _mod = sidecar_client
    dest = root / "out.bin"
    res = client.post(
        "/api/tools/write_file",
        json={
            "file_path": str(dest),
            "content_b64": base64.b64encode(b"\x00\x01\x02").decode("ascii"),
        },
    )
    assert res.status_code == 200, res.text
    assert dest.read_bytes() == b"\x00\x01\x02"


def test_download_tool_alias_route(sidecar_client):
    client, img, _root, _mod = sidecar_client
    res = client.post("/api/tools/download_file", json={"path": str(img)})
    assert res.status_code == 200
    assert res.headers.get("content-type", "").startswith("image/png")


def test_download_post_relative_to_workspace(sidecar_client):
    client, img, _root, _mod = sidecar_client
    res = client.post("/api/files/download", json={"path": str(img)})
    assert res.status_code == 200
    assert len(res.content) > 0
