# -*- coding: utf-8 -*-
"""Phase 1 HITL 与用户数据库挂载配置单元测试。"""
from __future__ import annotations

import json
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

from gibh_agent.core.user_database_settings_store import (
    get_database_mount_config,
    save_database_mount_config,
)
from gibh_agent.tools.hitl_tools import (
    HITL_SCENARIO_WHITELIST,
    Trigger_Expert_Annotation,
    is_hitl_scenario_allowed,
    resolve_ls_accessible_image_url,
)
from gibh_agent.utils.ls_client import (
    LabelStudioClient,
    LabelStudioClientError,
)


def test_hitl_whitelist_rejects_unknown_scenario():
    result = Trigger_Expert_Annotation(
        scenario_type="random_qc_report",
        project_title="test",
        image_path="https://example.com/a.jpg",
    )
    assert result["status"] == "error"
    assert "白名单" in result["message"]
    assert "allowed_scenarios" in result


def test_hitl_whitelist_allows_and_returns_hitl_required(tmp_path: Path, monkeypatch):
    monkeypatch.setenv("LABEL_STUDIO_API_KEY", "test-token")
    mock_project = {"id": 42, "title": "Expert Review"}
    with patch.object(LabelStudioClient, "create_project", return_value=mock_project):
        with patch.object(LabelStudioClient, "import_task", return_value={"task_count": 1}):
            with patch.object(LabelStudioClient, "project_url", return_value="http://127.0.0.1:8082/projects/42/data"):
                result = Trigger_Expert_Annotation(
                    scenario_type="scrna_cell_type_annotation",
                    project_title="scRNA Expert Review",
                    image_path="https://example.com/cell.png",
                )
    assert result["status"] == "hitl_required"
    assert result["ls_project_id"] == 42
    assert result["ls_project_url"] in (
        "http://127.0.0.1:8082/projects/42/data",
        "/label-studio/projects/42/data",
    )
    assert result["scenario_type"] == "scrna_cell_type_annotation"


def test_is_hitl_scenario_allowed():
    assert is_hitl_scenario_allowed("spatial_microenvironment") is True
    assert is_hitl_scenario_allowed("not_allowed") is False
    assert len(HITL_SCENARIO_WHITELIST) == 6
    assert is_hitl_scenario_allowed("generic_corpus_processing") is True


def test_ls_client_create_project_validation():
    client = LabelStudioClient(base_url="http://127.0.0.1:8082", api_key="test")
    with pytest.raises(LabelStudioClientError, match="title"):
        client.create_project("")


def test_ls_client_http_error():
    client = LabelStudioClient(base_url="http://127.0.0.1:8082", api_key="test")
    mock_resp = MagicMock()
    mock_resp.status_code = 401
    mock_resp.content = b'{"detail":"auth"}'
    mock_resp.json.return_value = {"detail": "auth"}
    with patch.object(client._session, "request", return_value=mock_resp):
        with pytest.raises(LabelStudioClientError, match="HTTP 401"):
            client.create_project("Demo")


def test_database_settings_local_volume_roundtrip(tmp_path, monkeypatch):
    monkeypatch.setenv("GIBH_USER_DATABASE_SETTINGS_DIR", str(tmp_path))
    from gibh_agent.core import user_database_settings_store as store

    store._CACHE.clear()
    mount = tmp_path / "data_vol"
    mount.mkdir()
    cfg = {
        "mount_type": "local_volume",
        "is_auto_ingestion_enabled": True,
        "local_volume": {"mount_path": str(mount)},
    }
    saved = save_database_mount_config("owner_test", cfg)
    assert saved["is_auto_ingestion_enabled"] is True
    loaded = get_database_mount_config("owner_test")
    assert loaded["local_volume"]["mount_path"] == str(mount)


def test_database_settings_test_connection_local(tmp_path):
    from gibh_agent.core.user_database_settings_store import probe_database_mount_connection

    mount = tmp_path / "vol"
    mount.mkdir()
    ok = probe_database_mount_connection(
        {
            "mount_type": "local_volume",
            "local_volume": {"mount_path": str(mount)},
        }
    )
    assert ok["status"] == "success"
    bad = probe_database_mount_connection(
        {
            "mount_type": "local_volume",
            "local_volume": {"mount_path": str(tmp_path / "missing")},
        }
    )
    assert bad["status"] == "error"


def test_trigger_expert_annotation_tasks_from_file(tmp_path: Path, monkeypatch):
    monkeypatch.setenv("LABEL_STUDIO_API_KEY", "test-token")
    tasks_file = tmp_path / "tasks.json"
    tasks_file.write_text(json.dumps([{"image": "https://example.com/x.jpg"}]), encoding="utf-8")
    mock_project = {"id": 7, "title": "T"}
    with patch.object(LabelStudioClient, "create_project", return_value=mock_project):
        with patch.object(LabelStudioClient, "import_task", return_value={}) as imp:
            with patch.object(LabelStudioClient, "project_url", return_value="http://127.0.0.1:8082/projects/7/data"):
                result = Trigger_Expert_Annotation(
                    scenario_type="radiomics_roi_validation",
                    project_title="ROI Review",
                    file_path=str(tasks_file),
                )
    assert result["status"] == "hitl_required"
    imp.assert_called_once()


def test_resolve_ls_accessible_image_url_uses_fetch_base(monkeypatch):
    monkeypatch.setenv("LS_IMAGE_FETCH_BASE_URL", "http://nginx")
    url = resolve_ls_accessible_image_url("/results/run1/umap_annotated.png")
    assert url == "http://nginx/results/run1/umap_annotated.png"
    assert url.startswith("http://")
    assert "192.168." not in url


def test_normalize_frontend_media_path_strips_toxic_https():
    from gibh_agent.tools.hitl_tools import normalize_frontend_media_path

    rel = normalize_frontend_media_path(
        "https://192.168.32.31:8028/results/run_20240602/umap_annotated.png"
    )
    assert rel == "/results/run_20240602/umap_annotated.png"


def test_resolve_ls_import_rejects_https_base(monkeypatch):
    monkeypatch.setenv("LS_IMAGE_FETCH_BASE_URL", "https://nginx")
    from gibh_agent.tools.hitl_tools import resolve_ls_import_image_url

    url = resolve_ls_import_image_url("/results/a.png")
    assert url.startswith("http://nginx/")
    assert not url.startswith("https://")


def test_resolve_ls_import_rewrites_toxic_absolute_url(monkeypatch):
    monkeypatch.setenv("LS_IMAGE_FETCH_BASE_URL", "http://nginx")
    from gibh_agent.tools.hitl_tools import resolve_ls_import_image_url

    toxic = "https://192.168.32.31:8028/results/run_20240602/umap_annotated.png"
    url = resolve_ls_import_image_url(toxic)
    assert url == "http://nginx/results/run_20240602/umap_annotated.png"


def test_auto_bootstrap_session_login(monkeypatch):
    monkeypatch.delenv("LABEL_STUDIO_API_KEY", raising=False)
    from gibh_agent.utils import ls_client as mod

    mod._BOOTSTRAP_SESSION_CACHE.clear()
    mock_session = MagicMock()
    mock_session.cookies.get.return_value = "sid"
    with patch.object(mod, "_wait_ls_health", return_value=True):
        with patch.object(mod, "_session_login", return_value=mock_session):
            session = mod._auto_bootstrap_session(
                "http://127.0.0.1:8082",
                username="system_admin@omics.local",
                password="test-pass",
                timeout_sec=5,
                force_refresh=True,
            )
    assert session is mock_session
    assert mod._BOOTSTRAP_SESSION_CACHE["http://127.0.0.1:8082"] is mock_session


def test_label_studio_client_lazy_bootstrap_on_request(monkeypatch):
    monkeypatch.delenv("LABEL_STUDIO_API_KEY", raising=False)
    mock_session = MagicMock()
    with patch(
        "gibh_agent.utils.ls_client._auto_bootstrap_session",
        return_value=mock_session,
    ) as boot:
        client = LabelStudioClient(
            base_url="http://127.0.0.1:8082",
            api_key="",
            auto_bootstrap=True,
        )
        assert client.api_key == ""
        client._ensure_authenticated()
    assert client._auth_mode == "session"
    assert client._session is mock_session
    boot.assert_called_once()


def test_trigger_expert_annotation_uses_auto_bootstrap(monkeypatch):
    monkeypatch.delenv("LABEL_STUDIO_API_KEY", raising=False)
    mock_project = {"id": 99, "title": "T"}
    with patch.object(LabelStudioClient, "create_project", return_value=mock_project):
        with patch.object(LabelStudioClient, "import_task", return_value={"task_count": 1}):
            with patch.object(
                LabelStudioClient,
                "__init__",
                lambda self, *a, **k: setattr(self, "api_key", "bootstrapped"),
            ):
                result = Trigger_Expert_Annotation(
                    scenario_type="scrna_cell_type_annotation",
                    project_title="Auto Bootstrap Test",
                    image_path="https://example.com/cell.png",
                )
    assert result["status"] == "hitl_required"
    assert result["ls_project_id"] == 99
