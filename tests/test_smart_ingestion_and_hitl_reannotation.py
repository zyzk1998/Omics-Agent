# -*- coding: utf-8 -*-
"""Smart Mount Discovery 与 HITL 重新标注放行单元测试。"""
from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace

from gibh_agent.api.routers import hitl_api
from gibh_agent.core.hitl_session_registry import register_hitl_session
from gibh_agent.core.ingestion_deploy import deploy_artifacts_to_mount
from gibh_agent.core.ingestion_mount_discovery import (
    discover_ingestion_mount_paths,
    get_default_ingestion_mount_path,
    probe_mount_path,
    validate_ingestion_mount_path,
)


def test_probe_mount_path_writable_tmp(tmp_path):
    info = probe_mount_path(str(tmp_path))
    assert info is not None
    assert info["writable"] is True
    assert Path(info["path"]).resolve() == tmp_path.resolve()


def test_validate_rejects_missing_path():
    assert validate_ingestion_mount_path("/nonexistent/gibh_mount_test") is None


def test_discover_includes_tmp_when_env_set(tmp_path, monkeypatch):
    monkeypatch.setenv("DEFAULT_INGESTION_DIR", str(tmp_path))
    paths = discover_ingestion_mount_paths()
    assert paths
    assert paths[0]["path"] == str(tmp_path.resolve())
    assert get_default_ingestion_mount_path() == str(tmp_path.resolve())


def test_deploy_artifacts_copies_archive_and_bundle(tmp_path):
    archive = tmp_path / "bundle_test.tar.gz"
    archive.write_bytes(b"fake-tar")
    bundle = tmp_path / "bundle_test"
    corpus = bundle / "corpus_archive"
    corpus.mkdir(parents=True)
    (corpus / "dataset.json").write_text('{"ok": true}', encoding="utf-8")

    mount = tmp_path / "mount_vol"
    mount.mkdir()

    result = deploy_artifacts_to_mount(
        archive_path=str(archive),
        mount_path=str(mount),
        session_id="sess-1",
        owner_id="owner-abc",
        bundle_dir=str(bundle),
    )
    dest_dir = Path(result["destination_dir"])
    assert (dest_dir / archive.name).is_file()
    assert (dest_dir / bundle.name / "corpus_archive" / "dataset.json").is_file()
    assert result["bundle_destination"]


def test_reannotation_resume_allowed_when_completed(monkeypatch, tmp_path):
    monkeypatch.setenv("GIBH_HITL_REGISTRY_DIR", str(tmp_path))
    from gibh_agent.core import hitl_session_registry as reg

    reg._REGISTRY_DIR = tmp_path  # noqa: SLF001

    snap = {
        "hitl_resumed": True,
        "hitl_final": True,
        "execution_snapshot": {"hitl_final": True, "hitl_annotations": [{"id": 1}]},
        "hitl": {"project_id": 42, "ls_url": "/label-studio/projects/42/data"},
    }
    register_hitl_session("sess-x", {"project_id": 42})

    class _FakeMsg:
        content = {"state_snapshot": snap}

    class _FakeQuery:
        def filter(self, *args, **kwargs):
            return self

        def order_by(self, *args, **kwargs):
            return self

        def first(self):
            return _FakeMsg()

    class _FakeDb:
        def query(self, *args, **kwargs):
            return _FakeQuery()

    monkeypatch.setattr(
        "gibh_agent.core.hitl_resume._load_latest_agent_snapshot",
        lambda db, sid, oid: (_FakeMsg(), snap),
    )

    assert hitl_api._snapshot_allows_reannotation_resume(_FakeDb(), "sess-x", "owner") is True
    session = SimpleNamespace(id="sess-x", status="completed")
    assert hitl_api._session_allows_hitl_resume(_FakeDb(), session, "owner") is True

    assert hitl_api._resolve_hitl_project_id("sess-x", snap) == 42
