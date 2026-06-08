# -*- coding: utf-8 -*-
"""Phase 3 缺口补救：打包器 / 注册表 / 唤醒逻辑单元测试。"""
from __future__ import annotations

import json
from pathlib import Path

from gibh_agent.core.data_packager import build_artifacts_archive
from gibh_agent.core.hitl_session_registry import (
    get_hitl_session,
    register_hitl_session,
    resolve_session_by_project,
)
from gibh_agent.core.orchestrator import AgentOrchestrator


def test_hitl_session_registry_roundtrip(tmp_path, monkeypatch):
    monkeypatch.setenv("GIBH_HITL_REGISTRY_DIR", str(tmp_path))
    from gibh_agent.core import hitl_session_registry as reg

    reg._REGISTRY_DIR = tmp_path  # noqa: SLF001
    register_hitl_session("sess-abc", {"project_id": 77, "ls_url": "http://127.0.0.1:8082/projects/77/data"})
    assert resolve_session_by_project(77) == "sess-abc"
    loaded = get_hitl_session("sess-abc")
    assert loaded and loaded.get("project_id") == 77


def test_data_packager_builds_archive(tmp_path, monkeypatch):
    monkeypatch.setenv("GIBH_ARTIFACTS_ARCHIVE_DIR", str(tmp_path / "archives"))
    report = "## 专家分析报告（最终版）\n\n测试内容"
    result = build_artifacts_archive(
        session_id="s1",
        owner_id="owner1",
        expert_report_markdown=report,
        steps_details=[{"step_id": "step1", "status": "success"}],
        hitl_annotations=[{"task": 1, "annotations": []}],
        hitl_meta={"scenario_type": "scrna_cell_type_annotation"},
    )
    assert result["status"] == "success"
    archive = Path(result["archive_path"])
    assert archive.is_file()
    manifest = json.loads((Path(result["bundle_dir"]) / "manifest.json").read_text(encoding="utf-8"))
    assert "expert_report_final.md" in manifest["files"]
    assert "hitl/annotations.json" in manifest["files"]


def test_extract_hitl_payload_unchanged():
    hitl = {"ls_project_url": "http://x", "ls_project_id": 3, "status": "hitl_required"}
    payload = AgentOrchestrator._extract_hitl_payload({"status": "hitl_required", "hitl": hitl})
    assert payload and payload.get("ls_project_id") == 3
