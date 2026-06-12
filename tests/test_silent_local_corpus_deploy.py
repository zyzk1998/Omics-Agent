# -*- coding: utf-8 -*-
"""静默本地落盘与入库检测。"""
from __future__ import annotations

import json
from pathlib import Path

import pytest

from gibh_agent.core.silent_local_corpus_deploy import (
    check_local_session_result_exists,
    silent_deploy_corpus_bundle_to_local_mount,
)


def test_silent_deploy_copies_corpus_to_session_result(tmp_path, monkeypatch):
    mount = tmp_path / "project"
    mount.mkdir()
    corpus_src = tmp_path / "corpus_hitl" / "sess1" / "sft_corpus_silver_test.json"
    corpus_src.parent.mkdir(parents=True)
    corpus_src.write_text(json.dumps([{"instruction": "hi"}]), encoding="utf-8")

    bundle = {
        "status": "success",
        "corpus_hitl_path": str(corpus_src),
        "archive_dir": "",
    }
    delivery = silent_deploy_corpus_bundle_to_local_mount(
        session_id="sess1",
        session_title="My Session",
        bundle=bundle,
        workspace_context={"workspace_path": str(mount)},
        local_workspace_mounted=True,
    )
    assert delivery is not None
    result_dir = Path(delivery["destination_dir"])
    assert (result_dir / "sft_corpus.json").is_file()
    check = check_local_session_result_exists(
        mount_path=str(mount),
        session_id="sess1",
        session_title="My Session",
    )
    assert check["has_local_result"] is True
    assert check["sft_corpus_paths"]


def test_silent_deploy_client_relay_for_windows_mount(tmp_path):
    corpus_src = tmp_path / "corpus_hitl" / "sess-win" / "sft_corpus_silver_test.json"
    corpus_src.parent.mkdir(parents=True)
    corpus_src.write_text(json.dumps([{"instruction": "hi"}]), encoding="utf-8")
    bundle = {
        "status": "success",
        "corpus_hitl_path": str(corpus_src),
        "archive_dir": "",
    }
    delivery = silent_deploy_corpus_bundle_to_local_mount(
        session_id="sess-win",
        session_title="Win Session",
        bundle=bundle,
        workspace_context={"host_workspace_path": r"C:\Users\dev\rawdata"},
        local_workspace_mounted=True,
    )
    assert delivery is not None
    assert delivery.get("strategy") == "client_sidecar_relay"
    assert delivery.get("client_deploy_package")
    pkg = delivery["client_deploy_package"]
    assert pkg["host_mount_path"] == r"C:\Users\dev\rawdata"
    assert pkg["session_id"] == "sess-win"

    out = silent_deploy_corpus_bundle_to_local_mount(
        session_id="x",
        session_title="t",
        bundle={"corpus_hitl_path": "/nope.json"},
        workspace_context={},
        local_workspace_mounted=False,
    )
    assert out is None
