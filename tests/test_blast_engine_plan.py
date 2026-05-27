# -*- coding: utf-8 -*-
"""BLAST 双引擎路由单元测试（不调用真实 blastn）。"""
from __future__ import annotations

from pathlib import Path

import pytest

from gibh_agent.skills.blast_engine import (
    MSG_REMOTE_TIMEOUT,
    plan_blastn_execution,
    user_facing_blastn_error,
)


def test_remote_plan_when_no_local_db(monkeypatch, tmp_path):
    monkeypatch.setenv("LOCAL_BLASTDB_PATH", str(tmp_path))
    plan = plan_blastn_execution(54, "nt", "", use_remote=True)
    assert plan["ok"] is True
    assert plan["engine"] == "remote"
    assert plan["use_remote"] is True
    assert plan["timeout_s"] == 600


def test_local_plan_when_db_files_exist(monkeypatch, tmp_path):
    db = tmp_path / "nt"
    (tmp_path / "nt.nhr").write_bytes(b"x")
    (tmp_path / "nt.nin").write_bytes(b"x")
    (tmp_path / "nt.nsq").write_bytes(b"x")
    monkeypatch.setenv("LOCAL_BLASTDB_PATH", str(tmp_path))
    plan = plan_blastn_execution(54, "nt", "", use_remote=True)
    assert plan["ok"] is True
    assert plan["engine"] == "local"
    assert plan["use_remote"] is False
    assert Path(plan["db_arg"]).name == "nt"


def test_user_facing_timeout_message():
    msg = user_facing_blastn_error(
        code=-1,
        stderr="[timeout] x",
        stdout="",
        engine="remote",
        timeout_s=600,
    )
    assert msg == MSG_REMOTE_TIMEOUT
    assert "Traceback" not in msg
    assert "TimeoutExpired" not in msg
