# -*- coding: utf-8 -*-
"""临时 RESULTS_DIR 扫描必须按 session_id 隔离，禁止穿透其它 run_* 目录。"""
from __future__ import annotations

from pathlib import Path

from gibh_agent.core.storage.dual_path import list_session_runtime_files


def test_list_session_runtime_files_ignores_unrelated_run_dirs(tmp_path, monkeypatch):
    results_root = tmp_path / "results"
    upload_root = tmp_path / "uploads"
    results_root.mkdir()
    upload_root.mkdir()

    sid_a = "session-aaa"
    sid_b = "session-bbb"

    run_a = results_root / "run_20260101_120000"
    run_b = results_root / "run_20260101_130000"
    run_a.mkdir()
    run_b.mkdir()
    (run_a / "only_a.txt").write_text("a", encoding="utf-8")
    (run_b / "only_b.txt").write_text("b", encoding="utf-8")

    monkeypatch.setenv("RESULTS_DIR", str(results_root))
    monkeypatch.setenv("UPLOAD_DIR", str(upload_root))

    inv_a = list_session_runtime_files(sid_a, referenced_paths=[str(run_a / "only_a.txt")])
    paths_a = {row["path"] for row in inv_a["results"]}
    assert str((run_a / "only_a.txt").resolve()) in paths_a
    assert str((run_b / "only_b.txt").resolve()) not in paths_a

    inv_empty = list_session_runtime_files(sid_a, referenced_paths=[])
    assert inv_empty["results"] == []


def test_list_session_runtime_files_scans_session_subdirs(tmp_path, monkeypatch):
    results_root = tmp_path / "results"
    upload_root = tmp_path / "uploads"
    results_root.mkdir()
    upload_root.mkdir()

    sid = "sess-corpus"
    corpus_dir = results_root / "corpus_archive" / sid
    corpus_dir.mkdir(parents=True)
    sample = corpus_dir / "dataset.json"
    sample.write_text("{}", encoding="utf-8")

    monkeypatch.setenv("RESULTS_DIR", str(results_root))
    monkeypatch.setenv("UPLOAD_DIR", str(upload_root))

    inv = list_session_runtime_files(sid)
    paths = {row["path"] for row in inv["results"]}
    assert str(sample.resolve()) in paths
