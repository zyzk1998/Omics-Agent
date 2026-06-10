# -*- coding: utf-8 -*-
"""Data Packager 须将 corpus_archive（解耦路径 + 图像目录）纳入一键入库归档。"""
from __future__ import annotations

import json
import tarfile
from pathlib import Path

from gibh_agent.core.data_packager import build_artifacts_archive


def test_build_artifacts_archive_includes_corpus_archive_tree(tmp_path, monkeypatch):
    monkeypatch.setenv("GIBH_ARTIFACTS_ARCHIVE_DIR", str(tmp_path / "archives"))
    results = tmp_path / "results"
    archive_src = results / "corpus_archive" / "sess-01"
    images_dir = archive_src / "images"
    reports_dir = archive_src / "reports"
    images_dir.mkdir(parents=True)
    reports_dir.mkdir(parents=True)
    (images_dir / "umap_01.png").write_bytes(b"\x89PNG\r\n\x1a\n")
    (reports_dir / "expert_report.md").write_text("# 初稿", encoding="utf-8")
    dataset = [
        {
            "id": "corpus_task_01",
            "image": "images/umap_01.png",
            "conversations": [
                {"from": "human", "value": "<image>\n描述"},
                {"from": "gpt", "value": "全局描述"},
            ],
        }
    ]
    (archive_src / "dataset.json").write_text(json.dumps(dataset, ensure_ascii=False), encoding="utf-8")
    monkeypatch.setenv("RESULTS_DIR", str(results))

    pack = build_artifacts_archive(
        session_id="sess-01",
        owner_id="owner-1",
        expert_report_markdown="# 初稿报告",
        steps_details=[{"step_name": "UMAP", "status": "success"}],
        skip_hitl=True,
        corpus_archive_dir=str(archive_src),
        corpus_modality="vlm",
    )

    assert pack["status"] == "success"
    manifest = pack["manifest"]
    files = manifest.get("files") or []
    assert "corpus_archive/dataset.json" in files
    assert "corpus_archive/images/umap_01.png" in files
    assert "corpus_archive/reports/expert_report.md" in files

    archive_path = Path(pack["archive_path"])
    with tarfile.open(archive_path, "r:gz") as tar:
        names = tar.getnames()
    assert any(n.endswith("corpus_archive/dataset.json") for n in names)
    assert any(n.endswith("umap_01.png") for n in names)
