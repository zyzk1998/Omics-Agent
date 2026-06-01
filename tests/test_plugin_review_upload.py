# -*- coding: utf-8 -*-
"""用户技能审核提交：暂存区落盘 + status=pending。"""
from __future__ import annotations

import json
import tempfile
from pathlib import Path

import pytest

from gibh_agent.core.skills_assets_layout import get_skills_review_staging_dir
from gibh_agent.plugin_system.router import (
    _infer_metadata_from_submission,
    _normalize_driver_type,
    _stage_submission_files,
)


@pytest.fixture()
def upload_tmp(monkeypatch):
    root = Path(tempfile.mkdtemp(prefix="gibh_review_"))
    monkeypatch.setenv("UPLOAD_DIR", str(root))
    return root


def test_normalize_driver_type_aliases():
    assert _normalize_driver_type("prompt") == "prompt"
    assert _normalize_driver_type("单脚本驱动") == "script"
    assert _normalize_driver_type("model") == "model"


def test_stage_submission_writes_prompt_and_manifest(upload_tmp):
    uid, extract_dir, manifest = _stage_submission_files(
        user_id="alice",
        driver_type="prompt",
        prompt_text="## Name\nDemo\n## Description\nHello",
        file_bytes=None,
        original_filename=None,
    )
    dest = Path(extract_dir)
    assert dest.is_dir()
    assert (dest / "prompt.md").is_file()
    assert (dest / "manifest.json").is_file()
    assert manifest["review_status"] == "pending"
    assert "prompt.md" in manifest["saved_files"]
    assert get_skills_review_staging_dir().name == "skills_review_staging"


def test_infer_metadata_from_prompt_text(upload_tmp):
    meta = _infer_metadata_from_submission(
        submission_uid="abc-123",
        driver_type="prompt",
        prompt_text="## Name\nMy Skill\n## Description\nDesc",
        file_bytes=None,
        original_filename=None,
    )
    assert meta["name"] == "my_skill"
    assert meta["display_name"] == "My Skill"
