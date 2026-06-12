# -*- coding: utf-8 -*-
"""Corpus 终点站：快照持久化字段与独立技能 Skip 媒体提取。"""
from __future__ import annotations

import json

from gibh_agent.core.execution_snapshot import update_corpus_asset_in_state
from gibh_agent.core.hitl_resume import (
    _extract_standalone_corpus_media_paths,
    _is_standalone_corpus_skill,
)


def test_update_corpus_asset_in_state_writes_execution_snapshot():
    snap = {"execution_snapshot": {"steps_details": []}}
    records = [{"instruction": "test", "input": "a", "output": "b"}]
    payload = {
        "standard": "silver",
        "modality": "nlp",
        "count": 1,
        "records": records,
        "corpus_json": json.dumps(records, ensure_ascii=False, indent=2),
        "archive_dir": "/tmp/corpus_archive/s1",
        "corpus_hitl_path": "/tmp/corpus_hitl/s1/sft.json",
    }
    update_corpus_asset_in_state(snap, payload)
    ex = snap["execution_snapshot"]
    assert ex["corpus_standard"] == "silver"
    assert ex["corpus_modality"] == "nlp"
    assert ex["corpus_count"] == 1
    assert ex["corpus_sft_records"][0]["instruction"] == "test"
    assert "instruction" in ex["corpus_sft_json"]
    assert snap["corpus_sft_json"] == ex["corpus_sft_json"]


def test_extract_standalone_corpus_media_from_hitl_registry():
    snap = {"workflow": {"skill_id": "skill_corpus_data_processing"}}
    reg = {
        "agent_key": "corpus_processing_agent",
        "image_paths": ["/uploads/demo/cluster.png"],
        "uploaded_file_paths": ["/uploads/demo/cluster.png"],
    }
    assert _is_standalone_corpus_skill(snap, reg, {})
    paths = _extract_standalone_corpus_media_paths(snap, reg, {})
    assert "/uploads/demo/cluster.png" in paths
