# -*- coding: utf-8 -*-
"""科学语料数据加工技能 · 注册与 HITL 触发冒烟。"""
from __future__ import annotations

import asyncio
import json
from pathlib import Path
from unittest.mock import patch

import pytest

from gibh_agent.agents.corpus_processing_agent import CorpusProcessingAgent
from gibh_agent.core.prompt_manager import create_default_prompt_manager
from gibh_agent.core.tool_registry import registry
from gibh_agent.tools.hitl_tools import (
    HITL_SCENARIO_WHITELIST,
    SCENARIO_TO_LS_XML,
    Trigger_Expert_Annotation,
    _normalize_tasks_for_scenario,
    is_hitl_scenario_allowed,
)


def test_skill_registered_in_tool_registry():
    from gibh_agent.skills import bootstrap_skills

    bootstrap_skills()
    tool = registry.get_tool("skill_corpus_data_processing")
    assert tool is not None


def test_generic_corpus_processing_in_hitl_whitelist():
    assert "generic_corpus_processing" in HITL_SCENARIO_WHITELIST
    assert set(SCENARIO_TO_LS_XML.keys()) == set(HITL_SCENARIO_WHITELIST)
    xml = SCENARIO_TO_LS_XML["generic_corpus_processing"]
    assert "<Image" in xml and "<RectangleLabels" in xml and "<Choices" in xml


def test_normalize_corpus_task_pads_text_for_image_only():
    tasks = _normalize_tasks_for_scenario(
        "generic_corpus_processing",
        [{"data": {"image": "http://nginx/uploads/test_corpus_umap.png"}}],
    )
    assert len(tasks) == 1
    data = tasks[0]["data"]
    assert "text" in data
    assert str(data["text"]).strip()
    assert data["image"].startswith("http://")


def test_corpus_agent_resolve_upload_image():
    agent = CorpusProcessingAgent(llm_client=None, prompt_manager=create_default_prompt_manager())
    media = agent.resolve_upload_media(["/uploads/demo/cell.png"], "")
    assert media["image_path"] == "/uploads/demo/cell.png"
    assert media["image_paths"] == ["/uploads/demo/cell.png"]


def test_corpus_agent_resolve_multiple_upload_images():
    agent = CorpusProcessingAgent(llm_client=None, prompt_manager=create_default_prompt_manager())
    media = agent.resolve_upload_media(
        ["/uploads/a.png", "/uploads/b.jpg", "/uploads/tasks.json"],
        "",
    )
    assert media["image_paths"] == ["/uploads/a.png", "/uploads/b.jpg"]
    assert media["file_path"] == "/uploads/tasks.json"


def test_corpus_agent_prefers_upload_over_demo_template():
    agent = CorpusProcessingAgent(llm_client=None, prompt_manager=create_default_prompt_manager())
    query = (
        '```json\n'
        '{"image_path": "/uploads/test_corpus_umap.png", "project_title": "演示"}\n'
        '```'
    )
    media = agent.resolve_upload_media(
        ["/uploads/user@example.com/20260608_171418/real_upload.png"],
        query,
    )
    assert media["image_path"] == "/uploads/user@example.com/20260608_171418/real_upload.png"


def test_corpus_agent_trigger_hitl_mock(tmp_path, monkeypatch):
    monkeypatch.setenv("LABEL_STUDIO_API_KEY", "test-token")
    monkeypatch.setenv("UPLOAD_DIR", str(tmp_path))
    (tmp_path / "x.png").write_bytes(b"\x89PNG\r\n\x1a\n")
    mock_project = {"id": 101, "title": "Corpus"}
    with patch(
        "gibh_agent.agents.corpus_processing_agent.Trigger_Expert_Annotation",
        wraps=Trigger_Expert_Annotation,
    ) as wrapped:
        with patch(
            "gibh_agent.tools.hitl_tools.LabelStudioClient.create_project",
            return_value=mock_project,
        ):
            with patch(
                "gibh_agent.tools.hitl_tools.LabelStudioClient.import_task",
                return_value={"task_count": 1},
            ):
                agent = CorpusProcessingAgent(
                    llm_client=None, prompt_manager=create_default_prompt_manager()
                )
                result = agent.trigger_corpus_hitl(image_path="/uploads/x.png")
    assert result["status"] == "hitl_required"
    assert result["scenario_type"] == "generic_corpus_processing"
    wrapped.assert_called_once()
    call_kw = wrapped.call_args.kwargs
    assert call_kw["scenario_type"] == "generic_corpus_processing"


def test_fallback_ai_only_vlm_global_description():
    agent = CorpusProcessingAgent(llm_client=None, prompt_manager=create_default_prompt_manager())
    structured = [
        {
            "task_id": "corpus_task_01",
            "image_rel_path": "images/corpus_task_01.png",
            "regions": [],
            "text_fields": {
                "analysis_draft": "Cluster 0 高表达 S100A8；Cluster 1 为 CXCL8 炎症相关群。",
            },
        }
    ]
    out = agent._fallback_ai_only_vlm_from_structured(structured)
    assert len(out) == 1
    assert out[0]["conversations"][0]["value"].startswith("<image>")
    gpt = out[0]["conversations"][1]["value"]
    assert "全局描述" in gpt
    assert "S100A8" in gpt
    assert "bbox" not in gpt.lower()


def test_build_ai_only_structured_tasks_from_draft_and_image():
    agent = CorpusProcessingAgent(llm_client=None, prompt_manager=create_default_prompt_manager())
    structured = agent._build_ai_only_structured_tasks(
        "## 专家分析报告（初稿）\n\nUMAP 显示 3 个主要细胞群。",
        ["data:image/png;base64,xyz"],
        [],
    )
    assert len(structured) == 1
    assert structured[0]["regions"] == []
    assert "UMAP" in structured[0]["text_fields"]["analysis_draft"]


def test_generate_corpus_archive_vlm_slim_image_path(tmp_path, monkeypatch):
    monkeypatch.setenv("RESULTS_DIR", str(tmp_path))
    agent = CorpusProcessingAgent(llm_client=None, prompt_manager=create_default_prompt_manager())
    archive_dir = tmp_path / "corpus_archive" / "test-sess"

    async def _run():
        return await agent.generate_corpus_archive_bundle(
            summary_context={
                "hitl_skipped": True,
                "expert_report_draft": "AI 初稿：T 细胞与 B 细胞群清晰可分。",
                "hitl_image_paths": ["data:image/png;base64,xyz"],
            },
            steps_results=[],
            steps_details=[],
            archive_dir=archive_dir,
            expert_report_md="AI 初稿",
        )

    bundle = asyncio.run(_run())
    assert bundle["status"] == "success"
    assert bundle["count"] == 1
    assert bundle["ai_only"] is True
    payload = json.loads((archive_dir / "dataset.json").read_text(encoding="utf-8"))
    assert payload[0]["id"] == "corpus_task_01"
    assert payload[0]["image"].startswith("images/")
    assert payload[0]["conversations"][0]["from"] == "human"


def test_fallback_vlm_corpus_bbox_norm_1000():
    agent = CorpusProcessingAgent(llm_client=None, prompt_manager=create_default_prompt_manager())
    structured = [
        {
            "task_id": "corpus_task_01",
            "image_rel_path": "images/corpus_task_01.png",
            "regions": [
                {
                    "label": "T-Cell",
                    "bbox_norm_1000": [100, 200, 300, 400],
                    "xywh_percent": [10.0, 20.0, 10.0, 20.0],
                }
            ],
            "text_fields": {},
        }
    ]
    out = agent._fallback_vlm_corpus_from_structured(structured)
    assert len(out) == 1
    assert out[0]["id"] == "corpus_task_01"
    assert out[0]["image"] == "images/corpus_task_01.png"
    gpt = out[0]["conversations"][1]["value"]
    assert "[100, 200, 300, 400]" in gpt
    assert "T-Cell" in gpt
    assert out[0]["conversations"][0]["value"].startswith("<image>")


def test_extract_structured_ls_tasks_from_export():
    agent = CorpusProcessingAgent(llm_client=None, prompt_manager=create_default_prompt_manager())
    export = [
        {
            "id": 8,
            "data": {"image": "data:image/png;base64,xyz", "text": "ctx"},
            "annotations": [
                {
                    "result": [
                        {
                            "type": "rectanglelabels",
                            "value": {
                                "x": 10,
                                "y": 20,
                                "width": 30,
                                "height": 40,
                                "rectanglelabels": ["B-Cell"],
                            },
                        }
                    ]
                }
            ],
        }
    ]
    structured = agent._extract_structured_ls_tasks(export, "")
    assert len(structured) == 1
    assert structured[0]["task_id"] == "8"
    assert structured[0]["regions"][0]["bbox_norm_1000"] == [200, 100, 600, 400]
    assert structured[0]["regions"][0]["label"] == "B-Cell"


def test_parse_sft_json_array():
    agent = CorpusProcessingAgent(llm_client=None, prompt_manager=create_default_prompt_manager())
    raw = '[{"instruction":"a","input":"b","output":"c"}]'
    parsed = agent._parse_sft_json(raw)
    assert len(parsed) == 1
    assert parsed[0]["instruction"] == "a"


def test_corpus_skill_execute_requires_upload():
    from gibh_agent.skills import bootstrap_skills
    from gibh_agent.skills.skill_corpus_data_processing import CorpusDataProcessingSkill

    bootstrap_skills()
    out = CorpusDataProcessingSkill().execute()
    assert out["status"] == "error"
