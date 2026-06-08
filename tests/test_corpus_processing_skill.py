# -*- coding: utf-8 -*-
"""科学语料数据加工技能 · 注册与 HITL 触发冒烟。"""
from __future__ import annotations

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


def test_corpus_agent_trigger_hitl_mock(monkeypatch):
    monkeypatch.setenv("LABEL_STUDIO_API_KEY", "test-token")
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
