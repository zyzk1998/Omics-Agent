# -*- coding: utf-8 -*-
"""Prompt-Engineered 软技能：注册与结构化输出单元测试（不调用真实 LLM）。"""
from __future__ import annotations

import json
from unittest.mock import patch

import pytest

from gibh_agent.core.skill_registry import discover_and_register_skills
from gibh_agent.core.tool_registry import registry
from gibh_agent.skills._prompt_skill_llm import normalize_mermaid_code
from gibh_agent.skills._prompt_skill_schemas import PPTOutline, Slide
from gibh_agent.skills.skill_mindmap_gen import MindmapGenSkill
from gibh_agent.skills.skill_ppt_outline import PptOutlineSkill


@pytest.fixture(scope="module", autouse=True)
def _register_skills():
    discover_and_register_skills(force_reload=True)


def test_prompt_skills_not_delegated_to_launch_worker():
    """Prompt 软技能须在 api-server 本地执行（LLM + prompt_specs），勿误委托 launch-skills。"""
    from gibh_agent.skills.launch_skill_demos import LAUNCH_ISOLATED_TOOL_IDS
    from gibh_agent.skills.launch_skill_isolated import should_delegate_to_launch_worker

    for tid in (
        "ppt_outline",
        "mindmap_gen",
        "weekly_report_writer",
        "tailored_resume",
        "academic_poster_generator",
        "academic_abstract_refiner",
        "email_manager",
        "pdf_extractor",
        "deep_research",
        "blueprint_drafter",
    ):
        assert tid not in LAUNCH_ISOLATED_TOOL_IDS
        assert should_delegate_to_launch_worker(tid) is False


def test_tools_registered():
    assert registry.get_tool("ppt_outline") is not None
    assert registry.get_tool("mindmap_gen") is not None
    assert registry.get_tool("weekly_report_writer") is not None
    assert registry.get_tool("email_manager") is not None
    assert registry.get_tool("deep_research") is not None
    assert registry.get_tool("blueprint_drafter") is not None


def test_prompt_specs_load():
    from gibh_agent.skills._prompt_skill_spec import load_prompt_spec

    assert len(load_prompt_spec("weekly_report_writer")) > 200
    assert len(load_prompt_spec("email_manager")) > 100


def test_ppt_outline_pydantic_schema():
    model = PPTOutline.model_validate(
        {
            "theme": "Demo",
            "slides": [
                {"title": "背景", "content_bullets": ["a", "b"]},
                {"title": "结论", "bullets": ["c"]},
            ],
        }
    )
    out = model.to_frontend_dict()
    assert out["total_pages"] == 2
    assert out["slides"][1]["bullets"] == ["c"]


def test_normalize_mermaid_strips_fence():
    code = normalize_mermaid_code("```mermaid\ngraph TD; A-->B\n```")
    assert code.startswith("graph TD")


@patch("gibh_agent.skills.skill_ppt_outline.llm_chat_structured")
def test_ppt_outline_execute_mock(mock_structured):
    mock_structured.return_value = PPTOutline(
        theme="T",
        slides=[Slide(title="P1", content_bullets=["x"])],
    )
    res = PptOutlineSkill().execute(user_request="test topic")
    assert res["status"] == "success"
    assert res["data"]["ppt_outline"]["slides"][0]["title"] == "P1"
    assert res["ppt_outline"]["slides"][0]["title"] == "P1"


@patch("gibh_agent.skills._prompt_skill_terminal.llm_chat_json")
def test_weekly_report_terminal_deliver_mock(mock_json):
    mock_json.return_value = {
        "action": "deliver",
        "markdown": "# 周报\n\n## 本周完成\n\n项目 A 里程碑\n\n## 下周计划\n\n推进 B",
        "missing_params": [],
        "summary": "周报终稿",
    }
    from gibh_agent.skills.skill_weekly_report_writer import WeeklyReportWriterSkill

    res = WeeklyReportWriterSkill().execute(
        user_request="周期 2026-05-19 至 2026-05-25；本周完成项目 A；下周推进 B"
    )
    assert res["status"] == "success"
    assert res["phase"] == "deliver"
    assert "markdown" in res
    assert "待确认项" not in res["markdown"]


def test_blueprint_drafter_html_mock():
    import gibh_agent.skills.skill_blueprint_drafter as bp_mod

    html = (
        "<!DOCTYPE html><html lang=\"zh-CN\"><head><meta charset=\"UTF-8\" />"
        "<title>Demo</title></head><body><div class=\"diagram-canvas\">X</div></body></html>"
    )
    with patch.object(bp_mod, "llm_chat_text", return_value=html):
        res = bp_mod.BlueprintDrafterSkill().execute(user_request="绘制系统架构图")
    assert res["status"] == "success"
    assert "html_content" in res
    assert "<!DOCTYPE html>" in res["html_content"]


@patch("gibh_agent.skills.skill_mindmap_gen.llm_chat_text")
def test_mindmap_gen_execute_mock(mock_text):
    mock_text.return_value = "graph TD; A[根] --> B[叶]"
    res = MindmapGenSkill().execute(user_request="logic tree")
    assert res["status"] == "success"
    assert "mermaid_code" in res
    assert res["mermaid_code"].startswith("graph TD")
