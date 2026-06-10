# -*- coding: utf-8 -*-
"""时光机 execution_snapshot 持久化与 HITL 快照保留测试。"""
from __future__ import annotations

import copy

from gibh_agent.core.execution_snapshot import (
    apply_execution_snapshot_to_state,
    mirror_execution_snapshot_to_message_content,
    update_diagnosis_in_state,
    update_expert_report_in_state,
)
from gibh_agent.core.hitl_resume import _seed_state_snapshot_from_prior


def test_update_expert_report_final_preserves_draft():
    state: dict = {}
    update_expert_report_in_state(state, "## 专家分析报告（初稿）\n\n初稿正文", is_final=False)
    update_expert_report_in_state(state, "## 专家分析报告（最终版）\n\n最终正文", is_final=True)
    ex = state["execution_snapshot"]
    assert "初稿" in ex["expert_report_draft_markdown"]
    assert "最终版" in ex["expert_report_markdown"]
    assert ex.get("hitl_final") is True


def test_apply_execution_snapshot_preserves_diagnosis_and_draft():
    state: dict = {
        "report": {"diagnosis_report": "规划阶段诊断正文"},
        "execution_snapshot": {
            "data_diagnosis_html": "规划阶段诊断正文",
            "expert_report_markdown": "## 专家分析报告（初稿）\n\nA",
            "expert_report_draft_markdown": "## 专家分析报告（初稿）\n\nA",
            "hitl_final": False,
        },
    }
    steps = [{"step_name": "PCA", "status": "success", "step_result": {"status": "success"}}]
    apply_execution_snapshot_to_state(state, steps, workflow_name="单细胞分析工作流")
    ex = state["execution_snapshot"]
    assert ex["data_diagnosis_html"] == "规划阶段诊断正文"
    assert "初稿" in ex["expert_report_draft_markdown"]
    assert len(ex["steps_details"]) == 1
    assert ex["workflow_name"] == "单细胞分析工作流"


def test_seed_state_snapshot_from_prior_does_not_wipe_execution_snapshot():
    prior = {
        "text": "气泡文本",
        "execution_snapshot": {
            "data_diagnosis_html": "诊断 HTML",
            "expert_report_markdown": "## 专家分析报告（初稿）",
            "expert_report_draft_markdown": "## 专家分析报告（初稿）",
            "steps_details": [{"step_name": "Leiden", "status": "success"}],
        },
        "process_log": [{"content": "完成", "state": "completed"}],
    }
    seeded = _seed_state_snapshot_from_prior(prior)
    assert seeded["execution_snapshot"]["data_diagnosis_html"] == "诊断 HTML"
    assert seeded["execution_snapshot"]["steps_details"]
    assert seeded["process_log"]


def test_mirror_execution_snapshot_to_message_content():
    state = {
        "execution_snapshot": {
            "version": 5,
            "data_diagnosis_html": "D",
            "expert_report_markdown": "F",
            "expert_report_draft_markdown": "Draft",
            "steps_details": [],
            "hitl_final": True,
        }
    }
    content: dict = {"state_snapshot": copy.deepcopy(state)}
    mirror_execution_snapshot_to_message_content(content, state)
    assert content["execution_snapshot"]["expert_report_draft_markdown"] == "Draft"
    assert content["execution_snapshot"]["hitl_final"] is True


def test_update_diagnosis_then_expert_chain():
    state: dict = {}
    update_diagnosis_in_state(state, "数据诊断：样本质量良好")
    update_expert_report_in_state(state, "## 专家分析报告（初稿）\n\n结论", is_final=False)
    apply_execution_snapshot_to_state(
        state,
        [{"step_name": "UMAP", "status": "success"}],
        workflow_name="RNA-seq",
    )
    ex = state["execution_snapshot"]
    assert "样本质量" in ex["data_diagnosis_html"]
    assert "结论" in ex["expert_report_markdown"]
    assert ex["steps_details"][0]["step_name"] == "UMAP"
