# -*- coding: utf-8 -*-
"""Executor 程序化软 HITL（rna_cell_annotation 后）单元测试。"""
from __future__ import annotations

from unittest.mock import patch

from gibh_agent.core.executor import WorkflowExecutor


def test_provision_scrna_hitl_after_cell_annotation_success():
    ex = WorkflowExecutor()
    step_detail = {
        "status": "success",
        "plot": "/results/run_1/umap_clusters.png",
        "step_result": {"data": {"status": "success"}},
    }
    mock_hitl = {
        "status": "hitl_required",
        "ls_project_id": 99,
        "ls_project_url": "http://127.0.0.1:8082/projects/99/data",
        "scenario_type": "scrna_cell_type_annotation",
    }
    with patch.object(WorkflowExecutor, "_provision_scrna_hitl_after_cell_annotation", return_value=mock_hitl):
        ex._soft_hitl_payload = None
        ex._scrna_hitl_provisioned = False
        hitl = ex._provision_scrna_hitl_after_cell_annotation(
            step_id="rna_cell_annotation",
            tool_id="rna_cell_annotation",
            step_detail=step_detail,
            steps_details=[step_detail],
            workflow_data={"_session_id": "sess-test", "workflow_name": "转录组"},
            workflow_name="转录组",
        )
    assert hitl and hitl.get("status") == "hitl_required"
    assert hitl.get("ls_project_id") == 99


def test_provision_scrna_hitl_skips_non_annotation_step():
    ex = WorkflowExecutor()
    step_detail = {"status": "success", "plot": "/results/x.png", "step_result": {}}
    hitl = ex._provision_scrna_hitl_after_cell_annotation(
        step_id="rna_umap",
        tool_id="rna_umap",
        step_detail=step_detail,
        steps_details=[step_detail],
        workflow_data={},
        workflow_name="RNA",
    )
    assert hitl is None
