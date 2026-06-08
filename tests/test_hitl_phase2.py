# -*- coding: utf-8 -*-
"""Phase 2 HITL 挂起与 XML 映射测试。"""
from __future__ import annotations

from gibh_agent.core.orchestrator import AgentOrchestrator
from gibh_agent.tools.hitl_tools import SCENARIO_TO_LS_XML, HITL_SCENARIO_WHITELIST


def test_scenario_to_ls_xml_covers_whitelist():
    assert set(SCENARIO_TO_LS_XML.keys()) == set(HITL_SCENARIO_WHITELIST)
    for xml in SCENARIO_TO_LS_XML.values():
        assert "<View>" in xml


def test_extract_hitl_payload_from_workflow_results():
    results = {
        "status": "hitl_required",
        "hitl": {
            "status": "hitl_required",
            "ls_project_url": "http://127.0.0.1:8082/projects/1/data",
            "ls_project_id": 1,
            "scenario_type": "scrna_cell_type_annotation",
        },
        "steps_details": [],
    }
    payload = AgentOrchestrator._extract_hitl_payload(results)
    assert payload is not None
    assert payload["ls_project_id"] == 1


def test_extract_hitl_payload_from_step_detail():
    results = {
        "steps_details": [
            {
                "status": "hitl_required",
                "step_result": {
                    "data": {
                        "status": "hitl_required",
                        "ls_project_url": "http://127.0.0.1:8082/projects/9/data",
                        "ls_project_id": 9,
                    }
                },
            }
        ]
    }
    payload = AgentOrchestrator._extract_hitl_payload(results)
    assert payload and payload.get("ls_project_id") == 9
