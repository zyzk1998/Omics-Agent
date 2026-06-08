# -*- coding: utf-8 -*-
"""visualize_differential_clustermap 工具烟测。"""
from __future__ import annotations

import csv
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
RUN = REPO / "data/results/run_20260602_154908"
CSV = RUN / "human_cachexia_differential_results.csv"
MATRIX = RUN / "human_cachexia_preprocessed.csv"


def _load_diff_results() -> dict:
    with CSV.open() as f:
        rows = list(csv.DictReader(f))
    results = []
    for r in rows:
        results.append(
            {
                "metabolite": r["metabolite"],
                "p_value": float(r["p_value"]),
                "log2fc": float(r.get("log2fc") or r["log2_fold_change"]),
                "fdr": float(r.get("fdr") or r["fdr_corrected_pvalue"]),
                "significant": str(r.get("significant", "")).lower() == "true",
                "vip": float(r.get("vip") or 0),
            }
        )
    return {"status": "success", "results": results}


def test_visualize_differential_clustermap_registered_and_runs(tmp_path):
    from gibh_agent.tools import load_all_tools
    from gibh_agent.core.tool_registry import registry

    load_all_tools()
    fn = registry.get_tool("visualize_differential_clustermap")
    assert fn is not None

    if not MATRIX.is_file() or not CSV.is_file():
        return

    out = tmp_path / "clustermap_test.png"
    res = fn(
        file_path=str(MATRIX),
        diff_results=_load_diff_results(),
        output_path=str(out),
        top_n=15,
    )
    assert res.get("status") == "success", res
    assert out.is_file()
    assert res.get("clustermap_path")
    assert res.get("n_metabolites", 0) >= 5


def test_hydrate_diff_results_from_csv_only_payload():
    from gibh_agent.tools.metabolomics.viz_helpers import hydrate_diff_results_payload

    if not CSV.is_file():
        return
    payload = {
        "status": "success",
        "output_path": str(CSV),
        "results": [],
    }
    hydrated = hydrate_diff_results_payload(payload)
    assert hydrated is not None
    assert len(hydrated.get("results") or []) > 0


def test_metabolomics_workflow_includes_clustermap_step():
    from gibh_agent.core.workflows.registry import WorkflowRegistry

    wf = WorkflowRegistry().get_workflow("Metabolomics")
    dag = wf.get_steps_dag()
    assert "visualize_differential_clustermap" in dag
    assert "differential_analysis" in dag["visualize_differential_clustermap"]
    meta = wf.get_step_metadata("visualize_differential_clustermap")
    assert meta["tool_id"] == "visualize_differential_clustermap"
