# -*- coding: utf-8 -*-
"""化学 RDKit 子进程与注册表：与 tests/test_drt_analysis_smoke.py 同级烟雾验收。"""
from __future__ import annotations

import importlib
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
RUNNER = ROOT / "gibh_agent" / "assets" / "chem_rdkit" / "run_chem_tool.py"


def test_runner_script_exists() -> None:
    assert RUNNER.is_file()


def test_chem_registry_tool_ids() -> None:
    importlib.import_module("gibh_agent.tools.chem_rdkit_tools")
    from gibh_agent.core.tool_registry import registry

    for tid in (
        "chem_bbb_assessment",
        "chem_pains_filter",
        "chem_brenk_filter",
        "chem_morgan_fingerprint",
        "chem_molecular_weight",
        "chem_tanimoto_similarity",
    ):
        assert registry.get_tool(tid) is not None, tid


def test_run_chem_tool_mol_weight_smoke(tmp_path: Path) -> None:
    pytest.importorskip("rdkit")
    out = tmp_path / "o"
    out.mkdir()
    r = subprocess.run(
        [
            sys.executable,
            str(RUNNER),
            "--tool",
            "mol_weight",
            "--smiles",
            "CCO",
            "--output-dir",
            str(out),
        ],
        capture_output=True,
        text=True,
        timeout=120,
    )
    assert r.returncode == 0, r.stderr + r.stdout
    assert (out / "chem_summary.json").is_file()
    assert (out / "structure.png").is_file()


def test_chem_tool_returns_workspace_fields(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    pytest.importorskip("rdkit")
    monkeypatch.setenv("RESULTS_DIR", str(tmp_path / "results"))

    ct = importlib.import_module("gibh_agent.tools.chem_rdkit_tools")

    inp = tmp_path / "smi.txt"
    inp.parent.mkdir(parents=True, exist_ok=True)
    inp.write_text("CCO\n", encoding="utf-8")

    out = ct.chem_molecular_weight(file_path=str(inp))
    assert out["status"] == "success", out
    assert out.get("markdown")
    assert out.get("image_urls") and len(out["image_urls"]) >= 1
    assert all(u.startswith("/results/") for u in out["image_urls"])
    assert (out.get("json_url") or "").startswith("/results/")
