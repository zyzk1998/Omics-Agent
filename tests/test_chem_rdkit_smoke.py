# -*- coding: utf-8 -*-
"""化学 RDKit 子进程与注册表：与 tests/test_drt_analysis_smoke.py 同级烟雾验收。"""
from __future__ import annotations

import importlib
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
RUNNER = ROOT / "gibh_agent" / "assets" / "chem_rdkit" / "run_chem_tool.py"
GI_SCRIPT = ROOT / "gibh_agent" / "assets" / "gi_absorption" / "gi_absorption.py"
MOLMASS_CLI = ROOT / "gibh_agent" / "assets" / "molmass" / "molmass_cli.py"
ELEMENT_QUERY = ROOT / "gibh_agent" / "assets" / "chemical_element_query" / "element_query.py"


def test_runner_script_exists() -> None:
    assert RUNNER.is_file()
    assert GI_SCRIPT.is_file()
    assert MOLMASS_CLI.is_file()
    assert ELEMENT_QUERY.is_file()


def test_chem_registry_tool_ids() -> None:
    importlib.import_module("gibh_agent.tools.chem_rdkit_tools")
    importlib.import_module("gibh_agent.tools.chem_gi_absorption_tools")
    importlib.import_module("gibh_agent.tools.chem_misc_tools")
    importlib.import_module("gibh_agent.tools.chem_rdkit_batch2_tools")
    from gibh_agent.core.tool_registry import registry

    for tid in (
        "chem_bbb_assessment",
        "chem_pains_filter",
        "chem_brenk_filter",
        "chem_morgan_fingerprint",
        "chem_molecular_weight",
        "chem_tanimoto_similarity",
        "chem_gi_absorption",
        "chem_molmass",
        "chem_element_query",
        "chem_molecule_image",
        "chem_openbabel",
        "chem_aromaticity_perception",
        "chem_lipinski_five",
        "chem_functional_groups",
        "chem_kekulization",
        "chem_pattern_fingerprint",
        "chem_tanimoto_matrix",
    ):
        assert registry.get_tool(tid) is not None, tid


def test_gi_absorption_script_smoke(tmp_path: Path) -> None:
    pytest.importorskip("rdkit")
    out = tmp_path / "gi_out"
    out.mkdir()
    r = subprocess.run(
        [
            sys.executable,
            str(GI_SCRIPT),
            "--smiles",
            "CCO",
            "--output-dir",
            str(out),
            "--format",
            "json",
        ],
        capture_output=True,
        text=True,
        timeout=120,
    )
    assert r.returncode == 0, r.stderr + r.stdout
    assert (out / "chem_summary.json").is_file()
    assert (out / "structure.png").is_file()


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


def test_molmass_cli_json_stdout(tmp_path: Path) -> None:
    r = subprocess.run(
        [sys.executable, str(MOLMASS_CLI), "--formula", "H2O", "--format", "json"],
        capture_output=True,
        text=True,
        timeout=60,
        cwd=str(tmp_path),
    )
    assert r.returncode == 0, r.stderr + r.stdout
    assert '"Status"' in r.stdout or '"status"' in r.stdout.lower()


def test_element_query_json_stdout(tmp_path: Path) -> None:
    r = subprocess.run(
        [sys.executable, str(ELEMENT_QUERY), "--query", "8", "--format", "json"],
        capture_output=True,
        text=True,
        timeout=60,
        cwd=str(tmp_path),
    )
    assert r.returncode == 0, r.stderr + r.stdout


def test_chem_molmass_tool_workspace(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setenv("RESULTS_DIR", str(tmp_path / "results"))
    misc = importlib.import_module("gibh_agent.tools.chem_misc_tools")
    out = misc.chem_molmass(formula_text="H2O")
    assert out["status"] == "success", out
    assert out.get("markdown")
    assert (out.get("json_url") or "").startswith("/results/")


def test_chem_gi_absorption_tool_workspace(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    pytest.importorskip("rdkit")
    monkeypatch.setenv("RESULTS_DIR", str(tmp_path / "results"))

    gi = importlib.import_module("gibh_agent.tools.chem_gi_absorption_tools")

    out = gi.chem_gi_absorption(smiles_text="CCO")
    assert out["status"] == "success", out
    assert out.get("markdown")
    assert out.get("image_urls") and len(out["image_urls"]) >= 1
    assert all(u.startswith("/results/") for u in out["image_urls"])
    assert (out.get("json_url") or "").startswith("/results/")


@pytest.mark.skipif(not shutil.which("obabel"), reason="openbabel CLI not installed")
def test_chem_openbabel_properties_smoke(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setenv("RESULTS_DIR", str(tmp_path / "results"))
    misc = importlib.import_module("gibh_agent.tools.chem_misc_tools")
    out = misc.chem_openbabel(smiles_text="CCO", compute_properties=True)
    assert out["status"] == "success", out
    assert out.get("markdown")
    flat = (out.get("data") or {}).get("result") or {}
    props = flat.get("properties") or {}
    if props:
        mw = props.get("molwt")
        assert mw not in ("??", "", None), props


@pytest.mark.skipif(not shutil.which("obabel"), reason="openbabel CLI not installed")
def test_chem_openbabel_convert_smoke(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setenv("RESULTS_DIR", str(tmp_path / "results"))
    misc = importlib.import_module("gibh_agent.tools.chem_misc_tools")
    out = misc.chem_openbabel(smiles_text="CCO", out_format="smi")
    assert out["status"] == "success", out


def test_obabel_append_line_parser() -> None:
    """多列制表输出须与 --append 顺序对齐（回归：勿只解析最后一列）。"""
    import importlib.util

    p = ROOT / "gibh_agent" / "assets" / "openbabel_converter" / "openbabel_converter.py"
    spec = importlib.util.spec_from_file_location("_obabel_skill_parser", p)
    assert spec and spec.loader
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    fn = mod._parse_obabel_append_line
    assert fn("CCO\t46.069\t-0.0014\t20.23\t1\t1\t0", 6) == [
        "46.069",
        "-0.0014",
        "20.23",
        "1",
        "1",
        "0",
    ]
    assert fn("CCO\t46.069 -0.0014 20.23 1 1 0", 6) == [
        "46.069",
        "-0.0014",
        "20.23",
        "1",
        "1",
        "0",
    ]
