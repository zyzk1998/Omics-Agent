"""
将 run_chem_tool.py 输出与 tests/fixtures/chem_2026_04_20/expected 下黄金 JSON 对齐（需 RDKit）。
"""
from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
FIX = ROOT / "tests/fixtures/chem_2026_04_20"
EXPECTED = FIX / "expected"
RUNNER = ROOT / "gibh_agent/assets/chem_rdkit/run_chem_tool.py"
SINGLE = FIX / "test_smiles_single.txt"
PAIR = FIX / "test_smiles_pair.txt"


def _run_tool(tool: str, out: Path, file_arg: Path) -> dict:
    r = subprocess.run(
        [sys.executable, str(RUNNER), "--tool", tool, "--file", str(file_arg), "--output-dir", str(out)],
        capture_output=True,
        text=True,
        timeout=120,
    )
    assert r.returncode == 0, (r.stderr, r.stdout)
    js = out / "chem_summary.json"
    assert js.is_file(), f"missing {js}"
    return json.loads(js.read_text(encoding="utf-8"))


@pytest.mark.parametrize(
    ("tool", "golden_name"),
    [
        ("bbb", "bbb.chem_summary.json"),
        ("pains", "pains.chem_summary.json"),
        ("brenk", "brenk.chem_summary.json"),
        ("morgan_fp", "morgan_fp.chem_summary.json"),
        ("mol_weight", "mol_weight.chem_summary.json"),
    ],
)
def test_chem_golden_single(tool: str, golden_name: str, tmp_path: Path) -> None:
    pytest.importorskip("rdkit")
    assert SINGLE.is_file(), SINGLE
    exp_path = EXPECTED / golden_name
    assert exp_path.is_file(), exp_path
    expected = json.loads(exp_path.read_text(encoding="utf-8"))
    got = _run_tool(tool, tmp_path / tool, SINGLE)
    assert got == expected


def test_chem_golden_tanimoto(tmp_path: Path) -> None:
    pytest.importorskip("rdkit")
    assert PAIR.is_file(), PAIR
    exp_path = EXPECTED / "tanimoto.chem_summary.json"
    assert exp_path.is_file(), exp_path
    expected = json.loads(exp_path.read_text(encoding="utf-8"))
    got = _run_tool("tanimoto", tmp_path / "tanimoto", PAIR)
    assert got == expected


def test_fixture_png_sidecars_when_rdkit(tmp_path: Path) -> None:
    pytest.importorskip("rdkit")
    got = _run_tool("mol_weight", tmp_path / "mw", SINGLE)
    assert (tmp_path / "mw" / "structure.png").is_file()
    assert got["tool"] == "mol_weight"
