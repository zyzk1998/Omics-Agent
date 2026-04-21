# -*- coding: utf-8 -*-
"""DRT 脚本与工具：CSV 解析 + 子进程退出码 + 成功载荷含 image_urls/markdown。"""
from __future__ import annotations

import json
import os
import subprocess
import sys
import tempfile
import textwrap
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "gibh_agent" / "assets" / "drt" / "drt_analysis.py"

CSV_SAMPLE = textwrap.dedent(
    """\
    Freq/Hz,Z'/ohm,Z''/ohm
    0.05,0.01,-0.001
    0.5,0.012,-0.005
    5.0,0.015,-0.025
    50.0,0.018,-0.08
    500.0,0.021,-0.15
    5000.0,0.024,-0.2
    """
)


def test_drt_script_runs(tmp_path: Path) -> None:
    assert SCRIPT.is_file()
    inp = tmp_path / "eis.csv"
    inp.write_text(CSV_SAMPLE.strip() + "\n", encoding="utf-8")
    out = tmp_path / "out"
    out.mkdir()
    r = subprocess.run(
        [sys.executable, str(SCRIPT), "--input", str(inp), "--output", str(out), "--lambda", "0.1"],
        capture_output=True,
        text=True,
        timeout=120,
        cwd=str(tmp_path),
    )
    assert r.returncode == 0, r.stderr + r.stdout
    assert (out / "drt_summary.json").is_file()
    assert (out / "drt_distribution.png").is_file()
    with open(out / "drt_summary.json", encoding="utf-8") as f:
        body = json.load(f)
    assert "rmse_complex_fit" in body


def test_eis_drt_tool_returns_image_urls(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
    monkeypatch.setenv("RESULTS_DIR", str(tmp_path / "results"))
    import importlib

    drt_tool = importlib.import_module("gibh_agent.tools.drt_tool")

    inp = tmp_path / "uploads" / "e.csv"
    inp.parent.mkdir(parents=True)
    inp.write_text(CSV_SAMPLE.strip() + "\n", encoding="utf-8")

    out = drt_tool.eis_drt_analysis(file_path=str(inp), regularization_lambda=0.1)
    assert out["status"] == "success", out
    assert out.get("markdown")
    assert out.get("image_urls") and len(out["image_urls"]) >= 1
    assert all(u.startswith("/results/") for u in out["image_urls"])
