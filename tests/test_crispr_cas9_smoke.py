# -*- coding: utf-8 -*-
"""CRISPR-Cas9 仿真工具冒烟：子进程 + 资产脚本存在性。"""
from __future__ import annotations

import asyncio
from pathlib import Path

from gibh_agent.tools.crispr_cas9_tool import crispr_cas9_simulation, _crispr_script


def test_crispr_script_present():
    p = _crispr_script()
    assert p.is_file(), f"missing {p}"


def test_crispr_cas9_simulation_runs():
    r = asyncio.run(
        crispr_cas9_simulation(
            guides_text="GACGTCAGTCTAGCTAGCTA",
            target_sequence="ATCGGACGTCAGTCTAGCTAGCTAGGCTAGCTAGCTAAGG",
            cell_line="HEK293",
            result_format="markdown",
            random_seed=42,
        )
    )
    assert r.get("status") == "success", r
    md = (r.get("markdown") or "").strip()
    assert "CRISPR" in md or "Cas9" in md or "gRNA" in md or "Guide" in md, md[:500]
    od = r.get("output_dir")
    assert od and Path(od).is_dir()
    assert (Path(od) / "original.txt").is_file()
