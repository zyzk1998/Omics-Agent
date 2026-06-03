#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
将已跑通组学管线的真实任务产出复制到前端静态目录，供技能详情 demo_visualization 引用。

默认来源（宿主机）：
  转录组：data/results/run_20260602_160748/
  空间组：data/results/run_20260602_164157/ + visium_loaded_spatial_scatter.png

用法（仓库根目录）：
  python3 scripts/build_omics_pipeline_demo_assets.py
"""
from __future__ import annotations

import shutil
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
OUT = REPO / "services/nginx/html/assets/images/demos/pipelines"

COPY_MAP = {
    "transcriptomics/qc_violin.png": REPO / "data/results/run_20260602_160748/qc_violin_1780387669.png",
    "transcriptomics/hvg_mean_variance.png": REPO / "data/results/run_20260602_160748/hvg_1780387671.png",
    "transcriptomics/umap_leiden.png": REPO / "data/results/run_20260602_160748/umap_1780387684.png",
    "transcriptomics/umap_multires.png": REPO / "data/results/run_20260602_160748/multires_leiden_umap.png",
    "transcriptomics/marker_dotplot.png": REPO / "data/results/run_20260602_160748/dotplot_celltype_1780387691.png",
    "spatial/spatial_scatter.png": REPO / "data/results/visium_loaded_spatial_scatter.png",
    "spatial/spatial_multires.png": REPO / "data/results/run_20260602_164157/spatial_multires_comparison.png",
    "spatial/spatial_autocorr.png": REPO / "data/results/spatial_autocorr_spatial_scatter.png",
}


def main() -> int:
    missing = []
    for rel, src in COPY_MAP.items():
        dest = OUT / rel
        dest.parent.mkdir(parents=True, exist_ok=True)
        if not src.is_file():
            missing.append(str(src))
            continue
        shutil.copy2(src, dest)
        print(f"OK {rel} <- {src.name}")
    if missing:
        print("WARN missing sources:", *missing, sep="\n  ")
        return 1
    print(f"Done. Assets under {OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
