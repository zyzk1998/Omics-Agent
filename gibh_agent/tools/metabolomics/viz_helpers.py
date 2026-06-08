# -*- coding: utf-8 -*-
"""代谢组可视化共享逻辑（Top 代谢物筛选 + 样本×代谢物聚类热图）。"""
from __future__ import annotations

import logging
import csv
from pathlib import Path
from typing import Any, Dict, List, Optional

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import seaborn as sns  # noqa: E402

logger = logging.getLogger(__name__)


def _resolve_fs_path(path: str, results_dir: str = "/app/results") -> str:
    """将 /results/... URL 或相对路径解析为容器内可读的绝对路径。"""
    import os

    p = str(path or "").strip()
    if not p:
        return p
    if p.startswith("/results/"):
        base = os.getenv("RESULTS_DIR", results_dir).rstrip("/")
        return str(Path(base) / p[len("/results/") :])
    if p.startswith("results/"):
        base = os.getenv("RESULTS_DIR", results_dir).rstrip("/")
        return str(Path(base) / p[len("results/") :])
    return p


def load_differential_results_from_csv(csv_path: str) -> List[Dict[str, Any]]:
    """从 differential_analysis 输出的 CSV 重建 results 列表。"""
    fs_path = _resolve_fs_path(csv_path)
    if not Path(fs_path).is_file():
        return []
    rows: List[Dict[str, Any]] = []
    with Path(fs_path).open(newline="", encoding="utf-8") as f:
        for r in csv.DictReader(f):
            rows.append(
                {
                    "metabolite": r.get("metabolite", ""),
                    "p_value": float(r.get("p_value") or 1.0),
                    "log2fc": float(r.get("log2fc") or r.get("log2_fold_change") or 0.0),
                    "log2_fold_change": float(r.get("log2_fold_change") or r.get("log2fc") or 0.0),
                    "fdr": float(r.get("fdr") or r.get("fdr_corrected_pvalue") or r.get("p_value") or 1.0),
                    "fdr_corrected_pvalue": float(
                        r.get("fdr_corrected_pvalue") or r.get("fdr") or r.get("p_value") or 1.0
                    ),
                    "significant": str(r.get("significant", "")).lower() == "true",
                    "vip": float(r.get("vip") or 0.0),
                }
            )
    return rows


def hydrate_diff_results_payload(payload: Any) -> Optional[Dict[str, Any]]:
    """
    将 differential_analysis 步骤产出规范化为含 results 列表的字典。
    优先内存 results；缺失时从 output_path / output_file CSV 回退加载。
    """
    if payload is None:
        return None
    if isinstance(payload, str):
        s = payload.strip()
        if s.startswith("<") and s.endswith(">"):
            return None
        if s.lower().endswith(".csv"):
            payload = {"output_path": s}
        else:
            return None
    if not isinstance(payload, dict):
        return None

    inner = payload.get("result") if isinstance(payload.get("result"), dict) else payload
    if not isinstance(inner, dict):
        return None

    results = inner.get("results")
    if isinstance(results, list) and len(results) > 0:
        return inner

    csv_path = inner.get("output_path") or inner.get("output_file") or inner.get("file_path")
    if csv_path and str(csv_path).lower().endswith(".csv"):
        loaded = load_differential_results_from_csv(str(csv_path))
        if loaded:
            out = dict(inner)
            out["results"] = loaded
            return out

    if isinstance(results, list):
        return inner
    return None


def select_top_differential_metabolites(
    results: List[Dict[str, Any]],
    *,
    top_n: int = 20,
    min_vip: float = 1.0,
) -> List[Dict[str, Any]]:
    """与 differential_analysis 一致的 Top 代谢物候选排序（VIP 优先，不足则回退 p 值）。"""
    top_candidates = [
        r for r in results if float(r.get("vip", 0) or 0) > min_vip and r.get("significant", False)
    ]
    top_candidates.sort(key=lambda x: (-float(x.get("vip", 0) or 0), float(x.get("p_value", 1) or 1)))
    top_list = top_candidates[:top_n]
    if len(top_list) < 5:
        sig_list = [r for r in results if r.get("significant", False)]
        sig_list.sort(key=lambda x: float(x.get("p_value", 1) or 1))
        top_list = (top_list + sig_list)[:top_n]
    if len(top_list) < 5:
        top_list = sorted(results, key=lambda x: float(x.get("p_value", 1.0) or 1.0))[:top_n]
    return top_list


def render_metabolite_clustermap(
    file_path: str,
    metabolite_names: List[str],
    output_path: str,
) -> Optional[str]:
    """
    绘制样本 × Top 差异代谢物 Z-score 聚类热图（代谢物为行、样本为列）。

    Returns:
        成功时返回 output_path，无可绘列时返回 None。
    """
    if not metabolite_names:
        return None
    df = pd.read_csv(_resolve_fs_path(file_path), index_col=0)
    plot_cols = [c for c in metabolite_names if c in df.columns]
    if not plot_cols:
        logger.warning("聚类热图：无匹配代谢物列 file=%s names=%s", file_path, metabolite_names[:5])
        return None
    mat = df[plot_cols].T
    mat_z = (mat - mat.mean(axis=1).values.reshape(-1, 1)) / (mat.std(axis=1).values.reshape(-1, 1) + 1e-8)
    out = Path(output_path)
    out.parent.mkdir(parents=True, exist_ok=True)
    g = sns.clustermap(
        mat_z,
        cmap="vlag",
        figsize=(10, max(6, len(plot_cols) * 0.3)),
        cbar_kws={"label": "Z-score"},
    )
    g.fig.suptitle("Top 差异代谢物表达 (Z-score)", y=1.02)
    g.savefig(out, bbox_inches="tight", dpi=150)
    plt.close()
    return str(out)
