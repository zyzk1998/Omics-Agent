"""
Spatial visualization (scatter plots).
"""
import logging
from pathlib import Path
from typing import Dict, Any, Optional

from ...core.tool_registry import registry

logger = logging.getLogger(__name__)


def _scatter_spatial_no_image(adata: Any, color_by: str, out: Path) -> None:
    """当无 tissue 图时用 obsm['spatial'] 做散点图，避免依赖 uns['spatial'][*]['images']。"""
    import numpy as np
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(6, 5))
    xy = adata.obsm["spatial"]
    valid = ~np.any(np.isnan(xy), axis=1)
    xy = xy[valid]
    if color_by in adata.obs.columns:
        raw = adata.obs[color_by].values[valid]
        # leiden 等为字符串类别，matplotlib.scatter 的 c 需要数值
        if hasattr(raw, "dtype") and (raw.dtype == object or (hasattr(raw.dtype, "kind") and raw.dtype.kind in ("U", "O", "S"))):
            import pandas as pd
            cat = pd.Categorical(raw)
            vals = cat.codes.astype(float)
            uniq = cat.categories.tolist()
            n_cat = len(uniq)
            cmap = plt.get_cmap("tab10" if n_cat <= 10 else "tab20")
            sc = ax.scatter(xy[:, 0], xy[:, 1], c=vals, s=5, cmap=cmap, vmin=-0.5, vmax=n_cat - 0.5)
            cbar = fig.colorbar(sc, ax=ax, ticks=np.arange(n_cat), label=color_by)
            cbar.ax.set_yticklabels(uniq)
        else:
            vals = np.asarray(raw, dtype=float)
            sc = ax.scatter(xy[:, 0], xy[:, 1], c=vals, s=5, cmap="viridis")
            fig.colorbar(sc, ax=ax, label=color_by)
    else:
        idx = list(adata.var_names).index(color_by)
        v = adata.X[valid, idx]
        vals = np.ravel(v.toarray() if hasattr(v, "toarray") else v)
        sc = ax.scatter(xy[:, 0], xy[:, 1], c=vals, s=5, cmap="viridis")
        fig.colorbar(sc, ax=ax, label=color_by)
    ax.set_title(f"Spatial: {color_by}")
    ax.set_xlabel("x"); ax.set_ylabel("y")
    fig.tight_layout()
    fig.savefig(str(out), dpi=150)
    plt.close(fig)


@registry.register(
    name="spatial_plot_scatter",
    description="Plot spatial scatter of spots/cells colored by a gene or observation (e.g. cluster). Requires obsm['spatial'].",
    category="Spatial",
    output_type="file_path",
)
def plot_spatial_scatter(
    h5ad_path: str,
    color_by: str = "total_counts",
    output_path: Optional[str] = None,
) -> Dict[str, Any]:
    """
    Generate spatial scatter plot colored by a feature.

    Args:
        h5ad_path: Path to h5ad file (AnnData with obsm['spatial']).
        color_by: Gene name or adata.obs column to color by (default: total_counts).
        output_path: Optional path to save figure; if not set, writes next to h5ad.

    Returns:
        Dict with status, plot_path, and optional error.
    """
    try:
        import anndata as ad
        import squidpy as sq
        import matplotlib
        matplotlib.use("Agg")
    except ImportError as e:
        logger.warning("squidpy/anndata not installed: %s", e)
        return {"status": "error", "error": "squidpy and anndata are required. Install with: pip install squidpy>=1.2.0 anndata"}
    p = Path(h5ad_path)
    if not p.is_file():
        return {"status": "error", "error": f"h5ad file not found: {h5ad_path}"}
    # 空字符串视为未提供，保存到 h5ad 同目录（通常在 results/run_xxx/），避免保存到 cwd 导致 404
    if output_path is None or (isinstance(output_path, str) and not output_path.strip()):
        output_path = str(p.parent / (p.stem + "_spatial_scatter.png"))
    out = Path(output_path)
    if out.suffix.lower() != ".png":
        out = out.parent / (out.stem + ".png" if out.stem else "spatial_scatter.png")
    out.parent.mkdir(parents=True, exist_ok=True)
    try:
        adata = ad.read_h5ad(p)
        if "spatial" not in adata.obsm:
            return {"status": "error", "error": "No obsm['spatial'] found in AnnData."}
        # 解析 color_by：scanpy 常用 n_counts，默认 total_counts；取首个存在的
        if color_by not in adata.var_names and color_by not in adata.obs.columns:
            for fallback in ("n_counts", "total_counts", "n_genes_by_counts"):
                if fallback in adata.obs.columns:
                    color_by = fallback
                    break
            else:
                if len(adata.obs.columns) > 0:
                    color_by = str(adata.obs.columns[0])
                else:
                    return {"status": "error", "error": "'total_counts' / 'n_counts' not in obs; no obs column to color by."}
        if color_by not in adata.var_names and color_by not in adata.obs.columns:
            return {"status": "error", "error": f"'{color_by}' not in var_names or obs.columns."}
        # 无组织图时用 img=False；若报 images 或 color/c 参数错误则用自绘散点图（支持 leiden 等分类列）
        try:
            sq.pl.spatial_scatter(adata, color=color_by, save=str(out), show=False, img=False)
        except Exception as e:
            err = str(e).lower()
            if "images" in err or "image" in err or "img" in err or "'c' argument" in err or "sequence of numbers" in err:
                _scatter_spatial_no_image(adata, color_by, out)
            else:
                try:
                    _scatter_spatial_no_image(adata, color_by, out)
                except Exception:
                    raise
        # 返回前端可识别的路径：优先 results/ 相对路径，便于 <img src="/results/...">
        abs_path = str(out.resolve())
        import os
        results_dir = os.getenv("RESULTS_DIR", "/app/results").rstrip("/")
        if results_dir and abs_path.startswith(results_dir):
            rel = abs_path[len(results_dir):].lstrip("/")
            plot_path_return = f"results/{rel}"
        else:
            plot_path_return = abs_path
        return {
            "status": "success",
            "plot_path": plot_path_return,
            "color_by": color_by,
            "h5ad_path": str(p.resolve()),
        }
    except Exception as e:
        logger.exception("plot_spatial_scatter failed: %s", e)
        return {"status": "error", "error": str(e), "h5ad_path": h5ad_path}
