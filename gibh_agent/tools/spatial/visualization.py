"""
Spatial visualization (scatter plots).
"""
import logging
from pathlib import Path
from typing import Dict, Any, Optional

from ...core.tool_registry import registry

logger = logging.getLogger(__name__)


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
    if output_path is None:
        output_path = str(p.parent / (p.stem + "_spatial_scatter.png"))
    out = Path(output_path)
    out.parent.mkdir(parents=True, exist_ok=True)
    try:
        adata = ad.read_h5ad(p)
        if "spatial" not in adata.obsm:
            return {"status": "error", "error": "No obsm['spatial'] found in AnnData."}
        if color_by not in adata.var_names and color_by not in adata.obs.columns:
            return {"status": "error", "error": f"'{color_by}' not in var_names or obs.columns."}
        sq.pl.spatial_scatter(adata, color=color_by, save=str(out), show=False)
        return {
            "status": "success",
            "plot_path": str(out.resolve()),
            "color_by": color_by,
            "h5ad_path": str(p.resolve()),
        }
    except Exception as e:
        logger.exception("plot_spatial_scatter failed: %s", e)
        return {"status": "error", "error": str(e), "h5ad_path": h5ad_path}
