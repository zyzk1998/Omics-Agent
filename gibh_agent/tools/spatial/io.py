"""
Spatial data loading (Visium / Space Ranger).
"""
import logging
from pathlib import Path
from typing import Dict, Any, Optional

from ...core.tool_registry import registry

logger = logging.getLogger(__name__)


@registry.register(
    name="spatial_load_visium_data",
    description="Load 10x Visium (Space Ranger) data from a directory. Returns path to saved AnnData or in-memory summary. Requires filtered_feature_bc_matrix.h5 and spatial folder.",
    category="Spatial",
    output_type="json",
)
def load_visium_data(
    data_dir: str,
    output_path: Optional[str] = None,
    counts_file: str = "filtered_feature_bc_matrix.h5",
) -> Dict[str, Any]:
    """
    Load Visium data from Space Ranger output directory.

    Args:
        data_dir: Path to the Space Ranger output root (contains 'spatial' and matrix).
        output_path: Optional path to save AnnData as .h5ad; if not set, does not save.
        counts_file: Which counts file to use (default: filtered_feature_bc_matrix.h5).

    Returns:
        Dict with status, adata_path (if saved), shape, and optional error.
    """
    try:
        import squidpy as sq
    except ImportError as e:
        logger.warning("squidpy not installed: %s", e)
        return {"status": "error", "error": "squidpy is required. Install with: pip install squidpy>=1.2.0"}
    p = Path(data_dir)
    if not p.is_dir():
        return {"status": "error", "error": f"Data directory not found: {data_dir}"}
    try:
        adata = sq.read.visium(path=str(p), counts_file=counts_file)
        n_obs, n_vars = adata.n_obs, adata.n_vars
        out_path = None
        if output_path:
            out = Path(output_path)
            out.parent.mkdir(parents=True, exist_ok=True)
            if not out.suffix or out.suffix.lower() != ".h5ad":
                out = out.parent / (out.name + ".h5ad")
            adata.write_h5ad(str(out))
            out_path = str(out.resolve())
        return {
            "status": "success",
            "adata_path": out_path,
            "data_dir": str(p.resolve()),
            "shape": {"n_obs": n_obs, "n_vars": n_vars},
        }
    except Exception as e:
        logger.exception("load_visium_data failed: %s", e)
        return {"status": "error", "error": str(e), "data_dir": data_dir}
