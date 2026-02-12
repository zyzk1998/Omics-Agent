"""
Spatial data loading (The Consumer) — Visium / Space Ranger.

Expects standard layout (after Assembler):
  data_dir/
    ├── <matrix>.h5   (name detected dynamically; do not hardcode)
    └── spatial/
         ├── tissue_hires_image.png
         └── scalefactors_json.json

Uses squidpy.read.visium(path, counts_file=...). The counts_file is either:
- Passed from the Sensor via workflow (file_metadata.matrix_file), or
- Resolved dynamically: prefer filtered_feature_bc_matrix.h5, else any *_feature_bc_matrix.h5 in data_dir.
"""
import logging
from pathlib import Path
from typing import Dict, Any, Optional

from ...core.tool_registry import registry

logger = logging.getLogger(__name__)

PREFERRED_COUNTS = "filtered_feature_bc_matrix.h5"
COUNTS_SUFFIX = "_feature_bc_matrix.h5"


def _fill_obsm_spatial_from_tissue_positions(adata: Any, data_dir: Path) -> None:
    """
    当 squidpy 未写入 obsm['spatial'] 时，从 spatial/tissue_positions_list.csv 补全。
    10x 格式：无表头或表头 barcode,in_tissue,array_row,array_col,px_row,px_col。
    """
    import numpy as np
    spatial_dir = data_dir / "spatial"
    if not spatial_dir.is_dir():
        logger.warning("spatial/ not found, cannot fill obsm['spatial']")
        return
    candidates = [
        spatial_dir / "tissue_positions_list.csv",
        spatial_dir / "tissue_positions.csv",
    ]
    tpath = None
    for c in candidates:
        if c.is_file():
            tpath = c
            break
    if not tpath:
        logger.warning("tissue_positions_list.csv not found in spatial/")
        return
    try:
        import pandas as pd
        df = pd.read_csv(tpath, header=None)
        if df.shape[1] < 2:
            logger.warning("tissue_positions has too few columns")
            return
        names = ["barcode", "in_tissue", "array_row", "array_col", "px_row", "px_col"]
        df.columns = names[: df.shape[1]]
        if str(df.iloc[0, 0]).strip().lower() == "barcode":
            df = df.iloc[1:].reset_index(drop=True)
        barcode_col = "barcode"
        if "array_row" in df.columns and "array_col" in df.columns:
            row_col = df[["array_row", "array_col"]].values.astype(float)
        elif "px_row" in df.columns and "px_col" in df.columns:
            row_col = df[["px_row", "px_col"]].values.astype(float)
        else:
            logger.warning("tissue_positions has no array_row/array_col or px_row/px_col")
            return
        barcodes = df[barcode_col].astype(str).values
        order = []
        for name in adata.obs_names:
            idx = np.where(barcodes == name)[0]
            if len(idx) > 0:
                order.append(idx[0])
            else:
                order.append(-1)
        coords = np.full((adata.n_obs, 2), np.nan, dtype=float)
        for i, j in enumerate(order):
            if j >= 0:
                coords[i] = row_col[j]
        if np.any(np.isnan(coords)):
            logger.warning("Some obs not in tissue_positions; spatial scatter may be incomplete")
        adata.obsm["spatial"] = coords
        logger.info("Filled adata.obsm['spatial'] from %s", tpath.name)
    except Exception as e:
        logger.warning("Failed to fill obsm['spatial'] from tissue_positions: %s", e)


def _derive_visium_root_from_path(path: Path) -> Path:
    """
    从文件或目录路径推导 Visium 根目录（含 spatial/ 与 matrix .h5 的目录）。
    - 若为文件（如 tissue_positions_list.csv）：若父目录名为 spatial 则返回 spatial 的父目录，否则返回文件所在目录。
    - 若为目录：直接返回。
    """
    if not path.exists():
        return path
    if path.is_dir():
        return path.resolve()
    # 文件路径：标准为 .../sample/spatial/tissue_positions_list.csv
    parent = path.parent
    if path.name == "tissue_positions_list.csv" or "tissue_positions" in path.name.lower():
        if parent.name == "spatial":
            return parent.parent.resolve()
        return parent.resolve()
    if "spatial" in path.parts:
        for i, part in enumerate(path.parts):
            if part == "spatial" and i > 0:
                return Path(*path.parts[:i]).resolve()
    return parent.resolve()


def _resolve_counts_file(data_dir: Path, counts_file: str) -> str:
    """
    Resolve which matrix .h5 to use (dynamic H5 loading).
    If counts_file is given and exists, use it (e.g. from Sensor's matrix_file).
    Else prefer filtered_feature_bc_matrix.h5; if missing, use first *_feature_bc_matrix.h5 in data_dir.
    """
    preferred = data_dir / counts_file
    if preferred.is_file():
        return counts_file
    for f in data_dir.iterdir():
        if f.is_file() and f.name.lower().endswith(COUNTS_SUFFIX):
            logger.info("Using matrix file %s (standard %s not found)", f.name, counts_file)
            return f.name
    return counts_file


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
        data_dir: Path to the Visium root (contains 'spatial/' and matrix .h5).
        output_path: Optional path to save AnnData as .h5ad; if not set, does not save.
        counts_file: Which counts file to use (default: filtered_feature_bc_matrix.h5).
                    Resolved dynamically if missing: falls back to any *_feature_bc_matrix.h5 in data_dir.

    Returns:
        Dict with status, adata_path (if saved), shape, and optional error.
    """
    try:
        import squidpy as sq
    except ImportError as e:
        logger.warning("squidpy not installed: %s", e)
        return {"status": "error", "error": "squidpy is required. Install with: pip install squidpy>=1.2.0"}
    p = Path(data_dir)
    if p.is_file():
        p = _derive_visium_root_from_path(p)
        logger.info("load_visium_data: 从文件路径推导 data_dir -> %s", p)
    if not p.exists() or not p.is_dir():
        return {"status": "error", "error": f"Data directory not found: {data_dir}"}
    resolved_counts = _resolve_counts_file(p, counts_file)
    try:
        # load_images=False 避免缺少 tissue_hires_image.png 时报错，仅加载表达与坐标
        adata = sq.read.visium(path=str(p), counts_file=resolved_counts, load_images=False)
        n_obs, n_vars = adata.n_obs, adata.n_vars
        # 若 squidpy 未写入 obsm['spatial']（部分版本/load_images=False 时），从 tissue_positions 补全
        if "spatial" not in adata.obsm:
            _fill_obsm_spatial_from_tissue_positions(adata, p)
        out_path = None
        if output_path:
            out = Path(output_path)
        else:
            # 未传 output_path 时写入默认路径，便于后续步骤链使用 adata_path
            import os
            results_dir = Path(os.getenv("RESULTS_DIR", "/app/results"))
            results_dir.mkdir(parents=True, exist_ok=True)
            out = results_dir / "visium_loaded.h5ad"
        out.parent.mkdir(parents=True, exist_ok=True)
        if not out.suffix or out.suffix.lower() != ".h5ad":
            out = out.parent / (out.name + ".h5ad" if not out.suffix else out.name)
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
