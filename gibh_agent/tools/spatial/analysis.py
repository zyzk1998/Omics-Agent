"""
Spatial analysis: QC, PCA, clustering, UMAP, neighbors graph, autocorrelation (Moran's I).
"""
import logging
from pathlib import Path
from typing import Dict, Any, Optional

from ...core.tool_registry import registry

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Preprocessing & dimensionality reduction (Scanpy pipeline)
# ---------------------------------------------------------------------------

@registry.register(
    name="spatial_preprocess_qc",
    description="Spatial QC and preprocessing: filter spots (min_genes), filter genes (min_cells), normalize total, log1p. Writes updated AnnData.",
    category="Spatial",
    output_type="json",
)
def spatial_preprocess_qc(
    h5ad_path: str,
    output_path: Optional[str] = None,
    min_genes: int = 200,
    min_cells: int = 3,
    target_sum: float = 1e4,
) -> Dict[str, Any]:
    """
    Filter spots and genes, normalize total counts, log1p transform.

    Args:
        h5ad_path: Path to h5ad file.
        output_path: Optional path to write updated h5ad; if None, overwrites input.
        min_genes: Minimum genes per spot to keep spot (default 200).
        min_cells: Minimum cells per gene to keep gene (default 3).
        target_sum: Target total counts per cell for normalization (default 1e4).

    Returns:
        Dict with status, output_path (h5ad written), n_obs, n_vars, and optional error.
    """
    try:
        import anndata as ad
        import scanpy as sc
    except ImportError as e:
        logger.warning("scanpy/anndata not installed: %s", e)
        return {"status": "error", "error": "scanpy and anndata are required. Install with: pip install scanpy anndata"}
    p = Path(h5ad_path)
    if not p.is_file():
        return {"status": "error", "error": f"h5ad file not found: {h5ad_path}"}
    try:
        adata = ad.read_h5ad(p)
        n_before_obs, n_before_vars = adata.n_obs, adata.n_vars
        sc.pp.filter_cells(adata, min_genes=min_genes)
        sc.pp.filter_genes(adata, min_cells=min_cells)
        sc.pp.normalize_total(adata, target_sum=target_sum)
        sc.pp.log1p(adata)
        out = Path(output_path) if output_path else p
        out = out.parent / (out.stem + ".h5ad") if out.suffix.lower() != ".h5ad" else out
        out.parent.mkdir(parents=True, exist_ok=True)
        adata.write_h5ad(str(out))
        return {
            "status": "success",
            "output_path": str(out.resolve()),
            "h5ad_path": str(p.resolve()),
            "n_obs": adata.n_obs,
            "n_vars": adata.n_vars,
            "n_obs_filtered": n_before_obs - adata.n_obs,
            "n_vars_filtered": n_before_vars - adata.n_vars,
        }
    except Exception as e:
        logger.exception("spatial_preprocess_qc failed: %s", e)
        return {"status": "error", "error": str(e), "h5ad_path": h5ad_path}


@registry.register(
    name="spatial_pca_reduction",
    description="Run PCA on spatial AnnData (scanpy). Requires log-normalized data (run spatial_preprocess_qc first).",
    category="Spatial",
    output_type="json",
)
def spatial_pca_reduction(
    h5ad_path: str,
    output_path: Optional[str] = None,
    n_comps: int = 50,
    svd_solver: str = "arpack",
) -> Dict[str, Any]:
    """
    Compute PCA and optionally save updated h5ad.

    Args:
        h5ad_path: Path to h5ad (e.g. after spatial_preprocess_qc).
        output_path: Optional path to write h5ad with obsm['X_pca'] and uns['pca'].
        n_comps: Number of components (default 50).
        svd_solver: 'arpack' or 'auto' (default arpack).

    Returns:
        Dict with status, output_path, n_comps, h5ad_path, and optional error.
    """
    try:
        import anndata as ad
        import scanpy as sc
    except ImportError as e:
        logger.warning("scanpy/anndata not installed: %s", e)
        return {"status": "error", "error": "scanpy and anndata are required. Install with: pip install scanpy anndata"}
    p = Path(h5ad_path)
    if not p.is_file():
        return {"status": "error", "error": f"h5ad file not found: {h5ad_path}"}
    try:
        adata = ad.read_h5ad(p)
        n_comps_actual = max(1, min(n_comps, adata.n_obs - 1, adata.n_vars - 1))
        sc.tl.pca(adata, n_comps=n_comps_actual, svd_solver=svd_solver)
        out = Path(output_path) if output_path else p
        out = out.parent / (out.stem + ".h5ad") if out.suffix.lower() != ".h5ad" else out
        out.parent.mkdir(parents=True, exist_ok=True)
        adata.write_h5ad(str(out))
        return {
            "status": "success",
            "output_path": str(out.resolve()),
            "h5ad_path": str(p.resolve()),
            "n_comps": adata.obsm["X_pca"].shape[1],
            "n_obs": adata.n_obs,
        }
    except Exception as e:
        logger.exception("spatial_pca_reduction failed: %s", e)
        return {"status": "error", "error": str(e), "h5ad_path": h5ad_path}


@registry.register(
    name="spatial_clustering",
    description="Compute kNN graph (from PCA) and Leiden clustering. Adds obs['leiden'] to AnnData.",
    category="Spatial",
    output_type="json",
)
def spatial_clustering(
    h5ad_path: str,
    output_path: Optional[str] = None,
    resolution: float = 0.5,
    n_neighbors: int = 15,
    use_rep: str = "X_pca",
    key_added: str = "leiden",
) -> Dict[str, Any]:
    """
    Build neighbors from PCA (or X), then run Leiden clustering.

    Args:
        h5ad_path: Path to h5ad (should have obsm['X_pca'] from spatial_pca_reduction).
        output_path: Optional path to write h5ad with obs[key_added].
        resolution: Leiden resolution (default 0.5).
        n_neighbors: Number of neighbors for graph (default 15).
        use_rep: Embedding to use ('X_pca', 'X', etc.) (default X_pca).
        key_added: Key in obs for cluster labels (default 'leiden').

    Returns:
        Dict with status, output_path, n_clusters, h5ad_path, and optional error.
    """
    try:
        import anndata as ad
        import scanpy as sc
    except ImportError as e:
        logger.warning("scanpy/anndata not installed: %s", e)
        return {"status": "error", "error": "scanpy and anndata are required. Install with: pip install scanpy anndata"}
    p = Path(h5ad_path)
    if not p.is_file():
        return {"status": "error", "error": f"h5ad file not found: {h5ad_path}"}
    try:
        adata = ad.read_h5ad(p)
        if use_rep == "X_pca" and "X_pca" not in adata.obsm:
            n_comps = max(1, min(50, adata.n_obs - 1, adata.n_vars - 1))
            sc.tl.pca(adata, n_comps=n_comps, svd_solver="arpack")
        rep = use_rep if use_rep in adata.obsm else "X"
        sc.pp.neighbors(adata, n_neighbors=n_neighbors, use_rep=rep)
        sc.tl.leiden(adata, resolution=resolution, key_added=key_added)
        n_clusters = int(adata.obs[key_added].astype(str).nunique())
        out = Path(output_path) if output_path else p
        out = out.parent / (out.stem + ".h5ad") if out.suffix.lower() != ".h5ad" else out
        out.parent.mkdir(parents=True, exist_ok=True)
        adata.write_h5ad(str(out))
        return {
            "status": "success",
            "output_path": str(out.resolve()),
            "h5ad_path": str(p.resolve()),
            "n_clusters": n_clusters,
            "key_added": key_added,
        }
    except Exception as e:
        logger.exception("spatial_clustering failed: %s", e)
        return {"status": "error", "error": str(e), "h5ad_path": h5ad_path}


@registry.register(
    name="spatial_umap",
    description="Compute UMAP embedding for visualization (scanpy). Uses existing neighbor graph or PCA.",
    category="Spatial",
    output_type="json",
)
def spatial_umap(
    h5ad_path: str,
    output_path: Optional[str] = None,
    min_dist: float = 0.5,
    spread: float = 1.0,
) -> Dict[str, Any]:
    """
    Run UMAP and add obsm['X_umap'].

    Args:
        h5ad_path: Path to h5ad (should have neighbors from spatial_clustering or run neighbors first).
        output_path: Optional path to write h5ad with obsm['X_umap'].
        min_dist: UMAP min_dist (default 0.5).
        spread: UMAP spread (default 1.0).

    Returns:
        Dict with status, output_path, h5ad_path, and optional error.
    """
    try:
        import anndata as ad
        import scanpy as sc
    except ImportError as e:
        logger.warning("scanpy/anndata not installed: %s", e)
        return {"status": "error", "error": "scanpy and anndata are required. Install with: pip install scanpy anndata"}
    p = Path(h5ad_path)
    if not p.is_file():
        return {"status": "error", "error": f"h5ad file not found: {h5ad_path}"}
    try:
        adata = ad.read_h5ad(p)
        if "neighbors" not in adata.uns:
            if "X_pca" in adata.obsm:
                sc.pp.neighbors(adata, n_neighbors=15, use_rep="X_pca")
            else:
                n_comps = max(1, min(50, adata.n_obs - 1, adata.n_vars - 1))
                sc.tl.pca(adata, n_comps=n_comps, svd_solver="arpack")
                sc.pp.neighbors(adata, n_neighbors=15, use_rep="X_pca")
        sc.tl.umap(adata, min_dist=min_dist, spread=spread)
        out = Path(output_path) if output_path else p
        out = out.parent / (out.stem + ".h5ad") if out.suffix.lower() != ".h5ad" else out
        out.parent.mkdir(parents=True, exist_ok=True)
        adata.write_h5ad(str(out))
        return {
            "status": "success",
            "output_path": str(out.resolve()),
            "h5ad_path": str(p.resolve()),
            "n_obs": adata.n_obs,
        }
    except Exception as e:
        logger.exception("spatial_umap failed: %s", e)
        return {"status": "error", "error": str(e), "h5ad_path": h5ad_path}


# ---------------------------------------------------------------------------
# Spatial graph & autocorrelation (existing)
# ---------------------------------------------------------------------------


@registry.register(
    name="spatial_calculate_neighbors",
    description="Compute spatial neighborhood graph from coordinates (obsm['spatial']). Saves connectivities to the AnnData and optionally writes updated h5ad.",
    category="Spatial",
    output_type="json",
)
def calculate_spatial_neighbors(
    h5ad_path: str,
    output_path: Optional[str] = None,
    coord_type: str = "generic",
    delaunay: bool = True,
) -> Dict[str, Any]:
    """
    Compute spatial graph and optionally save updated AnnData.

    Args:
        h5ad_path: Path to h5ad file (must have obsm['spatial']).
        output_path: Optional path to write AnnData with obsp updated.
        coord_type: 'generic' or 'grid' (default: generic).
        delaunay: Use Delaunay triangulation for graph (default: True).

    Returns:
        Dict with status, output_path (if saved), n_obs, and optional error.
    """
    try:
        import anndata as ad
        import squidpy as sq
    except ImportError as e:
        logger.warning("squidpy/anndata not installed: %s", e)
        return {"status": "error", "error": "squidpy and anndata are required. Install with: pip install squidpy>=1.2.0 anndata"}
    p = Path(h5ad_path)
    if not p.is_file():
        return {"status": "error", "error": f"h5ad file not found: {h5ad_path}"}
    try:
        adata = ad.read_h5ad(p)
        if "spatial" not in adata.obsm:
            return {"status": "error", "error": "No obsm['spatial'] found in AnnData."}
        sq.gr.spatial_neighbors(adata, coord_type=coord_type, delaunay=delaunay)
        out_path = None
        if output_path:
            out = Path(output_path)
            out.parent.mkdir(parents=True, exist_ok=True)
            if out.suffix.lower() != ".h5ad":
                out = out.parent / (out.stem + ".h5ad")
            adata.write_h5ad(str(out))
            out_path = str(out.resolve())
        return {
            "status": "success",
            "output_path": out_path,
            "n_obs": adata.n_obs,
            "h5ad_path": str(p.resolve()),
        }
    except Exception as e:
        logger.exception("calculate_spatial_neighbors failed: %s", e)
        return {"status": "error", "error": str(e), "h5ad_path": h5ad_path}


@registry.register(
    name="spatial_detect_autocorr",
    description="Calculate spatial autocorrelation (e.g. Moran's I) for genes to detect spatially variable genes (SVGs).",
    category="Spatial",
    output_type="json",
)
def detect_spatial_autocorr(
    h5ad_path: str,
    method: str = "moran",
    genes: Optional[str] = None,
) -> Dict[str, Any]:
    """
    Compute Moran's I (or other method) for genes.

    Args:
        h5ad_path: Path to h5ad file (with spatial graph or obsm['spatial']).
        method: 'moran' (default) or 'geary'.
        genes: Optional comma-separated gene names; if None, use first 500 vars.

    Returns:
        Dict with status, results (gene -> I, pval), n_genes, and optional error.
    """
    try:
        import anndata as ad
        import squidpy as sq
        import pandas as pd
    except ImportError as e:
        logger.warning("squidpy/anndata not installed: %s", e)
        return {"status": "error", "error": "squidpy and anndata are required. Install with: pip install squidpy>=1.2.0 anndata"}
    p = Path(h5ad_path)
    if not p.is_file():
        return {"status": "error", "error": f"h5ad file not found: {h5ad_path}"}
    try:
        adata = ad.read_h5ad(p)
        adata.obs_names_make_unique()
        adata.var_names_make_unique()
        if genes:
            gene_list = [g.strip() for g in genes.split(",") if g.strip() in adata.var_names]
            if not gene_list:
                return {"status": "error", "error": f"No valid genes in list; check var_names."}
            adata = adata[:, gene_list]
        if adata.n_vars > 500:
            adata = adata[:, adata.var_names[:500]]
        if adata.n_vars == 0:
            return {"status": "error", "error": "No genes to compute."}
        if "spatial_connectivities" not in adata.obsp and "spatial" in adata.obsm:
            sq.gr.spatial_neighbors(adata, coord_type="generic", delaunay=True)
        if "spatial_connectivities" not in adata.obsp:
            return {"status": "error", "error": "No spatial connectivities; run spatial_calculate_neighbors first or ensure obsm['spatial']."}
        sq.gr.spatial_autocorr(adata, mode=method, genes=adata.var_names.tolist())
        key = "moranI" if method == "moran" else "gearyI"
        if key not in adata.uns:
            return {"status": "error", "error": f"squidpy did not produce {key} in adata.uns."}
        res = adata.uns[key]
        if isinstance(res, pd.DataFrame):
            out = res.to_dict(orient="index")
        else:
            out = dict(res)
        return {
            "status": "success",
            "results": out,
            "method": method,
            "n_genes": len(out),
            "h5ad_path": str(p.resolve()),
        }
    except Exception as e:
        logger.exception("detect_spatial_autocorr failed: %s", e)
        return {"status": "error", "error": str(e), "h5ad_path": h5ad_path}
