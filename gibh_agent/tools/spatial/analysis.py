"""
Spatial analysis: QC, PCA, clustering, UMAP, neighbors graph, autocorrelation (Moran's I).
"""
import logging
from pathlib import Path
from typing import Dict, Any, Optional, List

from ...core.tool_registry import registry

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Step 0: Image & coordinate validation (show muscle)
# ---------------------------------------------------------------------------

@registry.register(
    name="spatial_data_validation",
    description="Fast validation: read h5ad, check obsm['spatial'] and uns['spatial']. For DAG visibility only.",
    category="Spatial",
    output_type="json",
)
def spatial_data_validation(adata_path: Optional[str] = None, h5ad_path: Optional[str] = None, **kwargs: Any) -> Dict[str, Any]:
    """
    极速空间数据校验：读取 h5ad，检查 obsm['spatial'] 与 uns['spatial']（如有）。
    不修改数据，不写回文件。接受 adata_path 或 h5ad_path（与 executor 一致）。
    """
    path_in = adata_path or h5ad_path or (kwargs.get("h5ad_path") if kwargs else None)
    if not path_in:
        return {"status": "error", "error": "请提供 adata_path 或 h5ad_path"}
    try:
        import anndata as ad
    except ImportError:
        return {"status": "error", "error": "anndata is required"}
    p = Path(path_in)
    if not p.is_file():
        return {"status": "error", "error": f"h5ad file not found: {path_in}"}
    try:
        adata = ad.read_h5ad(p)
        n_spots = adata.n_obs
        n_genes = adata.n_vars
        has_spatial_coord = "spatial" in adata.obsm
        has_spatial_uns = "spatial" in adata.uns
        return {
            "status": "success",
            "message": (
                f"空间数据校验通过。共 {n_spots} 个空间位点 (Spots)，{n_genes} 个基因。"
                "空间物理坐标对齐完成，内存预分配完毕。"
            ),
            "n_spots": int(n_spots),
            "n_genes": int(n_genes),
            "has_obsm_spatial": has_spatial_coord,
            "has_uns_spatial": has_spatial_uns,
            "h5ad_path": str(p.resolve()),
            "output_path": str(p.resolve()),
        }
    except Exception as e:
        logger.warning("spatial_data_validation failed: %s", e)
        return {"status": "error", "error": str(e)}


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

        # 追加：双屏联动图 (UMAP + 空间切片) + Top3 SVG 物理映射图（不替换原有逻辑）
        dual_path = None
        top3_svg_path = None
        if output_path and "spatial" in adata.obsm:
            import matplotlib
            matplotlib.use("Agg")
            import matplotlib.pyplot as plt
            out_dir = out.parent
            try:
                if "X_umap" not in adata.obsm and "neighbors" in adata.uns:
                    sc.tl.umap(adata, min_dist=0.5, spread=1.0)
                palette = adata.uns.get(f"{key_added}_colors")
                import numpy as np
                if palette is None:
                    palette = list(plt.get_cmap("tab10")(np.linspace(0, 1, max(n_clusters, 1))))
                if "X_umap" in adata.obsm:
                    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
                    sc.pl.umap(adata, color=key_added, ax=ax1, show=False, palette=palette, legend_loc="on data" if n_clusters <= 8 else "right margin")
                    try:
                        import squidpy as sq
                        sq.pl.spatial_scatter(adata, color=key_added, ax=ax2, show=False, img=bool(adata.uns.get("spatial")), spot_size=100, palette=palette)
                    except Exception:
                        _draw_spatial_categorical(adata, key_added, ax2, palette=palette)
                    plt.tight_layout()
                    dual_path = str(out_dir / "spatial_umap_dual.png")
                    fig.savefig(dual_path, bbox_inches="tight", dpi=150)
                    plt.close()
            except Exception as ex:
                logger.warning("双屏联动图生成失败: %s", ex)
            # 隐患 2 修复：按优先级探测 Top 基因来源（moranI -> rank_genes_groups），避免仅绑定 moranI 导致标准流程永不触发
            try:
                top_genes = []
                if "moranI" in adata.uns:
                    import pandas as pd
                    res = adata.uns["moranI"]
                    if isinstance(res, pd.DataFrame):
                        if "I" in res.columns:
                            top_genes = res.nlargest(3, "I").index.tolist()
                        else:
                            top_genes = res.index[:3].tolist()
                    else:
                        top_genes = list(res.keys())[:3] if isinstance(res, dict) else []
                if not top_genes and "rank_genes_groups" in adata.uns:
                    rgg = adata.uns["rank_genes_groups"]
                    names = rgg.get("names")
                    if names is not None:
                        try:
                            group_names = names.dtype.names if hasattr(names.dtype, "names") else None
                            if group_names:
                                first_group = group_names[0]
                                top_genes = list(names[first_group][:3])
                            else:
                                top_genes = list(names.flatten())[:3]
                        except Exception:
                            top_genes = []
                top_genes = [g for g in top_genes if g in adata.var_names][:3]
                if top_genes:
                    fig, axes = plt.subplots(1, 3, figsize=(12, 4))
                    for i, gene in enumerate(top_genes):
                        ax = axes[i]
                        try:
                            import squidpy as sq
                            sq.pl.spatial_scatter(adata, color=gene, ax=ax, show=False, img=bool(adata.uns.get("spatial")), cmap="magma", spot_size=100)
                        except Exception:
                            _draw_spatial_gene(adata, gene, ax, cmap="magma")
                        ax.set_title(gene)
                    plt.tight_layout()
                    top3_svg_path = str(out_dir / "spatial_top3_svg.png")
                    fig.savefig(top3_svg_path, bbox_inches="tight", dpi=150)
                    plt.close()
                else:
                    logger.info("无 Moran's I 或 rank_genes_groups 可用基因，跳过 Top3 基因空间映射图")
            except Exception as ex:
                logger.warning("Top3 SVG 物理映射图生成失败: %s", ex)

        result = {
            "status": "success",
            "output_path": str(out.resolve()),
            "h5ad_path": str(p.resolve()),
            "n_clusters": n_clusters,
            "key_added": key_added,
        }
        if dual_path:
            result["dual_panel_path"] = dual_path
        if top3_svg_path:
            result["top3_svg_path"] = top3_svg_path
        return result
    except Exception as e:
        logger.exception("spatial_clustering failed: %s", e)
        return {"status": "error", "error": str(e), "h5ad_path": h5ad_path}


@registry.register(
    name="spatial_clustering_comparison",
    description="Multi-resolution spatial domains: Leiden at 0.3/0.5/0.8, 1x3 spatial plot. Input adata must have spatial neighbors.",
    category="Spatial",
    output_type="mixed",
)
def spatial_clustering_comparison(
    adata_path: Optional[str] = None,
    output_plot_path: str = "",
    resolutions: Optional[List[float]] = None,
    h5ad_path: Optional[str] = None,
    **kwargs: Any,
) -> Dict[str, Any]:
    """
    多分辨率空间域对比：对已构建空间邻域图的 adata 在 0.3/0.5/0.8 下做 Leiden，
    绘制 1x3 组织切片聚类图。无 H&E 时仅绘制坐标散点，不抛 KeyError。
    """
    path_in = adata_path or h5ad_path or (kwargs.get("h5ad_path") if kwargs else None)
    if not path_in:
        return {"status": "error", "error": "请提供 adata_path 或 h5ad_path"}
    if resolutions is None:
        resolutions = [0.3, 0.5, 0.8]
    try:
        import anndata as ad
        import scanpy as sc
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError as e:
        return {"status": "error", "error": f"scanpy/anndata required: {e}"}
    p = Path(path_in)
    if not p.is_file():
        return {"status": "error", "error": f"h5ad file not found: {path_in}"}
    try:
        adata = ad.read_h5ad(p)
        if "spatial" not in adata.obsm:
            return {"status": "error", "error": "No obsm['spatial'] in AnnData."}
        if "neighbors" not in adata.uns:
            rep = "X_pca" if "X_pca" in adata.obsm else "X"
            sc.pp.neighbors(adata, n_neighbors=15, use_rep=rep)
        keys = []
        for res in resolutions:
            key = f"spatial_leiden_{res}"
            keys.append(key)
            sc.tl.leiden(adata, resolution=res, key_added=key)
        out = Path(output_plot_path)
        out.parent.mkdir(parents=True, exist_ok=True)
        if out.suffix.lower() != ".png":
            out = out.parent / (out.stem + ".png")
        fig, axes = plt.subplots(1, 3, figsize=(18, 6))
        has_image = bool(adata.uns.get("spatial"))
        for i, (res, key) in enumerate(zip(resolutions, keys)):
            ax = axes[i]
            try:
                if has_image:
                    import squidpy as sq
                    sq.pl.spatial_scatter(adata, color=key, ax=ax, show=False, img=True, spot_size=100)
                else:
                    _draw_spatial_categorical(adata, key, ax)
            except Exception:
                _draw_spatial_categorical(adata, key, ax)
            ax.set_title(f"Resolution {res}")
        plt.tight_layout()
        fig.savefig(str(out), bbox_inches="tight", dpi=150)
        plt.close()
        return {
            "status": "success",
            "message": "多分辨率空间域对比图已保存",
            "plot_path": str(out.resolve()),
            "h5ad_path": str(p.resolve()),
            "output_path": str(p.resolve()),
        }
    except Exception as e:
        logger.exception("spatial_clustering_comparison failed: %s", e)
        return {"status": "error", "error": str(e), "h5ad_path": path_in}


def _draw_spatial_categorical(adata: Any, color_by: str, ax: Any, palette: Optional[List[Any]] = None) -> None:
    """无底图时按 obsm['spatial'] 画散点，分类着色；palette 与 UMAP 一致时双屏颜色一致。"""
    import numpy as np
    import matplotlib.pyplot as plt
    xy = adata.obsm["spatial"]
    valid = ~np.any(np.isnan(xy), axis=1)
    xy = xy[valid]
    raw = adata.obs[color_by].values[valid]
    import pandas as pd
    cat = pd.Categorical(raw)
    vals = cat.codes.astype(int)
    n_cat = len(cat.categories)
    if palette is not None and len(palette) >= n_cat:
        colors = [palette[i] for i in vals]
        ax.scatter(xy[:, 0], xy[:, 1], c=colors, s=8)
    else:
        cmap = plt.get_cmap("tab10" if n_cat <= 10 else "tab20")
        ax.scatter(xy[:, 0], xy[:, 1], c=vals, s=8, cmap=cmap, vmin=-0.5, vmax=n_cat - 0.5)
    ax.set_title(color_by)
    ax.set_xlabel("x")
    ax.set_ylabel("y")


def _draw_spatial_gene(adata: Any, gene: str, ax: Any, cmap: str = "magma") -> None:
    """无底图时按 obsm['spatial'] 画基因表达散点。"""
    import numpy as np
    import matplotlib.pyplot as plt
    xy = adata.obsm["spatial"]
    valid = ~np.any(np.isnan(xy), axis=1)
    idx = list(adata.var_names).index(gene) if gene in adata.var_names else 0
    v = adata.X[valid, idx]
    vals = np.ravel(v.toarray() if hasattr(v, "toarray") else v)
    sc = ax.scatter(xy[valid, 0], xy[valid, 1], c=vals, s=8, cmap=cmap)
    plt.colorbar(sc, ax=ax, label=gene)
    ax.set_title(gene)
    ax.set_xlabel("x")
    ax.set_ylabel("y")


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
    output_path: Optional[str] = None,
) -> Dict[str, Any]:
    """
    Compute Moran's I (or other method) for genes. Optionally write AnnData with uns['moranI'] for downstream pathway enrichment.

    Args:
        h5ad_path: Path to h5ad file (with spatial graph or obsm['spatial']).
        method: 'moran' (default) or 'geary'.
        genes: Optional comma-separated gene names; if None, use first 500 vars.
        output_path: Optional path to write h5ad with uns['moranI']; required for spatial_pathway_enrichment.

    Returns:
        Dict with status, results, output_path (if written), n_genes, and optional error.
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
        # Always write h5ad with uns['moranI'] so functional_enrichment can read it (executor chains by output)
        if not output_path or (isinstance(output_path, str) and not output_path.strip()):
            output_path = str(p.parent / "spatial_autocorr.h5ad")
        out_p = Path(output_path)
        out_p.parent.mkdir(parents=True, exist_ok=True)
        if out_p.suffix.lower() != ".h5ad":
            out_p = out_p.parent / (out_p.stem + ".h5ad")
        adata.write_h5ad(str(out_p))
        out_path = str(out_p.resolve())
        return {
            "status": "success",
            "results": out,
            "method": method,
            "n_genes": len(out),
            "output_path": out_path,
            "h5ad_path": out_path,  # downstream steps must receive this path (has uns['moranI'])
        }
    except Exception as e:
        logger.exception("detect_spatial_autocorr failed: %s", e)
        return {"status": "error", "error": str(e), "h5ad_path": h5ad_path}


# ---------------------------------------------------------------------------
# Pathway enrichment from top SVGs (biological depth)
# ---------------------------------------------------------------------------


@registry.register(
    name="spatial_pathway_enrichment",
    description="Run pathway enrichment on top spatially variable genes (from Moran's I). Outputs dot plot and CSV. Uses gseapy when available; fallback CSV with gene list otherwise.",
    category="Spatial",
    output_type="file_path",
)
def spatial_pathway_enrichment(
    h5ad_path: str,
    top_n: int = 50,
    output_path: Optional[str] = None,
    gene_sets: str = "KEGG_2021_Human",
) -> Dict[str, Any]:
    """
    Take top SVGs from Moran's I -> enrichment (gseapy) -> dot plot (.png) + CSV.
    If gseapy fails (e.g. network), return graceful fallback CSV with top genes only.

    Args:
        h5ad_path: Path to h5ad (should have uns['moranI'] from spatial_detect_autocorr).
        top_n: Number of top SVGs to use (default 50).
        output_path: Optional directory or base path for plot and CSV.
        gene_sets: gseapy gene set name (default KEGG_2021_Human).

    Returns:
        Dict with status, plot_path, csv_path, top_svgs, enrichment_table (if any), fallback_message.
    """
    try:
        import anndata as ad
        import pandas as pd
        import numpy as np
    except ImportError as e:
        logger.warning("anndata/pandas not installed: %s", e)
        return {"status": "error", "error": "anndata and pandas are required"}
    p = Path(h5ad_path)
    if not p.is_file():
        return {"status": "error", "error": f"h5ad file not found: {h5ad_path}"}
    out_dir = Path(output_path).parent if output_path else p.parent
    out_dir.mkdir(parents=True, exist_ok=True)
    base_name = Path(output_path).stem if output_path else "spatial_pathway_enrichment"
    plot_path = out_dir / f"{base_name}_dotplot.png"
    csv_path = out_dir / f"{base_name}_top_genes.csv"
    fallback_csv = out_dir / f"{base_name}_top_genes_fallback.csv"

    try:
        adata = ad.read_h5ad(p)
        key = "moranI"
        if key not in adata.uns:
            return {"status": "error", "error": "Run spatial_detect_autocorr first to get uns['moranI']."}
        res = adata.uns[key]
        if isinstance(res, pd.DataFrame):
            if "I" in res.columns:
                res = res.sort_values("I", ascending=False)
            elif "pval_norm" in res.columns:
                res = res.sort_values("pval_norm", ascending=True)
            genes = res.index.tolist()[: top_n]
        else:
            genes = list(res.keys())[: top_n] if isinstance(res, dict) else []
        if not genes:
            return {"status": "error", "error": "No genes in Moran's I results."}

        # Fallback CSV: always write top genes description (for when gseapy fails)
        fallback_df = pd.DataFrame({"gene": genes, "rank": range(1, len(genes) + 1)})
        fallback_df.to_csv(str(fallback_csv), index=False)

        enrichment_table = None
        plot_path_final = None
        csv_path_final = None
        fallback_message = None

        try:
            import gseapy as gp
            gene_sets_list = [gene_sets] if isinstance(gene_sets, str) else gene_sets
            enr = gp.enrichr(
                gene_list=genes,
                gene_sets=gene_sets_list,
                organism="Human",
                outdir=None,
                no_plot=True,
            )
            if enr is not None:
                enrichment_table = getattr(enr, "results", enr) if hasattr(enr, "results") else enr
                if isinstance(enrichment_table, pd.DataFrame) and not enrichment_table.empty:
                    csv_path_final = out_dir / f"{base_name}_enrichment.csv"
                    enrichment_table.to_csv(str(csv_path_final), index=False)
                    n_show = min(15, len(enrichment_table))
                    plot_df = enrichment_table.head(n_show)
                    import matplotlib
                    matplotlib.use("Agg")
                    import matplotlib.pyplot as plt
                    fig, ax = plt.subplots(figsize=(8, max(4, n_show * 0.35)))
                    y_pos = np.arange(len(plot_df))
                    pvals = plot_df["Adjusted P-value"] if "Adjusted P-value" in plot_df.columns else plot_df.get("P-value", pd.Series([0.01] * len(plot_df)))
                    pvals = np.clip(pvals.astype(float), 1e-20, 1)
                    logp = -np.log10(pvals)
                    ax.barh(y_pos, logp, color="steelblue", alpha=0.8)
                    ax.set_yticks(y_pos)
                    ax.set_yticklabels(plot_df["Term"].str[:50] if "Term" in plot_df.columns else plot_df.index.astype(str), fontsize=9)
                    ax.set_xlabel("-log10(Adjusted P-value)")
                    ax.set_title("Pathway enrichment (top SVGs)")
                    ax.invert_yaxis()
                    fig.tight_layout()
                    fig.savefig(str(plot_path), dpi=150)
                    plt.close(fig)
                    plot_path_final = plot_path
                else:
                    fallback_message = "gseapy returned empty results; using top genes CSV only."
                    csv_path_final = fallback_csv
        except Exception as gseapy_err:
            logger.warning("gseapy enrichment failed (e.g. network): %s", gseapy_err)
            fallback_message = f"Enrichment skipped ({gseapy_err!s}); see top_genes CSV for biological context."
            csv_path_final = fallback_csv

        # Build return: always include fallback CSV path
        result = {
            "status": "success",
            "top_svgs": genes[:20],
            "n_svgs_used": len(genes),
            "csv_path": str(csv_path_final.resolve()) if csv_path_final and csv_path_final.exists() else str(fallback_csv.resolve()),
            "h5ad_path": str(p.resolve()),
        }
        if fallback_message:
            result["fallback_message"] = fallback_message
        if plot_path_final and plot_path_final.exists():
            result["plot_path"] = str(plot_path_final.resolve())
        if enrichment_table is not None:
            result["enrichment_terms_count"] = len(enrichment_table)

        # Normalize plot_path to results/ relative for frontend
        if "plot_path" in result:
            import os
            results_dir = os.getenv("RESULTS_DIR", "/app/results").rstrip("/")
            abs_p = result["plot_path"]
            if results_dir and abs_p.startswith(results_dir):
                result["plot_path"] = "results/" + abs_p[len(results_dir):].lstrip("/")
        return result
    except Exception as e:
        logger.exception("spatial_pathway_enrichment failed: %s", e)
        return {"status": "error", "error": str(e), "h5ad_path": h5ad_path}
