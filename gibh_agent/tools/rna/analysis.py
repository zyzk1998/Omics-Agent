"""
单细胞数据分析工具 - 标准化、降维、聚类
"""
import os
import time
import logging
from typing import Dict, Any, Optional
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from ...core.tool_registry import registry

logger = logging.getLogger(__name__)


@registry.register(
    name="rna_normalize",
    description="Normalizes single-cell RNA-seq data using total count normalization and log transformation. This is a standard preprocessing step before downstream analysis.",
    category="scRNA-seq",
    output_type="json"
)
def run_normalize(
    adata_path: str,
    target_sum: float = 1e4,
    output_dir: Optional[str] = None
) -> Dict[str, Any]:
    """
    执行数据标准化
    
    Args:
        adata_path: AnnData 文件路径（.h5ad）
        target_sum: 标准化目标总和（默认 10000）
        output_dir: 输出目录（可选）
    
    Returns:
        标准化结果字典
    """
    try:
        import scanpy as sc
        
        # 加载数据
        adata = sc.read_h5ad(adata_path)
        
        # 标准化
        sc.pp.normalize_total(adata, target_sum=target_sum)
        sc.pp.log1p(adata)
        
        # 保存结果
        output_h5ad = None
        if output_dir:
            os.makedirs(output_dir, exist_ok=True)
            output_h5ad = os.path.join(output_dir, "normalized.h5ad")
            adata.write(output_h5ad)
        
        return {
            "status": "success",
            "output_h5ad": output_h5ad,
            "summary": "LogNormalize 完成"
        }
    
    except ImportError:
        return {
            "status": "error",
            "error": "scanpy not installed. Please install: pip install scanpy"
        }
    except Exception as e:
        logger.error(f"❌ 标准化失败: {e}", exc_info=True)
        return {
            "status": "error",
            "error": str(e)
        }


@registry.register(
    name="rna_hvg",
    description="Identifies highly variable genes (HVG) in single-cell RNA-seq data. These genes show high variability across cells and are used for downstream dimensionality reduction.",
    category="scRNA-seq",
    output_type="mixed"
)
def run_hvg(
    adata_path: str,
    n_top_genes: int = 2000,
    output_dir: Optional[str] = None
) -> Dict[str, Any]:
    """
    寻找高变基因
    
    Args:
        adata_path: AnnData 文件路径（.h5ad）
        n_top_genes: 选择的高变基因数量
        output_dir: 输出目录（可选）
    
    Returns:
        高变基因结果字典
    """
    try:
        import scanpy as sc
        
        # 加载数据
        adata = sc.read_h5ad(adata_path)
        
        # 寻找高变基因
        sc.pp.highly_variable_genes(adata, n_top_genes=n_top_genes)
        
        # 生成图
        plot_path = None
        if output_dir:
            os.makedirs(output_dir, exist_ok=True)
            timestamp = int(time.time())
            plot_path = os.path.join(output_dir, f"hvg_{timestamp}.png")
            
            sc.pl.highly_variable_genes(adata, show=False)
            plt.savefig(plot_path, bbox_inches='tight', dpi=300)
            plt.close()
        
        # 过滤高变基因
        adata._inplace_subset_var(adata.var['highly_variable'])
        
        # 保存结果
        output_h5ad = None
        if output_dir:
            output_h5ad = os.path.join(output_dir, "hvg_filtered.h5ad")
            adata.write(output_h5ad)
        
        return {
            "status": "success",
            "n_hvg": n_top_genes,
            "plot_path": plot_path,
            "output_h5ad": output_h5ad,
            "summary": f"筛选 {n_top_genes} 个高变基因"
        }
    
    except ImportError:
        return {
            "status": "error",
            "error": "scanpy not installed. Please install: pip install scanpy"
        }
    except Exception as e:
        logger.error(f"❌ 高变基因筛选失败: {e}", exc_info=True)
        return {
            "status": "error",
            "error": str(e)
        }


@registry.register(
    name="rna_scale",
    description="Scales single-cell RNA-seq data to unit variance and zero mean. This is typically performed after normalization and HVG selection, before PCA.",
    category="scRNA-seq",
    output_type="json"
)
def run_scale(
    adata_path: str,
    max_value: float = 10.0,
    output_dir: Optional[str] = None
) -> Dict[str, Any]:
    """
    数据缩放
    
    Args:
        adata_path: AnnData 文件路径（.h5ad）
        max_value: 最大缩放值（用于裁剪）
        output_dir: 输出目录（可选）
    
    Returns:
        缩放结果字典
    """
    try:
        import scanpy as sc
        
        # 加载数据
        adata = sc.read_h5ad(adata_path)
        
        # 缩放
        sc.pp.scale(adata, max_value=max_value)
        
        # 保存结果
        output_h5ad = None
        if output_dir:
            os.makedirs(output_dir, exist_ok=True)
            output_h5ad = os.path.join(output_dir, "scaled.h5ad")
            adata.write(output_h5ad)
        
        return {
            "status": "success",
            "output_h5ad": output_h5ad,
            "summary": "数据缩放完成"
        }
    
    except ImportError:
        return {
            "status": "error",
            "error": "scanpy not installed. Please install: pip install scanpy"
        }
    except Exception as e:
        logger.error(f"❌ 数据缩放失败: {e}", exc_info=True)
        return {
            "status": "error",
            "error": str(e)
        }


@registry.register(
    name="rna_pca",
    description="Performs Principal Component Analysis (PCA) on single-cell RNA-seq data for dimensionality reduction. This is a key step before building neighborhood graphs and clustering.",
    category="scRNA-seq",
    output_type="mixed"
)
def run_pca(
    adata_path: str,
    n_comps: int = 50,
    svd_solver: str = "arpack",
    output_dir: Optional[str] = None
) -> Dict[str, Any]:
    """
    PCA 降维
    
    Args:
        adata_path: AnnData 文件路径（.h5ad）
        n_comps: 主成分数量
        svd_solver: SVD 求解器（"arpack" 或 "auto"）
        output_dir: 输出目录（可选）
    
    Returns:
        PCA 结果字典
    """
    try:
        import scanpy as sc
        
        # 加载数据
        adata = sc.read_h5ad(adata_path)
        
        # PCA
        sc.tl.pca(adata, n_comps=n_comps, svd_solver=svd_solver)
        
        # 生成方差解释图
        plot_path = None
        if output_dir:
            os.makedirs(output_dir, exist_ok=True)
            timestamp = int(time.time())
            plot_path = os.path.join(output_dir, f"pca_variance_{timestamp}.png")
            
            sc.pl.pca_variance_ratio(adata, log=True, show=False)
            plt.savefig(plot_path, bbox_inches='tight', dpi=300)
            plt.close()
        
        # 保存结果
        output_h5ad = None
        if output_dir:
            output_h5ad = os.path.join(output_dir, "pca.h5ad")
            adata.write(output_h5ad)
        
        # 提取解释方差
        explained_variance = {}
        if 'pca' in adata.uns:
            variance_ratio = adata.uns['pca']['variance_ratio']
            for i in range(min(10, len(variance_ratio))):
                explained_variance[f"PC{i+1}"] = float(variance_ratio[i])
        
        return {
            "status": "success",
            "n_comps": n_comps,
            "explained_variance": explained_variance,
            "plot_path": plot_path,
            "output_h5ad": output_h5ad,
            "summary": "PCA 降维完成"
        }
    
    except ImportError:
        return {
            "status": "error",
            "error": "scanpy not installed. Please install: pip install scanpy"
        }
    except Exception as e:
        logger.error(f"❌ PCA 失败: {e}", exc_info=True)
        return {
            "status": "error",
            "error": str(e)
        }


@registry.register(
    name="rna_neighbors",
    description="Computes neighborhood graph for single-cell data using PCA space. This graph is used for clustering and UMAP visualization.",
    category="scRNA-seq",
    output_type="json"
)
def run_neighbors(
    adata_path: str,
    n_neighbors: int = 10,
    n_pcs: int = 40,
    output_dir: Optional[str] = None
) -> Dict[str, Any]:
    """
    计算邻居图
    
    Args:
        adata_path: AnnData 文件路径（.h5ad）
        n_neighbors: 邻居数量
        n_pcs: 使用的 PC 数量
        output_dir: 输出目录（可选）
    
    Returns:
        邻居图结果字典
    """
    try:
        import scanpy as sc
        
        # 加载数据
        adata = sc.read_h5ad(adata_path)
        
        # 计算邻居
        sc.pp.neighbors(adata, n_neighbors=n_neighbors, n_pcs=n_pcs)
        
        # 保存结果
        output_h5ad = None
        if output_dir:
            os.makedirs(output_dir, exist_ok=True)
            output_h5ad = os.path.join(output_dir, "neighbors.h5ad")
            adata.write(output_h5ad)
        
        return {
            "status": "success",
            "n_neighbors": n_neighbors,
            "n_pcs": n_pcs,
            "output_h5ad": output_h5ad,
            "summary": "邻接图构建完成"
        }
    
    except ImportError:
        return {
            "status": "error",
            "error": "scanpy not installed. Please install: pip install scanpy"
        }
    except Exception as e:
        logger.error(f"❌ 邻居图计算失败: {e}", exc_info=True)
        return {
            "status": "error",
            "error": str(e)
        }


@registry.register(
    name="rna_clustering",
    description="Performs Leiden clustering on single-cell data using the neighborhood graph. Leiden clustering is a widely used method for identifying cell populations.",
    category="scRNA-seq",
    output_type="json"
)
def run_clustering(
    adata_path: str,
    resolution: float = 0.5,
    algorithm: str = "leiden",
    output_dir: Optional[str] = None
) -> Dict[str, Any]:
    """
    Leiden 聚类
    
    Args:
        adata_path: AnnData 文件路径（.h5ad）
        resolution: 聚类分辨率（越高，簇越多）
        algorithm: 聚类算法（"leiden" 或 "louvain"）
        output_dir: 输出目录（可选）
    
    Returns:
        聚类结果字典
    """
    try:
        import scanpy as sc
        
        # 加载数据
        adata = sc.read_h5ad(adata_path)
        
        # 聚类
        if algorithm == "leiden":
            sc.tl.leiden(adata, resolution=resolution)
            cluster_key = "leiden"
        elif algorithm == "louvain":
            sc.tl.louvain(adata, resolution=resolution)
            cluster_key = "louvain"
        else:
            return {
                "status": "error",
                "error": f"Unknown clustering algorithm: {algorithm}"
            }
        
        n_clusters = len(adata.obs[cluster_key].unique())
        
        # 保存结果
        output_h5ad = None
        if output_dir:
            os.makedirs(output_dir, exist_ok=True)
            output_h5ad = os.path.join(output_dir, f"{algorithm}_clustered.h5ad")
            adata.write(output_h5ad)
        
        return {
            "status": "success",
            "algorithm": algorithm,
            "resolution": resolution,
            "n_clusters": n_clusters,
            "cluster_key": cluster_key,
            "output_h5ad": output_h5ad,
            "summary": f"{algorithm.capitalize()} 聚类 (Res={resolution}): {n_clusters} 个簇"
        }
    
    except ImportError:
        return {
            "status": "error",
            "error": "scanpy not installed. Please install: pip install scanpy"
        }
    except Exception as e:
        logger.error(f"❌ 聚类失败: {e}", exc_info=True)
        return {
            "status": "error",
            "error": str(e)
        }


@registry.register(
    name="rna_umap",
    description="Generates UMAP (Uniform Manifold Approximation and Projection) visualization for single-cell data. UMAP is a popular dimensionality reduction technique for visualizing cell populations.",
    category="scRNA-seq",
    output_type="mixed"
)
def run_umap(
    adata_path: str,
    color_by: Optional[str] = None,
    output_dir: Optional[str] = None
) -> Dict[str, Any]:
    """
    UMAP 可视化
    
    Args:
        adata_path: AnnData 文件路径（.h5ad）
        color_by: 着色依据（如 "leiden", "total_counts" 等）
        output_dir: 输出目录（可选）
    
    Returns:
        UMAP 结果字典
    """
    try:
        import scanpy as sc
        
        # 加载数据
        adata = sc.read_h5ad(adata_path)
        
        # 计算 UMAP
        sc.tl.umap(adata)
        
        # 生成图
        plot_path = None
        if output_dir:
            os.makedirs(output_dir, exist_ok=True)
            timestamp = int(time.time())
            
            # 确定着色依据
            if color_by is None:
                # 自动选择：优先使用聚类结果
                if 'leiden' in adata.obs.columns:
                    color_by = 'leiden'
                elif 'louvain' in adata.obs.columns:
                    color_by = 'louvain'
                else:
                    color_by = 'total_counts'
            
            plot_path = os.path.join(output_dir, f"umap_{timestamp}.png")
            
            fig, ax = plt.subplots(figsize=(8, 6))
            sc.pl.umap(
                adata, 
                color=[color_by], 
                ax=ax, 
                show=False, 
                title="UMAP",
                legend_loc='on data' if color_by in ['leiden', 'louvain'] else 'right margin',
                frameon=False
            )
            plt.savefig(plot_path, bbox_inches='tight', dpi=300)
            plt.close()
        
        # 保存结果
        output_h5ad = None
        if output_dir:
            output_h5ad = os.path.join(output_dir, "umap.h5ad")
            adata.write(output_h5ad)
        
        return {
            "status": "success",
            "plot_path": plot_path,
            "output_h5ad": output_h5ad,
            "summary": "UMAP 生成完毕"
        }
    
    except ImportError:
        return {
            "status": "error",
            "error": "scanpy not installed. Please install: pip install scanpy"
        }
    except Exception as e:
        logger.error(f"❌ UMAP 失败: {e}", exc_info=True)
        return {
            "status": "error",
            "error": str(e)
        }


@registry.register(
    name="rna_tsne",
    description="Generates t-SNE (t-Distributed Stochastic Neighbor Embedding) visualization for single-cell data. t-SNE is an alternative dimensionality reduction method to UMAP.",
    category="scRNA-seq",
    output_type="mixed"
)
def run_tsne(
    adata_path: str,
    color_by: Optional[str] = None,
    output_dir: Optional[str] = None
) -> Dict[str, Any]:
    """
    t-SNE 可视化
    
    Args:
        adata_path: AnnData 文件路径（.h5ad）
        color_by: 着色依据（如 "leiden" 等）
        output_dir: 输出目录（可选）
    
    Returns:
        t-SNE 结果字典
    """
    try:
        import scanpy as sc
        
        # 加载数据
        adata = sc.read_h5ad(adata_path)
        
        # 检查细胞数（t-SNE 对大数据集较慢）
        if adata.n_obs > 5000:
            return {
                "status": "warning",
                "message": "细胞数过多（>5000），跳过 t-SNE（建议使用 UMAP）",
                "n_obs": adata.n_obs
            }
        
        # 计算 t-SNE
        sc.tl.tsne(adata)
        
        # 生成图
        plot_path = None
        if output_dir:
            os.makedirs(output_dir, exist_ok=True)
            timestamp = int(time.time())
            
            # 确定着色依据
            if color_by is None:
                if 'leiden' in adata.obs.columns:
                    color_by = 'leiden'
                elif 'louvain' in adata.obs.columns:
                    color_by = 'louvain'
                else:
                    color_by = 'total_counts'
            
            plot_path = os.path.join(output_dir, f"tsne_{timestamp}.png")
            
            fig, ax = plt.subplots(figsize=(8, 6))
            sc.pl.tsne(
                adata, 
                color=[color_by], 
                ax=ax, 
                show=False, 
                title="t-SNE",
                frameon=False
            )
            plt.savefig(plot_path, bbox_inches='tight', dpi=300)
            plt.close()
        
        # 保存结果
        output_h5ad = None
        if output_dir:
            output_h5ad = os.path.join(output_dir, "tsne.h5ad")
            adata.write(output_h5ad)
        
        return {
            "status": "success",
            "plot_path": plot_path,
            "output_h5ad": output_h5ad,
            "summary": "t-SNE 生成完毕"
        }
    
    except ImportError:
        return {
            "status": "error",
            "error": "scanpy not installed. Please install: pip install scanpy"
        }
    except Exception as e:
        logger.error(f"❌ t-SNE 失败: {e}", exc_info=True)
        return {
            "status": "error",
            "error": str(e)
        }


@registry.register(
    name="rna_find_markers",
    description="Identifies marker genes for each cluster using statistical tests (t-test, Wilcoxon, etc.). Marker genes help characterize cell populations.",
    category="scRNA-seq",
    output_type="json"
)
def run_find_markers(
    adata_path: str,
    cluster_key: str = "leiden",
    method: str = "t-test",
    n_genes: int = 5,
    output_dir: Optional[str] = None
) -> Dict[str, Any]:
    """
    寻找 Marker 基因
    
    Args:
        adata_path: AnnData 文件路径（.h5ad）
        cluster_key: 聚类列名（"leiden" 或 "louvain"）
        method: 统计方法（"t-test", "wilcoxon", "logreg"）
        n_genes: 每个簇返回的基因数量
        output_dir: 输出目录（可选）
    
    Returns:
        Marker 基因结果字典
    """
    try:
        import scanpy as sc
        import pandas as pd
        
        # 加载数据
        adata = sc.read_h5ad(adata_path)
        
        # 检查是否有聚类结果
        if cluster_key not in adata.obs.columns:
            return {
                "status": "error",
                "error": f"Cluster key '{cluster_key}' not found in data. Please run clustering first."
            }
        
        # 寻找 Marker 基因
        sc.tl.rank_genes_groups(adata, cluster_key, method=method)
        
        # 提取结果
        result = adata.uns['rank_genes_groups']
        groups = result['names'].dtype.names
        
        # 构建 Marker 基因表格
        markers_data = {}
        for group in groups:
            markers_data[f"{group}_names"] = result['names'][group][:n_genes].tolist()
            markers_data[f"{group}_pvals"] = result['pvals'][group][:n_genes].tolist()
            if 'logfoldchanges' in result:
                markers_data[f"{group}_logfc"] = result['logfoldchanges'][group][:n_genes].tolist()
        
        markers_df = pd.DataFrame(markers_data)
        
        # 保存结果
        output_csv = None
        if output_dir:
            os.makedirs(output_dir, exist_ok=True)
            output_csv = os.path.join(output_dir, "markers.csv")
            markers_df.to_csv(output_csv, index=False)
        
        return {
            "status": "success",
            "method": method,
            "n_clusters": len(groups),
            "n_genes_per_cluster": n_genes,
            "markers_table": markers_df.to_dict(orient='records'),
            "output_csv": output_csv,
            "summary": "Marker 基因鉴定完成"
        }
    
    except ImportError:
        return {
            "status": "error",
            "error": "scanpy not installed. Please install: pip install scanpy"
        }
    except Exception as e:
        logger.error(f"❌ Marker 基因鉴定失败: {e}", exc_info=True)
        return {
            "status": "error",
            "error": str(e)
        }

