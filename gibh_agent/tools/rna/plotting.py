"""
单细胞数据可视化工具
"""
import os
import time
import logging
from typing import Dict, Any, Optional, List
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from ...core.tool_registry import registry

logger = logging.getLogger(__name__)


@registry.register(
    name="rna_visualize_qc",
    description="Generates violin plots for QC metrics (n_genes_by_counts, total_counts, pct_counts_mt) to visualize data quality before and after filtering.",
    category="scRNA-seq",
    output_type="file_path"
)
def visualize_qc(
    adata_path: str,
    output_dir: str,
    metrics: Optional[List[str]] = None
) -> Dict[str, Any]:
    """
    生成 QC 指标的小提琴图
    
    Args:
        adata_path: AnnData 文件路径（.h5ad）
        output_dir: 输出目录
        metrics: 要可视化的指标列表（默认：n_genes_by_counts, total_counts, pct_counts_mt）
    
    Returns:
        包含图片路径的字典
    """
    try:
        import scanpy as sc
        
        # 加载数据
        adata = sc.read_h5ad(adata_path)
        
        # 默认指标
        if metrics is None:
            metrics = ['n_genes_by_counts', 'total_counts', 'pct_counts_mt']
        
        # 确保输出目录存在
        os.makedirs(output_dir, exist_ok=True)
        timestamp = int(time.time())
        plot_path = os.path.join(output_dir, f"qc_violin_{timestamp}.png")
        
        # 生成小提琴图
        sc.pl.violin(
            adata,
            metrics,
            jitter=0.4,
            multi_panel=True,
            show=False
        )
        plt.savefig(plot_path, bbox_inches='tight', dpi=300)
        plt.close()
        
        return {
            "status": "success",
            "plot_path": plot_path,
            "metrics": metrics,
            "summary": f"QC 小提琴图已保存: {plot_path}"
        }
    
    except ImportError:
        return {
            "status": "error",
            "error": "scanpy not installed. Please install: pip install scanpy"
        }
    except Exception as e:
        logger.error(f"❌ QC 可视化失败: {e}", exc_info=True)
        return {
            "status": "error",
            "error": str(e)
        }


@registry.register(
    name="rna_visualize_clustering",
    description="Generates UMAP/t-SNE plots colored by clustering results, cell type annotations, or other metadata. Essential for visualizing cell populations and their relationships.",
    category="scRNA-seq",
    output_type="file_path"
)
def visualize_clustering(
    adata_path: str,
    output_dir: str,
    color_by: Optional[List[str]] = None,
    basis: str = "umap",
    ncols: int = 2
) -> Dict[str, Any]:
    """
    生成聚类可视化图（UMAP/t-SNE）
    
    Args:
        adata_path: AnnData 文件路径（.h5ad）
        output_dir: 输出目录
        color_by: 着色依据的列（如 ['leiden', 'cell_type']）
        basis: 降维基础（'umap' 或 'tsne'）
        ncols: 图片列数
    
    Returns:
        包含图片路径的字典
    """
    try:
        import scanpy as sc
        
        # 加载数据
        adata = sc.read_h5ad(adata_path)
        
        # 检查降维结果是否存在
        if basis not in adata.obsm:
            return {
                "status": "error",
                "error": f"降维结果 '{basis}' 不存在。请先运行 rna_umap 或 rna_tsne。"
            }
        
        # 默认着色依据
        if color_by is None:
            # 自动检测可用的着色列
            color_by = []
            if 'leiden' in adata.obs.columns:
                color_by.append('leiden')
            if 'cell_type' in adata.obs.columns or 'celltype' in adata.obs.columns:
                color_by.append('cell_type' if 'cell_type' in adata.obs.columns else 'celltype')
            if 'n_genes_by_counts' in adata.obs.columns:
                color_by.append('n_genes_by_counts')
            
            if not color_by:
                color_by = ['total_counts']  # 默认使用总计数
        
        # 确保输出目录存在
        os.makedirs(output_dir, exist_ok=True)
        timestamp = int(time.time())
        plot_path = os.path.join(output_dir, f"clustering_{basis}_{timestamp}.png")
        
        # 生成 UMAP/t-SNE 图
        sc.pl.embedding(
            adata,
            basis=basis,
            color=color_by,
            ncols=ncols,
            show=False,
            save=f"_{timestamp}.png"
        )
        
        # Scanpy 保存到 figures 目录，我们需要移动它
        scanpy_fig_path = f"figures/{basis}_{timestamp}.png"
        if os.path.exists(scanpy_fig_path):
            import shutil
            shutil.move(scanpy_fig_path, plot_path)
        else:
            # 如果 Scanpy 没有保存，我们自己生成
            fig, axes = plt.subplots(1, len(color_by), figsize=(6 * len(color_by), 5))
            if len(color_by) == 1:
                axes = [axes]
            
            for i, color in enumerate(color_by):
                sc.pl.embedding(
                    adata,
                    basis=basis,
                    color=color,
                    ax=axes[i],
                    show=False
                )
            
            plt.savefig(plot_path, bbox_inches='tight', dpi=300)
            plt.close()
        
        return {
            "status": "success",
            "plot_path": plot_path,
            "basis": basis,
            "color_by": color_by,
            "summary": f"聚类可视化图已保存: {plot_path}"
        }
    
    except ImportError:
        return {
            "status": "error",
            "error": "scanpy not installed. Please install: pip install scanpy"
        }
    except Exception as e:
        logger.error(f"❌ 聚类可视化失败: {e}", exc_info=True)
        return {
            "status": "error",
            "error": str(e)
        }


@registry.register(
    name="rna_visualize_markers",
    description="Generates dot plots or heatmaps for marker genes to visualize their expression across cell types or clusters. Essential for validating cell type annotations.",
    category="scRNA-seq",
    output_type="file_path"
)
def visualize_markers(
    adata_path: str,
    output_dir: str,
    marker_genes: Optional[List[str]] = None,
    groupby: str = "leiden",
    plot_type: str = "dotplot",
    n_top_markers: int = 5
) -> Dict[str, Any]:
    """
    生成 Marker 基因可视化图
    
    Args:
        adata_path: AnnData 文件路径（.h5ad）
        output_dir: 输出目录
        marker_genes: Marker 基因列表（如果为 None，则使用前 n_top_markers 个）
        groupby: 分组依据（默认 'leiden'）
        plot_type: 图表类型（'dotplot' 或 'heatmap'）
        n_top_markers: 如果未指定 marker_genes，使用前 n 个 marker
    
    Returns:
        包含图片路径的字典
    """
    try:
        import scanpy as sc
        
        # 加载数据
        adata = sc.read_h5ad(adata_path)
        
        # 检查分组列是否存在
        if groupby not in adata.obs.columns:
            return {
                "status": "error",
                "error": f"分组列 '{groupby}' 不存在于数据中。请先运行聚类分析。"
            }
        
        # 如果没有指定 marker 基因，尝试从 uns 中获取
        if marker_genes is None:
            if 'rank_genes_groups' in adata.uns:
                # 从 rank_genes_groups 中提取 top markers
                marker_genes = []
                for group in adata.uns['rank_genes_groups']['names'].dtype.names:
                    top_genes = adata.uns['rank_genes_groups']['names'][group][:n_top_markers]
                    marker_genes.extend([g for g in top_genes if g not in marker_genes])
            else:
                return {
                    "status": "error",
                    "error": "未找到 marker 基因。请先运行 rna_find_markers 或手动指定 marker_genes。"
                }
        
        # 过滤不存在的基因
        marker_genes = [g for g in marker_genes if g in adata.var_names]
        
        if not marker_genes:
            return {
                "status": "error",
                "error": "指定的 marker 基因都不存在于数据中。"
            }
        
        # 确保输出目录存在
        os.makedirs(output_dir, exist_ok=True)
        timestamp = int(time.time())
        plot_path = os.path.join(output_dir, f"markers_{plot_type}_{timestamp}.png")
        
        # 生成可视化图
        if plot_type == "dotplot":
            sc.pl.dotplot(
                adata,
                marker_genes,
                groupby=groupby,
                show=False,
                save=f"_{timestamp}.png"
            )
        elif plot_type == "heatmap":
            sc.pl.heatmap(
                adata,
                marker_genes,
                groupby=groupby,
                show=False,
                save=f"_{timestamp}.png"
            )
        else:
            return {
                "status": "error",
                "error": f"不支持的图表类型: {plot_type}。支持的类型: 'dotplot', 'heatmap'"
            }
        
        # Scanpy 保存到 figures 目录，我们需要移动它
        scanpy_fig_path = f"figures/{plot_type}_{timestamp}.png"
        if os.path.exists(scanpy_fig_path):
            import shutil
            shutil.move(scanpy_fig_path, plot_path)
        else:
            # 如果 Scanpy 没有保存，我们自己生成
            if plot_type == "dotplot":
                sc.pl.dotplot(adata, marker_genes, groupby=groupby, show=False)
            else:
                sc.pl.heatmap(adata, marker_genes, groupby=groupby, show=False)
            
            plt.savefig(plot_path, bbox_inches='tight', dpi=300)
            plt.close()
        
        return {
            "status": "success",
            "plot_path": plot_path,
            "plot_type": plot_type,
            "marker_genes": marker_genes,
            "groupby": groupby,
            "summary": f"Marker 基因 {plot_type} 已保存: {plot_path}"
        }
    
    except ImportError:
        return {
            "status": "error",
            "error": "scanpy not installed. Please install: pip install scanpy"
        }
    except Exception as e:
        logger.error(f"❌ Marker 可视化失败: {e}", exc_info=True)
        return {
            "status": "error",
            "error": str(e)
        }

