"""
单细胞数据质量控制工具
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
    name="rna_qc_filter",
    description="Performs quality control filtering on single-cell RNA-seq data. Filters cells based on gene counts, total counts, and mitochondrial percentage. Calculates QC metrics and generates violin plots.",
    category="scRNA-seq",
    output_type="mixed"
)
def run_qc_filter(
    adata_path: str,
    min_genes: int = 200,
    max_mt: float = 20.0,
    min_cells: int = 3,
    output_dir: Optional[str] = None
) -> Dict[str, Any]:
    """
    执行质量控制过滤
    
    Args:
        adata_path: AnnData 文件路径（.h5ad）或 10x 目录路径
        min_genes: 每个细胞的最小基因数
        max_mt: 线粒体基因的最大百分比
        min_cells: 每个基因的最小细胞数
        output_dir: 输出目录（用于保存图片）
    
    Returns:
        包含以下键的字典:
        - status: "success" 或 "error"
        - n_obs_before: 过滤前细胞数
        - n_obs_after: 过滤后细胞数
        - n_vars_before: 过滤前基因数
        - n_vars_after: 过滤后基因数
        - plot_path: QC 小提琴图路径（如果生成）
        - error: 错误信息（如果失败）
    """
    try:
        import scanpy as sc
        
        # 加载数据
        if os.path.isdir(adata_path):
            adata = sc.read_10x_mtx(adata_path, var_names='gene_symbols', cache=False)
            adata.var_names_make_unique()
        elif adata_path.endswith('.h5ad'):
            adata = sc.read_h5ad(adata_path)
        else:
            adata = sc.read(adata_path)
        
        n_obs_before = adata.n_obs
        n_vars_before = adata.n_vars
        
        # 计算线粒体基因（严格处理 MT- 基因）
        # 支持多种命名约定：MT-, mt-, Mt-, 以及某些物种的 mt-
        adata.var['mt'] = (
            adata.var_names.str.startswith('MT-') |
            adata.var_names.str.startswith('mt-') |
            adata.var_names.str.startswith('Mt-') |
            adata.var_names.str.match('^MT-', case=False) |
            adata.var_names.str.match('^mt-', case=False)
        )
        sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, inplace=True)
        
        # 生成 QC 图
        plot_path = None
        if output_dir:
            os.makedirs(output_dir, exist_ok=True)
            timestamp = int(time.time())
            plot_path = os.path.join(output_dir, f"qc_violin_{timestamp}.png")
            
            sc.pl.violin(
                adata, 
                ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'], 
                jitter=0.4, 
                multi_panel=True, 
                show=False
            )
            plt.savefig(plot_path, bbox_inches='tight', dpi=300)
            plt.close()
        
        # 过滤
        sc.pp.filter_cells(adata, min_genes=min_genes)
        adata = adata[adata.obs.pct_counts_mt < max_mt, :]
        sc.pp.filter_genes(adata, min_cells=min_cells)
        
        n_obs_after = adata.n_obs
        n_vars_after = adata.n_vars
        
        # 保存过滤后的数据（如果指定了输出目录）
        output_h5ad = None
        if output_dir:
            output_h5ad = os.path.join(output_dir, "filtered.h5ad")
            adata.write(output_h5ad)
        
        return {
            "status": "success",
            "n_obs_before": n_obs_before,
            "n_obs_after": n_obs_after,
            "n_vars_before": n_vars_before,
            "n_vars_after": n_vars_after,
            "plot_path": plot_path,
            "output_h5ad": output_h5ad,
            "summary": f"过滤后剩余 {n_obs_after} 个细胞，{n_vars_after} 个基因"
        }
    
    except ImportError:
        return {
            "status": "error",
            "error": "scanpy not installed. Please install: pip install scanpy"
        }
    except Exception as e:
        logger.error(f"❌ QC 过滤失败: {e}", exc_info=True)
        return {
            "status": "error",
            "error": str(e)
        }


@registry.register(
    name="rna_doublet_detection",
    description="Detects doublets (multiple cells in one droplet) in single-cell RNA-seq data using Scrublet or similar methods.",
    category="scRNA-seq",
    output_type="json"
)
def run_doublet_detection(
    adata_path: str,
    method: str = "scrublet",
    expected_doublet_rate: float = 0.1,
    output_dir: Optional[str] = None
) -> Dict[str, Any]:
    """
    执行双联体检测
    
    Args:
        adata_path: AnnData 文件路径（.h5ad）
        method: 检测方法（"scrublet" 或其他）
        expected_doublet_rate: 预期双联体率
        output_dir: 输出目录（可选）
    
    Returns:
        检测结果字典
    """
    try:
        import scanpy as sc
        
        # 加载数据
        adata = sc.read_h5ad(adata_path)
        
        if method == "scrublet":
            try:
                import scrublet as scr
                
                # 运行 Scrublet
                scrub = scr.Scrublet(adata.X, expected_doublet_rate=expected_doublet_rate)
                doublet_scores, predicted_doublets = scrub.scrub_doublets()
                
                # 添加到 obs
                adata.obs['doublet_score'] = doublet_scores
                adata.obs['predicted_doublet'] = predicted_doublets
                
                n_doublets = predicted_doublets.sum()
                doublet_rate = (n_doublets / len(predicted_doublets)) * 100
                
                # 保存结果
                output_h5ad = None
                if output_dir:
                    os.makedirs(output_dir, exist_ok=True)
                    output_h5ad = os.path.join(output_dir, "doublet_detected.h5ad")
                    adata.write(output_h5ad)
                
                return {
                    "status": "success",
                    "n_doublets": int(n_doublets),
                    "doublet_rate": round(doublet_rate, 2),
                    "output_h5ad": output_h5ad,
                    "summary": f"检测到 {n_doublets} 个双联体（{doublet_rate:.2f}%）"
                }
            except ImportError:
                return {
                    "status": "error",
                    "error": "scrublet not installed. Please install: pip install scrublet"
                }
        else:
            return {
                "status": "error",
                "error": f"Unknown doublet detection method: {method}"
            }
    
    except ImportError:
        return {
            "status": "error",
            "error": "scanpy not installed. Please install: pip install scanpy"
        }
    except Exception as e:
        logger.error(f"❌ 双联体检测失败: {e}", exc_info=True)
        return {
            "status": "error",
            "error": str(e)
        }

