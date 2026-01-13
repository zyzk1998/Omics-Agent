"""
代谢组学可视化工具
"""
import logging
from typing import Dict, Any, Optional
from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns

from ...core.tool_registry import registry

logger = logging.getLogger(__name__)


@registry.register(
    name="metabolomics_volcano",
    description="Generates a volcano plot for differential analysis results. Shows log2 fold change vs -log10(p-value) with significance thresholds.",
    category="Metabolomics",
    output_type="image"
)
def plot_volcano(
    diff_results: Dict[str, Any],
    output_path: Optional[str] = None,
    fdr_threshold: float = 0.05,
    log2fc_threshold: float = 1.0
) -> Dict[str, Any]:
    """
    生成火山图
    
    Args:
        diff_results: 差异分析结果（包含 results 列表）
        output_path: 输出文件路径（如果为 None，自动生成）
        fdr_threshold: FDR 阈值
        log2fc_threshold: Log2FC 阈值
    
    Returns:
        包含图片路径的字典
    """
    try:
        results = diff_results.get("results", [])
        
        if not results:
            return {
                "status": "error",
                "error": "差异分析结果为空"
            }
        
        # 准备数据
        df = pd.DataFrame(results)
        df['-log10_fdr'] = -np.log10(df['fdr'].replace(0, 1e-10))
        df['significant'] = (
            (df['fdr'] < fdr_threshold) & 
            (np.abs(df['log2fc']) > log2fc_threshold)
        )
        
        # 生成图片
        if output_path is None:
            output_path = "./results/volcano_plot.png"
        
        output_path_obj = Path(output_path)
        output_path_obj.parent.mkdir(parents=True, exist_ok=True)
        
        plt.figure(figsize=(10, 8))
        
        # 非显著点
        non_sig = df[~df['significant']]
        plt.scatter(non_sig['log2fc'], non_sig['-log10_fdr'], 
                   alpha=0.5, color='gray', s=30, label='Not significant')
        
        # 显著点
        sig = df[df['significant']]
        if len(sig) > 0:
            up = sig[sig['log2fc'] > 0]
            down = sig[sig['log2fc'] < 0]
            
            if len(up) > 0:
                plt.scatter(up['log2fc'], up['-log10_fdr'], 
                           alpha=0.7, color='red', s=50, label='Up-regulated')
            if len(down) > 0:
                plt.scatter(down['log2fc'], down['-log10_fdr'], 
                           alpha=0.7, color='blue', s=50, label='Down-regulated')
        
        # 添加阈值线
        plt.axhline(y=-np.log10(fdr_threshold), color='black', linestyle='--', alpha=0.5)
        plt.axvline(x=log2fc_threshold, color='black', linestyle='--', alpha=0.5)
        plt.axvline(x=-log2fc_threshold, color='black', linestyle='--', alpha=0.5)
        
        plt.xlabel('Log2 Fold Change')
        plt.ylabel('-Log10 FDR')
        plt.title('Volcano Plot')
        plt.legend()
        plt.grid(True, alpha=0.3)
        plt.savefig(output_path, dpi=150, bbox_inches='tight')
        plt.close()
        
        return {
            "status": "success",
            "plot_path": str(output_path)
        }
    
    except Exception as e:
        logger.error(f"❌ 火山图生成失败: {e}", exc_info=True)
        return {
            "status": "error",
            "error": str(e)
        }


@registry.register(
    name="metabolomics_heatmap",
    description="Generates a heatmap for metabolite abundance data. Shows clustering patterns and can highlight group differences.",
    category="Metabolomics",
    output_type="image"
)
def plot_heatmap(
    file_path: str,
    output_path: Optional[str] = None,
    top_n: int = 50,
    group_column: Optional[str] = None
) -> Dict[str, Any]:
    """
    生成热图
    
    Args:
        file_path: 输入数据文件路径（CSV）
        output_path: 输出文件路径（如果为 None，自动生成）
        top_n: 显示前 N 个代谢物（按方差排序）
        group_column: 可选的分组列名（用于添加分组注释）
    
    Returns:
        包含图片路径的字典
    """
    try:
        # 读取数据
        df = pd.read_csv(file_path, index_col=0)
        
        # 提取数值列
        numeric_cols = df.select_dtypes(include=[np.number]).columns
        data = df[numeric_cols]
        
        # 选择前 N 个高方差代谢物
        if len(data.columns) > top_n:
            variances = data.var().sort_values(ascending=False)
            top_metabolites = variances.head(top_n).index
            data = data[top_metabolites]
        
        # 生成图片
        if output_path is None:
            output_path = "./results/heatmap.png"
        
        output_path_obj = Path(output_path)
        output_path_obj.parent.mkdir(parents=True, exist_ok=True)
        
        plt.figure(figsize=(12, 10))
        
        # 如果有分组信息，添加行注释
        row_colors = None
        if group_column and group_column in df.columns:
            groups = df[group_column]
            unique_groups = groups.unique()
            colors = plt.cm.Set3(np.linspace(0, 1, len(unique_groups)))
            group_color_map = dict(zip(unique_groups, colors))
            row_colors = [group_color_map[g] for g in groups]
        
        # 绘制热图
        sns.clustermap(
            data.T,  # 转置：行为代谢物，列为样本
            cmap='viridis',
            figsize=(12, 10),
            row_cluster=True,
            col_cluster=True,
            cbar_kws={"label": "Abundance"},
            row_colors=row_colors if row_colors else None
        )
        
        plt.savefig(output_path, dpi=150, bbox_inches='tight')
        plt.close()
        
        return {
            "status": "success",
            "plot_path": str(output_path)
        }
    
    except Exception as e:
        logger.error(f"❌ 热图生成失败: {e}", exc_info=True)
        return {
            "status": "error",
            "error": str(e)
        }

