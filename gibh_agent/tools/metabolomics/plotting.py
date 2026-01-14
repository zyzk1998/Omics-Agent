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
    name="visualize_volcano",
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
        
        # 转换为 DataFrame
        df = pd.DataFrame(results)
        
        # 提取必要的列
        if "log2fc" not in df.columns and "log2_fold_change" in df.columns:
            df["log2fc"] = df["log2_fold_change"]
        if "p_value" not in df.columns and "pvalue" in df.columns:
            df["p_value"] = df["pvalue"]
        if "fdr" not in df.columns and "fdr_corrected_pvalue" in df.columns:
            df["fdr"] = df["fdr_corrected_pvalue"]
        
        # 计算 -log10(p-value)
        df["neg_log10_p"] = -np.log10(df["p_value"] + 1e-300)  # 避免 log(0)
        
        # 生成图片
        if output_path is None:
            output_path = "./results/volcano_plot.png"
        
        output_path_obj = Path(output_path)
        output_path_obj.parent.mkdir(parents=True, exist_ok=True)
        
        plt.figure(figsize=(10, 8))
        
        # 分类点
        df["significance"] = "Not Significant"
        df.loc[(df["fdr"] < fdr_threshold) & (df["log2fc"].abs() < log2fc_threshold), "significance"] = "FDR Significant"
        df.loc[(df["fdr"] >= fdr_threshold) & (df["log2fc"].abs() >= log2fc_threshold), "significance"] = "FC Significant"
        df.loc[(df["fdr"] < fdr_threshold) & (df["log2fc"].abs() >= log2fc_threshold), "significance"] = "Both Significant"
        
        # 绘制散点图
        colors = {"Not Significant": "gray", "FDR Significant": "blue", "FC Significant": "orange", "Both Significant": "red"}
        for sig_type, color in colors.items():
            subset = df[df["significance"] == sig_type]
            plt.scatter(subset["log2fc"], subset["neg_log10_p"], c=color, alpha=0.6, label=sig_type, s=30)
        
        # 添加阈值线
        plt.axhline(y=-np.log10(fdr_threshold), color='r', linestyle='--', alpha=0.5, label=f'FDR = {fdr_threshold}')
        plt.axvline(x=log2fc_threshold, color='r', linestyle='--', alpha=0.5, label=f'Log2FC = {log2fc_threshold}')
        plt.axvline(x=-log2fc_threshold, color='r', linestyle='--', alpha=0.5)
        
        plt.xlabel("Log2 Fold Change")
        plt.ylabel("-Log10 P-value")
        plt.title("Volcano Plot")
        plt.legend()
        plt.grid(True, alpha=0.3)
        plt.tight_layout()
        plt.savefig(output_path_obj, dpi=150, bbox_inches='tight')
        plt.close()
        
        return {
            "status": "success",
            "plot_path": str(output_path_obj),
            "n_significant": len(df[df["significance"] == "Both Significant"])
        }
    
    except Exception as e:
        logger.error(f"❌ 火山图生成失败: {e}", exc_info=True)
        return {
            "status": "error",
            "error": str(e)
        }


@registry.register(
    name="visualize_pca",
    description="Generates an enhanced PCA visualization plot. Can use existing PCA results or re-run PCA with grouping information.",
    category="Metabolomics",
    output_type="image"
)
def plot_pca(
    pca_results: Optional[Dict[str, Any]] = None,
    file_path: Optional[str] = None,
    group_column: Optional[str] = None,
    output_path: Optional[str] = None,
    pc1: int = 1,
    pc2: int = 2
) -> Dict[str, Any]:
    """
    生成 PCA 可视化图
    
    Args:
        pca_results: PCA 分析结果（包含 pca_coordinates 和 explained_variance）
        file_path: 如果 pca_results 为 None，使用此文件路径重新运行 PCA
        group_column: 分组列名（用于着色）
        output_path: 输出文件路径
        pc1: 第一个主成分（1-indexed）
        pc2: 第二个主成分（1-indexed）
    
    Returns:
        包含图片路径的字典
    """
    try:
        # 如果提供了 pca_results，直接使用
        if pca_results and "pca_coordinates" in pca_results:
            coords_dict = pca_results["pca_coordinates"]
            coords_df = pd.DataFrame(coords_dict).T
            
            # 获取解释方差
            explained_var = pca_results.get("explained_variance", {})
            pc1_var = explained_var.get(f"PC{pc1}", 0) * 100
            pc2_var = explained_var.get(f"PC{pc2}", 0) * 100
            
            # 如果有 plot_path，直接返回
            if pca_results.get("plot_path"):
                return {
                    "status": "success",
                    "plot_path": pca_results["plot_path"],
                    "message": "使用已有的 PCA 图"
                }
        else:
            # 需要重新运行 PCA（简化版本）
            if not file_path:
                return {
                    "status": "error",
                    "error": "需要提供 pca_results 或 file_path"
                }
            
            # 这里可以调用 pca_analysis，但为了避免循环依赖，我们简化处理
            return {
                "status": "error",
                "error": "请先运行 pca_analysis，然后使用其返回的 plot_path"
            }
        
        # 生成图片
        if output_path is None:
            output_path = "./results/pca_plot.png"
        
        output_path_obj = Path(output_path)
        output_path_obj.parent.mkdir(parents=True, exist_ok=True)
        
        plt.figure(figsize=(10, 8))
        
        # 如果有分组信息，按组着色
        if group_column and group_column in coords_df.index:
            groups = coords_df[group_column]
            unique_groups = groups.unique()
            colors = plt.cm.Set3(np.linspace(0, 1, len(unique_groups)))
            for i, group in enumerate(unique_groups):
                mask = groups == group
                plt.scatter(
                    coords_df.loc[mask, f"PC{pc1}"],
                    coords_df.loc[mask, f"PC{pc2}"],
                    c=[colors[i]],
                    label=group,
                    alpha=0.6,
                    s=50
                )
            plt.legend()
        else:
            plt.scatter(coords_df[f"PC{pc1}"], coords_df[f"PC{pc2}"], alpha=0.6, s=50)
        
        plt.xlabel(f"PC{pc1} ({pc1_var:.2f}%)")
        plt.ylabel(f"PC{pc2} ({pc2_var:.2f}%)")
        plt.title("PCA Plot")
        plt.grid(True, alpha=0.3)
        plt.tight_layout()
        plt.savefig(output_path_obj, dpi=150, bbox_inches='tight')
        plt.close()
        
        return {
            "status": "success",
            "plot_path": str(output_path_obj)
        }
    
    except Exception as e:
        logger.error(f"❌ PCA 可视化失败: {e}", exc_info=True)
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
            data.T,
            cmap='viridis',
            row_cluster=True,
            col_cluster=True,
            figsize=(12, 10),
            row_colors=row_colors if row_colors else None
        )
        plt.savefig(output_path_obj, dpi=150, bbox_inches='tight')
        plt.close()
        
        return {
            "status": "success",
            "plot_path": str(output_path_obj)
        }
    
    except Exception as e:
        logger.error(f"❌ 热图生成失败: {e}", exc_info=True)
        return {
            "status": "error",
            "error": str(e)
        }
