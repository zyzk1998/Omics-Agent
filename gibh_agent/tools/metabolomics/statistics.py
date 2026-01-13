"""
代谢组学统计分析工具
"""
import logging
from typing import Dict, Any, Optional
from pathlib import Path
import pandas as pd
import numpy as np
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy import stats
from statsmodels.stats.multitest import multipletests

from ...core.tool_registry import registry

logger = logging.getLogger(__name__)


@registry.register(
    name="metabolomics_pca",
    description="Performs Principal Component Analysis (PCA) on metabolite abundance data. Returns PCA coordinates, explained variance, and optionally a PCA plot.",
    category="Metabolomics",
    output_type="mixed"  # 返回 JSON + 图片路径
)
def run_pca(
    file_path: str,
    n_components: int = 2,
    scale: bool = True,
    output_dir: Optional[str] = None
) -> Dict[str, Any]:
    """
    执行 PCA 分析
    
    Args:
        file_path: 输入数据文件路径（CSV）
        n_components: 主成分数量（默认 2）
        scale: 是否标准化数据（默认 True）
        output_dir: 输出目录（可选）
    
    Returns:
        包含以下键的字典:
        - status: "success" 或 "error"
        - pca_coordinates: PCA 坐标 (DataFrame 的 JSON 表示)
        - explained_variance: 解释方差比例
        - plot_path: PCA 图路径（如果生成）
        - error: 错误信息（如果失败）
    """
    try:
        # 读取数据
        df = pd.read_csv(file_path, index_col=0)
        
        # 提取数值列（排除非数值列）
        numeric_cols = df.select_dtypes(include=[np.number]).columns
        data = df[numeric_cols]
        
        # 数据预处理
        if scale:
            scaler = StandardScaler()
            data_scaled = scaler.fit_transform(data)
        else:
            data_scaled = data.values
        
        # 执行 PCA
        pca = PCA(n_components=n_components)
        pca_coords = pca.fit_transform(data_scaled)
        
        # 创建结果 DataFrame
        coords_df = pd.DataFrame(
            pca_coords,
            index=data.index,
            columns=[f"PC{i+1}" for i in range(n_components)]
        )
        
        # 生成图片（如果指定了输出目录）
        plot_path = None
        if output_dir:
            output_path = Path(output_dir)
            output_path.mkdir(parents=True, exist_ok=True)
            plot_path = str(output_path / "pca_plot.png")
            
            plt.figure(figsize=(10, 8))
            plt.scatter(coords_df.iloc[:, 0], coords_df.iloc[:, 1], alpha=0.6)
            plt.xlabel(f"PC1 ({pca.explained_variance_ratio_[0]:.2%})")
            plt.ylabel(f"PC2 ({pca.explained_variance_ratio_[1]:.2%})")
            plt.title("PCA Plot")
            plt.grid(True, alpha=0.3)
            plt.savefig(plot_path, dpi=150, bbox_inches='tight')
            plt.close()
        
        return {
            "status": "success",
            "pca_coordinates": coords_df.to_dict(orient='index'),
            "explained_variance": {
                f"PC{i+1}": float(ratio) 
                for i, ratio in enumerate(pca.explained_variance_ratio_)
            },
            "plot_path": plot_path
        }
    
    except Exception as e:
        logger.error(f"❌ PCA 分析失败: {e}", exc_info=True)
        return {
            "status": "error",
            "error": str(e)
        }


@registry.register(
    name="metabolomics_differential_analysis",
    description="Performs differential analysis between two groups in metabolite data. Uses t-test to identify significantly different metabolites. Returns p-values, FDR-corrected p-values, and log2 fold changes.",
    category="Metabolomics",
    output_type="json"
)
def run_differential_analysis(
    file_path: str,
    group_column: str,
    case_group: str,
    control_group: str,
    fdr_method: str = "fdr_bh"
) -> Dict[str, Any]:
    """
    执行差异代谢物分析
    
    Args:
        file_path: 输入数据文件路径（CSV，包含分组信息）
        group_column: 分组列名
        case_group: 实验组名称
        control_group: 对照组名称
        fdr_method: FDR 校正方法（默认 "fdr_bh"）
    
    Returns:
        包含以下键的字典:
        - status: "success" 或 "error"
        - results: 差异分析结果列表（每个代谢物一行）
        - error: 错误信息（如果失败）
    """
    try:
        # 读取数据
        df = pd.read_csv(file_path, index_col=0)
        
        # 检查分组列是否存在
        if group_column not in df.columns:
            return {
                "status": "error",
                "error": f"分组列 '{group_column}' 不存在于数据中"
            }
        
        # 分离分组
        groups = df[group_column]
        case_mask = groups == case_group
        control_mask = groups == control_group
        
        if not case_mask.any():
            return {
                "status": "error",
                "error": f"实验组 '{case_group}' 不存在"
            }
        if not control_mask.any():
            return {
                "status": "error",
                "error": f"对照组 '{control_group}' 不存在"
            }
        
        # 提取代谢物列（数值列，排除分组列）
        metabolite_cols = [
            col for col in df.columns 
            if col != group_column and pd.api.types.is_numeric_dtype(df[col])
        ]
        
        results = []
        p_values = []
        
        for metabolite in metabolite_cols:
            case_values = df.loc[case_mask, metabolite].dropna()
            control_values = df.loc[control_mask, metabolite].dropna()
            
            if len(case_values) < 2 or len(control_values) < 2:
                continue
            
            # T-test
            t_stat, p_val = stats.ttest_ind(case_values, control_values)
            
            # 计算 log2 fold change
            case_mean = case_values.mean()
            control_mean = control_values.mean()
            
            if control_mean > 0:
                log2fc = np.log2(case_mean / control_mean)
            else:
                log2fc = 0.0
            
            results.append({
                "metabolite": metabolite,
                "p_value": float(p_val),
                "log2fc": float(log2fc),
                "case_mean": float(case_mean),
                "control_mean": float(control_mean)
            })
            p_values.append(p_val)
        
        # FDR 校正
        if p_values:
            _, p_adjusted, _, _ = multipletests(p_values, method=fdr_method)
            
            # 添加 FDR 校正后的 p 值
            for i, result in enumerate(results):
                result["fdr"] = float(p_adjusted[i])
                result["significant"] = p_adjusted[i] < 0.05
        
        return {
            "status": "success",
            "results": results,
            "summary": {
                "total_metabolites": len(results),
                "significant_count": sum(1 for r in results if r.get("significant", False))
            }
        }
    
    except Exception as e:
        logger.error(f"❌ 差异分析失败: {e}", exc_info=True)
        return {
            "status": "error",
            "error": str(e)
        }

