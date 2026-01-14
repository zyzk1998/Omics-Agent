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
    name="pca_analysis",
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
    name="differential_analysis",
    description="Performs differential analysis between two groups in metabolite data. Supports t-test and Wilcoxon rank-sum test. Returns p-values, FDR-corrected p-values, and log2 fold changes. Automatically detects groups if not specified.",
    category="Metabolomics",
    output_type="json"
)
def run_differential_analysis(
    file_path: str,
    group_column: str,
    case_group: Optional[str] = None,
    control_group: Optional[str] = None,
    group1: Optional[str] = None,  # 别名，兼容 planner 发送的参数
    group2: Optional[str] = None,   # 别名，兼容 planner 发送的参数
    method: str = "t-test",
    p_value_threshold: float = 0.05,
    fold_change_threshold: float = 1.5,
    fdr_method: str = "fdr_bh",
    is_logged: bool = True,  # 🔥 CRITICAL FIX: 数据是否已Log2转换（默认True，因为SOP强制Log2转换）
    output_dir: Optional[str] = None,
    **kwargs  # 安全网：接受其他意外参数
) -> Dict[str, Any]:
    """
    执行差异代谢物分析
    
    Args:
        file_path: 输入数据文件路径（CSV，包含分组信息）
        group_column: 分组列名
        case_group: 实验组名称（可选，如果为 None 则自动检测）
        control_group: 对照组名称（可选，如果为 None 则自动检测）
        group1: 第一组名称（别名，兼容 planner 参数）
        group2: 第二组名称（别名，兼容 planner 参数）
        method: 统计方法（"t-test" 或 "wilcoxon"，默认 "t-test"）
        p_value_threshold: P 值阈值（默认 0.05）
        fold_change_threshold: 倍数变化阈值（默认 1.5）
        fdr_method: FDR 校正方法（默认 "fdr_bh"）
        is_logged: 数据是否已Log2转换（默认 True，因为SOP强制Log2转换）
        output_dir: 输出目录（如果提供，将保存结果到 CSV）
        **kwargs: 其他参数（安全网）
    
    Returns:
        包含以下键的字典:
        - status: "success" 或 "error"
        - results: 差异分析结果列表（每个代谢物一行）
        - output_path: 保存的结果文件路径（如果保存）
        - output_file: 保存的结果文件路径（别名，用于数据流传递）
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
        
        # 🔥 自动检测分组（如果未指定）
        # 优先使用 case_group/control_group，如果没有则使用 group1/group2
        if case_group is None and group1 is not None:
            case_group = group1
        if control_group is None and group2 is not None:
            control_group = group2
        
        # 如果仍然为 None，自动检测前两个唯一值
        if case_group is None or control_group is None:
            unique_groups = sorted(df[group_column].unique().tolist())  # 排序以确保一致性
            if len(unique_groups) < 2:
                return {
                    "status": "error",
                    "error": f"分组列 '{group_column}' 中只有 {len(unique_groups)} 个唯一值，需要至少 2 个组"
                }
            # 使用前两个唯一值
            case_group = unique_groups[0] if case_group is None else case_group
            control_group = unique_groups[1] if control_group is None else control_group
            logger.info(f"🔄 自动检测分组: case_group={case_group}, control_group={control_group}")
        
        # 🔥 将检测到的分组信息添加到返回结果中，供后续步骤使用
        detected_groups = {
            "case_group": case_group,
            "control_group": control_group,
            "group1": case_group,  # 别名
            "group2": control_group  # 别名
        }
        
        # 分离分组
        groups = df[group_column]
        case_mask = groups == case_group
        control_mask = groups == control_group
        
        if not case_mask.any():
            return {
                "status": "error",
                "error": f"实验组 '{case_group}' 不存在于分组列中。可用组: {list(unique_groups)}"
            }
        if not control_mask.any():
            return {
                "status": "error",
                "error": f"对照组 '{control_group}' 不存在于分组列中。可用组: {list(unique_groups)}"
            }
        
        # 提取代谢物列（数值列，排除分组列）
        metabolite_cols = [
            col for col in df.columns 
            if col != group_column and pd.api.types.is_numeric_dtype(df[col])
        ]
        
        results = []
        p_values = []
        
        # 🔥 根据 method 参数选择统计方法
        use_wilcoxon = method.lower() in ["wilcoxon", "wilcox", "ranksum", "mann-whitney"]
        
        for metabolite in metabolite_cols:
            case_values = df.loc[case_mask, metabolite].dropna()
            control_values = df.loc[control_mask, metabolite].dropna()
            
            if len(case_values) < 2 or len(control_values) < 2:
                continue
            
            # 🔥 根据方法选择统计检验
            if use_wilcoxon:
                # Wilcoxon rank-sum test (Mann-Whitney U test)
                try:
                    u_stat, p_val = stats.ranksums(case_values, control_values)
                except Exception as e:
                    logger.warning(f"⚠️ Wilcoxon 检验失败 ({metabolite}): {e}，跳过")
                    continue
            else:
                # T-test (默认)
                try:
                    t_stat, p_val = stats.ttest_ind(case_values, control_values)
                except Exception as e:
                    logger.warning(f"⚠️ T-test 失败 ({metabolite}): {e}，跳过")
                    continue
            
            # 🔥 CRITICAL FIX: 计算 log2 fold change
            # 对于已Log2转换的数据，使用减法：Log2FC = Mean_A - Mean_B
            # 对于原始数据，使用除法：Log2FC = log2(Mean_A / Mean_B)
            case_mean = case_values.mean()
            control_mean = control_values.mean()
            
            if is_logged:
                # 数据已Log2转换，使用减法
                log2fc = case_mean - control_mean
                logger.debug(f"✅ [Log2FC] 使用减法（数据已Log2转换）: {case_mean} - {control_mean} = {log2fc}")
            else:
                # 原始数据，使用除法
                if control_mean > 0:
                    log2fc = np.log2(case_mean / control_mean)
                else:
                    # 避免除零，使用小的epsilon
                    log2fc = np.log2(case_mean / (control_mean + 1e-9))
                logger.debug(f"✅ [Log2FC] 使用除法（原始数据）: log2({case_mean} / {control_mean}) = {log2fc}")
            
            results.append({
                "metabolite": metabolite,
                "p_value": float(p_val),
                "log2fc": float(log2fc),
                "log2_fold_change": float(log2fc),  # 别名，兼容 visualize_volcano
                "case_mean": float(case_mean),
                "control_mean": float(control_mean),
                "case_group": case_group,
                "control_group": control_group
            })
            p_values.append(p_val)
        
        # FDR 校正
        if p_values:
            _, p_adjusted, _, _ = multipletests(p_values, method=fdr_method)
            
            # 添加 FDR 校正后的 p 值
            for i, result in enumerate(results):
                result["fdr"] = float(p_adjusted[i])
                result["fdr_corrected_pvalue"] = float(p_adjusted[i])  # 别名，兼容 visualize_volcano
                # 🔥 使用用户指定的阈值判断显著性（确保返回 Python 原生 bool）
                result["significant"] = bool(
                    p_adjusted[i] < p_value_threshold and 
                    abs(result["log2fc"]) >= np.log2(fold_change_threshold)
                )
        
        # 🔥 保存结果到 CSV 文件（用于数据流传递和可视化）
        output_path = None
        if output_dir:
            output_path_obj = Path(output_dir)
            output_path_obj.mkdir(parents=True, exist_ok=True)
            
            # 生成输出文件路径
            input_filename = Path(file_path).stem
            output_path = str(output_path_obj / f"{input_filename}_differential_results.csv")
            
            # 转换为 DataFrame 并保存
            results_df = pd.DataFrame(results)
            results_df.to_csv(output_path, index=False)
            logger.info(f"💾 差异分析结果已保存: {output_path}")
        else:
            # 如果没有指定输出目录，尝试使用输入文件所在目录
            input_dir = Path(file_path).parent
            output_path = str(input_dir / "differential_results.csv")
            results_df = pd.DataFrame(results)
            results_df.to_csv(output_path, index=False)
            logger.info(f"💾 差异分析结果已保存: {output_path}")
        
        # 统计摘要
        significant_count = sum(1 for r in results if r.get("significant", False))
        
        return {
            "status": "success",
            "results": results,
            "output_path": output_path,
            "output_file": output_path,  # 别名，用于数据流传递
            "file_path": output_path,    # 另一个别名，确保兼容性
            "summary": {
                "total_metabolites": len(results),
                "significant_count": significant_count,
                "method": method,
                "case_group": case_group,
                "control_group": control_group,
                "p_value_threshold": p_value_threshold,
                "fold_change_threshold": fold_change_threshold
            }
        }
    
    except Exception as e:
        logger.error(f"❌ 差异分析失败: {e}", exc_info=True)
        return {
            "status": "error",
            "error": str(e)
        }

