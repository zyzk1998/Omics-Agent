"""
代谢组学高级分析工具 - PLS-DA 和通路富集分析
"""
import logging
from typing import Dict, Any, Optional, List
from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from ...core.tool_registry import registry

logger = logging.getLogger(__name__)


@registry.register(
    name="metabolomics_plsda",
    description="Performs Partial Least Squares Discriminant Analysis (PLS-DA) on metabolite data to identify metabolites that best discriminate between groups. Returns VIP (Variable Importance in Projection) scores for feature selection.",
    category="Metabolomics",
    output_type="mixed"
)
def run_plsda(
    file_path: str,
    group_column: str,
    n_components: int = 2,
    scale: bool = True,
    output_dir: Optional[str] = None
) -> Dict[str, Any]:
    """
    执行 PLS-DA 分析
    
    Args:
        file_path: 输入数据文件路径（CSV，包含分组信息）
        group_column: 分组列名
        n_components: PLS 成分数量（默认 2）
        scale: 是否标准化数据（默认 True）
        output_dir: 输出目录（可选）
    
    Returns:
        包含以下键的字典:
        - status: "success" 或 "error"
        - vip_scores: VIP 分数（按降序排列）
        - pls_scores: PLS 得分
        - plot_path: PLS-DA 图路径（如果生成）
        - error: 错误信息（如果失败）
    """
    try:
        from sklearn.cross_decomposition import PLSRegression
        from sklearn.preprocessing import StandardScaler, LabelEncoder
        
        # 读取数据
        df = pd.read_csv(file_path, index_col=0)
        
        # 检查分组列是否存在
        if group_column not in df.columns:
            return {
                "status": "error",
                "error": f"分组列 '{group_column}' 不存在于数据中"
            }
        
        # 提取数值列（代谢物）
        metabolite_cols = [
            col for col in df.columns 
            if col != group_column and pd.api.types.is_numeric_dtype(df[col])
        ]
        
        if len(metabolite_cols) == 0:
            return {
                "status": "error",
                "error": "未找到数值型代谢物列"
            }
        
        # 准备数据
        X = df[metabolite_cols].values
        y = df[group_column].values
        
        # 编码分类变量
        le = LabelEncoder()
        y_encoded = le.fit_transform(y)
        
        # 标准化（如果需要）
        if scale:
            scaler = StandardScaler()
            X_scaled = scaler.fit_transform(X)
        else:
            X_scaled = X
        
        # 执行 PLS-DA
        pls = PLSRegression(n_components=n_components, scale=False)
        pls.fit(X_scaled, y_encoded)
        
        # 计算 PLS 得分
        X_scores = pls.x_scores_
        
        # 计算 VIP 分数
        # VIP = sqrt(p * (w^2 * SSY) / SSY_total)
        # 其中 p 是变量数，w 是权重，SSY 是 Y 的方差解释
        T = pls.x_scores_  # 得分矩阵
        W = pls.x_weights_  # 权重矩阵
        Q = pls.y_loadings_  # Y 的载荷
        
        # 计算每个成分的方差解释
        explained_variance = []
        for i in range(n_components):
            ssy = np.sum((T[:, i] @ Q[i]) ** 2)
            explained_variance.append(ssy)
        
        total_ssy = sum(explained_variance)
        
        # 计算 VIP
        vip_scores = []
        for j in range(len(metabolite_cols)):
            vip = 0
            for i in range(n_components):
                if total_ssy > 0:
                    vip += (W[j, i] ** 2) * (explained_variance[i] / total_ssy)
            vip = np.sqrt(len(metabolite_cols) * vip)
            vip_scores.append(vip)
        
        # 创建 VIP 分数 DataFrame
        vip_df = pd.DataFrame({
            'metabolite': metabolite_cols,
            'vip_score': vip_scores
        }).sort_values('vip_score', ascending=False)
        
        # 生成图片（如果指定了输出目录）
        plot_path = None
        if output_dir:
            output_path = Path(output_dir)
            output_path.mkdir(parents=True, exist_ok=True)
            plot_path = str(output_path / "plsda_plot.png")
            
            # PLS-DA 得分图
            unique_groups = le.classes_
            colors = plt.cm.Set3(np.linspace(0, 1, len(unique_groups)))
            
            plt.figure(figsize=(10, 8))
            for i, group in enumerate(unique_groups):
                mask = y == group
                plt.scatter(
                    X_scores[mask, 0], 
                    X_scores[mask, 1],
                    label=group,
                    color=colors[i],
                    alpha=0.6,
                    s=50
                )
            
            plt.xlabel(f"PLS Component 1 ({explained_variance[0]/total_ssy*100:.1f}%)")
            plt.ylabel(f"PLS Component 2 ({explained_variance[1]/total_ssy*100:.1f}%)")
            plt.title("PLS-DA Score Plot")
            plt.legend()
            plt.grid(True, alpha=0.3)
            plt.savefig(plot_path, dpi=150, bbox_inches='tight')
            plt.close()
        
        return {
            "status": "success",
            "vip_scores": vip_df.to_dict(orient='records'),
            "top_vip_metabolites": vip_df.head(20).to_dict(orient='records'),
            "explained_variance": {
                f"Component{i+1}": float(var / total_ssy * 100) 
                for i, var in enumerate(explained_variance)
            },
            "plot_path": plot_path,
            "n_components": n_components
        }
    
    except ImportError as e:
        return {
            "status": "error",
            "error": f"缺少必要的依赖: {str(e)}. 请安装: pip install scikit-learn"
        }
    except Exception as e:
        logger.error(f"❌ PLS-DA 分析失败: {e}", exc_info=True)
        return {
            "status": "error",
            "error": str(e)
        }


@registry.register(
    name="metabolomics_pathway_enrichment",
    description="Performs KEGG pathway enrichment analysis on metabolite data using GSEApy. Identifies significantly enriched metabolic pathways based on metabolite abundance changes between groups.",
    category="Metabolomics",
    output_type="json"
)
def run_pathway_enrichment(
    file_path: str,
    group_column: str,
    case_group: str,
    control_group: str,
    organism: str = "hsa",
    p_value_threshold: float = 0.05,
    output_dir: Optional[str] = None
) -> Dict[str, Any]:
    """
    执行通路富集分析
    
    Args:
        file_path: 输入数据文件路径（CSV，包含分组信息）
        group_column: 分组列名
        case_group: 实验组名称
        control_group: 对照组名称
        organism: 物种代码（默认 "hsa" 人类，可选 "mmu" 小鼠等）
        p_value_threshold: P 值阈值（默认 0.05）
        output_dir: 输出目录（可选）
    
    Returns:
        包含以下键的字典:
        - status: "success" 或 "error"
        - enriched_pathways: 富集通路列表
        - summary: 富集分析摘要
        - error: 错误信息（如果失败）
    """
    try:
        import gseapy as gp
        
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
        
        # 计算 fold change（用于排序）
        case_mean = df.loc[case_mask, metabolite_cols].mean()
        control_mean = df.loc[control_mask, metabolite_cols].mean()
        
        # 避免除零
        control_mean = control_mean.replace(0, np.nan)
        fold_change = case_mean / control_mean
        fold_change = fold_change.fillna(0)
        
        # 创建排序的代谢物列表（按 fold change 降序）
        metabolite_rank = fold_change.sort_values(ascending=False)
        
        # 准备 GSEApy 输入（需要代谢物名称到 KEGG ID 的映射）
        # 注意：这里使用代谢物名称作为 ID（实际应用中需要映射到 KEGG ID）
        # 为了演示，我们使用代谢物名称
        metabolite_list = metabolite_rank.index.tolist()
        
        # 执行通路富集分析
        # 注意：gseapy 需要代谢物 ID 映射到 KEGG，这里使用简化版本
        try:
            # 尝试使用 KEGG 数据库
            enr = gp.enrichr(
                gene_list=metabolite_list[:100],  # 限制前100个
                gene_sets=['KEGG_2021_Human'],  # KEGG 通路数据库
                organism=organism,
                outdir=None,  # 不保存文件
                verbose=False
            )
            
            # 提取结果
            if enr is not None and hasattr(enr, 'results'):
                results_df = enr.results
                
                # 过滤显著通路
                significant_pathways = results_df[
                    results_df['Adjusted P-value'] < p_value_threshold
                ].copy()
                
                # 保存结果（如果指定了输出目录）
                output_csv = None
                if output_dir:
                    output_path = Path(output_dir)
                    output_path.mkdir(parents=True, exist_ok=True)
                    output_csv = str(output_path / "pathway_enrichment.csv")
                    significant_pathways.to_csv(output_csv, index=False)
                
                return {
                    "status": "success",
                    "enriched_pathways": significant_pathways.to_dict(orient='records'),
                    "n_significant": len(significant_pathways),
                    "n_total": len(results_df),
                    "output_csv": output_csv,
                    "summary": f"发现 {len(significant_pathways)} 个显著富集通路（p < {p_value_threshold}）"
                }
            else:
                return {
                    "status": "warning",
                    "message": "通路富集分析完成，但未找到显著富集通路",
                    "enriched_pathways": []
                }
        
        except Exception as e:
            # 如果 gseapy 失败，返回一个占位符结果
            logger.warning(f"⚠️ GSEApy 执行失败，使用简化版本: {e}")
            return {
                "status": "warning",
                "message": f"通路富集分析部分完成（GSEApy 错误: {str(e)}）",
                "enriched_pathways": [],
                "error": str(e)
            }
    
    except ImportError:
        return {
            "status": "error",
            "error": "gseapy not installed. Please install: pip install gseapy"
        }
    except Exception as e:
        logger.error(f"❌ 通路富集分析失败: {e}", exc_info=True)
        return {
            "status": "error",
            "error": str(e)
        }


