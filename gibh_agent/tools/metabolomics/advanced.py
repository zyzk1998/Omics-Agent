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
    output_dir: Optional[str] = None,
    **kwargs  # 🔥 CRITICAL FIX: 安全网，接受其他意外参数
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
        
        # 🔥 修复：检查分组列是否存在，如果不存在则尝试模糊匹配
        if group_column not in df.columns:
            # 尝试模糊匹配：忽略大小写、空格、下划线、连字符
            group_column_normalized = group_column.lower().replace(' ', '').replace('_', '').replace('-', '')
            matched_column = None
            
            for col in df.columns:
                col_normalized = col.lower().replace(' ', '').replace('_', '').replace('-', '')
                if col_normalized == group_column_normalized:
                    matched_column = col
                    logger.info(f"🔄 [PLS-DA] 模糊匹配分组列: '{group_column}' -> '{col}'")
                    break
            
            if matched_column:
                group_column = matched_column
            else:
                # 🔥 改进：区分显示元数据列（可能的分组列）和特征列（代谢物）
                metadata_cols = [col for col in df.columns if not pd.api.types.is_numeric_dtype(df[col])]
                numeric_cols = [col for col in df.columns if pd.api.types.is_numeric_dtype(df[col])]
                
                # 检查是否有唯一值较少的列（可能是分组列）
                potential_group_cols = []
                for col in df.columns:
                    unique_count = df[col].nunique()
                    if 2 <= unique_count <= 10:
                        potential_group_cols.append(f"{col} ({unique_count}个唯一值)")
                
                error_msg = f"分组列 '{group_column}' 不存在于数据中。\n\n"
                
                # 检查索引列是否可能包含分组信息
                if df.index.nunique() <= 10:
                    index_values = df.index.unique().tolist()[:10]
                    error_msg += f"⚠️ 索引列（第一列）有 {df.index.nunique()} 个唯一值，可能包含分组信息: {index_values}\n"
                    error_msg += "💡 提示：如果分组信息在索引列中，请将索引列转换为数据列。\n\n"
                
                if metadata_cols:
                    error_msg += f"可能的元数据列（非数值列）: {', '.join(metadata_cols[:10])}\n"
                else:
                    error_msg += "❌ 未找到非数值列（所有列都是数值型）。\n"
                    error_msg += "💡 数据格式要求：CSV文件应包含一列分组信息（如 'Group', 'Condition', 'Treatment' 等）。\n"
                    error_msg += "   数据格式示例：\n"
                    error_msg += "   SampleID,Group,Metabolite1,Metabolite2,...\n"
                    error_msg += "   Sample1,Control,1.2,3.4,...\n"
                    error_msg += "   Sample2,Treatment,2.3,4.5,...\n\n"
                
                if potential_group_cols:
                    error_msg += f"可能的分组列（唯一值2-10）: {', '.join(potential_group_cols[:10])}\n"
                
                if numeric_cols:
                    error_msg += f"特征列（代谢物，前5个）: {', '.join(numeric_cols[:5])}, ..."
                
                return {
                    "status": "error",
                    "error": error_msg.strip()
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
        
        # 🔥 CRITICAL FIX: 再次检查 NaN（防止在提取后出现）
        if np.isnan(X).any():
            logger.warning(f"⚠️ [PLS-DA] 提取后的数据仍包含 NaN，使用中位数填充...")
            from sklearn.impute import SimpleImputer
            imputer = SimpleImputer(strategy='median')
            X = imputer.fit_transform(X)
        
        # 编码分类变量
        le = LabelEncoder()
        y_encoded = le.fit_transform(y)
        
        # 标准化（如果需要）
        if scale:
            scaler = StandardScaler()
            X_scaled = scaler.fit_transform(X)
        else:
            X_scaled = X
        
        # 🔥 CRITICAL FIX: 最终检查 NaN
        if np.isnan(X_scaled).any():
            logger.error(f"❌ [PLS-DA] 标准化后的数据仍包含 NaN")
            return {
                "status": "error",
                "error": "数据预处理失败：标准化后的数据包含 NaN 值。请检查输入数据质量。"
            }
        
        # 执行 PLS-DA
        pls = PLSRegression(n_components=n_components, scale=False)
        pls.fit(X_scaled, y_encoded)
        
        # 计算 PLS 得分
        X_scores = pls.x_scores_
        
        # 计算 VIP 分数
        # VIP = sqrt(p * (w^2 * SSY) / SSY_total)
        # 其中 p 是变量数，w 是权重，SSY 是 Y 的方差解释
        T = pls.x_scores_  # 得分矩阵 (n_samples, n_components)
        W = pls.x_weights_  # 权重矩阵 (n_features, n_components)
        Q = pls.y_loadings_  # Y 的载荷 (n_components,) for single output
        
        # 🔥 CRITICAL FIX: 处理 y_loadings_ 的形状
        # 对于单输出 PLS，Q 是 1D 数组 (n_components,)
        # 对于多输出 PLS，Q 是 2D 数组 (n_components, n_targets)
        if Q.ndim == 1:
            # 单输出：Q[i] 是标量
            Q_flat = Q
        else:
            # 多输出：取第一列（或平均）
            Q_flat = Q[:, 0] if Q.shape[1] > 0 else Q.flatten()
        
        # 🔥 CRITICAL FIX: 确保 Q_flat 的长度足够
        # 如果 Q_flat 的长度小于 n_components，使用实际长度
        actual_n_components = min(n_components, len(Q_flat), T.shape[1])
        logger.debug(f"🔍 [PLS-DA] n_components={n_components}, Q_flat.shape={Q_flat.shape}, T.shape={T.shape}, actual_n_components={actual_n_components}")
        
        # 计算每个成分的方差解释
        explained_variance = []
        for i in range(actual_n_components):
            # 🔥 CRITICAL FIX: Q_flat[i] 是标量，使用元素乘法
            # SSY = sum((T[:, i] * Q_flat[i])^2) = Q_flat[i]^2 * sum(T[:, i]^2)
            if i < len(Q_flat) and i < T.shape[1]:
                ssy = Q_flat[i] ** 2 * np.sum(T[:, i] ** 2)
                explained_variance.append(ssy)
            else:
                # 如果索引越界，使用默认值
                logger.warning(f"⚠️ [PLS-DA] 索引 {i} 越界，跳过该成分")
                explained_variance.append(0.0)
        
        total_ssy = sum(explained_variance)
        
        # 计算 VIP
        vip_scores = []
        for j in range(len(metabolite_cols)):
            vip = 0
            for i in range(actual_n_components):
                if total_ssy > 0 and i < len(explained_variance) and i < W.shape[1]:
                    vip += (W[j, i] ** 2) * (explained_variance[i] / total_ssy)
            vip = np.sqrt(len(metabolite_cols) * vip) if vip > 0 else 0.0
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
            
            # 🔥 CRITICAL FIX: 处理不同数量的成分
            if actual_n_components >= 2 and X_scores.shape[1] >= 2:
                # 2D 散点图（Component 1 vs Component 2）
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
                
                comp1_var = explained_variance[0]/total_ssy*100 if len(explained_variance) > 0 and total_ssy > 0 else 0.0
                comp2_var = explained_variance[1]/total_ssy*100 if len(explained_variance) > 1 and total_ssy > 0 else 0.0
                plt.xlabel(f"PLS Component 1 ({comp1_var:.1f}%)")
                plt.ylabel(f"PLS Component 2 ({comp2_var:.1f}%)")
            elif actual_n_components == 1 and X_scores.shape[1] >= 1:
                # 1D 散点图（Component 1 only）
                for i, group in enumerate(unique_groups):
                    mask = y == group
                    plt.scatter(
                        X_scores[mask, 0], 
                        np.zeros(np.sum(mask)),
                        label=group,
                        color=colors[i],
                        alpha=0.6,
                        s=50
                    )
                
                comp1_var = explained_variance[0]/total_ssy*100 if len(explained_variance) > 0 and total_ssy > 0 else 0.0
                plt.xlabel(f"PLS Component 1 ({comp1_var:.1f}%)")
                plt.ylabel("")
            else:
                # 降级：使用前两个特征（如果可用）
                logger.warning(f"⚠️ [PLS-DA] 无法绘制标准 PLS-DA 图，使用降级方案")
                if X_scores.shape[1] >= 2:
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
                    plt.xlabel("PLS Component 1")
                    plt.ylabel("PLS Component 2")
                elif X_scores.shape[1] >= 1:
                    for i, group in enumerate(unique_groups):
                        mask = y == group
                        plt.scatter(
                            X_scores[mask, 0], 
                            np.zeros(np.sum(mask)),
                            label=group,
                            color=colors[i],
                            alpha=0.6,
                            s=50
                        )
                    plt.xlabel("PLS Component 1")
                    plt.ylabel("")
            
            plt.title("PLS-DA Score Plot")
            plt.legend()
            plt.grid(True, alpha=0.3)
            plt.savefig(plot_path, dpi=150, bbox_inches='tight')
            plt.close()
        
        # 🔥 Phase 2: Extract top VIP markers for AI report
        top_vip_markers = vip_df.head(10).to_dict(orient='records')
        top_vip_markers_formatted = [
            {
                "name": m.get("metabolite", "Unknown"),
                "vip": float(m.get("vip_score", 0.0))
            }
            for m in top_vip_markers
        ]
        
        return {
            "status": "success",
            "vip_scores": vip_df.to_dict(orient='records'),
            "top_vip_metabolites": vip_df.head(20).to_dict(orient='records'),
            "explained_variance": {
                f"Component{i+1}": float(var / total_ssy * 100) 
                for i, var in enumerate(explained_variance)
            },
            "plot_path": plot_path,
            "n_components": n_components,
            "summary": {
                "top_vip_markers": top_vip_markers_formatted,
                "n_components": n_components,
                "total_metabolites": len(metabolite_cols),
                "comp1_variance": float(explained_variance[0] / total_ssy * 100) if len(explained_variance) > 0 and total_ssy > 0 else 0.0,
                "comp2_variance": float(explained_variance[1] / total_ssy * 100) if len(explained_variance) > 1 and total_ssy > 0 else 0.0
            }
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
    case_group: Optional[str] = None,
    control_group: Optional[str] = None,
    organism: str = "hsa",
    p_value_threshold: float = 0.05,
    output_dir: Optional[str] = None,
    **kwargs  # 🔥 CRITICAL FIX: 安全网，接受其他意外参数
) -> Dict[str, Any]:
    """
    执行通路富集分析
    
    Args:
        file_path: 输入数据文件路径（CSV，包含分组信息）
        group_column: 分组列名
        case_group: 实验组名称（可选，如果为 None 或占位符，将自动检测）
        control_group: 对照组名称（可选，如果为 None 或占位符，将自动检测）
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
        
        # 🔥 CRITICAL FIX: Map organism codes to gseapy format
        # gseapy expects full names like "human", "mouse", not codes like "hsa", "mmu"
        organism_mapping = {
            "hsa": "human",
            "mmu": "mouse",
            "rno": "rat",
            "bta": "cow",  # 牛
            "gga": "chicken",  # 鸡
            "dre": "zebrafish",  # 斑马鱼
            "cel": "worm",  # 线虫
            "dme": "fly",  # 果蝇
        }
        
        # Convert organism code to full name if needed
        gseapy_organism = organism_mapping.get(organism.lower(), organism.lower())
        
        # Check if organism is supported
        supported_organisms = ["human", "mouse", "rat"]
        if gseapy_organism not in supported_organisms:
            logger.warning(f"⚠️ [Pathway Enrichment] Organism '{organism}' (mapped to '{gseapy_organism}') is not supported by gseapy database")
            return {
                "status": "warning",
                "message": f"Organism '{organism}' is not supported by the pathway enrichment database. Step skipped.",
                "enriched_pathways": [],
                "skipped_reason": f"Unsupported organism: {organism}"
            }
        
        # 读取数据
        df = pd.read_csv(file_path, index_col=0)
        
        # 🔥 CRITICAL FIX: 自动检测分组值（如果 case_group 或 control_group 是 None 或占位符）
        if not case_group or not control_group or case_group.startswith("<") or control_group.startswith("<"):
            logger.info(f"🔄 [Pathway Enrichment] 自动检测分组值...")
            if group_column not in df.columns:
                return {
                    "status": "error",
                    "error": f"分组列 '{group_column}' 不存在于数据中"
                }
            
            unique_groups = sorted(df[group_column].unique().tolist())
            if len(unique_groups) < 2:
                return {
                    "status": "error",
                    "error": f"分组列 '{group_column}' 中只有 {len(unique_groups)} 个唯一值，需要至少 2 个组"
                }
            
            # 使用前两个唯一值
            if not case_group or case_group.startswith("<"):
                case_group = unique_groups[0]
            if not control_group or control_group.startswith("<"):
                control_group = unique_groups[1]
            
            logger.info(f"✅ [Pathway Enrichment] 自动检测分组: case_group={case_group}, control_group={control_group}")
        
        # 🔥 修复：检查分组列是否存在，如果不存在则尝试模糊匹配
        if group_column not in df.columns:
            # 尝试模糊匹配：忽略大小写、空格、下划线、连字符
            group_column_normalized = group_column.lower().replace(' ', '').replace('_', '').replace('-', '')
            matched_column = None
            
            for col in df.columns:
                col_normalized = col.lower().replace(' ', '').replace('_', '').replace('-', '')
                if col_normalized == group_column_normalized:
                    matched_column = col
                    logger.info(f"🔄 [Pathway Enrichment] 模糊匹配分组列: '{group_column}' -> '{col}'")
                    break
            
            if matched_column:
                group_column = matched_column
            else:
                # 🔥 改进：区分显示元数据列（可能的分组列）和特征列（代谢物）
                metadata_cols = [col for col in df.columns if not pd.api.types.is_numeric_dtype(df[col])]
                numeric_cols = [col for col in df.columns if pd.api.types.is_numeric_dtype(df[col])]
                
                # 检查是否有唯一值较少的列（可能是分组列）
                potential_group_cols = []
                for col in df.columns:
                    unique_count = df[col].nunique()
                    if 2 <= unique_count <= 10:
                        potential_group_cols.append(f"{col} ({unique_count}个唯一值)")
                
                error_msg = f"分组列 '{group_column}' 不存在于数据中。\n\n"
                
                # 检查索引列是否可能包含分组信息
                if df.index.nunique() <= 10:
                    index_values = df.index.unique().tolist()[:10]
                    error_msg += f"⚠️ 索引列（第一列）有 {df.index.nunique()} 个唯一值，可能包含分组信息: {index_values}\n"
                    error_msg += "💡 提示：如果分组信息在索引列中，请将索引列转换为数据列。\n\n"
                
                if metadata_cols:
                    error_msg += f"可能的元数据列（非数值列）: {', '.join(metadata_cols[:10])}\n"
                else:
                    error_msg += "❌ 未找到非数值列（所有列都是数值型）。\n"
                    error_msg += "💡 数据格式要求：CSV文件应包含一列分组信息（如 'Group', 'Condition', 'Treatment' 等）。\n"
                    error_msg += "   数据格式示例：\n"
                    error_msg += "   SampleID,Group,Metabolite1,Metabolite2,...\n"
                    error_msg += "   Sample1,Control,1.2,3.4,...\n"
                    error_msg += "   Sample2,Treatment,2.3,4.5,...\n\n"
                
                if potential_group_cols:
                    error_msg += f"可能的分组列（唯一值2-10）: {', '.join(potential_group_cols[:10])}\n"
                
                if numeric_cols:
                    error_msg += f"特征列（代谢物，前5个）: {', '.join(numeric_cols[:5])}, ..."
                
                return {
                    "status": "error",
                    "error": error_msg.strip()
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

        def _fallback_metabolite_rank_output(reason: str) -> Dict[str, Any]:
            """KEGG/GSEApy 无匹配时：用组间 FC 排序代谢物，保证步骤仍有图/表/说明。"""
            top_n = metabolite_rank.replace([np.inf, -np.inf], np.nan).dropna().head(15)
            if top_n.empty:
                return {
                    "status": "warning",
                    "message": reason,
                    "enriched_pathways": [],
                    "summary": reason + " 且无法从丰度矩阵生成代谢物排序图。",
                }
            plot_path = None
            images: List[Dict[str, str]] = []
            tables: list = []
            output_csv = None
            table_rows = [
                {
                    "metabolite": str(idx)[:80],
                    "fold_change": round(float(val), 4) if np.isfinite(val) else None,
                }
                for idx, val in top_n.items()
            ]
            tables = [table_rows]
            if output_dir:
                out_p = Path(output_dir)
                out_p.mkdir(parents=True, exist_ok=True)
                output_csv = str(out_p / "pathway_metabolite_rank_fallback.csv")
                pd.DataFrame(table_rows).to_csv(output_csv, index=False)
                try:
                    fig, ax = plt.subplots(figsize=(9, max(4, len(top_n) * 0.35)))
                    labels = [str(x)[:45] for x in top_n.index]
                    vals = top_n.values.astype(float)
                    colors = ["#e74c3c" if v >= 1 else "#3498db" for v in vals]
                    ax.barh(range(len(labels)), vals, color=colors, alpha=0.85)
                    ax.set_yticks(range(len(labels)))
                    ax.set_yticklabels(labels, fontsize=8)
                    ax.axvline(1.0, color="gray", linestyle="--", alpha=0.6)
                    ax.set_xlabel("Fold change (case/control)")
                    ax.set_title("Top metabolites by |FC| (KEGG fallback)")
                    plt.tight_layout()
                    plot_path = str(out_p / "pathway_metabolite_rank_fallback.png")
                    fig.savefig(plot_path, bbox_inches="tight", dpi=150)
                    plt.close()
                    images.append({"title": "代谢物 FC 排序（KEGG 未匹配回退）", "path": plot_path})
                except Exception as fe:
                    logger.warning("通路步骤回退图失败: %s", fe)
            summary_text = (
                f"{reason} 已改为展示组间倍数变化 Top {len(top_n)} 代谢物（非 KEGG 通路富集，供探索参考）。"
            )
            out_fb = {
                "status": "warning",
                "message": reason,
                "enriched_pathways": [],
                "summary": summary_text,
                "metrics": {"fallback_metabolites": len(top_n)},
                "skipped_reason": "kegg_mapping_empty",
            }
            if plot_path:
                out_fb["plot_path"] = plot_path
            if images:
                out_fb["images"] = images
            if tables:
                out_fb["tables"] = tables
            if output_csv:
                out_fb["download_links"] = [
                    {"title": "下载代谢物 FC 排序表", "path": output_csv, "kind": "data"},
                ]
            return out_fb
        
        # 准备 GSEApy 输入（需要代谢物名称到 KEGG ID 的映射）
        # 注意：这里使用代谢物名称作为 ID（实际应用中需要映射到 KEGG ID）
        # 为了演示，我们使用代谢物名称
        metabolite_list = metabolite_rank.index.tolist()
        
        # 执行通路富集分析
        # 注意：gseapy 需要代谢物 ID 映射到 KEGG，这里使用简化版本
        try:
            # 🔥 CRITICAL FIX: Use mapped organism name and appropriate gene sets
            gene_sets_map = {
                "human": ['KEGG_2021_Human'],
                "mouse": ['KEGG_2021_Mouse'],
                "rat": ['KEGG_2021_Rat']
            }
            gene_sets = gene_sets_map.get(gseapy_organism, ['KEGG_2021_Human'])  # Default to human
            
            # 尝试使用 KEGG 数据库
            enr = gp.enrichr(
                gene_list=metabolite_list[:100],  # 限制前100个
                gene_sets=gene_sets,
                organism=gseapy_organism,  # 🔥 CRITICAL FIX: Use mapped organism name
                outdir=None,  # 不保存文件
                verbose=False
            )
            
            # 提取结果
            if enr is not None and hasattr(enr, 'results'):
                results_df = enr.results
                if results_df is None or len(results_df) == 0:
                    return _fallback_metabolite_rank_output(
                        "代谢物名称未能匹配 KEGG/Enrichr 基因集（常见于代谢物名非 HGNC 符号）。"
                    )
                term_col = "Term" if "Term" in results_df.columns else (
                    "Gene_set" if "Gene_set" in results_df.columns else (
                        "Pathway" if "Pathway" in results_df.columns else results_df.columns[0]
                    )
                )
                p_col = "Adjusted P-value" if "Adjusted P-value" in results_df.columns else (
                    "P-value" if "P-value" in results_df.columns else None
                )

                significant_pathways = results_df.copy()
                if p_col:
                    significant_pathways = results_df[results_df[p_col] < p_value_threshold].copy()

                top_pathways = []
                if len(significant_pathways) > 0 and p_col:
                    top_pathways_df = significant_pathways.nsmallest(5, p_col)
                    top_pathways = top_pathways_df[term_col].astype(str).tolist() if term_col in top_pathways_df.columns else []

                plot_df = results_df.copy()
                if p_col:
                    plot_df = plot_df.sort_values(p_col, ascending=True)
                plot_df = plot_df.head(15)

                output_csv = None
                plot_path = None
                images: List[Dict[str, str]] = []
                tables: list = []
                if output_dir:
                    output_path = Path(output_dir)
                    output_path.mkdir(parents=True, exist_ok=True)
                    output_csv = str(output_path / "pathway_enrichment.csv")
                    results_df.to_csv(output_csv, index=False)
                    if len(plot_df) > 0 and term_col in plot_df.columns:
                        try:
                            fig, ax = plt.subplots(figsize=(9, max(4, len(plot_df) * 0.35)))
                            labels = plot_df[term_col].astype(str).str.slice(0, 50).tolist()
                            if p_col:
                                scores = -np.log10(plot_df[p_col].clip(lower=1e-300).astype(float))
                            else:
                                scores = np.arange(len(plot_df), 0, -1)
                            y_pos = np.arange(len(labels))
                            colors = ["#2ecc71" if (p_col and row[p_col] < p_value_threshold) else "#bdc3c7"
                                      for _, row in plot_df.iterrows()]
                            ax.barh(y_pos, scores, color=colors, alpha=0.85)
                            ax.set_yticks(y_pos)
                            ax.set_yticklabels(labels, fontsize=8)
                            ax.set_xlabel("-log10(adjusted P-value)" if p_col else "rank score")
                            ax.set_title("KEGG 通路富集 (Top 15)")
                            plt.tight_layout()
                            plot_path = str(output_path / "pathway_enrichment_bar.png")
                            fig.savefig(plot_path, bbox_inches="tight", dpi=150)
                            plt.close()
                            images.append({"title": "通路富集条形图", "path": plot_path})
                        except Exception as pe:
                            logger.warning("通路富集图生成失败: %s", pe)

                table_rows = []
                for _, row in plot_df.iterrows():
                    table_rows.append({
                        "pathway": str(row.get(term_col, ""))[:80],
                        "adjusted_p": round(float(row[p_col]), 6) if p_col else None,
                        "significant": bool(p_col and row[p_col] < p_value_threshold),
                    })
                if table_rows:
                    tables = [table_rows]

                n_sig = len(significant_pathways)
                if n_sig == 0:
                    summary_text = (
                        f"共检索 {len(results_df)} 条通路，在校正 P<{p_value_threshold} 下无显著富集。"
                        f"下图与表格为按 P 值排序的 Top 通路（探索性展示）。"
                    )
                    status = "warning"
                else:
                    summary_text = (
                        f"显著富集通路 {n_sig} 条（共 {len(results_df)} 条，P<{p_value_threshold}）。"
                        f"代表通路: {', '.join(top_pathways[:3]) or '见下表'}。"
                    )
                    status = "success"

                out_pe = {
                    "status": status,
                    "enriched_pathways": significant_pathways.to_dict(orient='records'),
                    "n_significant": n_sig,
                    "n_total": len(results_df),
                    "output_csv": output_csv,
                    "summary": summary_text,
                    "summary_stats": {
                        "n_significant": n_sig,
                        "n_total": len(results_df),
                        "top_pathways": top_pathways,
                        "p_value_threshold": p_value_threshold,
                    },
                    "metrics": {"significant_pathways": n_sig, "total_pathways": len(results_df)},
                }
                if plot_path:
                    out_pe["plot_path"] = plot_path
                if images:
                    out_pe["images"] = images
                if tables:
                    out_pe["tables"] = tables
                if output_csv:
                    out_pe["download_links"] = [
                        {"title": "下载通路富集结果 CSV", "path": output_csv, "kind": "data"},
                    ]
                return out_pe
            else:
                return _fallback_metabolite_rank_output("GSEApy 未返回可解析的富集结果表。")
        
        except Exception as e:
            # 🔥 CRITICAL FIX: If gseapy fails, return warning (not error) so pipeline continues
            error_msg = str(e)
            logger.warning(f"⚠️ [Pathway Enrichment] GSEApy 执行失败: {error_msg}")
            
            # Check if it's an organism-related error
            if "organism" in error_msg.lower() or "not supported" in error_msg.lower():
                return {
                    "status": "warning",
                    "message": f"Organism '{organism}' is not supported by the pathway enrichment database. Step skipped.",
                    "enriched_pathways": [],
                    "skipped_reason": f"Unsupported organism: {organism}"
                }
            
            # Other errors: return warning (not error) to allow pipeline to continue
            return {
                "status": "warning",
                "message": f"Pathway enrichment could not be performed due to database limitations. Step skipped.",
                "enriched_pathways": [],
                "skipped_reason": error_msg
            }
    
    except ImportError as e:
        # 🔥 CRITICAL FIX: ImportError should return warning, not error, so pipeline continues
        logger.warning(f"⚠️ [Pathway Enrichment] gseapy 未安装: {e}")
        return {
            "status": "warning",
            "message": "Pathway enrichment could not be performed because gseapy library is not installed. Step skipped.",
            "enriched_pathways": [],
            "skipped_reason": "gseapy not installed"
        }
    except Exception as e:
        # 🔥 CRITICAL FIX: All other exceptions should return warning, not error
        logger.warning(f"⚠️ [Pathway Enrichment] 执行失败: {e}", exc_info=True)
        return {
            "status": "warning",
            "message": f"Pathway enrichment could not be performed due to an error. Step skipped.",
            "enriched_pathways": [],
            "skipped_reason": str(e)
        }


