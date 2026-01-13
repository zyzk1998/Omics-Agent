"""
代谢组学分析工具
支持代谢组学数据的下载、预处理和分析

包含两个类：
1. MetabolomicsToolkit: 标准工作流工具包（仅使用标准库）
2. MetabolomicsTool: 原有工具（保持兼容性）
"""
import os
# 🔧 修复：设置 Matplotlib 配置目录（避免权限问题）
if 'MPLCONFIGDIR' not in os.environ:
    import tempfile
    os.environ['MPLCONFIGDIR'] = tempfile.mkdtemp(prefix='matplotlib_')

import requests
import pandas as pd
import numpy as np
from typing import Dict, Any, Optional, List, Tuple
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from statsmodels.stats.multitest import multipletests
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import logging
import gc

logger = logging.getLogger(__name__)


# ============================================================================
# 🔥 Step 2: MetabolomicsToolkit - 标准工作流工具包
# ============================================================================

class MetabolomicsToolkit:
    """
    标准代谢组学工作流工具包
    
    使用标准库实现：
    - pandas, numpy, scipy, statsmodels, sklearn, seaborn, matplotlib
    
    不依赖外部 API 或 Web 服务
    """
    
    def __init__(self, output_dir: Optional[str] = None):
        """
        初始化工具包
        
        Args:
            output_dir: 输出目录（用于保存图片和结果）
        """
        if output_dir:
            self.output_dir = Path(output_dir)
        else:
            self.output_dir = Path(os.getcwd()) / "results" / "metabolomics"
        self.output_dir.mkdir(parents=True, exist_ok=True)
    
    def preprocess_data(
        self,
        df: pd.DataFrame,
        method: str = 'log2_scale',
        missing_imputation: str = 'min'
    ) -> pd.DataFrame:
        """
        预处理代谢组学数据
        
        Args:
            df: 输入 DataFrame（只包含数值列，不包含元数据）
            method: 预处理方法 ('log2_scale', 'zscore', 'none')
            missing_imputation: 缺失值填充方法 ('min', 'median', 'mean')
        
        Returns:
            预处理后的 DataFrame
        """
        df_processed = df.copy()
        
        # 1. 处理缺失值
        if missing_imputation == 'min':
            df_processed = df_processed.fillna(df_processed.min())
        elif missing_imputation == 'median':
            df_processed = df_processed.fillna(df_processed.median())
        elif missing_imputation == 'mean':
            df_processed = df_processed.fillna(df_processed.mean())
        else:
            df_processed = df_processed.fillna(0)
        
        # 2. 标准化
        if method == 'log2_scale':
            # Log2 转换（处理零值和负值）
            df_processed = df_processed.apply(lambda x: np.log2(x + 1))
            # Z-score 标准化
            scaler = StandardScaler()
            df_processed = pd.DataFrame(
                scaler.fit_transform(df_processed),
                columns=df_processed.columns,
                index=df_processed.index
            )
        elif method == 'zscore':
            scaler = StandardScaler()
            df_processed = pd.DataFrame(
                scaler.fit_transform(df_processed),
                columns=df_processed.columns,
                index=df_processed.index
            )
        # method == 'none': 不做标准化
        
        return df_processed
    
    def differential_analysis(
        self,
        df: pd.DataFrame,
        group_col: str,
        case_group: str,
        control_group: str
    ) -> pd.DataFrame:
        """
        差异代谢物分析
        
        Args:
            df: 包含代谢物数据和分组信息的 DataFrame
            group_col: 分组列名
            case_group: 实验组名称
            control_group: 对照组名称
        
        Returns:
            包含以下列的 DataFrame:
            - metabolite: 代谢物名称
            - p_value: T-test p值
            - fdr: FDR 校正后的 p值
            - log2fc: Log2 倍数变化
            - regulation: 'Up' 或 'Down'
        """
        # 分离分组信息
        if group_col not in df.columns:
            raise ValueError(f"分组列 '{group_col}' 不存在于 DataFrame 中")
        
        groups = df[group_col]
        case_mask = groups == case_group
        control_mask = groups == control_group
        
        if not case_mask.any():
            raise ValueError(f"实验组 '{case_group}' 不存在")
        if not control_mask.any():
            raise ValueError(f"对照组 '{control_group}' 不存在")
        
        # 提取代谢物列（数值列，排除分组列）
        metabolite_cols = [col for col in df.columns 
                          if col != group_col and pd.api.types.is_numeric_dtype(df[col])]
        
        results = []
        
        for metabolite in metabolite_cols:
            case_values = df.loc[case_mask, metabolite].dropna()
            control_values = df.loc[control_mask, metabolite].dropna()
            
            if len(case_values) < 2 or len(control_values) < 2:
                continue  # 跳过样本数不足的代谢物
            
            # T-test
            t_stat, p_value = stats.ttest_ind(case_values, control_values)
            
            # 计算 Log2FC
            case_mean = case_values.mean()
            control_mean = control_values.mean()
            
            # 避免除零或对数域错误
            if control_mean <= 0:
                log2fc = np.nan
            else:
                log2fc = np.log2(case_mean / control_mean) if case_mean > 0 else np.nan
            
            results.append({
                'metabolite': metabolite,
                'p_value': p_value,
                'log2fc': log2fc
            })
        
        # 转换为 DataFrame
        diff_df = pd.DataFrame(results)
        
        # FDR 校正（Benjamini-Hochberg）
        if len(diff_df) > 0:
            _, fdr, _, _ = multipletests(
                diff_df['p_value'].fillna(1.0),
                method='fdr_bh',
                alpha=0.05
            )
            diff_df['fdr'] = fdr
            
            # 判断上调/下调
            diff_df['regulation'] = diff_df.apply(
                lambda row: 'Up' if row['log2fc'] > 0 else 'Down' if not np.isnan(row['log2fc']) else 'N/A',
                axis=1
            )
        else:
            diff_df['fdr'] = np.nan
            diff_df['regulation'] = 'N/A'
        
        return diff_df
    
    def run_pca(
        self,
        df: pd.DataFrame,
        group_col: Optional[str] = None,
        n_components: int = 10
    ) -> Dict[str, Any]:
        """
        执行主成分分析 (PCA)
        
        Args:
            df: 输入 DataFrame（只包含数值列）
            group_col: 可选的分组列名（用于可视化）
            n_components: 主成分数量
        
        Returns:
            包含以下键的字典:
            - coordinates: PCA 坐标 (DataFrame)
            - explained_variance: 解释方差比例 (array)
            - explained_variance_ratio: 解释方差比例 (array)
            - components: 主成分载荷 (DataFrame)
            - groups: 分组信息（如果提供了 group_col）
        """
        # 确保只使用数值列
        numeric_df = df.select_dtypes(include=[np.number])
        
        # 执行 PCA
        n_components = min(n_components, numeric_df.shape[0], numeric_df.shape[1])
        pca = PCA(n_components=n_components)
        pca_result = pca.fit_transform(numeric_df)
        
        # 构建结果
        result = {
            'coordinates': pd.DataFrame(
                pca_result,
                columns=[f'PC{i+1}' for i in range(n_components)],
                index=numeric_df.index
            ),
            'explained_variance': pca.explained_variance_,
            'explained_variance_ratio': pca.explained_variance_ratio_,
            'components': pd.DataFrame(
                pca.components_.T,
                columns=[f'PC{i+1}' for i in range(n_components)],
                index=numeric_df.columns
            )
        }
        
        # 如果有分组信息，添加到结果中
        if group_col and group_col in df.columns:
            result['groups'] = df[group_col]
        
        return result
    
    def plot_volcano(
        self,
        diff_df: pd.DataFrame,
        output_path: Optional[str] = None,
        fdr_threshold: float = 0.05,
        log2fc_threshold: float = 1.0
    ) -> str:
        """
        生成火山图 (Volcano Plot)
        
        Args:
            diff_df: 差异分析结果 DataFrame（必须包含 'log2fc', 'fdr' 列）
            output_path: 输出文件路径（如果为 None，自动生成）
            fdr_threshold: FDR 阈值
            log2fc_threshold: Log2FC 阈值
        
        Returns:
            保存的图片路径
        """
        if output_path is None:
            output_path = str(self.output_dir / "volcano_plot.png")
        
        # 准备数据
        diff_df = diff_df.copy()
        diff_df['-log10_fdr'] = -np.log10(diff_df['fdr'].replace(0, 1e-10))
        
        # 判断显著性
        diff_df['significant'] = (
            (diff_df['fdr'] < fdr_threshold) & 
            (np.abs(diff_df['log2fc']) > log2fc_threshold)
        )
        
        # 绘图
        plt.figure(figsize=(10, 8))
        
        # 非显著点
        non_sig = diff_df[~diff_df['significant']]
        plt.scatter(non_sig['log2fc'], non_sig['-log10_fdr'], 
                   alpha=0.5, color='gray', s=30, label='Not significant')
        
        # 显著点
        sig = diff_df[diff_df['significant']]
        if len(sig) > 0:
            up = sig[sig['log2fc'] > 0]
            down = sig[sig['log2fc'] < 0]
            
            if len(up) > 0:
                plt.scatter(up['log2fc'], up['-log10_fdr'], 
                           alpha=0.7, color='red', s=50, label=f'Up (n={len(up)})')
            if len(down) > 0:
                plt.scatter(down['log2fc'], down['-log10_fdr'], 
                           alpha=0.7, color='blue', s=50, label=f'Down (n={len(down)})')
        
        # 添加阈值线
        plt.axhline(y=-np.log10(fdr_threshold), color='black', linestyle='--', alpha=0.5)
        plt.axvline(x=log2fc_threshold, color='black', linestyle='--', alpha=0.5)
        plt.axvline(x=-log2fc_threshold, color='black', linestyle='--', alpha=0.5)
        
        plt.xlabel('Log2 Fold Change', fontsize=12)
        plt.ylabel('-Log10 FDR', fontsize=12)
        plt.title('Volcano Plot', fontsize=14, fontweight='bold')
        plt.legend()
        plt.grid(alpha=0.3)
        plt.tight_layout()
        
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        plt.close()
        
        return output_path
    
    def plot_heatmap(
        self,
        df: pd.DataFrame,
        top_n: int = 50,
        output_path: Optional[str] = None,
        group_col: Optional[str] = None
    ) -> str:
        """
        生成聚类热图
        
        Args:
            df: 输入 DataFrame（只包含数值列）
            top_n: 选择变异最大的 top_n 个代谢物
            output_path: 输出文件路径（如果为 None，自动生成）
            group_col: 可选的分组列名（用于添加分组注释）
        
        Returns:
            保存的图片路径
        """
        if output_path is None:
            output_path = str(self.output_dir / "heatmap.png")
        
        # 选择变异最大的 top_n 个代谢物
        numeric_df = df.select_dtypes(include=[np.number])
        variances = numeric_df.var().sort_values(ascending=False)
        top_metabolites = variances.head(top_n).index
        top_df = numeric_df[top_metabolites]
        
        # 如果有分组信息，添加分组注释
        row_colors = None
        if group_col and group_col in df.columns:
            # 创建分组颜色映射
            unique_groups = df[group_col].unique()
            colors = plt.cm.Set3(np.linspace(0, 1, len(unique_groups)))
            group_color_map = dict(zip(unique_groups, colors))
            row_colors = df[group_col].map(group_color_map)
        
        # 绘制热图（clustermap 返回 figure 对象）
        if row_colors is not None:
            g = sns.clustermap(
                top_df.T,
                cmap='RdYlBu_r',
                center=0,
                robust=True,
                row_cluster=True,
                col_cluster=True,
                figsize=(12, max(8, len(top_df) * 0.2)),
                cbar_kws={'label': 'Normalized Intensity'},
                row_colors=row_colors if row_colors is not None else None
            )
        else:
            g = sns.clustermap(
                top_df.T,
                cmap='RdYlBu_r',
                center=0,
                robust=True,
                row_cluster=True,
                col_cluster=True,
                figsize=(12, max(8, len(top_df) * 0.2)),
                cbar_kws={'label': 'Normalized Intensity'}
            )
        
        # 保存图片
        g.savefig(output_path, dpi=300, bbox_inches='tight')
        plt.close(g.fig)
        
        return output_path


# ============================================================================
# 原有 MetabolomicsTool 类（保持兼容性）
# ============================================================================


class MetabolomicsTool:
    """
    代谢组学分析工具
    
    核心功能：
    - 下载演示数据集
    - 数据检查和预处理
    - 主成分分析 (PCA)
    - 差异代谢物分析
    - 可视化
    """
    
    def __init__(self, config: Dict[str, Any] = None):
        """
        初始化代谢组学工具
        
        Args:
            config: 配置字典
        """
        self.config = config or {}
        self.output_dir = Path(self.config.get("output_dir", os.path.join(os.getcwd(), "results", "metabolomics"))).resolve()
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        # 存储当前加载的数据
        self.data = None
        self.metadata = None
        self.metabolites = None
        
        # 工具映射表
        self.tool_map = {
            "download_demo_data": self.download_demo_data,
            "inspect_data": self.inspect_data,
            "preprocess_data": self.preprocess_data,
            "pca_analysis": self.pca_analysis,
            "differential_analysis": self.differential_analysis,
            "visualize_pca": self.visualize_pca,
            "visualize_volcano": self.visualize_volcano,
        }
    
    def download_demo_data(
        self,
        output_dir: Optional[str] = None,
        filename: str = "human_cachexia.csv"
    ) -> Dict[str, Any]:
        """
        下载演示数据集（Human Cachexia）
        
        使用官方 API 端点下载数据：
        https://rest.xialab.ca/api/download/metaboanalyst/human_cachexia.csv
        
        Args:
            output_dir: 输出目录（如果为 None，使用 self.output_dir）
            filename: 保存的文件名（默认为 human_cachexia.csv）
        
        Returns:
            包含下载状态的字典：
            {
                "status": "success" | "error",
                "message": "描述信息",
                "file_path": "文件路径",
                "file_size": 文件大小（字节）
            }
        """
        try:
            # 确定输出目录
            if output_dir is None:
                output_dir = self.output_dir
            else:
                os.makedirs(output_dir, exist_ok=True)
            
            # 构建完整文件路径
            file_path = os.path.join(output_dir, filename)
            
            # API 端点 URL
            url = "https://rest.xialab.ca/api/download/metaboanalyst/human_cachexia.csv"
            
            print(f"📥 正在从官方 API 下载 Human Cachexia 数据集...")
            print(f"   URL: {url}")
            print(f"   保存到: {file_path}")
            
            # 使用 requests.get 下载文件
            response = requests.get(url, timeout=30)
            response.raise_for_status()  # 如果状态码不是 200，抛出异常
            
            # 保存文件
            with open(file_path, 'wb') as f:
                f.write(response.content)
            
            # 验证文件是否为空
            file_size = os.path.getsize(file_path)
            if file_size == 0:
                os.remove(file_path)  # 删除空文件
                return {
                    "status": "error",
                    "message": "下载的文件为空",
                    "file_path": file_path,
                    "file_size": 0
                }
            
            print(f"✅ 下载成功！")
            print(f"   文件大小: {file_size / 1024:.2f} KB")
            print(f"   文件路径: {file_path}")
            
            return {
                "status": "success",
                "message": "Human Cachexia 数据集下载成功",
                "file_path": file_path,
                "file_size": file_size
            }
            
        except requests.exceptions.RequestException as e:
            error_msg = f"下载失败: {str(e)}"
            print(f"❌ {error_msg}")
            return {
                "status": "error",
                "message": error_msg,
                "file_path": file_path if 'file_path' in locals() else None,
                "file_size": 0
            }
        except Exception as e:
            error_msg = f"发生错误: {str(e)}"
            print(f"❌ {error_msg}")
            return {
                "status": "error",
                "message": error_msg,
                "file_path": file_path if 'file_path' in locals() else None,
                "file_size": 0
            }
    
    def inspect_data(self, file_path: str) -> Dict[str, Any]:
        """
        检查代谢组学数据文件，返回数据摘要
        
        🔧 重构：委托给 FileInspector（Universal Eyes）
        不再重复实现检查逻辑，统一使用核心检查器
        
        Args:
            file_path: CSV 文件路径
        
        Returns:
            包含数据摘要的字典（兼容原有格式）
        """
        logger.info(f"🔍 [CHECKPOINT] inspect_data START (Delegating to FileInspector)")
        logger.info(f"   File path: {file_path}")
        
        try:
            # 🔧 委托给 FileInspector
            from ..core.file_inspector import FileInspector
            upload_dir = os.getenv("UPLOAD_DIR", "/app/uploads")
            inspector = FileInspector(upload_dir)
            
            # 使用通用检查器
            result = inspector.inspect_file(file_path)
            
            if result.get("status") == "success" and result.get("file_type") == "tabular":
                # 转换为代谢组学工具期望的格式（保持兼容性）
                summary = result.get("data", {}).get("summary", {})
                data_range = result.get("data_range", {})
                potential_groups = result.get("potential_groups", {})
                
                # 提取分组信息（如果有）
                group_info = {}
                if potential_groups:
                    # 使用第一个潜在分组列
                    first_group_col = list(potential_groups.keys())[0]
                    group_info = {
                        "column": first_group_col,
                        "groups": {str(v): 0 for v in potential_groups[first_group_col]["values"]},
                        "n_groups": potential_groups[first_group_col]["n_unique"]
                    }
                
                # 构建兼容格式的结果
                compatible_result = {
                    "status": "success",
                    "file_path": file_path,
                    "n_samples": summary.get("n_samples", "N/A"),
                    "n_metabolites": summary.get("n_features", 0),
                    "metadata_columns": result.get("metadata_columns", []),
                    "metabolite_columns": result.get("feature_columns", [])[:10],
                    "total_metabolite_columns": result.get("total_feature_columns", 0),
                    "missing_values": {
                        "total": 0,  # 不提供具体数值，只提供百分比
                        "percentage": summary.get("missing_rate", 0)
                    },
                    "group_info": group_info,
                    "data_statistics": {
                        "min": data_range.get("min", 0),
                        "max": data_range.get("max", 0),
                        "mean": data_range.get("mean", 0),
                        "median": data_range.get("median", 0)
                    },
                    # 前端可用的数据
                    "data": {
                        "summary": {
                            "n_samples": summary.get("n_samples", "N/A"),
                            "n_metabolites": summary.get("n_features", 0),
                            "missing_percentage": summary.get("missing_rate", 0),
                            "group_info": group_info,
                            "data_range": data_range,
                            "is_sampled": summary.get("is_sampled", False)
                        }
                    }
                }
                
                logger.info(f"✅ [CHECKPOINT] inspect_data SUCCESS (via FileInspector)")
                logger.info(f"   Samples: {summary.get('n_samples')}, Features: {summary.get('n_features')}, Missing: {summary.get('missing_rate', 0):.2f}%")
                return compatible_result
            else:
                # 检查失败或非表格文件
                logger.error(f"❌ [CHECKPOINT] inspect_data FAILED: {result.get('error', 'Unknown error')}")
                return result
                
        except Exception as e:
            import traceback
            error_traceback = traceback.format_exc()
            logger.error("=" * 80)
            logger.error(f"❌ [CHECKPOINT] inspect_data FAILED")
            logger.error(f"   File path: {file_path}")
            logger.error(f"   Error type: {type(e).__name__}")
            logger.error(f"   Error message: {str(e)}")
            logger.error(f"   Full traceback:")
            logger.error(error_traceback)
            logger.error("=" * 80)
            return {
                "status": "error",
                "error": str(e),
                "error_type": type(e).__name__,
                "file_path": file_path
            }
    
    def preprocess_data(
        self,
        file_path: str,
        missing_threshold: float = 0.5,
        normalization: str = "log2",
        scale: bool = True
    ) -> Dict[str, Any]:
        """
        预处理代谢组学数据
        
        Args:
            file_path: CSV 文件路径
            missing_threshold: 缺失值阈值（超过此比例的代谢物将被移除）
            normalization: 标准化方法 ("log2", "zscore", "none")
            scale: 是否进行缩放（StandardScaler）
        
        Returns:
            预处理结果字典
        """
        logger.info(f"🔍 [CHECKPOINT] preprocess_data START")
        logger.info(f"   File path: {file_path}")
        logger.info(f"   File exists? {os.path.exists(file_path)}")
        logger.info(f"   Parameters: missing_threshold={missing_threshold}, normalization={normalization}, scale={scale}")
        
        try:
            logger.info(f"🔧 开始预处理数据: {file_path}")
            
            # 🔥 修复：如果文件不存在，尝试智能路径解析
            if not os.path.exists(file_path):
                logger.warning(f"⚠️ 文件不存在: {file_path}，尝试智能路径解析...")
                from ..core.file_inspector import FileInspector
                upload_dir = os.getenv("UPLOAD_DIR", "/app/uploads")
                inspector = FileInspector(upload_dir)
                resolved_path, searched_paths = inspector._resolve_actual_path(file_path)
                if resolved_path:
                    file_path = resolved_path
                    logger.info(f"✅ 找到文件: {file_path}")
                else:
                    error_msg = f"文件未找到: {file_path}\n已搜索路径: {searched_paths[:5]}"
                    logger.error(f"❌ {error_msg}")
                    return {
                        "status": "error",
                        "message": error_msg,
                        "data": {}
                    }
            
            # 读取数据
            logger.info(f"   Attempting to read CSV file: {file_path}")
            df = pd.read_csv(file_path)
            logger.info(f"   ✅ CSV file read successfully: {len(df)} rows, {len(df.columns)} columns")
            
            # 分离元数据和代谢物数据
            metadata_cols = []
            metabolite_cols = []
            
            for col in df.columns:
                if pd.api.types.is_numeric_dtype(df[col]):
                    metabolite_cols.append(col)
                else:
                    metadata_cols.append(col)
            
            # 提取元数据
            self.metadata = df[metadata_cols].copy()
            self.metabolites = df[metabolite_cols].copy()
            
            # 1. 处理缺失值：移除缺失值比例过高的代谢物
            missing_ratio = self.metabolites.isnull().sum() / len(self.metabolites)
            valid_metabolites = missing_ratio[missing_ratio < missing_threshold].index
            removed_count = len(metabolite_cols) - len(valid_metabolites)
            
            self.metabolites = self.metabolites[valid_metabolites]
            
            # 2. 填充剩余缺失值（使用中位数）
            self.metabolites = self.metabolites.fillna(self.metabolites.median())
            
            # 3. 标准化
            if normalization == "log2":
                # Log2 转换（处理零值和负值）
                self.metabolites = self.metabolites.apply(lambda x: np.log2(x + 1))
                logger.info("✅ 已应用 Log2 转换")
            elif normalization == "zscore":
                from scipy.stats import zscore
                self.metabolites = self.metabolites.apply(zscore, axis=0)
                logger.info("✅ 已应用 Z-score 标准化")
            
            # 4. 缩放（可选）
            if scale:
                scaler = StandardScaler()
                metabolite_array = scaler.fit_transform(self.metabolites)
                self.metabolites = pd.DataFrame(
                    metabolite_array,
                    columns=self.metabolites.columns,
                    index=self.metabolites.index
                )
                logger.info("✅ 已应用 StandardScaler 缩放")
            
            # 保存预处理后的数据
            preprocessed_path = self.output_dir / "preprocessed_data.csv"
            preprocessed_df = pd.concat([self.metadata, self.metabolites], axis=1)
            preprocessed_df.to_csv(preprocessed_path, index=False)
            
            # 生成预处理后的数据预览（前5行）
            preprocessed_preview = preprocessed_df.head(5).to_dict('records')
            
            result = {
                "status": "success",
                "message": "数据预处理完成",
                "n_samples": len(self.metabolites),
                "n_metabolites": len(self.metabolites.columns),
                "removed_metabolites": removed_count,
                "normalization": normalization,
                "scaled": scale,
                "preprocessed_file": str(preprocessed_path),
                # 前端可用的数据
                "data": {
                    "preview": preprocessed_preview,  # 预处理后的前5行数据
                    "summary": {
                        "n_samples": len(self.metabolites),
                        "n_metabolites": len(self.metabolites.columns),
                        "removed_metabolites": removed_count,
                        "normalization": normalization,
                        "scaled": scale
                    }
                }
            }
            
            logger.info(f"✅ 预处理完成: {result['n_metabolites']} 个代谢物保留")
            logger.info(f"✅ [CHECKPOINT] preprocess_data SUCCESS")
            return result
            
        except Exception as e:
            import traceback
            error_traceback = traceback.format_exc()
            logger.error("=" * 80)
            logger.error(f"❌ [CHECKPOINT] preprocess_data FAILED")
            logger.error(f"   File path: {file_path}")
            logger.error(f"   Error type: {type(e).__name__}")
            logger.error(f"   Error message: {str(e)}")
            logger.error(f"   Full traceback:")
            logger.error(error_traceback)
            logger.error("=" * 80)
            return {
                "status": "error",
                "error": str(e),
                "error_type": type(e).__name__
            }
    
    def pca_analysis(
        self,
        n_components: int = 10,
        file_path: Optional[str] = None
    ) -> Dict[str, Any]:
        """
        执行主成分分析 (PCA)
        
        Args:
            n_components: 主成分数量
            file_path: 数据文件路径（如果未预处理，需要提供）
        
        Returns:
            PCA 分析结果
        """
        try:
            logger.info("=" * 80)
            logger.info("👉 [STEP 1] pca_analysis START")
            logger.info(f"   Parameters: n_components={n_components}, file_path={file_path}")
            logger.info("=" * 80)
            
            # 如果数据未加载，从文件读取
            logger.info("👉 [STEP 2] Checking if data is loaded...")
            if self.metabolites is None:
                logger.info("   Data not loaded, need to read from file")
                if file_path is None:
                    logger.error("❌ [STEP 2] 数据未加载且未提供文件路径")
                    return {
                        "status": "error",
                        "error": "数据未加载且未提供文件路径"
                    }
                
                # 🔥 修复：如果文件不存在，尝试智能路径解析
                if not os.path.exists(file_path):
                    logger.warning(f"⚠️ 文件不存在: {file_path}，尝试智能路径解析...")
                    from ..core.file_inspector import FileInspector
                    upload_dir = os.getenv("UPLOAD_DIR", "/app/uploads")
                    inspector = FileInspector(upload_dir)
                    resolved_path, searched_paths = inspector._resolve_actual_path(file_path)
                    if resolved_path:
                        file_path = resolved_path
                        logger.info(f"✅ 找到文件: {file_path}")
                    else:
                        error_msg = f"文件未找到: {file_path}\n已搜索路径: {searched_paths[:5]}"
                        logger.error(f"❌ {error_msg}")
                        return {
                            "status": "error",
                            "error": error_msg,
                            "data": {}
                        }
                
                # 读取预处理后的数据或原始数据
                logger.info(f"   Reading CSV from: {file_path}")
                df = pd.read_csv(file_path)
                logger.info(f"   Loaded CSV. Shape: {df.shape}")
                metadata_cols = [col for col in df.columns if not pd.api.types.is_numeric_dtype(df[col])]
                logger.info(f"   Found {len(metadata_cols)} metadata columns: {metadata_cols}")
                self.metabolites = df.drop(columns=metadata_cols)
                logger.info(f"   Metabolites shape: {self.metabolites.shape}")
                if len(metadata_cols) > 0:
                    self.metadata = df[metadata_cols]
                    logger.info(f"   Metadata shape: {self.metadata.shape}")
                logger.info("✅ [STEP 2] Data loaded")
                gc.collect()
            else:
                logger.info(f"✅ [STEP 2] Data already loaded. Metabolites shape: {self.metabolites.shape}")
            
            # 执行 PCA
            logger.info("👉 [STEP 3] Initializing PCA...")
            try:
                actual_n_components = min(n_components, len(self.metabolites.columns), len(self.metabolites))
                logger.info(f"   Actual n_components: {actual_n_components} (requested: {n_components}, features: {len(self.metabolites.columns)}, samples: {len(self.metabolites)})")
                pca = PCA(n_components=actual_n_components)
                logger.info("✅ [STEP 3] PCA initialized")
            except Exception as e:
                logger.error(f"❌ [STEP 3] Failed to initialize PCA: {e}", exc_info=True)
                raise
            
            logger.info("👉 [STEP 4] Fitting and transforming PCA...")
            try:
                pca_result = pca.fit_transform(self.metabolites)
                logger.info(f"✅ [STEP 4] PCA completed. Result shape: {pca_result.shape}")
            except Exception as e:
                logger.error(f"❌ [STEP 4] Failed to fit/transform PCA: {e}", exc_info=True)
                raise
            
            # 计算解释方差
            logger.info("👉 [STEP 5] Calculating explained variance...")
            try:
                explained_variance = pca.explained_variance_ratio_
                cumulative_variance = np.cumsum(explained_variance)
                logger.info(f"✅ [STEP 5] Explained variance calculated. PC1: {explained_variance[0]*100:.2f}%")
            except Exception as e:
                logger.error(f"❌ [STEP 5] Failed to calculate variance: {e}", exc_info=True)
                raise
            
            # 保存 PCA 结果
            logger.info("👉 [STEP 6] Creating PCA DataFrame...")
            try:
                pca_df = pd.DataFrame(
                    pca_result,
                    columns=[f"PC{i+1}" for i in range(pca_result.shape[1])],
                    index=self.metabolites.index
                )
                logger.info(f"✅ [STEP 6] PCA DataFrame created. Shape: {pca_df.shape}")
            except Exception as e:
                logger.error(f"❌ [STEP 6] Failed to create DataFrame: {e}", exc_info=True)
                raise
            
            # 如果有元数据，合并
            logger.info("👉 [STEP 7] Merging with metadata...")
            try:
                if self.metadata is not None:
                    pca_df = pd.concat([self.metadata, pca_df], axis=1)
                    logger.info(f"✅ [STEP 7] Merged with metadata. Final shape: {pca_df.shape}")
                else:
                    logger.info("✅ [STEP 7] No metadata to merge")
            except Exception as e:
                logger.error(f"❌ [STEP 7] Failed to merge metadata: {e}", exc_info=True)
                raise
            
            logger.info("👉 [STEP 8] Saving PCA results to CSV...")
            pca_path = self.output_dir / "pca_results.csv"
            try:
                pca_df.to_csv(pca_path, index=False)
                logger.info(f"✅ [STEP 8] Saved to {pca_path}")
                gc.collect()  # 强制垃圾回收
            except Exception as e:
                logger.error(f"❌ [STEP 8] Failed to save CSV: {e}", exc_info=True)
                raise
            
            # PCA 结果预览（前5行，包含中文列名）
            logger.info("👉 [STEP 9] Generating preview table...")
            try:
                pca_preview_df = pca_df.head(5).copy()
                pca_preview = []
                for _, row in pca_preview_df.iterrows():
                    preview_row = {}
                    for col in pca_preview_df.columns:
                        preview_row[col] = row[col]
                        # 添加中文列名映射（如果列名是 PC1, PC2 等）
                        if col.startswith("PC"):
                            preview_row[f"主成分{col}"] = row[col]
                    pca_preview.append(preview_row)
                logger.info(f"✅ [STEP 9] Preview table generated. Rows: {len(pca_preview)}")
            except Exception as e:
                logger.error(f"❌ [STEP 9] Failed to generate preview: {e}", exc_info=True)
                raise
            
            # 生成统计表格（前10个主成分的解释方差，中英文双语）
            logger.info("👉 [STEP 10] Generating variance table...")
            try:
                variance_table = []
                for i in range(min(10, len(explained_variance))):
                    variance_table.append({
                        "主成分": f"PC{i+1}",
                        "PC": f"PC{i+1}",
                        "解释方差": f"{explained_variance[i]*100:.2f}%",
                        "Explained Variance": f"{explained_variance[i]*100:.2f}%",
                        "累积方差": f"{cumulative_variance[i]*100:.2f}%",
                        "Cumulative Variance": f"{cumulative_variance[i]*100:.2f}%"
                    })
                logger.info(f"✅ [STEP 10] Variance table generated. Rows: {len(variance_table)}")
            except Exception as e:
                logger.error(f"❌ [STEP 10] Failed to generate variance table: {e}", exc_info=True)
                raise
            
            # 获取载荷表格（中英文双语）
            logger.info("👉 [STEP 11] Getting top loadings...")
            try:
                top_loadings_raw = self._get_top_loadings(pca.components_[0], self.metabolites.columns, 10)
                top_loadings = [
                    {
                        "代谢物": item["metabolite"],
                        "Metabolite": item["metabolite"],
                        "载荷值": round(item["loading"], 4),
                        "Loading": round(item["loading"], 4)
                    }
                    for item in top_loadings_raw
                ]
                logger.info(f"✅ [STEP 11] Top loadings retrieved. Count: {len(top_loadings)}")
            except Exception as e:
                logger.error(f"❌ [STEP 11] Failed to get loadings: {e}", exc_info=True)
                raise
            
            logger.info("👉 [STEP 12] Building result dictionary...")
            result = {
                "status": "success",
                "message": "PCA 分析完成",
                "n_components": pca_result.shape[1],
                "explained_variance": {
                    "PC1": float(explained_variance[0]),
                    "PC2": float(explained_variance[1]) if len(explained_variance) > 1 else 0.0,
                    "PC3": float(explained_variance[2]) if len(explained_variance) > 2 else 0.0,
                },
                "cumulative_variance_pc10": float(cumulative_variance[min(9, len(cumulative_variance)-1)]),
                "pca_file": str(pca_path),
                "loadings": {
                    "top_10_pc1": top_loadings_raw  # 保留原始格式以兼容
                },
                # 前端可用的数据
                "data": {
                    "preview": pca_preview,  # PCA 结果前5行（包含中文列名）
                    "tables": {
                        "variance_table": variance_table,  # 解释方差表格（中英文双语）
                        "top_loadings": top_loadings  # 载荷表格（中英文双语）
                    }
                }
            }
            
            logger.info("=" * 80)
            logger.info("✅ [STEP 13] pca_analysis SUCCESS")
            logger.info(f"   PC1 explains {result['explained_variance']['PC1']*100:.2f}% variance")
            logger.info("=" * 80)
            gc.collect()  # 最终垃圾回收
            return result
            
        except Exception as e:
            import traceback
            error_traceback = traceback.format_exc()
            logger.error("=" * 80)
            logger.error(f"❌ [STEP X] pca_analysis FAILED")
            logger.error(f"   Error type: {type(e).__name__}")
            logger.error(f"   Error message: {str(e)}")
            logger.error(f"   Full traceback:")
            logger.error(error_traceback)
            logger.error("=" * 80)
            gc.collect()
            return {
                "status": "error",
                "error": str(e)
            }
    
    def _get_top_loadings(self, loadings: np.ndarray, metabolite_names: List[str], n: int = 10) -> List[Dict[str, Any]]:
        """获取载荷最高的代谢物"""
        indices = np.argsort(np.abs(loadings))[-n:][::-1]
        return [
            {
                "metabolite": metabolite_names[i],
                "loading": float(loadings[i])
            }
            for i in indices
        ]
    
    def differential_analysis(
        self,
        group_column: str,
        file_path: Optional[str] = None,
        method: str = "t-test",
        p_value_threshold: float = 0.05,
        fold_change_threshold: float = 1.5,
        group1: Optional[str] = None,
        group2: Optional[str] = None
    ) -> Dict[str, Any]:
        """
        执行差异代谢物分析
        
        Args:
            group_column: 分组列名
            file_path: 数据文件路径（如果未预处理，需要提供）
            method: 统计方法 ("t-test", "mann-whitney")
            p_value_threshold: P 值阈值
            fold_change_threshold: 倍数变化阈值
        
        Returns:
            差异分析结果
        """
        try:
            logger.info("=" * 80)
            logger.info("👉 [STEP 1] differential_analysis START")
            logger.info(f"   Parameters: group_column={group_column}, method={method}, p_value_threshold={p_value_threshold}, fold_change_threshold={fold_change_threshold}")
            logger.info("=" * 80)
            
            # 如果数据未加载，从文件读取
            logger.info("👉 [STEP 2] Checking if data is loaded...")
            if self.metabolites is None or self.metadata is None:
                logger.info("   Data not loaded, need to read from file")
                if file_path is None:
                    logger.error("❌ [STEP 2] 数据未加载且未提供文件路径")
                    return {
                        "status": "error",
                        "error": "数据未加载且未提供文件路径"
                    }
                logger.info(f"   Reading CSV from: {file_path}")
                df = pd.read_csv(file_path)
                logger.info(f"   Loaded CSV. Shape: {df.shape}")
                metadata_cols = [col for col in df.columns if not pd.api.types.is_numeric_dtype(df[col])]
                logger.info(f"   Found {len(metadata_cols)} metadata columns: {metadata_cols}")
                self.metabolites = df.drop(columns=metadata_cols)
                self.metadata = df[metadata_cols]
                logger.info(f"   Metabolites shape: {self.metabolites.shape}, Metadata shape: {self.metadata.shape}")
                logger.info("✅ [STEP 2] Data loaded")
                gc.collect()
            else:
                logger.info(f"✅ [STEP 2] Data already loaded. Metabolites: {self.metabolites.shape}, Metadata: {self.metadata.shape}")
            
            # 检查分组列是否存在
            logger.info("👉 [STEP 3] Validating group column...")
            if group_column not in self.metadata.columns:
                logger.error(f"❌ [STEP 3] 分组列 '{group_column}' 不存在。可用列: {list(self.metadata.columns)}")
                # 🔥 修复 2: 优雅失败 - 返回结构化错误而不是抛出异常
                return {
                    "status": "error",
                    "error": f"分组列 '{group_column}' 不存在。可用列: {list(self.metadata.columns)}",
                    "message": f"分组列 '{group_column}' 不存在。可用列: {list(self.metadata.columns)}",
                    "available_columns": list(self.metadata.columns),
                    "data": {}  # 空数据，避免后续步骤崩溃
                }
            logger.info(f"✅ [STEP 3] Group column '{group_column}' found")
            
            # 获取分组
            logger.info("👉 [STEP 4] Getting unique groups...")
            groups = self.metadata[group_column].unique()
            logger.info(f"   Found {len(groups)} groups: {list(groups)}")
            
            if len(groups) < 2:
                logger.error(f"❌ [STEP 4] 需要至少2个分组，但找到 {len(groups)} 个: {groups}")
                return {
                    "status": "error",
                    "error": f"需要至少2个分组，但找到 {len(groups)} 个: {groups}"
                }
            elif len(groups) > 2:
                # 如果用户指定了 group1 和 group2，使用它们
                if group1 and group2:
                    if group1 not in groups or group2 not in groups:
                        logger.error(f"❌ [STEP 4] 指定的分组不存在。可用分组: {list(groups)}")
                        return {
                            "status": "error",
                            "error": f"指定的分组不存在。可用分组: {list(groups)}"
                        }
                    logger.info(f"✅ [STEP 4] Using specified groups: {group1} vs {group2}")
                    # 使用指定的分组
                    pass  # group1 和 group2 已经设置
                else:
                    logger.warning(f"⚠️ [STEP 4] 检测到 {len(groups)} 个分组，需要用户选择")
                    # 返回需要用户选择的信息
                    return {
                        "status": "need_selection",
                        "message": f"检测到 {len(groups)} 个分组: {list(groups)}",
                        "groups": list(groups),
                        "error": f"需要选择2个分组进行比较，但找到 {len(groups)} 个分组: {groups}"
                    }
            
            # 如果没有指定分组，使用前两个
            if not group1 or not group2:
                group1, group2 = groups[0], groups[1]
                logger.info(f"✅ [STEP 4] Using first two groups: {group1} vs {group2}")
            
            logger.info("👉 [STEP 5] Creating group masks...")
            group1_mask = self.metadata[group_column] == group1
            group2_mask = self.metadata[group_column] == group2
            logger.info(f"   Group1 ({group1}): {group1_mask.sum()} samples")
            logger.info(f"   Group2 ({group2}): {group2_mask.sum()} samples")
            logger.info("✅ [STEP 5] Group masks created")
            
            # 对每个代谢物执行统计检验
            logger.info("👉 [STEP 6] Running statistical tests for each metabolite...")
            logger.info(f"   Total metabolites: {len(self.metabolites.columns)}")
            results = []
            processed = 0
            for metabolite in self.metabolites.columns:
                group1_data = self.metabolites.loc[group1_mask, metabolite].values
                group2_data = self.metabolites.loc[group2_mask, metabolite].values
                
                # 计算均值
                mean1 = np.mean(group1_data)
                mean2 = np.mean(group2_data)
                
                # 计算倍数变化
                if mean2 != 0:
                    fold_change = mean1 / mean2
                    log2_fc = np.log2(fold_change) if fold_change > 0 else 0
                else:
                    fold_change = np.inf
                    log2_fc = 0
                
                # 统计检验
                if method == "t-test":
                    try:
                        stat, p_value = stats.ttest_ind(group1_data, group2_data)
                    except:
                        p_value = 1.0
                        stat = 0
                elif method == "mann-whitney":
                    try:
                        stat, p_value = stats.mannwhitneyu(group1_data, group2_data, alternative='two-sided')
                    except:
                        p_value = 1.0
                        stat = 0
                else:
                    p_value = 1.0
                    stat = 0
                
                # 判断是否显著
                is_significant = p_value < p_value_threshold and abs(log2_fc) > np.log2(fold_change_threshold)
                
                results.append({
                    "metabolite": metabolite,
                    "group1_mean": float(mean1),
                    "group2_mean": float(mean2),
                    "fold_change": float(fold_change),
                    "log2_fold_change": float(log2_fc),
                    "p_value": float(p_value),
                    "statistic": float(stat),
                    "significant": is_significant
                })
                processed += 1
                if processed % 100 == 0:
                    logger.info(f"   Processed {processed}/{len(self.metabolites.columns)} metabolites...")
            
            logger.info(f"✅ [STEP 6] Statistical tests completed. Processed {len(results)} metabolites")
            gc.collect()  # 强制垃圾回收
            
            # 转换为 DataFrame
            logger.info("👉 [STEP 7] Converting results to DataFrame...")
            try:
                results_df = pd.DataFrame(results)
                results_df = results_df.sort_values("p_value")
                logger.info(f"✅ [STEP 7] DataFrame created. Shape: {results_df.shape}")
            except Exception as e:
                logger.error(f"❌ [STEP 7] Failed to create DataFrame: {e}", exc_info=True)
                raise
            
            # 保存结果
            logger.info("👉 [STEP 8] Saving results to CSV...")
            diff_path = self.output_dir / "differential_analysis.csv"
            try:
                results_df.to_csv(diff_path, index=False)
                logger.info(f"✅ [STEP 8] Saved to {diff_path}")
                gc.collect()
            except Exception as e:
                logger.error(f"❌ [STEP 8] Failed to save CSV: {e}", exc_info=True)
                raise
            
            # 统计显著代谢物
            logger.info("👉 [STEP 9] Identifying significant metabolites...")
            significant = results_df[results_df["significant"] == True]
            logger.info(f"✅ [STEP 9] Found {len(significant)} significant metabolites out of {len(results_df)} total")
            
            # 生成差异分析结果表格（前20个显著代谢物，包含中文列名）
            top_significant_df = significant.head(20).copy()
            # 添加中文列名映射
            top_significant_table = []
            for _, row in top_significant_df.iterrows():
                top_significant_table.append({
                    "代谢物": row["metabolite"],
                    "Metabolite": row["metabolite"],  # 保留英文列名以兼容
                    "P值": round(row["p_value"], 6),
                    "P-value": round(row["p_value"], 6),
                    "倍数变化": round(row["fold_change"], 3),
                    "Fold Change": round(row["fold_change"], 3),
                    "Log2倍数变化": round(row["log2_fold_change"], 3),
                    "Log2 Fold Change": round(row["log2_fold_change"], 3),
                    "组1均值": round(row["group1_mean"], 3),
                    "Group1 Mean": round(row["group1_mean"], 3),
                    "组2均值": round(row["group2_mean"], 3),
                    "Group2 Mean": round(row["group2_mean"], 3),
                    "状态": "上调" if row["log2_fold_change"] > 0 else "下调",
                    "Status": "Up-regulated" if row["log2_fold_change"] > 0 else "Down-regulated"
                })
            
            # 生成统计摘要表格（中英文双语）
            summary_table = [
                {
                    "类别": "总代谢物数",
                    "Category": "Total Metabolites",
                    "数量": len(results_df),
                    "Count": len(results_df)
                },
                {
                    "类别": "显著代谢物数",
                    "Category": "Significant Metabolites",
                    "数量": len(significant),
                    "Count": len(significant)
                },
                {
                    "类别": "上调代谢物",
                    "Category": "Up-regulated",
                    "数量": len(significant[significant["log2_fold_change"] > 0]),
                    "Count": len(significant[significant["log2_fold_change"] > 0])
                },
                {
                    "类别": "下调代谢物",
                    "Category": "Down-regulated",
                    "数量": len(significant[significant["log2_fold_change"] < 0]),
                    "Count": len(significant[significant["log2_fold_change"] < 0])
                }
            ]
            
            result = {
                "status": "success",
                "message": "差异分析完成",
                "n_total": len(results_df),
                "n_significant": len(significant),
                "groups": {
                    "group1": str(group1),
                    "group2": str(group2)
                },
                "top_significant": top_significant_table,
                "results_file": str(diff_path),
                # 前端可用的数据
                "data": {
                    "tables": {
                        "top_significant": top_significant_table,  # 前20个显著代谢物表格
                        "summary": summary_table  # 统计摘要表格
                    },
                    "summary": {
                        "n_total": len(results_df),
                        "n_significant": len(significant),
                        "groups": {
                            "group1": str(group1),
                            "group2": str(group2)
                        }
                    }
                }
            }
            
            logger.info("=" * 80)
            logger.info("✅ [STEP 12] differential_analysis SUCCESS")
            logger.info(f"   Total metabolites: {result['n_total']}")
            logger.info(f"   Significant metabolites: {result['n_significant']}")
            logger.info("=" * 80)
            gc.collect()  # 最终垃圾回收
            return result
            
        except Exception as e:
            import traceback
            error_traceback = traceback.format_exc()
            logger.error("=" * 80)
            logger.error(f"❌ [STEP X] differential_analysis FAILED")
            logger.error(f"   Error type: {type(e).__name__}")
            logger.error(f"   Error message: {str(e)}")
            logger.error(f"   Full traceback:")
            logger.error(error_traceback)
            logger.error("=" * 80)
            gc.collect()
            return {
                "status": "error",
                "error": str(e)
            }
    
    def visualize_pca(
        self,
        group_column: Optional[str] = None,
        pca_file: Optional[str] = None,
        pc1: int = 1,
        pc2: int = 2
    ) -> Dict[str, Any]:
        """
        可视化 PCA 结果
        
        Args:
            group_column: 用于着色的分组列
            pca_file: PCA 结果文件路径
            pc1: 第一个主成分（1-based）
            pc2: 第二个主成分（1-based）
        
        Returns:
            可视化结果
        """
        try:
            logger.info("=" * 80)
            logger.info("👉 [STEP 1] visualize_pca START")
            logger.info(f"   Parameters: group_column={group_column}, pca_file={pca_file}, pc1={pc1}, pc2={pc2}")
            logger.info("=" * 80)
            
            # 读取 PCA 结果
            logger.info("👉 [STEP 2] Loading PCA results file...")
            if pca_file is None:
                pca_file = self.output_dir / "pca_results.csv"
            
            logger.info(f"   Checking file: {pca_file}")
            if not os.path.exists(pca_file):
                logger.error(f"❌ [STEP 2] PCA 结果文件不存在: {pca_file}")
                return {
                    "status": "error",
                    "error": f"PCA 结果文件不存在: {pca_file}"
                }
            
            logger.info(f"   File exists. Reading CSV...")
            df = pd.read_csv(pca_file)
            logger.info(f"✅ [STEP 2] Loaded PCA data. Shape: {df.shape}, Columns: {list(df.columns)[:5]}...")
            gc.collect()  # 强制垃圾回收
            
            # 查找主成分列
            logger.info("👉 [STEP 3] Validating principal components...")
            pc_cols = [col for col in df.columns if col.startswith("PC")]
            logger.info(f"   Found {len(pc_cols)} PC columns: {pc_cols[:5]}...")
            
            if len(pc_cols) < max(pc1, pc2):
                logger.error(f"❌ [STEP 3] 主成分数量不足: 需要 PC{pc1} 和 PC{pc2}, 但只有 {len(pc_cols)} 个")
                return {
                    "status": "error",
                    "error": f"主成分数量不足: 需要 PC{pc1} 和 PC{pc2}"
                }
            
            pc1_col = f"PC{pc1}"
            pc2_col = f"PC{pc2}"
            logger.info(f"✅ [STEP 3] Using {pc1_col} and {pc2_col}")
            
            # 创建图形
            logger.info("👉 [STEP 4] Creating matplotlib figure...")
            try:
                plt.figure(figsize=(10, 8))
                logger.info("✅ [STEP 4] Figure created")
            except Exception as e:
                logger.error(f"❌ [STEP 4] Failed to create figure: {e}", exc_info=True)
                raise
            
            # 绘制散点图
            logger.info("👉 [STEP 5] Plotting scatter points...")
            try:
                if group_column and group_column in df.columns:
                    # 按分组着色
                    groups = df[group_column].unique()
                    logger.info(f"   Found {len(groups)} groups: {list(groups)}")
                    colors = plt.cm.Set3(np.linspace(0, 1, len(groups)))
                    for i, group in enumerate(groups):
                        mask = df[group_column] == group
                        logger.info(f"   Plotting group {group}: {mask.sum()} points")
                        plt.scatter(
                            df.loc[mask, pc1_col],
                            df.loc[mask, pc2_col],
                            label=str(group),
                            alpha=0.7,
                            s=100,
                            c=[colors[i]]
                        )
                    plt.legend(title=group_column)
                else:
                    logger.info("   Plotting without grouping")
                    plt.scatter(df[pc1_col], df[pc2_col], alpha=0.7, s=100)
                logger.info("✅ [STEP 5] Scatter points plotted")
            except Exception as e:
                logger.error(f"❌ [STEP 5] Failed to plot scatter: {e}", exc_info=True)
                plt.close('all')
                gc.collect()
                raise
            
            # 设置标签和标题
            logger.info("👉 [STEP 6] Setting labels and title...")
            try:
                plt.xlabel(f"{pc1_col}", fontsize=12)
                plt.ylabel(f"{pc2_col}", fontsize=12)
                plt.title(f"PCA Plot: {pc1_col} vs {pc2_col}", fontsize=14)
                plt.grid(True, alpha=0.3)
                logger.info("✅ [STEP 6] Labels and title set")
            except Exception as e:
                logger.error(f"❌ [STEP 6] Failed to set labels: {e}", exc_info=True)
                plt.close('all')
                gc.collect()
                raise
            
            # 保存图片
            logger.info("👉 [STEP 7] Saving figure to file...")
            plot_path = self.output_dir / f"pca_plot_pc{pc1}_pc{pc2}.png"
            logger.info(f"   Output path: {plot_path}")
            try:
                plt.tight_layout()
                logger.info("   Layout adjusted, saving...")
                plt.savefig(plot_path, dpi=300, bbox_inches='tight')
                logger.info(f"✅ [STEP 7] Figure saved to {plot_path}")
            except Exception as e:
                logger.error(f"❌ [STEP 7] Failed to save figure: {e}", exc_info=True)
                plt.close('all')
                gc.collect()
                raise
            
            # CRITICAL: 关闭所有图形以释放内存
            logger.info("👉 [STEP 8] Closing matplotlib figures and freeing memory...")
            try:
                plt.close('all')
                gc.collect()
                logger.info("✅ [STEP 8] Figures closed and memory freed")
            except Exception as e:
                logger.warning(f"⚠️ [STEP 8] Error closing figures: {e}")
                gc.collect()  # 即使关闭失败也强制 GC
            
            # 转换为相对路径（相对于 output_dir）
            logger.info("👉 [STEP 9] Converting to relative path...")
            relative_path = str(plot_path)
            if os.path.isabs(relative_path):
                relative_path = os.path.relpath(relative_path, self.output_dir)
            relative_path = relative_path.replace("\\", "/")
            logger.info(f"✅ [STEP 9] Relative path: {relative_path}")
            
            result = {
                "status": "success",
                "message": "PCA 可视化完成",
                "plot_path": str(plot_path),
                "plot_file": str(plot_path),  # 兼容旧字段名
                # 前端可用的数据
                "data": {
                    "images": [relative_path]  # 相对路径，前端可以直接使用
                }
            }
            
            logger.info("=" * 80)
            logger.info("✅ [STEP 10] visualize_pca SUCCESS")
            logger.info(f"   Plot saved: {plot_path}")
            logger.info("=" * 80)
            return result
            
        except Exception as e:
            import traceback
            error_traceback = traceback.format_exc()
            logger.error("=" * 80)
            logger.error(f"❌ [STEP X] visualize_pca FAILED")
            logger.error(f"   Error type: {type(e).__name__}")
            logger.error(f"   Error message: {str(e)}")
            logger.error(f"   Full traceback:")
            logger.error(error_traceback)
            logger.error("=" * 80)
            
            # 确保清理资源
            try:
                plt.close('all')
            except:
                pass
            gc.collect()
            
            return {
                "status": "error",
                "error": str(e)
            }
    
    def visualize_volcano(
        self,
        diff_file: Optional[str] = None,
        p_value_threshold: float = 0.05,
        fold_change_threshold: float = 1.5
    ) -> Dict[str, Any]:
        """
        绘制火山图（Volcano Plot）
        
        Args:
            diff_file: 差异分析结果文件
            p_value_threshold: P 值阈值
            fold_change_threshold: 倍数变化阈值
        
        Returns:
            可视化结果
        """
        try:
            logger.info("=" * 80)
            logger.info("👉 [STEP 1] visualize_volcano START")
            logger.info(f"   Parameters: diff_file={diff_file}, p_value_threshold={p_value_threshold}, fold_change_threshold={fold_change_threshold}")
            logger.info("=" * 80)
            
            # 读取差异分析结果
            logger.info("👉 [STEP 2] Loading differential analysis results...")
            if diff_file is None:
                diff_file = self.output_dir / "differential_analysis.csv"
            
            logger.info(f"   Checking file: {diff_file}")
            if not os.path.exists(diff_file):
                logger.error(f"❌ [STEP 2] 差异分析结果文件不存在: {diff_file}")
                # 🔥 修复 3: 生成空占位图，避免 UI 崩溃
                return self._generate_empty_volcano_plot("差异分析结果文件不存在，无法生成火山图")
            
            logger.info(f"   File exists. Reading CSV...")
            df = pd.read_csv(diff_file)
            logger.info(f"✅ [STEP 2] Loaded differential data. Shape: {df.shape}, Columns: {list(df.columns)[:5]}...")
            
            # 🔥 修复 3: 检查数据是否有效（是否有 p_value 列）
            if "p_value" not in df.columns or len(df) == 0:
                logger.error(f"❌ [STEP 2] 差异分析数据无效：缺少 p_value 列或数据为空")
                return self._generate_empty_volcano_plot("差异分析数据无效，无法生成火山图")
            
            gc.collect()  # 强制垃圾回收
            
            # 计算 -log10(p_value)
            logger.info("👉 [STEP 3] Calculating -log10(p_value)...")
            try:
                df["neg_log10_p"] = -np.log10(df["p_value"] + 1e-10)  # 避免 log(0)
                logger.info(f"✅ [STEP 3] Calculated -log10(p). Range: [{df['neg_log10_p'].min():.2f}, {df['neg_log10_p'].max():.2f}]")
            except Exception as e:
                logger.error(f"❌ [STEP 3] Failed to calculate -log10(p): {e}", exc_info=True)
                raise
            
            # 分类：显著上调、显著下调、不显著
            logger.info("👉 [STEP 4] Categorizing metabolites...")
            try:
                df["category"] = "Not Significant"
                df.loc[
                    (df["p_value"] < p_value_threshold) & (df["log2_fold_change"] > np.log2(fold_change_threshold)),
                    "category"
                ] = "Up"
                df.loc[
                    (df["p_value"] < p_value_threshold) & (df["log2_fold_change"] < -np.log2(fold_change_threshold)),
                    "category"
                ] = "Down"
                
                category_counts = df["category"].value_counts().to_dict()
                logger.info(f"✅ [STEP 4] Categorized. Counts: {category_counts}")
            except Exception as e:
                logger.error(f"❌ [STEP 4] Failed to categorize: {e}", exc_info=True)
                raise
            
            # 创建火山图
            logger.info("👉 [STEP 5] Creating matplotlib figure...")
            try:
                plt.figure(figsize=(12, 8))
                logger.info("✅ [STEP 5] Figure created")
            except Exception as e:
                logger.error(f"❌ [STEP 5] Failed to create figure: {e}", exc_info=True)
                raise
            
            # 绘制散点图
            logger.info("👉 [STEP 6] Plotting scatter points by category...")
            try:
                colors = {"Up": "red", "Down": "blue", "Not Significant": "gray"}
                for cat in ["Not Significant", "Up", "Down"]:
                    mask = df["category"] == cat
                    count = mask.sum()
                    logger.info(f"   Plotting {cat}: {count} points")
                    if count > 0:
                        plt.scatter(
                            df.loc[mask, "log2_fold_change"],
                            df.loc[mask, "neg_log10_p"],
                            label=cat,
                            alpha=0.6,
                            s=50,
                            c=colors[cat]
                        )
                logger.info("✅ [STEP 6] Scatter points plotted")
            except Exception as e:
                logger.error(f"❌ [STEP 6] Failed to plot scatter: {e}", exc_info=True)
                plt.close('all')
                gc.collect()
                raise
            
            # 添加阈值线
            logger.info("👉 [STEP 7] Adding threshold lines...")
            try:
                plt.axhline(y=-np.log10(p_value_threshold), color='black', linestyle='--', alpha=0.5, label=f'p={p_value_threshold}')
                plt.axvline(x=np.log2(fold_change_threshold), color='black', linestyle='--', alpha=0.5)
                plt.axvline(x=-np.log2(fold_change_threshold), color='black', linestyle='--', alpha=0.5)
                logger.info("✅ [STEP 7] Threshold lines added")
            except Exception as e:
                logger.error(f"❌ [STEP 7] Failed to add threshold lines: {e}", exc_info=True)
                plt.close('all')
                gc.collect()
                raise
            
            # 设置标签和标题
            logger.info("👉 [STEP 8] Setting labels and title...")
            try:
                plt.xlabel("Log2 Fold Change", fontsize=12)
                plt.ylabel("-Log10 P-value", fontsize=12)
                plt.title("Volcano Plot: Differential Metabolites", fontsize=14)
                plt.legend()
                plt.grid(True, alpha=0.3)
                logger.info("✅ [STEP 8] Labels and title set")
            except Exception as e:
                logger.error(f"❌ [STEP 8] Failed to set labels: {e}", exc_info=True)
                plt.close('all')
                gc.collect()
                raise
            
            # 保存图片
            logger.info("👉 [STEP 9] Saving figure to file...")
            plot_path = self.output_dir / "volcano_plot.png"
            logger.info(f"   Output path: {plot_path}")
            try:
                plt.tight_layout()
                logger.info("   Layout adjusted, saving...")
                plt.savefig(plot_path, dpi=300, bbox_inches='tight')
                logger.info(f"✅ [STEP 9] Figure saved to {plot_path}")
            except Exception as e:
                logger.error(f"❌ [STEP 9] Failed to save figure: {e}", exc_info=True)
                plt.close('all')
                gc.collect()
                raise
            
            # CRITICAL: 关闭所有图形以释放内存
            logger.info("👉 [STEP 10] Closing matplotlib figures and freeing memory...")
            try:
                plt.close('all')
                gc.collect()
                logger.info("✅ [STEP 10] Figures closed and memory freed")
            except Exception as e:
                logger.warning(f"⚠️ [STEP 10] Error closing figures: {e}")
                gc.collect()  # 即使关闭失败也强制 GC
            
            # 转换为相对路径（相对于 output_dir）
            logger.info("👉 [STEP 11] Converting to relative path...")
            relative_path = str(plot_path)
            if os.path.isabs(relative_path):
                relative_path = os.path.relpath(relative_path, self.output_dir)
            relative_path = relative_path.replace("\\", "/")
            logger.info(f"✅ [STEP 11] Relative path: {relative_path}")
            
            result = {
                "status": "success",
                "message": "火山图生成完成",
                "plot_path": str(plot_path),
                "plot_file": str(plot_path),  # 兼容旧字段名
                # 前端可用的数据
                "data": {
                    "images": [relative_path]  # 相对路径，前端可以直接使用
                }
            }
            
            logger.info("=" * 80)
            logger.info("✅ [STEP 12] visualize_volcano SUCCESS")
            logger.info(f"   Plot saved: {plot_path}")
            logger.info("=" * 80)
            return result
            
        except Exception as e:
            import traceback
            error_traceback = traceback.format_exc()
            logger.error("=" * 80)
            logger.error(f"❌ [STEP X] visualize_volcano FAILED")
            logger.error(f"   Error type: {type(e).__name__}")
            logger.error(f"   Error message: {str(e)}")
            logger.error(f"   Full traceback:")
            logger.error(error_traceback)
            logger.error("=" * 80)
            
            # 确保清理资源
            try:
                plt.close('all')
            except:
                pass
            gc.collect()
            
            # 🔥 修复 3: 生成空占位图，避免 UI 崩溃
            return self._generate_empty_volcano_plot(f"火山图生成失败: {str(e)}")
    
    def _generate_empty_volcano_plot(self, error_message: str) -> Dict[str, Any]:
        """
        生成空占位火山图（当分析失败时）
        
        🔥 修复 3: 工具健壮性 - 生成占位图，避免 UI 崩溃
        
        Args:
            error_message: 错误消息
        
        Returns:
            包含占位图的字典
        """
        try:
            import matplotlib.pyplot as plt
            import numpy as np
            
            logger.info("👉 [Placeholder] Generating empty volcano plot placeholder...")
            
            # 创建空图
            plt.figure(figsize=(12, 8))
            plt.text(0.5, 0.5, f"Analysis Failed - No Data\n\n{error_message}", 
                    ha='center', va='center', fontsize=14, 
                    bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
            plt.xlabel("Log2 Fold Change", fontsize=12)
            plt.ylabel("-Log10 P-value", fontsize=12)
            plt.title("Volcano Plot: Analysis Failed", fontsize=14)
            plt.xlim(-5, 5)
            plt.ylim(0, 5)
            plt.grid(True, alpha=0.3)
            
            # 保存占位图
            plot_path = self.output_dir / "volcano_plot.png"
            plt.tight_layout()
            plt.savefig(plot_path, dpi=300, bbox_inches='tight')
            plt.close('all')
            
            relative_path = str(plot_path)
            if os.path.isabs(relative_path):
                relative_path = os.path.relpath(relative_path, self.output_dir)
            relative_path = relative_path.replace("\\", "/")
            
            logger.info(f"✅ [Placeholder] Empty plot saved: {plot_path}")
            
            return {
                "status": "error",
                "error": error_message,
                "message": error_message,
                "plot_path": str(plot_path),
                "plot_file": str(plot_path),
                "data": {
                    "images": [relative_path]
                }
            }
        except Exception as e:
            logger.error(f"❌ [Placeholder] Failed to generate empty plot: {e}", exc_info=True)
            return {
                "status": "error",
                "error": f"无法生成占位图: {str(e)}",
                "message": error_message
            }

