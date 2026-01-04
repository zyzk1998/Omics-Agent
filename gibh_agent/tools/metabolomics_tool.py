"""
代谢组学分析工具
支持代谢组学数据的下载、预处理和分析
"""
import os
import requests
import pandas as pd
import numpy as np
from typing import Dict, Any, Optional, List
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import logging

logger = logging.getLogger(__name__)


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
        
        Args:
            file_path: CSV 文件路径
        
        Returns:
            包含数据摘要的字典
        """
        try:
            logger.info(f"📖 正在检查数据文件: {file_path}")
            
            # 读取 CSV 文件
            df = pd.read_csv(file_path)
            
            # 识别元数据列（通常是前几列，如 Patient ID, Group 等）
            # 假设第一列是样本ID，第二列是分组信息
            metadata_cols = []
            metabolite_cols = []
            
            for col in df.columns:
                # 检查是否是数值列（代谢物数据）
                if pd.api.types.is_numeric_dtype(df[col]):
                    metabolite_cols.append(col)
                else:
                    metadata_cols.append(col)
            
            # 基本信息
            n_samples = len(df)
            n_metabolites = len(metabolite_cols)
            
            # 检查缺失值
            missing_counts = df[metabolite_cols].isnull().sum()
            total_missing = missing_counts.sum()
            missing_percentage = (total_missing / (n_samples * n_metabolites)) * 100 if n_metabolites > 0 else 0
            
            # 检查分组信息（如果有）
            group_info = {}
            if len(metadata_cols) > 1:
                group_col = metadata_cols[1]  # 假设第二列是分组
                if group_col in df.columns:
                    group_counts = df[group_col].value_counts().to_dict()
                    group_info = {
                        "column": group_col,
                        "groups": group_counts,
                        "n_groups": len(group_counts)
                    }
            
            # 数据范围
            metabolite_data = df[metabolite_cols]
            data_stats = {
                "min": float(metabolite_data.min().min()),
                "max": float(metabolite_data.max().max()),
                "mean": float(metabolite_data.mean().mean()),
                "median": float(metabolite_data.median().median()),
            }
            
            # 预览前几行
            preview = df.head(3).to_dict('records')
            
            result = {
                "status": "success",
                "file_path": file_path,
                "n_samples": n_samples,
                "n_metabolites": n_metabolites,
                "metadata_columns": metadata_cols,
                "metabolite_columns": metabolite_cols[:10],  # 只显示前10个
                "total_metabolite_columns": len(metabolite_cols),
                "missing_values": {
                    "total": int(total_missing),
                    "percentage": round(missing_percentage, 2)
                },
                "group_info": group_info,
                "data_statistics": data_stats,
                "preview": preview
            }
            
            logger.info(f"✅ 数据检查完成: {n_samples} 个样本, {n_metabolites} 个代谢物")
            return result
            
        except Exception as e:
            logger.error(f"❌ 数据检查失败: {e}", exc_info=True)
            return {
                "status": "error",
                "error": str(e),
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
        try:
            logger.info(f"🔧 开始预处理数据: {file_path}")
            
            # 读取数据
            df = pd.read_csv(file_path)
            
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
            
            result = {
                "status": "success",
                "message": "数据预处理完成",
                "n_samples": len(self.metabolites),
                "n_metabolites": len(self.metabolites.columns),
                "removed_metabolites": removed_count,
                "normalization": normalization,
                "scaled": scale,
                "preprocessed_file": str(preprocessed_path)
            }
            
            logger.info(f"✅ 预处理完成: {result['n_metabolites']} 个代谢物保留")
            return result
            
        except Exception as e:
            logger.error(f"❌ 预处理失败: {e}", exc_info=True)
            return {
                "status": "error",
                "error": str(e)
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
            logger.info(f"📊 开始 PCA 分析 (n_components={n_components})")
            
            # 如果数据未加载，从文件读取
            if self.metabolites is None:
                if file_path is None:
                    return {
                        "status": "error",
                        "error": "数据未加载且未提供文件路径"
                    }
                # 读取预处理后的数据或原始数据
                df = pd.read_csv(file_path)
                metadata_cols = [col for col in df.columns if not pd.api.types.is_numeric_dtype(df[col])]
                self.metabolites = df.drop(columns=metadata_cols)
                if len(metadata_cols) > 0:
                    self.metadata = df[metadata_cols]
            
            # 执行 PCA
            pca = PCA(n_components=min(n_components, len(self.metabolites.columns), len(self.metabolites)))
            pca_result = pca.fit_transform(self.metabolites)
            
            # 计算解释方差
            explained_variance = pca.explained_variance_ratio_
            cumulative_variance = np.cumsum(explained_variance)
            
            # 保存 PCA 结果
            pca_df = pd.DataFrame(
                pca_result,
                columns=[f"PC{i+1}" for i in range(pca_result.shape[1])],
                index=self.metabolites.index
            )
            
            # 如果有元数据，合并
            if self.metadata is not None:
                pca_df = pd.concat([self.metadata, pca_df], axis=1)
            
            pca_path = self.output_dir / "pca_results.csv"
            pca_df.to_csv(pca_path, index=False)
            
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
                    "top_10_pc1": self._get_top_loadings(pca.components_[0], self.metabolites.columns, 10)
                }
            }
            
            logger.info(f"✅ PCA 完成: PC1 解释 {result['explained_variance']['PC1']*100:.2f}% 方差")
            return result
            
        except Exception as e:
            logger.error(f"❌ PCA 分析失败: {e}", exc_info=True)
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
            logger.info(f"🔬 开始差异代谢物分析 (分组: {group_column})")
            
            # 如果数据未加载，从文件读取
            if self.metabolites is None or self.metadata is None:
                if file_path is None:
                    return {
                        "status": "error",
                        "error": "数据未加载且未提供文件路径"
                    }
                df = pd.read_csv(file_path)
                metadata_cols = [col for col in df.columns if not pd.api.types.is_numeric_dtype(df[col])]
                self.metabolites = df.drop(columns=metadata_cols)
                self.metadata = df[metadata_cols]
            
            # 检查分组列是否存在
            if group_column not in self.metadata.columns:
                return {
                    "status": "error",
                    "error": f"分组列 '{group_column}' 不存在。可用列: {list(self.metadata.columns)}"
                }
            
            # 获取分组
            groups = self.metadata[group_column].unique()
            if len(groups) < 2:
                return {
                    "status": "error",
                    "error": f"需要至少2个分组，但找到 {len(groups)} 个: {groups}"
                }
            elif len(groups) > 2:
                # 如果用户指定了 group1 和 group2，使用它们
                if group1 and group2:
                    if group1 not in groups or group2 not in groups:
                        return {
                            "status": "error",
                            "error": f"指定的分组不存在。可用分组: {list(groups)}"
                        }
                    # 使用指定的分组
                    pass  # group1 和 group2 已经设置
                else:
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
            group1_mask = self.metadata[group_column] == group1
            group2_mask = self.metadata[group_column] == group2
            
            # 对每个代谢物执行统计检验
            results = []
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
            
            # 转换为 DataFrame
            results_df = pd.DataFrame(results)
            results_df = results_df.sort_values("p_value")
            
            # 保存结果
            diff_path = self.output_dir / "differential_analysis.csv"
            results_df.to_csv(diff_path, index=False)
            
            # 统计显著代谢物
            significant = results_df[results_df["significant"] == True]
            
            result = {
                "status": "success",
                "message": "差异分析完成",
                "n_total": len(results_df),
                "n_significant": len(significant),
                "groups": {
                    "group1": str(group1),
                    "group2": str(group2)
                },
                "top_significant": significant.head(20).to_dict('records'),
                "results_file": str(diff_path)
            }
            
            logger.info(f"✅ 差异分析完成: {result['n_significant']} 个显著代谢物")
            return result
            
        except Exception as e:
            logger.error(f"❌ 差异分析失败: {e}", exc_info=True)
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
            logger.info(f"📈 生成 PCA 可视化图")
            
            # 读取 PCA 结果
            if pca_file is None:
                pca_file = self.output_dir / "pca_results.csv"
            
            if not os.path.exists(pca_file):
                return {
                    "status": "error",
                    "error": f"PCA 结果文件不存在: {pca_file}"
                }
            
            df = pd.read_csv(pca_file)
            
            # 查找主成分列
            pc_cols = [col for col in df.columns if col.startswith("PC")]
            if len(pc_cols) < max(pc1, pc2):
                return {
                    "status": "error",
                    "error": f"主成分数量不足: 需要 PC{pc1} 和 PC{pc2}"
                }
            
            pc1_col = f"PC{pc1}"
            pc2_col = f"PC{pc2}"
            
            # 创建图形
            plt.figure(figsize=(10, 8))
            
            if group_column and group_column in df.columns:
                # 按分组着色
                groups = df[group_column].unique()
                colors = plt.cm.Set3(np.linspace(0, 1, len(groups)))
                for i, group in enumerate(groups):
                    mask = df[group_column] == group
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
                plt.scatter(df[pc1_col], df[pc2_col], alpha=0.7, s=100)
            
            plt.xlabel(f"{pc1_col}", fontsize=12)
            plt.ylabel(f"{pc2_col}", fontsize=12)
            plt.title(f"PCA Plot: {pc1_col} vs {pc2_col}", fontsize=14)
            plt.grid(True, alpha=0.3)
            
            # 保存图片
            plot_path = self.output_dir / f"pca_plot_pc{pc1}_pc{pc2}.png"
            plt.tight_layout()
            plt.savefig(plot_path, dpi=300, bbox_inches='tight')
            plt.close()
            
            result = {
                "status": "success",
                "message": "PCA 可视化完成",
                "plot_path": str(plot_path),
                "plot_file": str(plot_path)  # 兼容旧字段名
            }
            
            logger.info(f"✅ PCA 图已保存: {plot_path}")
            return result
            
        except Exception as e:
            logger.error(f"❌ PCA 可视化失败: {e}", exc_info=True)
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
            logger.info(f"📈 生成火山图")
            
            # 读取差异分析结果
            if diff_file is None:
                diff_file = self.output_dir / "differential_analysis.csv"
            
            if not os.path.exists(diff_file):
                return {
                    "status": "error",
                    "error": f"差异分析结果文件不存在: {diff_file}"
                }
            
            df = pd.read_csv(diff_file)
            
            # 计算 -log10(p_value)
            df["neg_log10_p"] = -np.log10(df["p_value"] + 1e-10)  # 避免 log(0)
            
            # 分类：显著上调、显著下调、不显著
            df["category"] = "Not Significant"
            df.loc[
                (df["p_value"] < p_value_threshold) & (df["log2_fold_change"] > np.log2(fold_change_threshold)),
                "category"
            ] = "Up"
            df.loc[
                (df["p_value"] < p_value_threshold) & (df["log2_fold_change"] < -np.log2(fold_change_threshold)),
                "category"
            ] = "Down"
            
            # 创建火山图
            plt.figure(figsize=(12, 8))
            
            colors = {"Up": "red", "Down": "blue", "Not Significant": "gray"}
            for cat in ["Not Significant", "Up", "Down"]:
                mask = df["category"] == cat
                plt.scatter(
                    df.loc[mask, "log2_fold_change"],
                    df.loc[mask, "neg_log10_p"],
                    label=cat,
                    alpha=0.6,
                    s=50,
                    c=colors[cat]
                )
            
            # 添加阈值线
            plt.axhline(y=-np.log10(p_value_threshold), color='black', linestyle='--', alpha=0.5, label=f'p={p_value_threshold}')
            plt.axvline(x=np.log2(fold_change_threshold), color='black', linestyle='--', alpha=0.5)
            plt.axvline(x=-np.log2(fold_change_threshold), color='black', linestyle='--', alpha=0.5)
            
            plt.xlabel("Log2 Fold Change", fontsize=12)
            plt.ylabel("-Log10 P-value", fontsize=12)
            plt.title("Volcano Plot: Differential Metabolites", fontsize=14)
            plt.legend()
            plt.grid(True, alpha=0.3)
            
            # 保存图片
            plot_path = self.output_dir / "volcano_plot.png"
            plt.tight_layout()
            plt.savefig(plot_path, dpi=300, bbox_inches='tight')
            plt.close()
            
            result = {
                "status": "success",
                "message": "火山图生成完成",
                "plot_path": str(plot_path),
                "plot_file": str(plot_path)  # 兼容旧字段名
            }
            
            logger.info(f"✅ 火山图已保存: {plot_path}")
            return result
            
        except Exception as e:
            logger.error(f"❌ 火山图生成失败: {e}", exc_info=True)
            return {
                "status": "error",
                "error": str(e)
            }

