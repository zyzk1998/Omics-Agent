"""
数据诊断器 - 统一的数据质量评估和参数推荐模块
支持多种组学类型：scRNA-seq, Metabolomics, Bulk RNA-seq 等
"""
import logging
from typing import Dict, Any, Optional, List
import pandas as pd
import numpy as np

logger = logging.getLogger(__name__)


class DataDiagnostician:
    """
    数据诊断器
    
    职责：
    1. 基于文件元数据计算统计指标
    2. 根据组学类型应用不同的诊断策略
    3. 返回结构化的统计事实（供 LLM 生成报告）
    """
    
    def __init__(self):
        """初始化数据诊断器"""
        pass
    
    def analyze(
        self,
        file_metadata: Dict[str, Any],
        omics_type: str,
        dataframe: Optional[pd.DataFrame] = None
    ) -> Dict[str, Any]:
        """
        分析数据并返回统计事实
        
        Args:
            file_metadata: FileInspector 返回的文件元数据
            omics_type: 组学类型 ("scRNA", "Metabolomics", "BulkRNA", "default")
            dataframe: 可选的数据预览（前几行）
        
        Returns:
            包含统计事实的字典
        """
        if file_metadata.get("status") != "success":
            return {
                "status": "error",
                "error": file_metadata.get("error", "Unknown error"),
                "stats": {}
            }
        
        # 根据组学类型分发到不同的分析器
        if omics_type.lower() in ["scrna", "scrna-seq", "single_cell", "single-cell"]:
            return self._analyze_scRNA(file_metadata, dataframe)
        elif omics_type.lower() in ["metabolomics", "metabolomic", "metabonomics"]:
            return self._analyze_metabolomics(file_metadata, dataframe)
        elif omics_type.lower() in ["bulkrna", "bulk_rna", "bulk-rna", "rna-seq"]:
            return self._analyze_bulkRNA(file_metadata, dataframe)
        else:
            return self._analyze_default(file_metadata, dataframe)
    
    def _analyze_scRNA(
        self,
        file_metadata: Dict[str, Any],
        dataframe: Optional[pd.DataFrame] = None
    ) -> Dict[str, Any]:
        """
        分析单细胞转录组数据
        
        Args:
            file_metadata: 文件元数据（来自 FileInspector 或 ScanpyTool.inspect_file）
            dataframe: 可选的数据预览
        
        Returns:
            统计事实字典
        """
        stats = {}
        
        # 🔥 CRITICAL FIX: 从 file_metadata 提取信息，支持多种格式
        # 优先级：n_samples/n_features > n_obs/n_vars > shape.rows/cols
        
        # 方法1: 尝试从通用键获取（DataDiagnostician 期望的格式）
        n_samples = file_metadata.get("n_samples")
        n_features = file_metadata.get("n_features")
        
        # 方法2: 如果通用键缺失或为零，尝试从 scRNA-seq 格式获取
        if (n_samples is None or n_samples == 0) and "n_obs" in file_metadata:
            n_samples = file_metadata.get("n_obs", 0)
            logger.debug(f"DEBUG: Using n_obs -> n_samples: {n_samples}")
        
        if (n_features is None or n_features == 0) and "n_vars" in file_metadata:
            n_features = file_metadata.get("n_vars", 0)
            logger.debug(f"DEBUG: Using n_vars -> n_features: {n_features}")
        
        # 方法3: 如果仍然缺失，尝试从 shape 获取
        if (n_samples is None or n_samples == 0):
            shape = file_metadata.get("shape", {})
            n_samples = shape.get("rows", 0)
            logger.debug(f"DEBUG: Using shape.rows -> n_samples: {n_samples}")
        
        if (n_features is None or n_features == 0):
            shape = file_metadata.get("shape", {})
            n_features = shape.get("cols", 0)
            logger.debug(f"DEBUG: Using shape.cols -> n_features: {n_features}")
        
        # 设置统计值
        stats["n_cells"] = n_samples if n_samples else 0
        stats["n_genes"] = n_features if n_features else 0
        
        # 如果存在 n_obs/n_vars（ScanpyTool 格式），提取额外信息
        if "n_obs" in file_metadata:
            stats["has_qc_metrics"] = file_metadata.get("has_qc_metrics", False)
            stats["is_normalized"] = file_metadata.get("is_normalized", False)
            stats["max_value"] = file_metadata.get("max_value", 0)
            stats["min_value"] = file_metadata.get("min_value", 0)
            
            # 如果有 QC 指标，尝试提取
            if stats["has_qc_metrics"] and dataframe is not None:
                try:
                    if "n_genes_by_counts" in dataframe.columns:
                        stats["median_genes_per_cell"] = dataframe["n_genes_by_counts"].median()
                    if "total_counts" in dataframe.columns:
                        stats["median_counts_per_cell"] = dataframe["total_counts"].median()
                    if "pct_counts_mt" in dataframe.columns:
                        stats["median_mt_percent"] = dataframe["pct_counts_mt"].median()
                except Exception as e:
                    logger.warning(f"⚠️ 提取 QC 指标失败: {e}")
        
        logger.debug(f"DEBUG: Final stats - n_cells: {stats['n_cells']}, n_genes: {stats['n_genes']}")
        
        # 数据质量评估
        stats["data_quality"] = self._assess_scRNA_quality(stats)
        
        # 参数推荐逻辑
        recommendations = self._recommend_scRNA_params(stats)
        stats["recommendations"] = recommendations
        
        return {
            "status": "success",
            "omics_type": "scRNA",
            "stats": stats
        }
    
    def _analyze_metabolomics(
        self,
        file_metadata: Dict[str, Any],
        dataframe: Optional[pd.DataFrame] = None
    ) -> Dict[str, Any]:
        """
        分析代谢组学数据
        
        Args:
            file_metadata: 文件元数据（来自 FileInspector）
            dataframe: 可选的数据预览
        
        Returns:
            统计事实字典
        """
        stats = {}
        
        # 从 file_metadata 提取信息
        summary = file_metadata.get("data", {}).get("summary", {}) or file_metadata
        
        stats["n_samples"] = summary.get("n_samples", file_metadata.get("n_samples", 0))
        stats["n_metabolites"] = summary.get("n_features", file_metadata.get("n_features", 0))
        stats["missing_rate"] = summary.get("missing_rate", file_metadata.get("missing_rate", 0))
        
        # 数据范围
        data_range = file_metadata.get("data_range", {})
        stats["min_value"] = data_range.get("min", 0)
        stats["max_value"] = data_range.get("max", 0)
        stats["mean_value"] = data_range.get("mean", 0)
        stats["median_value"] = data_range.get("median", 0)
        
        # 零值率（如果 dataframe 可用）
        if dataframe is not None:
            try:
                # 只计算数值列
                numeric_cols = [col for col in dataframe.columns if pd.api.types.is_numeric_dtype(dataframe[col])]
                if numeric_cols:
                    zero_rate = (dataframe[numeric_cols] == 0).sum().sum() / (len(dataframe) * len(numeric_cols)) * 100
                    stats["zero_rate"] = round(zero_rate, 2)
            except Exception as e:
                logger.warning(f"⚠️ 计算零值率失败: {e}")
        
        # 分组信息
        potential_groups = file_metadata.get("potential_groups", {})
        # 🔥 修复：potential_groups 是字典，不是列表
        if isinstance(potential_groups, dict):
            stats["has_groups"] = len(potential_groups) > 0
            stats["group_columns"] = list(potential_groups.keys())[:5]  # 最多返回5个
        else:
            stats["has_groups"] = False
            stats["group_columns"] = []
        
        # 数据质量评估
        stats["data_quality"] = self._assess_metabolomics_quality(stats)
        
        # 参数推荐逻辑
        recommendations = self._recommend_metabolomics_params(stats)
        stats["recommendations"] = recommendations
        
        return {
            "status": "success",
            "omics_type": "Metabolomics",
            "stats": stats
        }
    
    def _analyze_bulkRNA(
        self,
        file_metadata: Dict[str, Any],
        dataframe: Optional[pd.DataFrame] = None
    ) -> Dict[str, Any]:
        """
        分析 Bulk RNA-seq 数据
        
        Args:
            file_metadata: 文件元数据
            dataframe: 可选的数据预览
        
        Returns:
            统计事实字典
        """
        stats = {}
        
        # 基本统计
        stats["n_samples"] = file_metadata.get("n_samples", 0)
        stats["n_genes"] = file_metadata.get("n_features", 0)
        stats["missing_rate"] = file_metadata.get("missing_rate", 0)
        
        # 数据范围
        data_range = file_metadata.get("data_range", {})
        stats["min_value"] = data_range.get("min", 0)
        stats["max_value"] = data_range.get("max", 0)
        
        # 数据质量评估
        stats["data_quality"] = self._assess_bulkRNA_quality(stats)
        
        # 参数推荐
        recommendations = self._recommend_bulkRNA_params(stats)
        stats["recommendations"] = recommendations
        
        return {
            "status": "success",
            "omics_type": "BulkRNA",
            "stats": stats
        }
    
    def _analyze_default(
        self,
        file_metadata: Dict[str, Any],
        dataframe: Optional[pd.DataFrame] = None
    ) -> Dict[str, Any]:
        """
        默认分析（通用统计）
        
        Args:
            file_metadata: 文件元数据
            dataframe: 可选的数据预览
        
        Returns:
            统计事实字典
        """
        stats = {}
        
        # 基本统计
        stats["n_rows"] = file_metadata.get("shape", {}).get("rows", 0)
        stats["n_cols"] = file_metadata.get("shape", {}).get("cols", 0)
        stats["missing_rate"] = file_metadata.get("missing_rate", 0)
        
        # 数据范围
        data_range = file_metadata.get("data_range", {})
        stats["min_value"] = data_range.get("min", None)
        stats["max_value"] = data_range.get("max", None)
        
        return {
            "status": "success",
            "omics_type": "default",
            "stats": stats
        }
    
    # ================= 质量评估方法 =================
    
    def _assess_scRNA_quality(self, stats: Dict[str, Any]) -> str:
        """评估单细胞数据质量"""
        n_cells = stats.get("n_cells", 0)
        n_genes = stats.get("n_genes", 0)
        has_qc = stats.get("has_qc_metrics", False)
        
        if n_cells == 0 or n_genes == 0:
            return "数据为空"
        
        if n_cells < 100:
            return "数据量较小，建议增加样本"
        elif n_cells > 10000:
            return "数据量较大，建议使用更严格的过滤参数"
        else:
            return "数据量适中"
    
    def _assess_metabolomics_quality(self, stats: Dict[str, Any]) -> str:
        """评估代谢组数据质量"""
        missing_rate = stats.get("missing_rate", 0)
        n_samples = stats.get("n_samples", 0)
        n_metabolites = stats.get("n_metabolites", 0)
        
        if n_samples == 0 or n_metabolites == 0:
            return "数据为空"
        
        quality_issues = []
        
        if missing_rate > 20:
            quality_issues.append(f"缺失值率较高 ({missing_rate:.1f}%)")
        elif missing_rate > 10:
            quality_issues.append(f"缺失值率中等 ({missing_rate:.1f}%)")
        
        if n_samples < 10:
            quality_issues.append("样本数较少")
        
        if quality_issues:
            return "数据质量一般：" + "；".join(quality_issues)
        else:
            return "数据质量良好"
    
    def _assess_bulkRNA_quality(self, stats: Dict[str, Any]) -> str:
        """评估 Bulk RNA 数据质量"""
        missing_rate = stats.get("missing_rate", 0)
        n_samples = stats.get("n_samples", 0)
        
        if missing_rate > 10:
            return f"缺失值率较高 ({missing_rate:.1f}%)"
        elif n_samples < 3:
            return "样本数过少，无法进行可靠的统计分析"
        else:
            return "数据质量良好"
    
    # ================= 参数推荐方法 =================
    
    def _recommend_scRNA_params(self, stats: Dict[str, Any]) -> Dict[str, Any]:
        """推荐单细胞分析参数"""
        recommendations = {}
        n_cells = stats.get("n_cells", 0)
        
        # min_genes 推荐
        if n_cells > 10000:
            recommendations["min_genes"] = {"default": 200, "recommended": 500, "reason": "数据量大，需更严格过滤"}
        elif n_cells > 5000:
            recommendations["min_genes"] = {"default": 200, "recommended": 300, "reason": "数据量较大，建议提高阈值"}
        else:
            recommendations["min_genes"] = {"default": 200, "recommended": 200, "reason": "数据量适中，使用默认值"}
        
        # resolution 推荐
        if n_cells > 10000:
            recommendations["resolution"] = {"default": 0.5, "recommended": 0.8, "reason": "细胞数多，提高分辨率以发现细分亚群"}
        elif n_cells > 5000:
            recommendations["resolution"] = {"default": 0.5, "recommended": 0.6, "reason": "细胞数较多，适度提高分辨率"}
        else:
            recommendations["resolution"] = {"default": 0.5, "recommended": 0.5, "reason": "使用默认分辨率"}
        
        # max_mt 推荐（如果有 QC 指标）
        if stats.get("median_mt_percent"):
            mt_percent = stats["median_mt_percent"]
            if mt_percent < 5:
                recommendations["max_mt"] = {"default": 20, "recommended": 5, "reason": "线粒体基因比例低，可降低阈值"}
            elif mt_percent > 15:
                recommendations["max_mt"] = {"default": 20, "recommended": 25, "reason": "线粒体基因比例较高，适当放宽阈值"}
            else:
                recommendations["max_mt"] = {"default": 20, "recommended": 20, "reason": "使用默认值"}
        
        return recommendations
    
    def _recommend_metabolomics_params(self, stats: Dict[str, Any]) -> Dict[str, Any]:
        """推荐代谢组分析参数"""
        recommendations = {}
        missing_rate = stats.get("missing_rate", 0)
        max_value = stats.get("max_value", 0)
        has_groups = stats.get("has_groups", False)
        
        # 缺失值处理推荐
        if missing_rate > 20:
            recommendations["missing_imputation"] = {
                "default": "None",
                "recommended": "KNN/Min",
                "reason": f"缺失值率较高 ({missing_rate:.1f}%)，建议进行插补"
            }
        elif missing_rate > 10:
            recommendations["missing_imputation"] = {
                "default": "None",
                "recommended": "Min (最小值插补)",
                "reason": f"缺失值率中等 ({missing_rate:.1f}%)，建议最小值插补"
            }
        else:
            recommendations["missing_imputation"] = {
                "default": "None",
                "recommended": "None",
                "reason": "缺失值率较低，可不进行插补"
            }
        
        # Log2 转换推荐
        if max_value > 1000:
            recommendations["normalization"] = {
                "default": "None",
                "recommended": "Log2",
                "reason": f"数据最大值较大 ({max_value:.0f})，建议进行 Log2 转换以降低数据偏态"
            }
        elif max_value > 100:
            recommendations["normalization"] = {
                "default": "None",
                "recommended": "Log2 (可选)",
                "reason": f"数据范围较大，可考虑 Log2 转换"
            }
        else:
            recommendations["normalization"] = {
                "default": "None",
                "recommended": "None",
                "reason": "数据范围适中，可不进行 Log2 转换"
            }
        
        # 差异分析推荐
        if has_groups:
            recommendations["differential_method"] = {
                "default": "T-test",
                "recommended": "T-test/ANOVA",
                "reason": "检测到分组信息，建议进行差异分析"
            }
        else:
            recommendations["differential_method"] = {
                "default": "None",
                "recommended": "None",
                "reason": "未检测到分组信息，无法进行差异分析"
            }
        
        return recommendations
    
    def _recommend_bulkRNA_params(self, stats: Dict[str, Any]) -> Dict[str, Any]:
        """推荐 Bulk RNA 分析参数"""
        recommendations = {}
        missing_rate = stats.get("missing_rate", 0)
        
        if missing_rate > 10:
            recommendations["filtering"] = {
                "default": "None",
                "recommended": "Filter low expression",
                "reason": f"缺失值率较高 ({missing_rate:.1f}%)，建议过滤低表达基因"
            }
        
        return recommendations


