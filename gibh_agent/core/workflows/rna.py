"""
scRNA-seq 工作流

定义单细胞转录组分析的10步标准流程和依赖关系。
"""
from typing import Dict, Any, List
from .base import BaseWorkflow


class RNAWorkflow(BaseWorkflow):
    """
    scRNA-seq 工作流
    
    标准10步流程：
    1. rna_cellranger_count - Cell Ranger 计数（如果输入是 FASTQ）
    2. rna_convert_cellranger_to_h5ad - 转换为 H5AD（如果输入是 FASTQ）
    3. rna_qc_filter - 质量控制过滤
    4. rna_doublet_detection - 双联体检测
    5. rna_normalize - 数据标准化
    6. rna_hvg - 高变基因筛选
    7. rna_scale - 数据缩放
    8. rna_pca - 主成分分析
    9. rna_neighbors - 构建邻接图
    10. rna_umap - UMAP 降维
    11. rna_clustering - Leiden 聚类
    12. rna_find_markers - Marker 基因检测
    13. rna_cell_annotation - 细胞类型注释
    14. rna_export_results - 结果导出
    """
    
    def get_name(self) -> str:
        """获取工作流名称"""
        return "RNA"
    
    def get_description(self) -> str:
        """获取工作流描述"""
        return "单细胞转录组标准分析流程：Cell Ranger（可选）-> QC -> 预处理 -> 降维 -> 聚类 -> Marker 检测 -> 注释 -> 导出"
    
    def get_steps_dag(self) -> Dict[str, List[str]]:
        """
        获取步骤依赖图
        
        Returns:
            依赖图字典
        """
        return {
            # 步骤1: Cell Ranger 计数（无依赖，仅当输入是 FASTQ 时）
            "rna_cellranger_count": [],
            
            # 步骤2: 转换为 H5AD（依赖：rna_cellranger_count）
            "rna_convert_cellranger_to_h5ad": ["rna_cellranger_count"],
            
            # 步骤3: QC 过滤（依赖：输入文件或 rna_convert_cellranger_to_h5ad）
            "rna_qc_filter": [],  # 可以从原始 H5AD 或 Cell Ranger 输出开始
            
            # 步骤4: 双联体检测（依赖：rna_qc_filter）
            "rna_doublet_detection": ["rna_qc_filter"],
            
            # 步骤5: 数据标准化（依赖：rna_doublet_detection）
            "rna_normalize": ["rna_doublet_detection"],
            
            # 步骤6: 高变基因筛选（依赖：rna_normalize）
            "rna_hvg": ["rna_normalize"],
            
            # 步骤7: 数据缩放（依赖：rna_hvg）
            "rna_scale": ["rna_hvg"],
            
            # 步骤8: PCA 分析（依赖：rna_scale）
            "rna_pca": ["rna_scale"],
            
            # 步骤9: 构建邻接图（依赖：rna_pca）
            "rna_neighbors": ["rna_pca"],
            
            # 步骤10: UMAP 降维（依赖：rna_neighbors）
            "rna_umap": ["rna_neighbors"],
            
            # 步骤11: Leiden 聚类（依赖：rna_neighbors）
            "rna_clustering": ["rna_neighbors"],
            
            # 步骤12: Marker 基因检测（依赖：rna_clustering）
            "rna_find_markers": ["rna_clustering"],
            
            # 步骤13: 细胞类型注释（依赖：rna_find_markers）
            "rna_cell_annotation": ["rna_find_markers"],
            
            # 步骤14: 结果导出（依赖：rna_cell_annotation）
            "rna_export_results": ["rna_cell_annotation"],
        }
    
    def get_step_metadata(self, step_id: str) -> Dict[str, Any]:
        """
        获取步骤元数据
        
        Args:
            step_id: 步骤ID
            
        Returns:
            步骤元数据字典
        """
        metadata_map = {
            "rna_cellranger_count": {
                "name": "Cell Ranger 计数（异步）",
                "description": "SOP规则：如果输入是 FASTQ 文件，必须首先使用 Cell Ranger 进行计数",
                "tool_id": "rna_cellranger_count",
                "default_params": {
                    "localcores": 8,
                    "localmem": 32,
                    "create_bam": False
                }
            },
            "rna_convert_cellranger_to_h5ad": {
                "name": "转换为 H5AD 格式",
                "description": "SOP规则：将 Cell Ranger 输出转换为 H5AD 格式",
                "tool_id": "rna_convert_cellranger_to_h5ad",
                "default_params": {}
            },
            "rna_qc_filter": {
                "name": "质量控制过滤",
                "description": "SOP规则：必须首先进行质量控制，过滤低质量细胞和线粒体基因",
                "tool_id": "rna_qc_filter",
                "default_params": {
                    "min_genes": 200,
                    "max_mt": 20.0
                }
            },
            "rna_doublet_detection": {
                "name": "双联体检测",
                "description": "SOP规则：必须进行双联体检测以去除双联体细胞",
                "tool_id": "rna_doublet_detection",
                "default_params": {}
            },
            "rna_normalize": {
                "name": "数据标准化",
                "description": "SOP规则：必须进行LogNormalize标准化，为后续分析做准备",
                "tool_id": "rna_normalize",
                "default_params": {
                    "target_sum": 10000
                }
            },
            "rna_hvg": {
                "name": "高变基因筛选",
                "description": "SOP规则：必须筛选高变基因以降低计算复杂度",
                "tool_id": "rna_hvg",
                "default_params": {
                    "n_top_genes": 2000
                }
            },
            "rna_scale": {
                "name": "数据缩放",
                "description": "SOP规则：必须缩放数据以准备 PCA",
                "tool_id": "rna_scale",
                "default_params": {}
            },
            "rna_pca": {
                "name": "主成分分析 (PCA)",
                "description": "SOP规则：必须进行PCA分析以降低维度",
                "tool_id": "rna_pca",
                "default_params": {
                    "n_comps": 50
                }
            },
            "rna_neighbors": {
                "name": "构建邻接图",
                "description": "SOP规则：必须构建邻接图以准备 UMAP 和聚类",
                "tool_id": "rna_neighbors",
                "default_params": {
                    "n_neighbors": 15
                }
            },
            "rna_umap": {
                "name": "UMAP 降维",
                "description": "SOP规则：必须进行UMAP降维以可视化数据",
                "tool_id": "rna_umap",
                "default_params": {}
            },
            "rna_clustering": {
                "name": "Leiden 聚类",
                "description": "SOP规则：必须进行Leiden聚类以识别细胞群体",
                "tool_id": "rna_clustering",
                "default_params": {
                    "resolution": 0.5
                }
            },
            "rna_find_markers": {
                "name": "Marker 基因检测",
                "description": "SOP规则：必须检测Marker基因以识别细胞类型",
                "tool_id": "rna_find_markers",
                "default_params": {}
            },
            "rna_cell_annotation": {
                "name": "细胞类型注释",
                "description": "SOP规则：必须进行细胞类型注释以理解生物学意义",
                "tool_id": "rna_cell_annotation",
                "default_params": {}
            },
            "rna_export_results": {
                "name": "结果导出",
                "description": "SOP规则：必须导出结果（H5AD、CSV、图表）",
                "tool_id": "rna_export_results",
                "default_params": {}
            }
        }
        
        if step_id not in metadata_map:
            raise ValueError(f"未知的步骤ID: {step_id}")
        
        return metadata_map[step_id]

