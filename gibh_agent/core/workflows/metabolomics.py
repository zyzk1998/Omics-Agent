"""
代谢组学工作流

定义代谢组学分析的7步标准流程和依赖关系。
"""
from typing import Dict, Any, List
from .base import BaseWorkflow


class MetabolomicsWorkflow(BaseWorkflow):
    """
    代谢组学工作流
    
    标准流程（含工作量外显与可视化升维）：
    1. metabo_data_validation - 数据校验与稀疏性检查（秀肌肉前置）
    2. inspect_data - 数据检查
    3. preprocess_data - 数据预处理
    4. pca_analysis - PCA 分析
    5. metabolomics_plsda - PLS-DA 分析（如果有分组）
    6. metabo_model_comparison - 多维模型对比（PCA/PLS-DA/VIP 1x3 图）
    7. differential_analysis - 差异分析（如果有分组）
    8. visualize_volcano - 火山图可视化（如果有分组）
    9. metabolomics_pathway_enrichment - 通路富集分析（如果有分组）
    """
    
    def get_name(self) -> str:
        """获取工作流名称"""
        return "Metabolomics"
    
    def get_description(self) -> str:
        """获取工作流描述"""
        return "代谢组学标准分析流程：数据检查 -> 预处理 -> PCA -> 监督分析（可选）-> 差异分析（可选）-> 可视化 -> 通路富集（可选）"
    
    def get_steps_dag(self) -> Dict[str, List[str]]:
        """
        获取步骤依赖图
        
        Returns:
            依赖图字典
        """
        return {
            # 步骤1: 数据校验与稀疏性检查（无依赖，秀肌肉前置）
            "metabo_data_validation": [],
            
            # 步骤2: 数据检查（依赖：metabo_data_validation）
            "inspect_data": ["metabo_data_validation"],
            
            # 步骤3: 数据预处理（依赖：inspect_data）
            "preprocess_data": ["inspect_data"],
            
            # 步骤4: PCA 分析（依赖：preprocess_data）
            "pca_analysis": ["preprocess_data"],
            
            # 步骤5: PLS-DA 分析（依赖：preprocess_data，可选）
            "metabolomics_plsda": ["preprocess_data"],
            
            # 步骤6: 多维模型对比（依赖：preprocess_data，输出 pass-through 给差异分析）
            "metabo_model_comparison": ["preprocess_data"],
            
            # 步骤7: 差异分析（依赖：metabo_model_comparison，使用其 pass-through 的 data_path）
            "differential_analysis": ["metabo_model_comparison"],
            
            # 步骤8: 火山图可视化（依赖：differential_analysis）
            "visualize_volcano": ["differential_analysis"],
            
            # 步骤9: 通路富集分析（依赖：differential_analysis）
            "metabolomics_pathway_enrichment": ["differential_analysis"],
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
            "metabo_data_validation": {
                "name": "数据校验与稀疏性检查",
                "description": "快速校验丰度矩阵与样本分组：shape、缺失率、零值比例，用于流程前置展示",
                "tool_id": "metabo_data_validation",
                "default_params": {}
            },
            "inspect_data": {
                "name": "数据检查",
                "description": "SOP规则：必须首先进行数据质量评估，检查缺失值、数据范围等",
                "tool_id": "inspect_data",
                "default_params": {}
            },
            "preprocess_data": {
                "name": "数据预处理",
                "description": "SOP规则：必须进行Log2转换和标准化，缺失值处理",
                "tool_id": "preprocess_data",
                "default_params": {
                    "log_transform": True,
                    "standardize": True,
                    "missing_imputation": "min"
                }
            },
            "pca_analysis": {
                "name": "主成分分析 (PCA)",
                "description": "SOP规则：必须进行PCA分析以探索数据结构和降维",
                "tool_id": "pca_analysis",
                "default_params": {
                    "n_components": 2,
                    "scale": True
                }
            },
            "metabolomics_plsda": {
                "name": "PLS-DA 分析",
                "description": "SOP规则：检测到分组列时，必须进行监督分析（PLS-DA）以识别组间差异",
                "tool_id": "metabolomics_plsda",
                "default_params": {
                    "n_components": 2
                }
            },
            "metabo_model_comparison": {
                "name": "多维模型对比",
                "description": "PCA + PLS-DA + VIP 1x3 对比图（无监督/有监督/特征重要性），输入需已插补与标准化",
                "tool_id": "metabo_model_comparison",
                "default_params": {
                    "data_path": "<preprocess_data_output>",
                    "meta_path": "<preprocess_data_output>",
                    "output_plot_path": "<output_dir>/metabo_model_comparison.png"
                }
            },
            "differential_analysis": {
                "name": "差异代谢物分析",
                "description": "SOP规则：必须进行差异分析以识别显著差异的代谢物",
                "tool_id": "differential_analysis",
                "default_params": {
                    "method": "t-test",
                    "p_value_threshold": 0.05,
                    "fold_change_threshold": 1.5
                }
            },
            "visualize_volcano": {
                "name": "火山图可视化",
                "description": "SOP规则：必须可视化差异分析结果，展示显著差异代谢物",
                "tool_id": "visualize_volcano",
                "default_params": {
                    "fdr_threshold": 0.05,
                    "log2fc_threshold": 1.0
                }
            },
            "metabolomics_pathway_enrichment": {
                "name": "通路富集分析",
                "description": "SOP规则：必须进行通路富集分析以理解差异代谢物的生物学意义",
                "tool_id": "metabolomics_pathway_enrichment",
                "default_params": {
                    "organism": "hsa",
                    "p_value_threshold": 0.05
                }
            }
        }
        
        if step_id not in metadata_map:
            raise ValueError(f"未知的步骤ID: {step_id}")
        
        return metadata_map[step_id]

