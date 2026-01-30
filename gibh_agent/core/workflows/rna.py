"""
scRNA-seq 工作流

定义单细胞转录组分析的10步标准流程和依赖关系。
"""
from typing import Dict, Any, List, Optional
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
    
    def generate_template(
        self,
        target_steps: Optional[List[str]] = None,
        file_metadata: Optional[Dict[str, Any]] = None
    ) -> Dict[str, Any]:
        """
        生成工作流模板（RNA特定实现）
        
        🔥 TASK 1 FIX: 根据文件类型调整步骤顺序
        - 如果输入是FASTQ，必须先做 cellranger_count -> convert_cellranger_to_h5ad -> qc_filter
        - 如果输入是H5AD/10x，直接从 qc_filter 开始，跳过 cellranger 步骤
        
        Args:
            target_steps: 用户请求的步骤列表（如果为 None，返回完整工作流）
            file_metadata: 文件元数据（可选，用于填充参数）
            
        Returns:
            符合前端格式的工作流配置字典
        """
        # 如果没有指定目标步骤，返回完整工作流
        if target_steps is None:
            target_steps = list(self.steps_dag.keys())
        
        # 🔥 TASK 1 FIX: 根据文件类型调整步骤顺序
        file_type = file_metadata.get("file_type", "") if file_metadata else ""
        file_path = file_metadata.get("file_path", "") if file_metadata else ""
        
        # 检测是否为FASTQ文件
        is_fastq = False
        if file_type == "fastq":
            is_fastq = True
        elif file_path:
            # 检查文件路径或目录名
            import os
            if os.path.isdir(file_path):
                # 检查目录中是否包含FASTQ文件
                try:
                    fastq_files = [f for f in os.listdir(file_path) if f.endswith(('.fastq', '.fq', '.fastq.gz', '.fq.gz'))]
                    if fastq_files:
                        is_fastq = True
                except (OSError, PermissionError):
                    pass
            elif file_path.lower().endswith(('.fastq', '.fq', '.fastq.gz', '.fq.gz')):
                is_fastq = True
        
        # 🔥 TASK 1 FIX: 如果是FASTQ文件，确保包含cellranger步骤
        if is_fastq:
            # 如果目标步骤中没有cellranger相关步骤，添加它们
            if "rna_cellranger_count" not in target_steps:
                # 检查是否需要添加cellranger步骤
                # 如果用户请求的是"全流程"或包含"cellranger"，添加这些步骤
                target_steps = ["rna_cellranger_count", "rna_convert_cellranger_to_h5ad"] + [s for s in target_steps if s not in ["rna_cellranger_count", "rna_convert_cellranger_to_h5ad"]]
                logger.info(f"✅ [RNAWorkflow] 检测到FASTQ文件，添加Cell Ranger步骤: {target_steps[:2]}")
        else:
            # 如果不是FASTQ文件，移除cellranger步骤（如果存在）
            target_steps = [s for s in target_steps if s not in ["rna_cellranger_count", "rna_convert_cellranger_to_h5ad"]]
            logger.info(f"✅ [RNAWorkflow] 非FASTQ文件，跳过Cell Ranger步骤")
        
        # 解析依赖
        resolved_steps = self.resolve_dependencies(target_steps)
        
        # 🔥 CRITICAL FIX: 确保 resolved_steps 不为空
        if not resolved_steps:
            logger.warning(f"⚠️ [RNAWorkflow] resolve_dependencies 返回空列表，使用完整工作流")
            resolved_steps = list(self.steps_dag.keys())
        
        # 🔥 CRITICAL FIX: 再次确保不为空
        if not resolved_steps:
            logger.error(f"❌ [RNAWorkflow] steps_dag 为空，无法生成模板")
            raise ValueError("工作流 DAG 为空，无法生成模板")
        
        # 🔥 TASK 1 FIX: 如果是FASTQ，确保步骤顺序正确
        if is_fastq:
            # 确保cellranger步骤在qc_filter之前
            if "rna_cellranger_count" in resolved_steps and "rna_qc_filter" in resolved_steps:
                cellranger_idx = resolved_steps.index("rna_cellranger_count")
                convert_idx = resolved_steps.index("rna_convert_cellranger_to_h5ad") if "rna_convert_cellranger_to_h5ad" in resolved_steps else -1
                qc_idx = resolved_steps.index("rna_qc_filter")
                
                # 如果qc_filter在cellranger之前，重新排序
                if qc_idx < cellranger_idx:
                    # 移除这些步骤
                    resolved_steps.remove("rna_cellranger_count")
                    if convert_idx >= 0:
                        resolved_steps.remove("rna_convert_cellranger_to_h5ad")
                    resolved_steps.remove("rna_qc_filter")
                    
                    # 重新插入到正确位置
                    # 找到第一个下游步骤的位置
                    downstream_steps = ["rna_doublet_detection", "rna_normalize", "rna_hvg"]
                    insert_pos = 0
                    for ds in downstream_steps:
                        if ds in resolved_steps:
                            insert_pos = resolved_steps.index(ds)
                            break
                    
                    # 在正确位置插入
                    resolved_steps.insert(insert_pos, "rna_cellranger_count")
                    if convert_idx >= 0:
                        resolved_steps.insert(insert_pos + 1, "rna_convert_cellranger_to_h5ad")
                    resolved_steps.insert(insert_pos + (2 if convert_idx >= 0 else 1), "rna_qc_filter")
                    
                    logger.info(f"✅ [RNAWorkflow] 重新排序步骤，确保FASTQ流程正确: {resolved_steps[:5]}")
        
        # 生成步骤配置
        steps = []
        for step_id in resolved_steps:
            step_meta = self.get_step_metadata(step_id)
            
            # 构建步骤配置
            step_config = {
                "id": step_id,
                "step_id": step_id,
                "tool_id": step_meta.get("tool_id", step_id),
                "name": step_meta.get("name", step_id),
                "step_name": step_meta.get("name", step_id),
                "description": step_meta.get("description", ""),
                "desc": step_meta.get("description", "")[:100],
                "selected": True,
                "params": step_meta.get("default_params", {}).copy()
            }
            
            # 🔥 TASK 1 FIX: 根据文件类型和步骤类型填充参数
            if file_metadata:
                file_path = file_metadata.get("file_path")
                if file_path:
                    # CellRanger步骤需要fastqs_path
                    if step_id == "rna_cellranger_count":
                        step_config["params"]["fastqs_path"] = file_path
                        # 设置其他CellRanger参数
                        if "sample_id" not in step_config["params"]:
                            import os
                            sample_id = os.path.basename(file_path).replace("_fastqs", "").replace("fastqs", "").replace("_", "-")
                            if not sample_id or sample_id == "":
                                sample_id = "sample"
                            step_config["params"]["sample_id"] = sample_id
                        if "transcriptome_path" not in step_config["params"]:
                            step_config["params"]["transcriptome_path"] = "/opt/refdata-gex-GRCh38-2020-A"  # 默认值
                        if "output_dir" not in step_config["params"]:
                            import os
                            output_dir = os.path.join(os.path.dirname(file_path), "cellranger_output")
                            step_config["params"]["output_dir"] = output_dir
                    # Convert步骤需要cellranger输出路径（占位符，执行时自动填充）
                    elif step_id == "rna_convert_cellranger_to_h5ad":
                        step_config["params"]["cellranger_matrix_dir"] = "<rna_cellranger_count_output>"
                    # QC和其他步骤需要adata_path
                    elif step_id in ["rna_qc_filter", "rna_doublet_detection", "rna_normalize", "rna_hvg", 
                                    "rna_scale", "rna_pca", "rna_neighbors", "rna_umap", "rna_clustering",
                                    "rna_find_markers", "rna_cell_annotation"]:
                        if is_fastq:
                            # FASTQ流程：qc_filter应该使用convert步骤的输出
                            step_config["params"]["adata_path"] = "<rna_convert_cellranger_to_h5ad_output>" if step_id == "rna_qc_filter" else "<previous_step_output>"
                        else:
                            # 非FASTQ流程：直接使用输入文件
                            step_config["params"]["adata_path"] = file_path
            else:
                # Plan-First模式：使用占位符
                if step_id == "rna_cellranger_count":
                    step_config["params"]["fastqs_path"] = "<PENDING_UPLOAD>"
                elif step_id == "rna_convert_cellranger_to_h5ad":
                    step_config["params"]["cellranger_matrix_dir"] = "<rna_cellranger_count_output>"
                elif step_id in ["rna_qc_filter", "rna_doublet_detection", "rna_normalize", "rna_hvg",
                                "rna_scale", "rna_pca", "rna_neighbors", "rna_umap", "rna_clustering",
                                "rna_find_markers", "rna_cell_annotation"]:
                    step_config["params"]["adata_path"] = "<PENDING_UPLOAD>"
            
            steps.append(step_config)
        
        # 构建工作流配置
        workflow_name = self._generate_workflow_name(target_steps, file_metadata)
        if is_fastq:
            workflow_name = "RNA 全流程分析（含Cell Ranger）"
        
        return {
            "type": "workflow_config",
            "workflow_data": {
                "workflow_name": workflow_name,
                "name": workflow_name,
                "steps": steps
            },
            "file_paths": [file_metadata.get("file_path")] if file_metadata and file_metadata.get("file_path") else []
        }

