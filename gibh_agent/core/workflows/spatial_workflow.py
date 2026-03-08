"""
Spatial transcriptomics workflow (10x Visium) — production-grade pipeline.

DAG: load_data -> spatial_data_validation -> qc_norm -> dimensionality_reduction
      -> clustering -> spatial_clustering_comparison -> spatial_neighbors
      -> spatial_autocorr -> functional_enrichment -> plot_clusters -> plot_genes
"""
import logging
from typing import Dict, Any, List, Optional

from .base import BaseWorkflow

logger = logging.getLogger(__name__)


class SpatialWorkflow(BaseWorkflow):
    """
    Spatial transcriptomics (Visium) workflow — full Scanpy/Squidpy pipeline.

    Steps:
    0. load_data               - Load Visium (spatial_load_visium_data)
    1. spatial_data_validation - Image & coordinate validation (秀肌肉)
    2. qc_norm                 - QC + normalize + log1p (spatial_preprocess_qc)
    3. dimensionality_reduction - PCA (spatial_pca_reduction)
    4. clustering              - Leiden clustering (spatial_clustering)
    5. spatial_clustering_comparison - Multi-res spatial domains 1x3 (spatial_clustering_comparison)
    6. spatial_neighbors       - Spatial graph from coordinates (spatial_calculate_neighbors)
    7. spatial_autocorr        - Moran's I for SVGs (spatial_detect_autocorr)
    8. functional_enrichment   - Pathway enrichment (spatial_pathway_enrichment)
    9. plot_clusters           - Spatial scatter colored by leiden (spatial_plot_scatter)
    10. plot_genes             - Spatial scatter colored by gene/SVG (spatial_plot_scatter)
    """

    def get_name(self) -> str:
        """Workflow name for registry routing (domain)."""
        return "Spatial"

    def get_description(self) -> str:
        """Short description of the workflow."""
        return (
            "空间转录组（Visium）标准流程：加载 -> QC/标准化 -> PCA -> 聚类 -> "
            "空间邻域 -> 空间自相关（Moran's I）-> 空间可视化（聚类/基因）"
        )

    def get_steps_dag(self) -> Dict[str, List[str]]:
        """Step dependency DAG. Higher steps depend on lower steps."""
        return {
            "load_data": [],
            "spatial_data_validation": ["load_data"],
            "qc_norm": ["load_data"],
            "dimensionality_reduction": ["qc_norm"],
            "clustering": ["dimensionality_reduction"],
            "spatial_clustering_comparison": ["clustering"],
            "spatial_neighbors": ["spatial_clustering_comparison"],
            "spatial_autocorr": ["spatial_neighbors"],
            "functional_enrichment": ["spatial_autocorr"],
            "plot_clusters": ["clustering"],
            "plot_genes": ["spatial_autocorr"],
        }

    def get_step_metadata(self, step_id: str) -> Dict[str, Any]:
        """Step display name, description, tool_id, default_params."""
        metadata_map = {
            "spatial_data_validation": {
                "name": "空间数据校验",
                "description": "快速校验 obsm['spatial'] 与 uns['spatial']，用于流程前置展示",
                "tool_id": "spatial_data_validation",
                "default_params": {"h5ad_path": ""},
            },
            "load_data": {
                "name": "加载 Visium 数据",
                "description": "从 Space Ranger 输出目录加载 10x Visium 数据，得到 AnnData（含 obsm['spatial']）",
                "tool_id": "spatial_load_visium_data",
                "default_params": {
                    "data_dir": "",
                    "output_path": "",
                    "counts_file": "filtered_feature_bc_matrix.h5",
                },
            },
            "qc_norm": {
                "name": "QC 与标准化",
                "description": "过滤低质量 spot/基因，归一化总计数，log1p 变换（spatial_preprocess_qc）",
                "tool_id": "spatial_preprocess_qc",
                "default_params": {
                    "h5ad_path": "",
                    "output_path": "",
                    "min_genes": 200,
                    "min_cells": 3,
                    "target_sum": 1e4,
                },
            },
            "dimensionality_reduction": {
                "name": "PCA 降维",
                "description": "对标准化后的表达矩阵做 PCA（spatial_pca_reduction）",
                "tool_id": "spatial_pca_reduction",
                "default_params": {
                    "h5ad_path": "",
                    "output_path": "",
                    "n_comps": 50,
                    "svd_solver": "arpack",
                },
            },
            "clustering": {
                "name": "Leiden 聚类",
                "description": "基于 PCA 构建 kNN 图并做 Leiden 聚类，得到 obs['leiden']",
                "tool_id": "spatial_clustering",
                "default_params": {
                    "h5ad_path": "",
                    "output_path": "",
                    "resolution": 0.5,
                    "n_neighbors": 15,
                    "use_rep": "X_pca",
                    "key_added": "leiden",
                },
            },
            "spatial_clustering_comparison": {
                "name": "多分辨率空间域对比",
                "description": "在 0.3/0.5/0.8 分辨率下 Leiden 空间聚类，1x3 组织切片图",
                "tool_id": "spatial_clustering_comparison",
                "default_params": {
                    "h5ad_path": "",
                    "output_plot_path": "<output_dir>/spatial_multires_comparison.png",
                    "resolutions": [0.3, 0.5, 0.8],
                },
            },
            "spatial_neighbors": {
                "name": "计算空间邻域",
                "description": "基于 obsm['spatial'] 构建空间邻接图（Delaunay），用于 Moran's I",
                "tool_id": "spatial_calculate_neighbors",
                "default_params": {
                    "h5ad_path": "",
                    "output_path": "",
                    "coord_type": "generic",
                    "delaunay": True,
                },
            },
            "spatial_autocorr": {
                "name": "空间自相关 (Moran's I)",
                "description": "计算基因的空间自相关，识别空间可变基因 (SVGs)",
                "tool_id": "spatial_detect_autocorr",
                "default_params": {
                    "h5ad_path": "",
                    "output_path": "",
                    "method": "moran",
                    "genes": None,
                },
            },
            "functional_enrichment": {
                "name": "通路富集 (SVGs)",
                "description": "对 top 空间可变基因做通路富集（gseapy），输出点图与 CSV。h5ad_path 由 spatial_autocorr 输出链式传入（含 uns['moranI']）。",
                "tool_id": "spatial_pathway_enrichment",
                "default_params": {
                    "h5ad_path": "",  # resolved from <spatial_autocorr> output by executor
                    "top_n": 50,
                    "output_path": "",
                    "gene_sets": "KEGG_2021_Human",
                },
            },
            "plot_clusters": {
                "name": "空间图（按聚类着色）",
                "description": "在空间坐标上按 Leiden 聚类着色绘制散点图（离散色板）",
                "tool_id": "spatial_plot_scatter",
                "default_params": {
                    "h5ad_path": "",
                    "color_by": "leiden",
                    "output_path": "",
                },
            },
            "plot_genes": {
                "name": "空间图（按基因/SVG 着色）",
                "description": "在空间坐标上按基因或总计数着色；可设为 spatial_autocorr 得到的 top SVG",
                "tool_id": "spatial_plot_scatter",
                "default_params": {
                    "h5ad_path": "",
                    "color_by": "total_counts",
                    "output_path": "",
                },
            },
        }
        if step_id not in metadata_map:
            raise ValueError(f"未知的步骤ID: {step_id}")
        return metadata_map[step_id]

    def generate_template(
        self,
        target_steps: Optional[List[str]] = None,
        file_metadata: Optional[Dict[str, Any]] = None,
    ) -> Dict[str, Any]:
        """
        Generate workflow template; fill data_dir / h5ad_path from file_metadata.
        Executor will chain h5ad_path from previous step output where applicable.
        """
        if target_steps is None:
            target_steps = list(self.steps_dag.keys())
        resolved_steps = self.resolve_dependencies(target_steps)
        if not resolved_steps:
            resolved_steps = list(self.steps_dag.keys())

        file_path = (file_metadata or {}).get("file_path")
        placeholder = "<PENDING_UPLOAD>" if not file_path else file_path
        # 步骤依赖：用于 h5ad_path 链（qc_norm 用 load_data 输出，以此类推）
        step_deps = self.get_steps_dag()

        steps = []
        for step_id in resolved_steps:
            step_meta = self.get_step_metadata(step_id)
            params = dict(step_meta.get("default_params", {}))
            tool_id = step_meta.get("tool_id", step_id)

            if step_id == "load_data":
                params["data_dir"] = placeholder
            else:
                deps = step_deps.get(step_id, [])
                # 取最后一个依赖作为 h5ad 来源（如 qc_norm -> load_data, clustering -> dimensionality_reduction）
                dep = deps[-1] if deps else "load_data"
                params["h5ad_path"] = f"<{dep}>"

            step_config = {
                "id": step_id,
                "step_id": step_id,
                "tool_id": tool_id,
                "name": step_meta.get("name", step_id),
                "step_name": step_meta.get("name", step_id),
                "description": step_meta.get("description", ""),
                "desc": (step_meta.get("description", ""))[:100],
                "selected": True,
                "params": params,
            }
            steps.append(step_config)

        workflow_name = self._generate_workflow_name(target_steps, file_metadata)
        is_template = not file_path
        return {
            "type": "workflow_config",
            "workflow_data": {
                "workflow_name": workflow_name,
                "name": workflow_name,
                "steps": steps,
            },
            "file_paths": [file_path] if file_path else [],
            "template_mode": is_template,
        }
