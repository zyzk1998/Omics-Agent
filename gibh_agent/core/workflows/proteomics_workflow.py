"""
蛋白组学标准 DAG（数据库搜库主线，13 步）

domain_name: proteomics

De novo、修饰组定位等为并列或修饰专线，不与常规 DDA/DIA 全蛋白定量主线串联。
"""
from __future__ import annotations

import logging
from typing import Any, Dict, List, Optional

import networkx as nx

from .base import BaseWorkflow

logger = logging.getLogger(__name__)


class ProteomicsWorkflow(BaseWorkflow):
    def __init__(self) -> None:
        super().__init__()
        self._validate_dag(self.steps_dag)

    def _validate_dag(self, dag: Dict[str, List[str]]) -> None:
        g = nx.DiGraph()
        for node, deps in dag.items():
            g.add_node(node)
            for d in deps:
                g.add_edge(d, node)
        if not nx.is_directed_acyclic_graph(g):
            raise ValueError("Proteomics workflow DAG must be acyclic")

    def get_name(self) -> str:
        return "proteomics"

    def get_description(self) -> str:
        return (
            "蛋白组学数据库搜库主线：格式转换与质控、谱预处理、搜库、FDR 重打分、"
            "蛋白推断、定量、插补与批次校正、标准化 QC、差异分析、标志物、"
            "富集、PPI、全景报告"
        )

    def get_steps_dag(self) -> Dict[str, List[str]]:
        return {
            "step_prot_raw_qc": [],
            "step_prot_spectrum_pre": ["step_prot_raw_qc"],
            "step_prot_db_search": ["step_prot_spectrum_pre"],
            "step_prot_fdr_rescore": ["step_prot_db_search"],
            "step_prot_inference": ["step_prot_fdr_rescore"],
            "step_prot_quant": ["step_prot_inference"],
            "step_prot_impute_batch": ["step_prot_quant"],
            "step_prot_norm_qc": ["step_prot_impute_batch"],
            "step_prot_dea": ["step_prot_norm_qc"],
            "step_prot_biomarker": ["step_prot_dea"],
            "step_prot_enrichment": ["step_prot_biomarker"],
            "step_prot_ppi": ["step_prot_enrichment"],
            "step_prot_report": ["step_prot_ppi"],
        }

    def get_step_metadata(self, step_id: str) -> Dict[str, Any]:
        m = {
            "step_prot_raw_qc": {
                "name": "原始数据转换与基础质控",
                "description": "RAW→mzML，TIC/BPC 与仪器级系统评价",
                "tool_id": "proteomics_raw_qc_conversion",
                "default_params": {"file_path": "", "mz_tolerance_ppm": 10},
            },
            "step_prot_spectrum_pre": {
                "name": "谱图预处理与特征提取",
                "description": "基线扣除、平滑、同位素解卷积与精确峰提取",
                "tool_id": "proteomics_spectrum_preprocessing",
                "default_params": {"file_path": "", "snr_threshold": 3.0},
            },
            "step_prot_db_search": {
                "name": "搜索引擎比对与肽段鉴定",
                "description": "DDA/DIA 与序列库或谱图库匹配（MaxQuant/DIA-NN 逻辑）",
                "tool_id": "proteomics_database_search",
                "default_params": {"file_path": "", "fragment_tol_da": 0.05},
            },
            "step_prot_fdr_rescore": {
                "name": "假阳性控制与重打分",
                "description": "Target-Decoy + Percolator 等多层 FDR",
                "tool_id": "proteomics_fdr_rescoring",
                "default_params": {"file_path": "", "target_fdr": 0.01},
            },
            "step_prot_inference": {
                "name": "蛋白质推断与定性分组",
                "description": "最大简约原则解析共有肽段归属",
                "tool_id": "proteomics_protein_inference",
                "default_params": {"file_path": "", "razor_min_peptides": 2},
            },
            "step_prot_quant": {
                "name": "蛋白质丰度定量计算",
                "description": "LFQ XIC 或 TMT/iTRAQ 报告离子强度",
                "tool_id": "proteomics_quantification",
                "default_params": {"file_path": "", "lfq_ratio_type": "Median"},
            },
            "step_prot_impute_batch": {
                "name": "缺失值插补与批次效应校正",
                "description": "KNN 等插补与 ComBat 等批次校正",
                "tool_id": "proteomics_imputation_batch_correction",
                "default_params": {"file_path": "", "knn_neighbors": 5},
            },
            "step_prot_norm_qc": {
                "name": "数据标准化与质量控制",
                "description": "中位数/分位数归一化与 PCA、相关性热图及离群剔除",
                "tool_id": "proteomics_normalization_qc",
                "default_params": {"file_path": "", "norm_method": "quantile"},
            },
            "step_prot_dea": {
                "name": "差异蛋白表达统计",
                "description": "Limma 等 FC 与校正 P 值",
                "tool_id": "proteomics_differential_analysis",
                "default_params": {"file_path": "", "group_column": "", "log2fc_cutoff": 1.0},
            },
            "step_prot_biomarker": {
                "name": "标志物机器学习筛选",
                "description": "RF/SVM/LASSO 与特征重要性",
                "tool_id": "proteomics_biomarker_discovery",
                "default_params": {"file_path": "", "n_top_features": 50},
            },
            "step_prot_enrichment": {
                "name": "功能富集与通路分析",
                "description": "GO/KEGG/Reactome 与 GSEA",
                "tool_id": "proteomics_functional_enrichment",
                "default_params": {"file_path": "", "enrich_padj": 0.05},
            },
            "step_prot_ppi": {
                "name": "蛋白互作网络与核心节点",
                "description": "STRING PPI 与 Degree/Betweenness Hub",
                "tool_id": "proteomics_ppi_network_analysis",
                "default_params": {"file_path": "", "string_score_cutoff": 400},
            },
            "step_prot_report": {
                "name": "蛋白组学多维全景报告",
                "description": "定性定量、标志物与通路的标准化深度报告",
                "tool_id": "proteomics_clinical_reporting",
                "default_params": {"file_path": "", "report_depth": "standard"},
            },
        }
        if step_id not in m:
            raise ValueError(f"未知的蛋白组学步骤: {step_id}")
        return m[step_id]

    def generate_template(
        self,
        target_steps: Optional[List[str]] = None,
        file_metadata: Optional[Dict[str, Any]] = None,
    ) -> Dict[str, Any]:
        if target_steps is None:
            target_steps = list(self.get_steps_dag().keys())
        resolved = self.resolve_dependencies(target_steps)
        if not resolved:
            resolved = list(self.get_steps_dag().keys())

        file_path = (file_metadata or {}).get("file_path")
        placeholder = "<PENDING_UPLOAD>" if not file_path else file_path
        dag = self.get_steps_dag()

        steps: List[Dict[str, Any]] = []
        for sid in resolved:
            meta = self.get_step_metadata(sid)
            params = dict(meta.get("default_params") or {})
            if sid == "step_prot_raw_qc":
                params["file_path"] = placeholder
            else:
                deps = dag.get(sid, [])
                dep = deps[-1] if deps else None
                params["file_path"] = f"<{dep}>" if dep else placeholder

            steps.append(
                {
                    "id": sid,
                    "step_id": sid,
                    "tool_id": meta["tool_id"],
                    "name": meta["name"],
                    "step_name": meta["name"],
                    "description": meta["description"],
                    "desc": meta["description"][:100],
                    "selected": True,
                    "params": params,
                }
            )

        wf_name = self._generate_workflow_name(target_steps, file_metadata)
        is_template = not file_path
        return {
            "type": "workflow_config",
            "workflow_data": {
                "workflow_name": wf_name,
                "name": wf_name,
                "domain_name": self.get_name(),
                "steps": steps,
            },
            "file_paths": [file_path] if file_path else [],
            "template_mode": is_template,
        }
