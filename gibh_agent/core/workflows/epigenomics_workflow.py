"""
表观遗传组学标准 DAG（ATAC-seq / ChIP-seq 开放与结合主线，13 步）

domain_name: epigenomics

WGBS/RRBS、Hi-C 等为独立实验技术路线，不与此主线串联；顺式调控与多组学整合在峰与足迹之后收尾。
"""
from __future__ import annotations

import logging
from typing import Any, Dict, List, Optional

import networkx as nx

from .base import BaseWorkflow

logger = logging.getLogger(__name__)


class EpigenomicsWorkflow(BaseWorkflow):
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
            raise ValueError("Epigenomics workflow DAG must be acyclic")

    def get_name(self) -> str:
        return "epigenomics"

    def get_description(self) -> str:
        return (
            "ATAC-seq/ChIP-seq（开放度与蛋白质-DNA 结合）主线：质控、比对、过滤去重、"
            "移位/片段、Peak、IDR、共有峰矩阵、注释、差异、Motif、足迹、"
            "顺式调控推断、表观-转录整合"
        )

    def get_steps_dag(self) -> Dict[str, List[str]]:
        return {
            "step_epi_raw_qc": [],
            "step_epi_align": ["step_epi_raw_qc"],
            "step_epi_post_filter": ["step_epi_align"],
            "step_epi_shift": ["step_epi_post_filter"],
            "step_epi_peak": ["step_epi_shift"],
            "step_epi_idr": ["step_epi_peak"],
            "step_epi_consensus": ["step_epi_idr"],
            "step_epi_peak_anno": ["step_epi_consensus"],
            "step_epi_diff": ["step_epi_peak_anno"],
            "step_epi_motif": ["step_epi_diff"],
            "step_epi_footprint": ["step_epi_motif"],
            "step_epi_cis": ["step_epi_footprint"],
            "step_epi_multi": ["step_epi_cis"],
        }

    def get_step_metadata(self, step_id: str) -> Dict[str, Any]:
        m = {
            "step_epi_raw_qc": {
                "name": "原始数据质控与接头清洗",
                "description": "测序质量评估并切除接头与低质量碱基（FastQC/TrimGalore）",
                "tool_id": "epigenomics_raw_qc_trimming",
                "default_params": {"file_path": "", "input_dir": "", "trim_quality_threshold": 20},
            },
            "step_epi_align": {
                "name": "参考基因组定向比对",
                "description": "ATAC/ChIP 等读段比对（Bowtie2/BWA；本主线不含亚硫酸氢盐/Hi-C）",
                "tool_id": "epigenomics_alignment",
                "default_params": {
                    "file_path": "",
                    "reference_id": "hg38",
                    "threads": 8,
                    "mismatch_penalty": 4,
                },
            },
            "step_epi_post_filter": {
                "name": "比对后过滤与去重",
                "description": "MAPQ 过滤、线粒体污染去除、PCR duplicate",
                "tool_id": "epigenomics_post_align_filtering",
                "default_params": {"file_path": "", "min_mapq": 30},
            },
            "step_epi_shift": {
                "name": "移位校正与片段特征分析",
                "description": "ATAC Tn5 移位校正与核小体片段分布",
                "tool_id": "epigenomics_shift_fragment_analysis",
                "default_params": {"file_path": "", "shift_correction_bp": 4},
            },
            "step_epi_peak": {
                "name": "染色质开放区/结合位点检测",
                "description": "ATAC/DNase 或 ChIP 组蛋白/TF（MACS2/Genrich）",
                "tool_id": "epigenomics_peak_calling",
                "default_params": {
                    "file_path": "",
                    "qvalue_threshold": 0.05,
                    "broad_peak": False,
                },
            },
            "step_epi_idr": {
                "name": "生物学重复一致性检验",
                "description": "FRiP 与 IDR 评估重复间 Peak 一致性",
                "tool_id": "epigenomics_reproducibility_idr",
                "default_params": {"file_path": "", "idr_threshold": 0.05},
            },
            "step_epi_consensus": {
                "name": "共有峰集合构建与矩阵生成",
                "description": "Consensus peaks 与读段计数丰度矩阵",
                "tool_id": "epigenomics_consensus_peak_counting",
                "default_params": {"file_path": "", "merge_distance_bp": 100},
            },
            "step_epi_peak_anno": {
                "name": "峰区域基因组特征深度注释",
                "description": "启动子、外显子、内含子、增强子、UTR 与靶基因",
                "tool_id": "epigenomics_peak_annotation",
                "default_params": {"file_path": "", "promoter_window_bp": 2000},
            },
            "step_epi_diff": {
                "name": "差异富集与结合状态分析",
                "description": "DESeq2/DiffBind 差异开放或差异结合",
                "tool_id": "epigenomics_diff_accessibility",
                "default_params": {"file_path": "", "padj_cutoff": 0.05},
            },
            "step_epi_motif": {
                "name": "核心转录因子基序发现",
                "description": "JASPAR 扫描与 de novo Motif",
                "tool_id": "epigenomics_motif_discovery",
                "default_params": {"file_path": "", "motif_e_value": 1e-4},
            },
            "step_epi_footprint": {
                "name": "转录因子物理足迹分析",
                "description": "核酸酶切割保护区足迹",
                "tool_id": "epigenomics_tf_footprinting",
                "default_params": {"file_path": "", "nuc_resolution_bp": 10},
            },
            "step_epi_cis": {
                "name": "顺式调控元件靶基因推断",
                "description": "共开放性与增强子-启动子交互网络",
                "tool_id": "epigenomics_cis_regulatory_interactions",
                "default_params": {"file_path": "", "max_distance_bp": 500000},
            },
            "step_epi_multi": {
                "name": "表观-转录多组学调控整合",
                "description": "联合 RNA-seq 构建 GRN",
                "tool_id": "epigenomics_multiomics_integration",
                "default_params": {"file_path": "", "data_path": "", "min_correlation": 0.3},
            },
        }
        if step_id not in m:
            raise ValueError(f"未知的表观遗传组学步骤: {step_id}")
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

            if sid == "step_epi_raw_qc":
                params["file_path"] = placeholder
                params["input_dir"] = ""
            elif sid == "step_epi_align":
                params["file_path"] = "<step_epi_raw_qc>"
                params.pop("input_dir", None)
            elif sid == "step_epi_multi":
                deps = dag.get(sid, [])
                dep = deps[-1] if deps else None
                params["file_path"] = f"<{dep}>" if dep else placeholder
                params["data_path"] = "<PENDING_UPLOAD>"
            else:
                deps = dag.get(sid, [])
                dep = deps[-1] if deps else None
                params["file_path"] = f"<{dep}>" if dep else placeholder
                params.pop("input_dir", None)
                params.pop("data_path", None)

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
