"""
基因组学标准 DAG（胚系主线，12 步）

domain_name: genomics

体细胞肿瘤检测、PGx 等为独立业务线，不与此主线串联；VQSR 与位点标准化合并为一步。
"""
from __future__ import annotations

import logging
from typing import Any, Dict, List, Optional

import networkx as nx

from .base import BaseWorkflow

logger = logging.getLogger(__name__)


class GenomicsWorkflow(BaseWorkflow):
    """WGS/WES 胚系变异主线：步骤 tool_id 与 `omics_genomics_pipeline_tools` / Registry 对齐。"""

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
            raise ValueError("Genomics workflow DAG must be acyclic")

    def get_name(self) -> str:
        return "genomics"

    def get_description(self) -> str:
        return (
            "基因组学胚系主线：质控、修剪、比对、去重、BQSR、胚系变异、"
            "CNV、SV、VQSR 与位点标准化、注释、ACMG、临床报告"
        )

    def get_steps_dag(self) -> Dict[str, List[str]]:
        return {
            "step_genomics_raw_qc": [],
            "step_genomics_read_trim": ["step_genomics_raw_qc"],
            "step_genomics_align": ["step_genomics_read_trim"],
            "step_genomics_mark_dup": ["step_genomics_align"],
            "step_genomics_bqsr": ["step_genomics_mark_dup"],
            "step_genomics_germline": ["step_genomics_bqsr"],
            "step_genomics_cnv": ["step_genomics_germline"],
            "step_genomics_sv": ["step_genomics_cnv"],
            "step_genomics_vqsr": ["step_genomics_sv"],
            "step_genomics_anno": ["step_genomics_vqsr"],
            "step_genomics_acmg": ["step_genomics_anno"],
            "step_genomics_report": ["step_genomics_acmg"],
        }

    def get_step_metadata(self, step_id: str) -> Dict[str, Any]:
        m = {
            "step_genomics_raw_qc": {
                "name": "原始数据质控",
                "description": "FASTQ 序列质量、GC 偏好及污染源评估（FastQC/MultiQC）",
                "tool_id": "genomics_raw_qc",
                "default_params": {"file_path": "", "input_dir": "", "min_read_length": 36},
            },
            "step_genomics_read_trim": {
                "name": "序列清洗与接头过滤",
                "description": "切除低质量末端、去接头、过滤极短读段（fastp）",
                "tool_id": "genomics_read_trimming",
                "default_params": {"file_path": "", "input_dir": "", "quality_cutoff": 20},
            },
            "step_genomics_align": {
                "name": "参考基因组比对",
                "description": "高精度比对生成 BAM（BWA-MEM / DRAGEN）",
                "tool_id": "genomics_alignment",
                "default_params": {
                    "file_path": "",
                    "reference_id": "hg38",
                    "threads": 8,
                    "mismatch_penalty": 4,
                    "gap_open_penalty": 6,
                },
            },
            "step_genomics_mark_dup": {
                "name": "排序与重复标记",
                "description": "坐标排序并标记 PCR 冗余重复（Picard）",
                "tool_id": "genomics_mark_duplicates",
                "default_params": {"file_path": "", "optical_duplicate_distance": 2500},
            },
            "step_genomics_bqsr": {
                "name": "碱基质量重校准",
                "description": "BQSR 校正测序仪系统性碱基质量误差",
                "tool_id": "genomics_bqsr",
                "default_params": {"file_path": "", "bqsr_max_cycles": 2},
            },
            "step_genomics_germline": {
                "name": "胚系变异检测",
                "description": "遗传性 SNP/Indel（HaplotypeCaller / DeepVariant）",
                "tool_id": "genomics_germline_calling",
                "default_params": {
                    "file_path": "",
                    "reference_id": "hg38",
                    "min_base_quality": 20,
                    "min_mapping_quality": 20,
                    "stand_call_conf": 30.0,
                },
            },
            "step_genomics_cnv": {
                "name": "拷贝数变异分析",
                "description": "深度与 HMM 检测 CNV（GATK CNV）",
                "tool_id": "genomics_cnv_calling",
                "default_params": {"file_path": "", "cnv_bin_width": 1000},
            },
            "step_genomics_sv": {
                "name": "结构变异分析",
                "description": "Split-read / discordant pairs（Manta）",
                "tool_id": "genomics_sv_calling",
                "default_params": {"file_path": "", "min_sv_len": 50},
            },
            "step_genomics_vqsr": {
                "name": "变异质量校准、过滤与位点标准化",
                "description": "VQSR/VQSLOD 过滤；常与 bcftools norm 左对齐、多等位拆分同管线执行",
                "tool_id": "genomics_vqsr_filtering",
                "default_params": {
                    "file_path": "",
                    "reference_id": "hg38",
                    "tranche_sensitivity": 99.0,
                },
            },
            "step_genomics_anno": {
                "name": "多维多态性注释",
                "description": "转录本映射与蛋白影响（VEP / SnpEff）",
                "tool_id": "genomics_variant_annotation",
                "default_params": {"file_path": "", "pick_allele": "canonical"},
            },
            "step_genomics_acmg": {
                "name": "临床致病性评级",
                "description": "dbSNP/ClinVar 与 ACMG/AMP 分级",
                "tool_id": "genomics_acmg_classification",
                "default_params": {"file_path": "", "pp2_ba1_threshold": 0.05},
            },
            "step_genomics_report": {
                "name": "全景变异报告生成",
                "description": "汇总变异、QC 与注释的标准化 PDF 报告",
                "tool_id": "genomics_clinical_reporting",
                "default_params": {"file_path": "", "report_style": "germline_v1"},
            },
        }
        if step_id not in m:
            raise ValueError(f"未知的基因组学步骤: {step_id}")
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

            if sid == "step_genomics_raw_qc":
                params["file_path"] = placeholder
                params["input_dir"] = ""
            elif sid == "step_genomics_read_trim":
                params["file_path"] = "<step_genomics_raw_qc>"
                params["input_dir"] = ""
            else:
                deps = dag.get(sid, [])
                dep = deps[-1] if deps else None
                params["file_path"] = f"<{dep}>" if dep else placeholder
                params.pop("input_dir", None)

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
