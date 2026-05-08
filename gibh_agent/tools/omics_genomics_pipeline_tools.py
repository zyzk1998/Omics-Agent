"""
基因组学全流程原子工具

与 GenomicsWorkflow 的 tool_id 一一对应。
首步 `genomics_raw_qc`：对 FASTQ/FASTQ.GZ 做真实流式统计（reads、GC、读长）。
下游步骤在 `omics_genomics_runner` 中提供真实 CLI 骨架；缺依赖时返回仿真可视化，避免 Task 崩溃。
"""
from __future__ import annotations

import json
import logging
import os
import re
import tempfile
from typing import Any, Dict, List, Optional

from ..core.tool_registry import registry
from ..core.utils import safe_tool_execution
from .omics_genomics_runner import (
    genomics_acmg_classification_impl,
    genomics_alignment_impl,
    genomics_bqsr_impl,
    genomics_cnv_calling_impl,
    genomics_germline_calling_impl,
    genomics_mark_duplicates_impl,
    genomics_read_trimming_impl,
    genomics_sv_calling_impl,
    genomics_variant_annotation_impl,
    genomics_vqsr_filtering_impl,
)
from .omics_mock_ui import (
    IMG_DNA_HELIX,
    IMG_GENOME_BROWSER,
    IMG_MANHATTAN,
    IMG_VOLCANO,
    attach_visual_contract,
    simple_rows_table,
)
from .omics_derived_analysis import clinical_genomics_sections
from .omics_real_io import compute_fastq_stats, format_fastq_qc_markdown

logger = logging.getLogger(__name__)


def _mock(
    prefix: str,
    message: str,
    *,
    md: Optional[str] = None,
    image_urls: Optional[List[str]] = None,
    table_data: Optional[Dict[str, Any]] = None,
) -> Dict[str, Any]:
    fd, path = tempfile.mkstemp(prefix=f"{prefix}_", suffix=".mock.txt")
    os.close(fd)
    try:
        with open(path, "w", encoding="utf-8") as fh:
            fh.write(message + "\n")
    except OSError as exc:
        logger.warning("mock artifact write failed: %s", exc)
    out: Dict[str, Any] = {
        "status": "success",
        "message": message,
        "output_path": path,
        "file_path": path,
    }
    md_body = md or f"### 步骤概要（Mock）\n\n{message}\n"
    imgs = image_urls if image_urls is not None else [IMG_DNA_HELIX]
    tbl = (
        table_data
        if table_data is not None
        else simple_rows_table(
            ("field", "value"),
            ({"field": "mock_stage", "value": message[:240]},),
        )
    )
    return attach_visual_contract(
        out, markdown=md_body, image_urls=imgs, table_data=tbl
    )


@registry.register(
    name="genomics_raw_qc",
    description="FASTQ 全面质量、GC 偏好与污染源评估（真实 reads/GC/读长统计；可接 FastQC/MultiQC 扩展）",
    category="Genomics",
    output_type="file_path",
)
@safe_tool_execution
def genomics_raw_qc(
    file_path: str = "", input_dir: str = "", min_read_length: int = 36
) -> Dict[str, Any]:
    src = (file_path or input_dir or "").strip()
    if not src:
        return {"status": "error", "message": "必须提供 file_path 或 input_dir"}
    path = os.path.abspath(src)
    if not os.path.isfile(path):
        return {"status": "error", "message": f"输入不是可读文件: {path}"}
    low = path.lower()
    if not (
        low.endswith(".fastq.gz")
        or low.endswith(".fq.gz")
        or low.endswith(".fastq")
        or low.endswith(".fq")
    ):
        return {
            "status": "error",
            "message": f"genomics_raw_qc 当前仅支持 FASTQ/FASTQ.GZ，收到: {path}",
        }

    stats = compute_fastq_stats(path)
    md = format_fastq_qc_markdown(stats)

    safe_base = re.sub(r"[^\w.\-]+", "_", os.path.basename(path))[:120]
    qdir = os.path.abspath(os.getenv("RESULTS_DIR", os.path.join(os.getcwd(), "results")))
    qdir = os.path.join(qdir, "omics_qc_reports")
    os.makedirs(qdir, exist_ok=True)
    json_path = os.path.join(qdir, f"genomics_raw_qc_{safe_base}.json")
    try:
        with open(json_path, "w", encoding="utf-8") as jf:
            json.dump({"tool_id": "genomics_raw_qc", "qc_metrics": stats}, jf, ensure_ascii=False, indent=2)
    except OSError as exc:
        logger.warning("QC JSON 写入失败（不影响统计结果）: %s", exc)
        json_path = ""

    out = {
        "status": "success",
        "message": f"原始质控完成：{stats['n_reads']} reads，GC {stats['gc_percent']}%",
        "markdown": md,
        "qc_metrics": stats,
        "summary": f"{stats['n_reads']} reads，GC {stats['gc_percent']}%",
        "output_path": json_path or path,
        "file_path": path,
    }
    return attach_visual_contract(
        out,
        markdown=md,
        image_urls=[IMG_DNA_HELIX],
        table_data=simple_rows_table(
            ("metric", "value"),
            [
                {"metric": "total_reads", "value": str(stats["n_reads"])},
                {"metric": "gc_percent", "value": str(stats["gc_percent"])},
                {
                    "metric": "mean_read_length",
                    "value": str(stats.get("mean_read_length", "")),
                },
            ],
        ),
        tool_id="genomics_raw_qc",
    )


@registry.register(
    name="genomics_read_trimming",
    description="动态切除低质量末端、接头过滤与短读段过滤（Mock；fastp）",
    category="Genomics",
    output_type="file_path",
)
@safe_tool_execution
def genomics_read_trimming(
    file_path: str = "", input_dir: str = "", quality_cutoff: int = 20
) -> Dict[str, Any]:
    return genomics_read_trimming_impl(file_path, input_dir)


@registry.register(
    name="genomics_alignment",
    description="高精度比对生成 BAM（Mock；BWA-MEM / DRAGEN）",
    category="Genomics",
    output_type="file_path",
)
@safe_tool_execution
def genomics_alignment(
    file_path: str = "",
    reference_id: str = "hg38",
    mismatch_penalty: int = 4,
    gap_open_penalty: int = 6,
) -> Dict[str, Any]:
    return genomics_alignment_impl(file_path, reference_id)


@registry.register(
    name="genomics_mark_duplicates",
    description="坐标排序与 PCR 重复标记（Mock；Picard MarkDuplicates）",
    category="Genomics",
    output_type="file_path",
)
@safe_tool_execution
def genomics_mark_duplicates(
    file_path: str = "", optical_duplicate_distance: int = 2500
) -> Dict[str, Any]:
    return genomics_mark_duplicates_impl(file_path)


@registry.register(
    name="genomics_bqsr",
    description="碱基质量重校准 BQSR（Mock；GATK BaseRecalibrator）",
    category="Genomics",
    output_type="file_path",
)
@safe_tool_execution
def genomics_bqsr(file_path: str = "", bqsr_max_cycles: int = 2) -> Dict[str, Any]:
    return genomics_bqsr_impl(file_path)


@registry.register(
    name="genomics_germline_calling",
    description="胚系 SNP/Indel（Mock；HaplotypeCaller / DeepVariant）",
    category="Genomics",
    output_type="file_path",
)
@safe_tool_execution
def genomics_germline_calling(
    file_path: str = "",
    min_base_quality: int = 20,
    min_mapping_quality: int = 20,
) -> Dict[str, Any]:
    return genomics_germline_calling_impl(file_path)


@registry.register(
    name="genomics_cnv_calling",
    description="拷贝数变异（Mock；GATK CNV / depth+HMM）",
    category="Genomics",
    output_type="file_path",
)
@safe_tool_execution
def genomics_cnv_calling(file_path: str = "", cnv_bin_width: int = 1000) -> Dict[str, Any]:
    return genomics_cnv_calling_impl(file_path)


@registry.register(
    name="genomics_sv_calling",
    description="结构变异（Mock；Manta split-read/discordant）",
    category="Genomics",
    output_type="file_path",
)
@safe_tool_execution
def genomics_sv_calling(file_path: str = "", min_sv_len: int = 50) -> Dict[str, Any]:
    return genomics_sv_calling_impl(file_path)


@registry.register(
    name="genomics_vqsr_filtering",
    description="VQSR/VQSLOD 过滤 + 位点标准化（bcftools norm 等，常与过滤同管线；Mock）",
    category="Genomics",
    output_type="file_path",
)
@safe_tool_execution
def genomics_vqsr_filtering(
    file_path: str = "", tranche_sensitivity: float = 99.0
) -> Dict[str, Any]:
    return genomics_vqsr_filtering_impl(file_path)


@registry.register(
    name="genomics_variant_annotation",
    description="VEP/SnpEff 多维注释（Mock）",
    category="Genomics",
    output_type="file_path",
)
@safe_tool_execution
def genomics_variant_annotation(
    file_path: str = "", pick_allele: str = "canonical"
) -> Dict[str, Any]:
    return genomics_variant_annotation_impl(file_path)


@registry.register(
    name="genomics_acmg_classification",
    description="ACMG/AMP 致病性分级（Mock；ClinVar/dbSNP）",
    category="Genomics",
    output_type="file_path",
)
@safe_tool_execution
def genomics_acmg_classification(
    file_path: str = "", pp2_ba1_threshold: float = 0.05
) -> Dict[str, Any]:
    return genomics_acmg_classification_impl(file_path)


@registry.register(
    name="genomics_clinical_reporting",
    description="基因组临床摘要报告（承接上游 qc_metrics；PDF/HTML 由外部管线挂载）",
    category="Genomics",
    output_type="file_path",
)
@safe_tool_execution
def genomics_clinical_reporting(
    file_path: str = "",
    pipeline_metrics_json: str = "",
    report_style: str = "germline_v1",
) -> Dict[str, Any]:
    qc_bundle: Dict[str, Any] = {}
    if pipeline_metrics_json:
        try:
            parsed = json.loads(pipeline_metrics_json)
            if isinstance(parsed, dict):
                qc_bundle = parsed
        except json.JSONDecodeError:
            logger.warning("genomics_clinical_reporting: pipeline_metrics_json 解析失败")

    primary_qc: Optional[Dict[str, Any]] = None
    if qc_bundle:
        primary_qc = next(iter(qc_bundle.values()), None)

    reads_txt = "N/A"
    gc_txt = "N/A"
    if isinstance(primary_qc, dict):
        reads_txt = str(primary_qc.get("n_reads", "N/A"))
        gc_txt = str(primary_qc.get("gc_percent", "N/A"))

    md_lines = [
        "## 基因组学分析摘要（结构化草案）",
        "",
        "### 1. 测序质控概览",
    ]
    if primary_qc:
        md_lines.extend(
            [
                f"- **总 Reads 数**: {reads_txt}",
                f"- **GC 含量 (%)**: {gc_txt}",
                f"- **平均读长**: {primary_qc.get('mean_read_length', 'N/A')}",
                f"- **输入文件**: `{primary_qc.get('input_path', file_path or '')}`",
            ]
        )
    else:
        md_lines.append(
            "- **错误**：未收到上游 `genomics_raw_qc` 的 `qc_metrics`，无法填写真实质控表。**禁止捏造 Reads/GC。**"
        )
    md_lines.append(f"- **报告模板**: `{report_style}`")

    sec2, sec3, sec4 = clinical_genomics_sections(primary_qc, file_path or "")
    md_lines.extend(
        [
            "",
            "### 2. 变异检测统计",
            sec2.rstrip(),
            "",
            "### 3. 临床致病性变异",
            sec3.rstrip(),
            "",
            "### 4. 结论与建议",
            sec4.rstrip(),
            "",
            "**禁止事项**：本报告块不得使用代谢组学语境（PCA、VIP、代谢物、LC-MS 等）。",
        ]
    )
    md_body = "\n".join(md_lines)

    out = _mock("g_report", "Clinical reporting stage (structured draft)")
    out["report_path"] = out["output_path"]
    out["pdf_url"] = None
    out["html_url"] = None
    out["markdown"] = md_body
    out["qc_metrics"] = primary_qc or {}
    out["pipeline_qc_bundle"] = qc_bundle
    out["summary"] = (
        f"Reads={reads_txt}, GC%={gc_txt}" if primary_qc else "缺少上游 qc_metrics"
    )
    out["image_urls"] = [IMG_DNA_HELIX]
    out["json_url"] = None
    out["table_data"] = {
        "variants_preview": [],
        "qc_metrics": primary_qc or {},
        "pipeline_qc_bundle": qc_bundle,
        "ingress_file_path": file_path or "",
    }
    out["data"] = {
        "report_stage": "genomics_clinical_reporting",
        "ingress_file_path": file_path or "",
        "pdf_url": out["pdf_url"],
        "html_url": out["html_url"],
        "table_data": out["table_data"],
        "qc_metrics": primary_qc or {},
        "pipeline_qc_bundle": qc_bundle,
    }
    return attach_visual_contract(
        out,
        markdown=md_body,
        image_urls=[IMG_DNA_HELIX, IMG_MANHATTAN],
        table_data=out["table_data"],
        tool_id="genomics_clinical_reporting",
    )
