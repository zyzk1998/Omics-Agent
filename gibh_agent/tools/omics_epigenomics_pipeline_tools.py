"""
表观遗传组学全流程原子工具

与 EpigenomicsWorkflow 的 tool_id 一一对应。
首步对 FASTQ/FASTQ.GZ 做与基因组相同的真实流式质控统计。
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
from .omics_epigenomics_runner import (
    epigenomics_alignment_impl,
    epigenomics_peak_calling_impl,
    epigenomics_post_align_filtering_impl,
    epigenomics_shift_fragment_analysis_impl,
)
from .omics_derived_analysis import (
    epigenomics_multiomics_markdown,
    fastq_bundle_from_path,
)
from .omics_pipeline_env import modal_tool_with_degradation, write_temp_mock_artifact
from .omics_mock_ui import (
    IMG_CHIP_WORKFLOW,
    IMG_DNA_HELIX,
    IMG_SEQUENCE_LOGO,
    IMG_VOLCANO,
    attach_visual_contract,
    simple_rows_table,
)
from .omics_real_io import compute_fastq_stats, format_fastq_qc_markdown

logger = logging.getLogger(__name__)


def _epigenomics_modal_derived(
    file_path: str,
    prefix: str,
    msg: str,
    step: str,
    *,
    image_urls: List[str],
    tool_human_label: str,
    candidate_exes: Optional[List[str]],
    markdown_sim: str,
    table_data: Dict[str, Any],
    tool_id: str = "",
) -> Dict[str, Any]:
    """禁止 FASTQ 读段代理冒充 Peak/Motif 等真结果；统一走 modal（缺 CLI 则 error）。"""
    return modal_tool_with_degradation(
        prefix,
        msg,
        markdown_sim=markdown_sim,
        image_urls=image_urls,
        table_data=table_data,
        tool_human_label=tool_human_label,
        candidate_exes=candidate_exes,
        real_runner=None,
        tool_id=tool_id or None,
    )


def _mock(
    prefix: str,
    msg: str,
    *,
    md: Optional[str] = None,
    image_urls: Optional[List[str]] = None,
    table_data: Optional[Dict[str, Any]] = None,
) -> Dict[str, Any]:
    fd, path = tempfile.mkstemp(prefix=f"{prefix}_", suffix=".mock.txt")
    os.close(fd)
    with open(path, "w", encoding="utf-8") as fh:
        fh.write(msg + "\n")
    out: Dict[str, Any] = {
        "status": "success",
        "message": msg,
        "output_path": path,
        "file_path": path,
    }
    md_body = md or f"### 表观组学步骤（Mock）\n\n{msg}\n"
    imgs = image_urls if image_urls is not None else [IMG_DNA_HELIX]
    tbl = (
        table_data
        if table_data is not None
        else simple_rows_table(
            ("stage", "detail"),
            ({"stage": prefix, "detail": msg[:200]},),
        )
    )
    return attach_visual_contract(
        out, markdown=md_body, image_urls=imgs, table_data=tbl
    )


@registry.register(
    name="epigenomics_raw_qc_trimming",
    description="染色质测序 FASTQ 原始质控统计（真实 reads/GC/读长）；接头修剪可后续接 TrimGalore/fastp",
    category="Epigenomics",
    output_type="file_path",
)
@safe_tool_execution
def epigenomics_raw_qc_trimming(
    file_path: str = "", input_dir: str = "", trim_quality_threshold: int = 20
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
            "message": f"表观首步当前仅支持 FASTQ/FASTQ.GZ，收到: {path}",
        }

    stats = compute_fastq_stats(path)
    md = format_fastq_qc_markdown(stats).replace(
        "原始测序数据质控（真实统计）", "ATAC/ChIP 原始 FASTQ 质控（真实统计）"
    )

    safe_base = re.sub(r"[^\w.\-]+", "_", os.path.basename(path))[:120]
    qdir = os.path.abspath(os.getenv("RESULTS_DIR", os.path.join(os.getcwd(), "results")))
    qdir = os.path.join(qdir, "omics_qc_reports")
    os.makedirs(qdir, exist_ok=True)
    json_path = os.path.join(qdir, f"epigenomics_raw_qc_{safe_base}.json")
    try:
        with open(json_path, "w", encoding="utf-8") as jf:
            json.dump({"tool_id": "epigenomics_raw_qc_trimming", "qc_metrics": stats}, jf, ensure_ascii=False, indent=2)
    except OSError as exc:
        logger.warning("QC JSON 写入失败: %s", exc)
        json_path = ""

    out = {
        "status": "success",
        "message": f"表观原始质控：{stats['n_reads']} reads，GC {stats['gc_percent']}%",
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
        tool_id="epigenomics_raw_qc_trimming",
    )


@registry.register(
    name="epigenomics_alignment",
    description="ATAC/ChIP 读段比对（Bowtie2：threads、--mp；缺索引时降级仿真）",
    category="Epigenomics",
    output_type="file_path",
)
@safe_tool_execution
def epigenomics_alignment(
    file_path: str = "",
    reference_id: str = "hg38",
    threads: int = 8,
    mismatch_penalty: int = 4,
) -> Dict[str, Any]:
    return epigenomics_alignment_impl(
        file_path,
        reference_id,
        threads=threads,
        mismatch_penalty=mismatch_penalty,
    )


@registry.register(
    name="epigenomics_post_align_filtering",
    description="MAPQ 过滤、chrM 去除、PCR duplicate（Mock）",
    category="Epigenomics",
    output_type="file_path",
)
@safe_tool_execution
def epigenomics_post_align_filtering(
    file_path: str = "", min_mapq: int = 30
) -> Dict[str, Any]:
    return epigenomics_post_align_filtering_impl(file_path)


@registry.register(
    name="epigenomics_shift_fragment_analysis",
    description="ATAC-seq Tn5 移位校正与片段分布（Mock）",
    category="Epigenomics",
    output_type="file_path",
)
@safe_tool_execution
def epigenomics_shift_fragment_analysis(
    file_path: str = "", shift_correction_bp: int = 4
) -> Dict[str, Any]:
    return epigenomics_shift_fragment_analysis_impl(file_path)


@registry.register(
    name="epigenomics_peak_calling",
    description="MACS2 callpeak（-q qvalue、broad peak；缺对照 BAM 时降级仿真）",
    category="Epigenomics",
    output_type="file_path",
)
@safe_tool_execution
def epigenomics_peak_calling(
    file_path: str = "",
    qvalue_threshold: float = 0.05,
    broad_peak: bool = False,
) -> Dict[str, Any]:
    return epigenomics_peak_calling_impl(
        file_path,
        qvalue_threshold=qvalue_threshold,
        broad_peak=broad_peak,
    )


@registry.register(
    name="epigenomics_reproducibility_idr",
    description="FRiP 与 IDR 生物学重复一致性（Mock）",
    category="Epigenomics",
    output_type="file_path",
)
@safe_tool_execution
def epigenomics_reproducibility_idr(
    file_path: str = "", idr_threshold: float = 0.05
) -> Dict[str, Any]:
    return _epigenomics_modal_derived(
        file_path,
        "e_idr",
        "Reproducibility / IDR",
        "idr",
        image_urls=[IMG_VOLCANO],
        tool_human_label="IDR / chip-seq-tools",
        candidate_exes=["idr"],
        markdown_sim=(
            "### FRiP / IDR（仿真）\n\n"
            "- **FRiP**：0.078\n"
            "- IDR 合并峰：**21,340**\n"
        ),
        table_data=simple_rows_table(
            ("replicate_pair", "idr_peaks"),
            [
                {"replicate_pair": "rep1_vs_rep2", "idr_peaks": "21340"},
            ],
        ),
        tool_id="epigenomics_reproducibility_idr",
    )


@registry.register(
    name="epigenomics_consensus_peak_counting",
    description="Consensus peaks 与计数矩阵（Mock）",
    category="Epigenomics",
    output_type="file_path",
)
@safe_tool_execution
def epigenomics_consensus_peak_counting(
    file_path: str = "", merge_distance_bp: int = 100
) -> Dict[str, Any]:
    return _epigenomics_modal_derived(
        file_path,
        "e_cons",
        "Consensus peaks & count matrix",
        "consensus",
        image_urls=[IMG_CHIP_WORKFLOW],
        tool_human_label="bedtools / featureCounts / DiffBind",
        candidate_exes=["bedtools", "featureCounts"],
        markdown_sim=(
            "### Consensus peaks 与计数矩阵（仿真）\n\n"
            "- 共识峰：**19,882**；矩阵维度：峰 × 样本 **19882 × 12**\n"
        ),
        table_data=simple_rows_table(
            ("sample", "fragments_in_peaks"),
            [
                {"sample": "S01", "fragments_in_peaks": "4.2e6"},
                {"sample": "S02", "fragments_in_peaks": "3.9e6"},
            ],
        ),
        tool_id="epigenomics_consensus_peak_counting",
    )


@registry.register(
    name="epigenomics_peak_annotation",
    description="Peak 映射至启动子/增强子等基因组元件（Mock）",
    category="Epigenomics",
    output_type="file_path",
)
@safe_tool_execution
def epigenomics_peak_annotation(
    file_path: str = "", promoter_window_bp: int = 2000
) -> Dict[str, Any]:
    return _epigenomics_modal_derived(
        file_path,
        "e_panno",
        "Peak annotation",
        "panno",
        image_urls=[IMG_CHIP_WORKFLOW],
        tool_human_label="HOMER / ChIPseeker",
        candidate_exes=["annotatePeaks.pl", "Rscript"],
        markdown_sim=(
            "### Peak 基因组注释（仿真）\n\n"
            "| 区域类型 | 占比 |\n"
            "|----------|------|\n"
            "| Promoter | 31% |\n"
            "| Distal | 54% |\n"
        ),
        table_data=simple_rows_table(
            ("annotation", "fraction"),
            [
                {"annotation": "promoter", "fraction": "0.31"},
                {"annotation": "distal_intergenic", "fraction": "0.54"},
            ],
        ),
        tool_id="epigenomics_peak_annotation",
    )


@registry.register(
    name="epigenomics_diff_accessibility",
    description="DESeq2/DiffBind 差异开放/结合（Mock）",
    category="Epigenomics",
    output_type="file_path",
)
@safe_tool_execution
def epigenomics_diff_accessibility(
    file_path: str = "", padj_cutoff: float = 0.05
) -> Dict[str, Any]:
    return _epigenomics_modal_derived(
        file_path,
        "e_diff",
        "Differential accessibility",
        "diff",
        image_urls=[IMG_VOLCANO],
        tool_human_label="DESeq2 / DiffBind（R）",
        candidate_exes=["Rscript"],
        markdown_sim=(
            "### 差异开放分析（仿真 / DESeq2）\n\n"
            "| 对比 | 上调峰 | 下调峰 |\n"
            "|------|--------|--------|\n"
            "| Trt vs Ctrl | 3,420 | 2,881 |\n"
        ),
        table_data=simple_rows_table(
            ("peak_id", "log2FC"),
            [
                {"peak_id": "chr1:12345-12890", "log2FC": "2.41"},
                {"peak_id": "chr17:41234-41890", "log2FC": "-1.88"},
            ],
        ),
        tool_id="epigenomics_diff_accessibility",
    )


@registry.register(
    name="epigenomics_motif_discovery",
    description="JASPAR 扫描与 de novo Motif（Mock）",
    category="Epigenomics",
    output_type="file_path",
)
@safe_tool_execution
def epigenomics_motif_discovery(
    file_path: str = "", motif_e_value: float = 1e-4
) -> Dict[str, Any]:
    return _epigenomics_modal_derived(
        file_path,
        "e_motif",
        "Motif discovery",
        "motif",
        image_urls=[IMG_SEQUENCE_LOGO, IMG_CHIP_WORKFLOW],
        tool_human_label="MEME / Homer motif tools",
        candidate_exes=["meme", "findMotifsGenome.pl"],
        markdown_sim=(
            "### Motif 发现（仿真 / MEME-JASPAR）\n\n"
            "| Motif | E-value |\n"
            "|-------|--------|\n"
            "| AP-1 | 1e-18 |\n"
        ),
        table_data=simple_rows_table(
            ("motif", "e_value"),
            [
                {"motif": "AP-1", "e_value": "1e-18"},
                {"motif": "CTCF", "e_value": "4e-12"},
            ],
        ),
        tool_id="epigenomics_motif_discovery",
    )


@registry.register(
    name="epigenomics_tf_footprinting",
    description="单碱基分辨率 TF 足迹（Mock）",
    category="Epigenomics",
    output_type="file_path",
)
@safe_tool_execution
def epigenomics_tf_footprinting(
    file_path: str = "", nuc_resolution_bp: int = 10
) -> Dict[str, Any]:
    return _epigenomics_modal_derived(
        file_path,
        "e_ft",
        "TF footprinting",
        "footprint",
        image_urls=[IMG_DNA_HELIX],
        tool_human_label="HINT-ATAC / TOBIAS",
        candidate_exes=["tobias", "rgt-hint"],
        markdown_sim=(
            "### TF 足迹（仿真）\n\n"
            "- CTCF 足迹深度：**0.82**（相对背景）\n"
        ),
        table_data=simple_rows_table(
            ("tf", "footprint_score"),
            [
                {"tf": "CTCF", "footprint_score": "0.82"},
                {"tf": "NRF1", "footprint_score": "0.76"},
            ],
        ),
        tool_id="epigenomics_tf_footprinting",
    )


@registry.register(
    name="epigenomics_cis_regulatory_interactions",
    description="共开放性与增强子-启动子互作推断（Mock）",
    category="Epigenomics",
    output_type="file_path",
)
@safe_tool_execution
def epigenomics_cis_regulatory_interactions(
    file_path: str = "", max_distance_bp: int = 500000
) -> Dict[str, Any]:
    return _epigenomics_modal_derived(
        file_path,
        "e_cis",
        "Cis-regulatory interactions",
        "cis",
        image_urls=[IMG_CHIP_WORKFLOW],
        tool_human_label="HiC-Pro / Peak2Gene / cicero",
        candidate_exes=["python3"],
        markdown_sim=(
            "### 顺式调控互作（仿真）\n\n"
            "| Enhancer | Target promoter | Score |\n"
            "|----------|-----------------|-------|\n"
            "| chr17:42kb | BRCA1 TSS | 0.91 |\n"
        ),
        table_data=simple_rows_table(
            ("enhancer", "target_gene", "score"),
            [
                {"enhancer": "chr17:42kb", "target_gene": "BRCA1", "score": "0.91"},
            ],
        ),
        tool_id="epigenomics_cis_regulatory_interactions",
    )


@registry.register(
    name="epigenomics_multiomics_integration",
    description="联合 RNA-seq 的 GRN 整合（Mock）",
    category="Epigenomics",
    output_type="file_path",
)
@safe_tool_execution
def epigenomics_multiomics_integration(
    file_path: str = "",
    data_path: str = "",
    min_correlation: float = 0.3,
) -> Dict[str, Any]:
    out = _mock("e_multi", "Epigenome-transcriptome GRN placeholder")
    if data_path:
        out["data_path"] = data_path
    b = fastq_bundle_from_path(file_path)
    md_body, tbl_extra = epigenomics_multiomics_markdown(file_path, data_path, b)
    out["markdown"] = md_body
    out["omics_analysis_mode"] = "fastq_stream_proxy" if b else "placeholder"
    if b:
        out["proxy_metrics"] = b
    td: Dict[str, Any] = {
        "peaks_preview": [],
        "accessibility_matrix_hint": [],
        "file_path": file_path or "",
        "data_path": data_path or "",
    }
    if b:
        td["proxy_context"] = tbl_extra
    out["image_urls"] = [IMG_CHIP_WORKFLOW]
    out["json_url"] = None
    out["pdf_url"] = None
    out["html_url"] = None
    out["table_data"] = td
    out["data"] = {
        "report_stage": "epigenomics_multiomics_integration",
        "file_path": file_path or "",
        "data_path": data_path or "",
        "pdf_url": out["pdf_url"],
        "html_url": out["html_url"],
        "table_data": out["table_data"],
    }
    return attach_visual_contract(
        out,
        markdown=out["markdown"],
        image_urls=out["image_urls"],
        table_data=out["table_data"],
        tool_id="epigenomics_multiomics_integration",
    )
