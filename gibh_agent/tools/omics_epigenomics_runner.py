"""表观组学下游：Bowtie2 / samtools / MACS2 等 CLI 骨架与仿真降级。"""
from __future__ import annotations

import logging
import os
import shutil
import subprocess
import tempfile
from typing import Any, Dict, List, Optional

from .omics_mock_ui import (
    IMG_CHIP_WORKFLOW,
    IMG_DNA_HELIX,
    IMG_PEAK_GENOME_PIE,
    attach_visual_contract,
    simple_rows_table,
)
from .omics_derived_analysis import (
    epigenomics_proxy_alignment,
    epigenomics_proxy_peak,
    epigenomics_proxy_post_filter,
    epigenomics_proxy_shift,
    fastq_bundle_from_path,
)
from .omics_pipeline_env import (
    degraded_banner_md,
    resolve_bowtie2_index_prefix,
    run_or_degrade,
    smoke_subprocess,
    write_temp_mock_artifact,
)

logger = logging.getLogger(__name__)


def _contract(
    prefix: str,
    msg: str,
    *,
    markdown: str,
    image_urls: List[str],
    table_data: Dict[str, Any],
    extra: Optional[Dict[str, Any]] = None,
    tool_id: str = "",
) -> Dict[str, Any]:
    path = write_temp_mock_artifact(prefix, msg)
    out: Dict[str, Any] = {
        "status": "success",
        "message": msg,
        "output_path": path,
        "file_path": path,
    }
    if extra:
        out.update(extra)
    return attach_visual_contract(
        out,
        markdown=markdown,
        image_urls=image_urls,
        table_data=table_data,
        tool_id=tool_id or None,
    )


def _degraded_alignment(reference_id: str, file_path: str = "") -> Dict[str, Any]:
    b = fastq_bundle_from_path(file_path)
    if b:
        md, tbl, msg = epigenomics_proxy_alignment(b, reference_id)
        out = _contract(
            "e_aln",
            msg,
            markdown=md,
            image_urls=[IMG_CHIP_WORKFLOW],
            table_data=tbl,
            extra={"omics_analysis_mode": "fastq_stream_proxy", "proxy_metrics": b},
            tool_id="epigenomics_alignment",
        )
        out["bam_path"] = out["output_path"]
        return out
    md = degraded_banner_md("Bowtie2 / BWA 与索引前缀（GIBH_BOWTIE2_*）") + (
        "### ATAC/ChIP 比对（仿真）\n\n"
        f"- 参考：`{reference_id}`\n"
        "- 比对率（仿真）：**96.8%**\n"
        "- 线粒体占比（仿真）：**2.1%**\n"
    )
    out = _contract(
        "e_aln",
        f"Epigenomics alignment simulated ref={reference_id}",
        markdown=md,
        image_urls=[IMG_CHIP_WORKFLOW],
        table_data=simple_rows_table(
            ("metric", "value"),
            [
                {"metric": "mapped_rate", "value": "96.8%"},
                {"metric": "mitochondrial_frac", "value": "2.1%"},
            ],
        ),
        tool_id="epigenomics_alignment",
    )
    out["bam_path"] = out["output_path"]
    return out


def _try_bowtie2_alignment(file_path: str, reference_id: str) -> Optional[Dict[str, Any]]:
    bt2 = shutil.which("bowtie2")
    samtools = shutil.which("samtools")
    idx = resolve_bowtie2_index_prefix(reference_id)
    fq = (file_path or "").strip()
    if not (bt2 and samtools and idx and fq and os.path.isfile(fq)):
        return None
    low = fq.lower()
    if not (
        low.endswith(".fastq.gz")
        or low.endswith(".fq.gz")
        or low.endswith(".fastq")
        or low.endswith(".fq")
    ):
        return None
    work = tempfile.mkdtemp(prefix="e_bt2_")
    sam_p = os.path.join(work, "aligned.sam")
    bam_p = os.path.join(work, "aligned.bam")
    try:
        with open(sam_p, "wb") as sf:
            subprocess.run(
                [bt2, "-x", idx, "-U", fq, "--threads", "4"],
                stdout=sf,
                stderr=subprocess.PIPE,
                check=True,
                timeout=7200,
            )
        subprocess.run(
            [samtools, "view", "-bS", sam_p, "-o", bam_p],
            check=True,
            timeout=3600,
        )
        md = (
            "### ATAC/ChIP 比对（bowtie2 + samtools）\n\n"
            f"- **BAM**: `{bam_p}`\n"
        )
        out = {
            "status": "success",
            "message": "比对完成（bowtie2）",
            "output_path": bam_p,
            "file_path": bam_p,
            "bam_path": bam_p,
        }
        return attach_visual_contract(
            out,
            markdown=md,
            image_urls=[IMG_CHIP_WORKFLOW],
            table_data=simple_rows_table(
                ("artifact", "path"),
                [
                    {"artifact": "SAM", "path": sam_p},
                    {"artifact": "BAM", "path": bam_p},
                ],
            ),
            tool_id="epigenomics_alignment",
        )
    except (subprocess.CalledProcessError, OSError, subprocess.TimeoutExpired) as exc:
        logger.warning("bowtie2 alignment failed: %s", exc)
        return None


def epigenomics_alignment_impl(file_path: str = "", reference_id: str = "hg38") -> Dict[str, Any]:
    return run_or_degrade(
        lambda: _try_bowtie2_alignment(file_path, reference_id),
        lambda: _degraded_alignment(reference_id, file_path),
        ctx="epigenomics_alignment",
    )


def _degraded_post_filter(file_path: str = "") -> Dict[str, Any]:
    b = fastq_bundle_from_path(file_path)
    if b:
        md, tbl, msg = epigenomics_proxy_post_filter(b)
        return _contract(
            "e_post_filt",
            msg,
            markdown=md,
            image_urls=[IMG_CHIP_WORKFLOW],
            table_data=tbl,
            extra={"omics_analysis_mode": "fastq_stream_proxy", "proxy_metrics": b},
            tool_id="epigenomics_post_align_filtering",
        )
    md = degraded_banner_md("samtools view/filter 与排序 BAM") + (
        "### MAPQ 过滤与去重（仿真）\n\n"
        "| 项目 | 剩余 reads |\n|------|------------|\n"
        "| 去 chrM | 98.2% |\n"
        "| MAPQ≥30 | 96.5% |\n"
    )
    return _contract(
        "e_post_filt",
        "Post-align filtering (simulated)",
        markdown=md,
        image_urls=[IMG_CHIP_WORKFLOW],
        table_data=simple_rows_table(
            ("filter_stage", "fraction_remaining"),
            [
                {"filter_stage": "MAPQ>=30", "fraction_remaining": "0.965"},
                {"filter_stage": "dedup", "fraction_remaining": "0.822"},
            ],
        ),
        tool_id="epigenomics_post_align_filtering",
    )


def _try_post_filter(file_path: str) -> Optional[Dict[str, Any]]:
    samtools = shutil.which("samtools")
    bam = (file_path or "").strip()
    if not (samtools and bam and os.path.isfile(bam) and bam.lower().endswith(".bam")):
        return None
    work = tempfile.mkdtemp(prefix="e_sfilt_")
    filt_bam = os.path.join(work, "filtered.bam")
    try:
        subprocess.run(
            [
                samtools,
                "view",
                "-b",
                "-q",
                "30",
                bam,
                "-o",
                filt_bam,
            ],
            check=True,
            timeout=7200,
        )
        md = "### MAPQ 过滤（samtools view -q 30）\n\n" f"- **BAM**: `{filt_bam}`\n"
        return _contract(
            "e_post_filt",
            "samtools filter ok",
            markdown=md,
            image_urls=[IMG_CHIP_WORKFLOW],
            table_data=simple_rows_table(
                ("artifact", "path"),
                [{"artifact": "filtered_bam", "path": filt_bam}],
            ),
            tool_id="epigenomics_post_align_filtering",
        )
    except (subprocess.CalledProcessError, OSError, subprocess.TimeoutExpired) as exc:
        logger.warning("post filter failed: %s", exc)
        return None


def epigenomics_post_align_filtering_impl(file_path: str = "") -> Dict[str, Any]:
    return run_or_degrade(
        lambda: _try_post_filter(file_path),
        lambda: _degraded_post_filter(file_path),
        ctx="epigenomics_post_align_filtering",
    )


def _degraded_peak(file_path: str = "") -> Dict[str, Any]:
    b = fastq_bundle_from_path(file_path)
    if b:
        md, tbl, msg = epigenomics_proxy_peak(b)
        out = _contract(
            "e_peak",
            msg,
            markdown=md,
            image_urls=[IMG_PEAK_GENOME_PIE, IMG_CHIP_WORKFLOW],
            table_data=tbl,
            extra={"omics_analysis_mode": "fastq_stream_proxy", "proxy_metrics": b},
            tool_id="epigenomics_peak_calling",
        )
        out["bed_path"] = out["output_path"]
        return out
    md = degraded_banner_md("MACS2 / Genrich 与配对对照 BAM") + (
        "### Peak calling（仿真）\n\n"
        "| 指标 | 仿真值 |\n|------|--------|\n"
        "| Peaks | 28,420 |\n"
    )
    out = _contract(
        "e_peak",
        "Peak calling (simulated)",
        markdown=md,
        image_urls=[IMG_PEAK_GENOME_PIE, IMG_CHIP_WORKFLOW],
        table_data=simple_rows_table(
            ("chrom", "peaks"),
            [
                {"chrom": "chr1", "peaks": "4120"},
                {"chrom": "chr17", "peaks": "2888"},
            ],
        ),
        tool_id="epigenomics_peak_calling",
    )
    out["bed_path"] = out["output_path"]
    return out


def _try_macs2_peak(file_path: str) -> Optional[Dict[str, Any]]:
    macs2 = shutil.which("macs2")
    if not macs2:
        return None
    if not smoke_subprocess([macs2, "--help"], timeout=120):
        return None
    logger.info("MACS2 可用；完整 peak calling 需 treatment/control BAM；仿真降级。")
    return None


def epigenomics_peak_calling_impl(file_path: str = "") -> Dict[str, Any]:
    return run_or_degrade(
        lambda: _try_macs2_peak(file_path),
        lambda: _degraded_peak(file_path),
        ctx="epigenomics_peak_calling",
    )


def _degraded_shift(file_path: str = "") -> Dict[str, Any]:
    b = fastq_bundle_from_path(file_path)
    if b:
        md, tbl, msg = epigenomics_proxy_shift(b)
        return _contract(
            "e_shift",
            msg,
            markdown=md,
            image_urls=[IMG_DNA_HELIX],
            table_data=tbl,
            extra={"omics_analysis_mode": "fastq_stream_proxy", "proxy_metrics": b},
            tool_id="epigenomics_shift_fragment_analysis",
        )
    md = degraded_banner_md("alignmentSieve/deepTools 等片段统计工具") + (
        "### Tn5 移位校正与片段分布（仿真）\n\n"
        "- 主峰间距：**~200 bp**\n"
    )
    return _contract(
        "e_shift",
        "Shift analysis (simulated)",
        markdown=md,
        image_urls=[IMG_DNA_HELIX],
        table_data=simple_rows_table(
            ("fragment_bp", "density_peak"),
            [
                {"fragment_bp": "120", "density_peak": "low"},
                {"fragment_bp": "200", "density_peak": "high"},
            ],
        ),
        tool_id="epigenomics_shift_fragment_analysis",
    )


def epigenomics_shift_fragment_analysis_impl(file_path: str = "") -> Dict[str, Any]:
    return run_or_degrade(
        lambda: None,
        lambda: _degraded_shift(file_path),
        ctx="epigenomics_shift_fragment_analysis",
    )
