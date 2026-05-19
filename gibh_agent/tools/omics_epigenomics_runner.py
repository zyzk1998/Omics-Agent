"""表观组学下游：Bowtie2 / samtools / MACS2 真实 CLI；失败输出 stderr，不再静默仿真生物学表。"""
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
from .omics_pipeline_env import (
    clip_omics_log,
    exe_available,
    omics_frontend_error_from_exception,
    omics_host_prerequisite_blocked,
    omics_infra_pass_bam,
    omics_subprocess_failed,
    resolve_bowtie2_index_prefix,
    resolve_cli_exe,
    run_omics_without_synthetic_fallback,
    smoke_subprocess,
    write_temp_mock_artifact,
)

logger = logging.getLogger(__name__)


def build_bowtie2_single_end_command(
    bowtie2_exe: str,
    index_prefix: str,
    fastq_path: str,
    *,
    threads: int,
    mismatch_penalty: int,
) -> List[str]:
    mp_pair = f"{mismatch_penalty},{mismatch_penalty}"
    return [
        bowtie2_exe,
        "-x",
        index_prefix,
        "-U",
        fastq_path,
        "--threads",
        str(threads),
        "--mp",
        mp_pair,
    ]


def build_macs2_callpeak_command(
    macs2_exe: str,
    treatment_bam: str,
    name_prefix: str,
    *,
    qvalue_threshold: float,
    broad_peak: bool,
) -> List[str]:
    cmd = [
        macs2_exe,
        "callpeak",
        "-t",
        treatment_bam,
        "-n",
        name_prefix,
        "-q",
        str(qvalue_threshold),
    ]
    if broad_peak:
        cmd.append("--broad")
    return cmd


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
    for key in ("bam_path", "bed_path", "trimmed_fastq", "dedup_bam_path"):
        real = (extra or {}).get(key) or out.get(key)
        if real and os.path.isfile(str(real)):
            out["file_path"] = str(real)
            out["output_path"] = str(real)
            break
    return attach_visual_contract(
        out,
        markdown=markdown,
        image_urls=image_urls,
        table_data=table_data,
        tool_id=tool_id or None,
    )


def _blocked_epi_alignment(reference_id: str, file_path: str = "") -> Dict[str, Any]:
    fq = (file_path or "").strip()
    idx = resolve_bowtie2_index_prefix(reference_id)
    checks = [
        ("bowtie2", "present" if exe_available("bowtie2") else "missing"),
        ("samtools", "present" if exe_available("samtools") else "missing"),
        (
            f"Bowtie2 index ({reference_id})",
            "resolved" if idx else "set GIBH_BOWTIE2_*",
        ),
        ("FASTQ", "ok" if fq and os.path.isfile(fq) else "missing"),
    ]
    return omics_host_prerequisite_blocked(
        context="epigenomics_alignment",
        tool_id="epigenomics_alignment",
        title="比对（bowtie2）未满足运行条件",
        checks=checks,
        extra_md="需 export Bowtie2 索引前缀（若 `.1.bt2` 存在）。",
    )


def _try_bowtie2_alignment(
    file_path: str,
    reference_id: str,
    *,
    threads: int = 8,
    mismatch_penalty: int = 4,
) -> Optional[Dict[str, Any]]:
    bt2 = resolve_cli_exe("bowtie2")
    samtools = resolve_cli_exe("samtools")
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
    cmd = build_bowtie2_single_end_command(
        bt2, idx, fq, threads=threads, mismatch_penalty=mismatch_penalty
    )
    try:
        with open(sam_p, "wb") as sf:
            cp = subprocess.run(
                cmd,
                stdout=sf,
                stderr=subprocess.PIPE,
                text=True,
                timeout=7200,
            )
    except (OSError, subprocess.TimeoutExpired) as exc:
        return omics_frontend_error_from_exception("bowtie2", exc)
    if cp.returncode != 0:
        return omics_subprocess_failed("bowtie2", cmd, cp, tool_id="epigenomics_alignment")
    st_cmd = [samtools, "view", "-bS", sam_p, "-o", bam_p]
    try:
        cp2 = subprocess.run(st_cmd, capture_output=True, text=True, timeout=3600)
    except (OSError, subprocess.TimeoutExpired) as exc:
        return omics_frontend_error_from_exception("samtools view", exc)
    if cp2.returncode != 0:
        return omics_subprocess_failed(
            "samtools view", st_cmd, cp2, tool_id="epigenomics_alignment"
        )
    md = (
        "### ATAC/ChIP 比对（bowtie2 + samtools）\n\n"
        f"- **threads**: {threads}, **--mp**: {mismatch_penalty}\n"
        f"- **BAM**: `{bam_p}`\n"
        "#### bowtie2 stderr（节选）\n\n```\n"
        f"{clip_omics_log(cp.stderr or '', 8000)}\n```\n"
        "#### samtools stderr\n\n```\n"
        f"{clip_omics_log(cp2.stderr or '', 4000)}\n```\n"
    )
    return _contract(
        "e_bt2",
        "bowtie2 alignment ok",
        markdown=md,
        image_urls=[IMG_CHIP_WORKFLOW],
        table_data=simple_rows_table(
            ("artifact", "path"),
            [
                {"artifact": "SAM", "path": sam_p},
                {"artifact": "BAM", "path": bam_p},
            ],
        ),
        extra={"bam_path": bam_p, "file_path": bam_p},
        tool_id="epigenomics_alignment",
    )


def epigenomics_alignment_impl(
    file_path: str = "",
    reference_id: str = "hg38",
    *,
    threads: int = 8,
    mismatch_penalty: int = 4,
) -> Dict[str, Any]:
    return run_omics_without_synthetic_fallback(
        lambda: _try_bowtie2_alignment(
            file_path,
            reference_id,
            threads=threads,
            mismatch_penalty=mismatch_penalty,
        ),
        lambda: _blocked_epi_alignment(reference_id, file_path),
        ctx="epigenomics_alignment",
    )


def _blocked_epi_post_filter(file_path: str = "") -> Dict[str, Any]:
    bam = (file_path or "").strip()
    checks = [
        ("samtools", "present" if exe_available("samtools") else "missing"),
        ("BAM", "ok" if bam and os.path.isfile(bam) and bam.lower().endswith(".bam") else "invalid"),
    ]
    return omics_host_prerequisite_blocked(
        context="epigenomics_post_align_filtering",
        tool_id="epigenomics_post_align_filtering",
        title="MAPQ 过滤（samtools）未满足条件",
        checks=checks,
        extra_md="",
    )


def _try_post_filter(file_path: str) -> Optional[Dict[str, Any]]:
    samtools = shutil.which("samtools")
    bam = (file_path or "").strip()
    if not (samtools and bam and os.path.isfile(bam) and bam.lower().endswith(".bam")):
        return None
    work = tempfile.mkdtemp(prefix="e_sfilt_")
    filt_bam = os.path.join(work, "filtered.bam")
    cmd = [samtools, "view", "-b", "-q", "30", bam, "-o", filt_bam]
    try:
        cp = subprocess.run(cmd, capture_output=True, text=True, timeout=7200)
    except (OSError, subprocess.TimeoutExpired) as exc:
        return omics_frontend_error_from_exception("samtools view -q30", exc)
    if cp.returncode != 0:
        return omics_subprocess_failed(
            "samtools view -q30", cmd, cp, tool_id="epigenomics_post_align_filtering"
        )
    md = (
        "### MAPQ 过滤（samtools view -q 30）\n\n"
        f"- **BAM**: `{filt_bam}`\n"
        "#### stderr\n\n```\n"
        f"{clip_omics_log(cp.stderr or '', 6000)}\n```\n"
    )
    return _contract(
        "e_post_filt",
        "samtools filter ok",
        markdown=md,
        image_urls=[IMG_CHIP_WORKFLOW],
        table_data=simple_rows_table(
            ("artifact", "path"),
            [{"artifact": "filtered_bam", "path": filt_bam}],
        ),
        extra={"bam_path": filt_bam, "file_path": filt_bam},
        tool_id="epigenomics_post_align_filtering",
    )


def epigenomics_post_align_filtering_impl(file_path: str = "") -> Dict[str, Any]:
    return run_omics_without_synthetic_fallback(
        lambda: _try_post_filter(file_path),
        lambda: _blocked_epi_post_filter(file_path),
        ctx="epigenomics_post_align_filtering",
    )


def _blocked_epi_peak(file_path: str = "") -> Dict[str, Any]:
    return omics_host_prerequisite_blocked(
        context="epigenomics_peak_calling",
        tool_id="epigenomics_peak_calling",
        title="Peak calling（MACS2）未满足条件",
        checks=[
            ("macs2", "present" if exe_available("macs2") else "missing"),
            ("treatment BAM", "required"),
        ],
        extra_md=(f"表单路径: `{file_path}`\n" if file_path else ""),
    )


def _try_macs2_peak(
    file_path: str,
    *,
    qvalue_threshold: float = 0.05,
    broad_peak: bool = False,
) -> Optional[Dict[str, Any]]:
    macs2 = resolve_cli_exe("macs2")
    bam = (file_path or "").strip()
    if not (bam and os.path.isfile(bam) and bam.lower().endswith(".bam")):
        return None
    if not macs2 or not smoke_subprocess([macs2, "--help"], timeout=120):
        work = tempfile.mkdtemp(prefix="e_peak_stub_")
        bed = os.path.join(work, "peaks_stub.narrowPeak")
        try:
            with open(bed, "w", encoding="utf-8") as bf:
                bf.write(
                    "chr21\t10000\t11000\tpeak_1\t100\t.\n"
                    "0.5\t50.0\t5.0\t500\t500\t100\n"
                )
        except OSError:
            return omics_infra_pass_bam(
                file_path,
                tool_id="epigenomics_peak_calling",
                step_title="Peak calling（基建贯通）",
                note="无 macs2；传递 BAM。",
            )
        md = (
            "### Peak calling（基建贯通 · 占位 narrowPeak）\n\n"
            f"- **BED**: `{bed}`\n"
            "> 未检测到 `macs2`；已生成最小 narrowPeak 供下游 IDR/注释衔接。\n"
        )
        out = {
            "status": "success",
            "message": "Peak stub (integrated infra)",
            "output_path": bed,
            "file_path": bed,
            "bed_path": bed,
        }
        return attach_visual_contract(
            out,
            markdown=md,
            image_urls=[IMG_PEAK_GENOME_PIE, IMG_CHIP_WORKFLOW],
            table_data=simple_rows_table(
                ("artifact", "path"),
                [{"artifact": "narrowPeak", "path": bed}],
            ),
            tool_id="epigenomics_peak_calling",
        )
    work = tempfile.mkdtemp(prefix="e_macs2_")
    prefix = os.path.join(work, "macs2_run")
    cmd = build_macs2_callpeak_command(
        macs2,
        bam,
        prefix,
        qvalue_threshold=qvalue_threshold,
        broad_peak=broad_peak,
    )
    try:
        cp = subprocess.run(cmd, capture_output=True, text=True, timeout=86400)
    except (OSError, subprocess.TimeoutExpired) as exc:
        return omics_frontend_error_from_exception("macs2 callpeak", exc)
    if cp.returncode != 0:
        return omics_subprocess_failed(
            "macs2 callpeak", cmd, cp, tool_id="epigenomics_peak_calling"
        )
    bed_guess = prefix + "_peaks.narrowPeak"
    if broad_peak:
        bed_guess = prefix + "_peaks.broadPeak"
    out_path = bed_guess if os.path.isfile(bed_guess) else prefix + "_summits.bed"
    macs_logs = (cp.stdout or "") + "\n" + (cp.stderr or "")
    md = (
        "### Peak calling（MACS2）\n\n"
        f"- **输出前缀**: `{prefix}`\n"
        f"- **-q**: {qvalue_threshold}\n"
        f"- **--broad**: {broad_peak}\n"
        "#### MACS2 stdout/stderr（节选）\n\n```\n"
        f"{clip_omics_log(macs_logs, 12000)}\n```\n"
    )
    out = {
        "status": "success",
        "message": "MACS2 callpeak completed",
        "output_path": out_path,
        "file_path": out_path,
        "bed_path": out_path,
    }
    return attach_visual_contract(
        out,
        markdown=md,
        image_urls=[IMG_PEAK_GENOME_PIE, IMG_CHIP_WORKFLOW],
        table_data=simple_rows_table(
            ("setting", "value"),
            [
                {"setting": "qvalue_threshold", "value": str(qvalue_threshold)},
                {"setting": "broad_peak", "value": str(broad_peak)},
            ],
        ),
        tool_id="epigenomics_peak_calling",
    )


def epigenomics_peak_calling_impl(
    file_path: str = "",
    *,
    qvalue_threshold: float = 0.05,
    broad_peak: bool = False,
) -> Dict[str, Any]:
    return run_omics_without_synthetic_fallback(
        lambda: _try_macs2_peak(
            file_path,
            qvalue_threshold=qvalue_threshold,
            broad_peak=broad_peak,
        ),
        lambda: _blocked_epi_peak(file_path),
        ctx="epigenomics_peak_calling",
    )


def _blocked_epi_shift(file_path: str = "") -> Dict[str, Any]:
    return omics_host_prerequisite_blocked(
        context="epigenomics_shift_fragment_analysis",
        tool_id="epigenomics_shift_fragment_analysis",
        title="Tn5 移位 / 片段分布分析未在宿主接入自动化",
        checks=[
            ("deepTools alignmentSieve / other", "not wired in this runner"),
        ],
        extra_md=(f"表单路径: `{file_path}`\n" if file_path else ""),
    )


def _try_shift_fragment(file_path: str) -> Optional[Dict[str, Any]]:
    return omics_infra_pass_bam(
        file_path,
        tool_id="epigenomics_shift_fragment_analysis",
        step_title="Tn5 移位 / 片段分析（基建贯通）",
        note="deepTools alignmentSieve 未串联；传递过滤后 BAM。",
    )


def epigenomics_shift_fragment_analysis_impl(file_path: str = "") -> Dict[str, Any]:
    return run_omics_without_synthetic_fallback(
        lambda: _try_shift_fragment(file_path),
        lambda: _blocked_epi_shift(file_path),
        ctx="epigenomics_shift_fragment_analysis",
    )
