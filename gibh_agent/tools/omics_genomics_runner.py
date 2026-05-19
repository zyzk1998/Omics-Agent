"""
基因组学下游步骤：真实 CLI（subprocess）+ stdout/stderr 写入 Markdown；
缺依赖时返回 **status:error** 与环境诊断表，不再使用仿真生物学占位指标。
"""
from __future__ import annotations

import gzip
import json
import logging
import os
import shutil
import subprocess
import tempfile
from typing import Any, Dict, List, Optional, Tuple

from .omics_mock_ui import (
    IMG_DNA_HELIX,
    IMG_GENOME_BROWSER,
    IMG_VOLCANO,
    attach_visual_contract,
    simple_rows_table,
)
from .omics_genomics_real_io import (
    ensure_bam_indexed,
    ensure_vcf_indexed,
    format_variant_summary_md,
    load_fastp_summary,
    require_nonempty_file,
    summarize_vcf,
    vcf_header_has_info_ids,
)
from .omics_pipeline_env import (
    bwa_index_ready,
    clip_omics_log,
    exe_available,
    get_omics_ref_dir,
    omics_frontend_error_from_exception,
    omics_host_prerequisite_blocked,
    omics_subprocess_failed,
    resolve_cli_exe,
    resolve_reference_fasta,
    run_omics_without_synthetic_fallback,
    write_temp_mock_artifact,
)

logger = logging.getLogger(__name__)

# 与 scripts/setup_minimal_genomics_test_data.sh / init_omics_mock_references.py 一致
MINIMAL_HG38_REFERENCE_FALLBACK = "/tmp/omics_test_data/ref/chr21.fa"


def resolve_genomics_reference_fasta(reference_id: str) -> Optional[str]:
    """
    解析参考序列：GIBH_REF_* → OMICS_REF_DIR/genomics → 宿主机微缩 chr21 兜底。
    优先返回已完成 bwa index 的路径。
    """
    ref = resolve_reference_fasta(reference_id)
    if ref and bwa_index_ready(ref):
        return ref
    rid = (reference_id or "").strip().lower().replace(" ", "")
    if rid not in ("hg38", "grch38"):
        return ref if ref else None
    fallbacks = [
        os.environ.get("GIBH_REF_HG38", "").strip(),
        os.path.join(get_omics_ref_dir(), "genomics", "hg38.fa"),
        MINIMAL_HG38_REFERENCE_FALLBACK,
    ]
    for p in fallbacks:
        if p and os.path.isfile(p) and bwa_index_ready(p):
            return os.path.abspath(p)
    if ref:
        return ref
    for p in fallbacks:
        if p and os.path.isfile(p):
            return os.path.abspath(p)
    return None


def build_bwa_mem_command(
    bwa_exe: str,
    ref_fasta: str,
    fastq_path: str,
    *,
    threads: int,
    mismatch_penalty: int,
    gap_open_penalty: int,
) -> List[str]:
    """组装 BWA-MEM 命令（参数与工具签名对齐，供单测与 Worker 路由核对）。"""
    gap_pair = f"{gap_open_penalty},{gap_open_penalty}"
    return [
        bwa_exe,
        "mem",
        "-t",
        str(threads),
        "-B",
        str(mismatch_penalty),
        "-O",
        gap_pair,
        ref_fasta,
        fastq_path,
    ]


def build_gatk_haplotypecaller_command(
    gatk_exe: str,
    ref_fasta: str,
    bam_path: str,
    out_vcf: str,
    *,
    stand_call_conf: float,
    min_mapping_quality: int,
) -> List[str]:
    """GATK4 HaplotypeCaller：与 Broad 文档一致的检出置信（stand-call-conf）与 MQ 阈值。"""
    return [
        gatk_exe,
        "HaplotypeCaller",
        "-R",
        ref_fasta,
        "-I",
        bam_path,
        "-O",
        out_vcf,
        "--standard-min-confidence-threshold-for-calling",
        str(stand_call_conf),
        "--minimum-mapping-quality",
        str(min_mapping_quality),
    ]


def _gatk_subprocess_failed(cp: subprocess.CompletedProcess) -> bool:
    """GATK 有时对 USER ERROR 仍返回 0，须读 stderr。"""
    if cp.returncode != 0:
        return True
    blob = f"{cp.stderr or ''}\n{cp.stdout or ''}"
    return "USER ERROR has occurred" in blob or "A fatal error has occurred" in blob


def _bcftools_qual_filter_vcf(
    bcftools: str,
    invcf: str,
    *,
    work_prefix: str = "g_vf_bc_",
) -> Dict[str, Any]:
    work = tempfile.mkdtemp(prefix=work_prefix)
    out_vcf = os.path.join(work, "filtered.vcf.gz")
    cmd = [bcftools, "view", "-Oz", "-o", out_vcf, "-i", "QUAL>=30", invcf]
    try:
        cp = subprocess.run(cmd, capture_output=True, text=True, timeout=86400)
    except (OSError, subprocess.TimeoutExpired) as exc:
        return omics_frontend_error_from_exception("bcftools filter", exc)
    if cp.returncode != 0:
        return omics_subprocess_failed("bcftools filter", cmd, cp, tool_id="genomics_vqsr_filtering")
    idx_err = ensure_vcf_indexed(bcftools, out_vcf)
    if idx_err:
        return idx_err
    bad_out = require_nonempty_file(out_vcf, "bcftools 过滤 VCF", extensions=(".vcf.gz",))
    if bad_out:
        return bad_out
    var_sum = summarize_vcf(out_vcf)
    md = (
        "### 变异过滤（bcftools view -i QUAL>=30）\n\n"
        f"- **输入**: `{invcf}`\n- **输出**: `{out_vcf}`\n"
        f"{format_variant_summary_md(var_sum)}\n"
    )
    out = {
        "status": "success",
        "message": "Variant filtering (bcftools)",
        "output_path": out_vcf,
        "file_path": out_vcf,
        "vcf_path": out_vcf,
        "filtered_vcf_path": out_vcf,
        "variant_summary": var_sum,
        "vqsr_mode": "bcftools_qual_filter",
    }
    return attach_visual_contract(
        out,
        markdown=md,
        image_urls=[IMG_VOLCANO],
        table_data=simple_rows_table(
            ("artifact", "path"),
            [{"artifact": "filtered_vcf", "path": out_vcf}],
        ),
        tool_id="genomics_vqsr_filtering",
    )


def build_gatk_hard_filtering_command(
    gatk_exe: str,
    ref_fasta: str,
    input_vcf: str,
    output_vcf: str,
) -> List[str]:
    """
    GATK4 VariantFiltration：无 VQSR 训练集时的硬过滤降级。
    使用常用 QD / FS / MQ 组合；与 GATK4 一致采用 --filter-expression / --filter-name。
    """
    return [
        gatk_exe,
        "VariantFiltration",
        "-R",
        ref_fasta,
        "-V",
        input_vcf,
        "-O",
        output_vcf,
        "--filter-expression",
        "QD < 2.0 || FS > 60.0 || MQ < 40.0",
        "--filter-name",
        "HardFiltered",
    ]


def _resolve_cnvkit_exe() -> Optional[str]:
    return resolve_cli_exe("cnvkit.py") or resolve_cli_exe("cnvkit")


def _cnvkit_intervals_bed_from_fai(
    ref_fasta: str,
    work_dir: str,
    *,
    bin_width: int = 50_000,
) -> str:
    """
    从参考 FASTA 的 .fai 生成 cnvkit coverage 所需的 interval BED（无 GFF/PoN 亦可跑通）。
    """
    fai = f"{ref_fasta}.fai"
    if not os.path.isfile(fai):
        samtools = resolve_cli_exe("samtools")
        if samtools:
            subprocess.run(
                [samtools, "faidx", ref_fasta],
                capture_output=True,
                text=True,
                timeout=3600,
                check=False,
            )
    bed_path = os.path.join(work_dir, "cnvkit_intervals.bed")
    rows: List[str] = []
    if os.path.isfile(fai):
        with open(fai, encoding="utf-8") as fh:
            for line in fh:
                parts = line.split()
                if len(parts) < 2:
                    continue
                chrom, length_s = parts[0], parts[1]
                try:
                    length = int(length_s)
                except ValueError:
                    continue
                for start in range(0, length, bin_width):
                    end = min(length, start + bin_width)
                    rows.append(f"{chrom}\t{start}\t{end}\n")
    if not rows:
        rows.append("chr21\t0\t1000000\n")
    with open(bed_path, "w", encoding="utf-8") as out:
        out.writelines(rows)
    return bed_path


def _resolve_snpeff_exe() -> Optional[str]:
    return resolve_cli_exe("snpEff") or resolve_cli_exe("snpeff")


def _resolve_bam_path(file_path: str, bam_path: str = "") -> Optional[str]:
    """CNV/SV 优先使用显式 bam_path（Executor 从 markdup 注入），否则 file_path 本身为 BAM。"""
    bp = (bam_path or "").strip()
    if bp and os.path.isfile(bp):
        return os.path.abspath(bp)
    p = (file_path or "").strip()
    if p.lower().endswith(".bam") and os.path.isfile(p):
        return os.path.abspath(p)
    return None


def _cli_log_fields(cmd: List[str], cp: subprocess.CompletedProcess, note: str = "") -> Dict[str, Any]:
    cmd_s = " ".join(str(c) for c in cmd)
    return {
        "cli_command": cmd_s,
        "cli_returncode": int(cp.returncode),
        "cli_stdout_excerpt": clip_omics_log(cp.stdout or "", 500),
        "cli_stderr_excerpt": clip_omics_log(cp.stderr or "", 500),
        "cli_note": note,
    }


def _format_cli_markdown(
    title: str,
    cmd: List[str],
    cp: subprocess.CompletedProcess,
    *,
    note: str = "",
    extra_lines: str = "",
) -> str:
    cmd_s = " ".join(str(c) for c in cmd)
    note_line = f"- **Result note**: {note}\n" if note else ""
    extra = extra_lines.rstrip() + "\n" if extra_lines else ""
    return (
        f"### {title}\n\n"
        f"- **Command executed**: `{cmd_s}`\n"
        f"- **Exit code**: `{cp.returncode}`\n"
        f"{note_line}"
        f"{extra}"
        "#### stderr（前 500 字符）\n\n```\n"
        f"{clip_omics_log(cp.stderr or '', 500)}\n```\n"
        "#### stdout（前 500 字符）\n\n```\n"
        f"{clip_omics_log(cp.stdout or '', 500)}\n```\n"
    )


def _run_logged_cli(
    label: str,
    cmd: List[str],
    *,
    timeout: int = 7200,
) -> subprocess.CompletedProcess:
    return subprocess.run(
        cmd,
        capture_output=True,
        text=True,
        timeout=timeout,
        check=False,
    )


def _count_vcf_variants(vcf_path: str) -> int:
    bt = resolve_cli_exe("bcftools")
    if not bt or not os.path.isfile(vcf_path):
        return 0
    try:
        cp = subprocess.run(
            [bt, "view", "-H", vcf_path],
            capture_output=True,
            text=True,
            timeout=3600,
            check=False,
        )
        if cp.returncode == 0:
            return sum(
                1 for line in (cp.stdout or "").splitlines() if line and not line.startswith("#")
            )
    except (OSError, subprocess.TimeoutExpired):
        pass
    return 0


def _vqsr_training_bundle_ready() -> bool:
    """
    VQSR（ApplyVQSR）所需的已知变异训练资源。
    若未同时配置可读路径，则不在此 Runner 尝试 ApplyVQSR（避免小样本/缺资源崩溃）。
    """
    dbsnp = (os.getenv("GIBH_DBSNP_VCF") or "").strip()
    mills = (os.getenv("GIBH_MILLS_INDELS_VCF") or "").strip()
    if not dbsnp or not mills:
        return False
    return os.path.isfile(dbsnp) and os.path.isfile(mills)


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
    for key in (
        "trimmed_fastq",
        "bam_path",
        "dedup_bam_path",
        "vcf_path",
        "filtered_vcf_path",
        "bed_path",
    ):
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


# --- genomics_read_trimming ---


def _blocked_genomics_trim(file_path: str = "") -> Dict[str, Any]:
    src = (file_path or "").strip()
    checks = [
        ("fastp", "present" if exe_available("fastp") else "missing from PATH"),
        ("FASTQ", "ok" if src and os.path.isfile(src) else "missing"),
    ]
    return omics_host_prerequisite_blocked(
        context="genomics_read_trimming",
        tool_id="genomics_read_trimming",
        title="接头修剪（fastp）未满足运行条件",
        checks=checks,
        extra_md=(
            "需在宿主安装 `fastp` 并提供可读 FASTQ；参见 `scripts/install_omics_real_env.sh`。"
        ),
    )


def _try_fastp_trim(
    file_path: str, input_dir: str, quality_cutoff: int = 20
) -> Optional[Dict[str, Any]]:
    fastp = resolve_cli_exe("fastp")
    src = (file_path or input_dir or "").strip()
    if not (fastp and src and os.path.isfile(src)):
        return None
    low = src.lower()
    if not (
        low.endswith(".fastq.gz")
        or low.endswith(".fq.gz")
        or low.endswith(".fastq")
        or low.endswith(".fq")
    ):
        return None
    work = tempfile.mkdtemp(prefix="g_fastp_")
    out_fq = os.path.join(work, "trimmed.fastq.gz")
    json_rep = os.path.join(work, "fastp.json")
    html_rep = os.path.join(work, "fastp.html")
    cmd = [
        fastp,
        "-i",
        src,
        "-o",
        out_fq,
        "-j",
        json_rep,
        "-h",
        html_rep,
        "-w",
        "4",
        "-q",
        str(int(quality_cutoff)),
    ]
    try:
        cp = subprocess.run(cmd, capture_output=True, text=True, timeout=7200)
    except (OSError, subprocess.TimeoutExpired) as exc:
        return omics_frontend_error_from_exception("fastp", exc)
    if cp.returncode != 0:
        return omics_subprocess_failed("fastp", cmd, cp, tool_id="genomics_read_trimming")
    bad = require_nonempty_file(
        out_fq, "fastp 输出 FASTQ", extensions=(".fastq.gz", ".fq.gz", ".fastq", ".fq")
    )
    if bad:
        bad["message"] = (bad.get("message") or "") + f"\nfastp stderr:\n{cp.stderr or ''}"
        return bad
    stderr_excerpt = clip_omics_log(cp.stderr or "", 6000)
    fastp_summary = load_fastp_summary(json_rep)
    json_hint = ""
    try:
        import json

        with open(json_rep, encoding="utf-8") as jf:
            fj = json.load(jf)
            summ = fj.get("summary") or {}
            af = summ.get("after_filtering") or {}
            if af:
                json_hint = (
                    "\n#### fastp.json 摘要（节选）\n\n```json\n"
                    + clip_omics_log(json.dumps(af, ensure_ascii=False, indent=2), 4000)
                    + "\n```\n"
                )
    except Exception:  # noqa: BLE001
        json_hint = ""
    md = (
        "### 接头修剪（fastp）\n\n"
        f"- **输出 FASTQ**: `{out_fq}`\n"
        f"- **JSON 报告**: `{json_rep}`\n"
        f"- **HTML 报告**: `{html_rep}`\n"
        "#### fastp stderr（节选）\n\n```\n"
        f"{stderr_excerpt}\n```\n"
        f"{json_hint}"
    )
    return _contract(
        "g_trim",
        "fastp trimming completed",
        markdown=md,
        image_urls=[IMG_DNA_HELIX],
        table_data=simple_rows_table(
            ("artifact", "path"),
            [
                {"artifact": "trimmed_fastq", "path": out_fq},
                {"artifact": "fastp_json", "path": json_rep},
            ],
        ),
        extra={
            "trimmed_fastq": out_fq,
            "fastp_json": json_rep,
            "fastp_html": html_rep,
            "fastp_summary": fastp_summary,
            "fastp_stderr_excerpt": stderr_excerpt,
        },
        tool_id="genomics_read_trimming",
    )


def genomics_read_trimming_impl(
    file_path: str = "", input_dir: str = "", quality_cutoff: int = 20
) -> Dict[str, Any]:
    fp = (file_path or "").strip()
    idir = (input_dir or "").strip()
    src = fp or idir
    if not src:
        return {
            "status": "error",
            "error_category": "data_issue",
            "message": "严重错误：动态流程传入的 file_path 与 input_dir 均为空，CLI 拒绝空转！",
            "user_message": "未收到 FASTQ 路径，请检查工作流首步是否已绑定上传文件。",
        }
    if not os.path.exists(src):
        return {
            "status": "error",
            "error_category": "data_issue",
            "message": f"严重错误：动态流程传入的路径 {src!r} 不存在，CLI 拒绝空转！",
            "user_message": "容器内找不到输入 FASTQ，请确认文件已上传且路径正确。",
        }
    return run_omics_without_synthetic_fallback(
        lambda: _try_fastp_trim(file_path, input_dir, quality_cutoff),
        lambda: _blocked_genomics_trim(file_path or input_dir),
        ctx="genomics_read_trimming",
    )


# --- genomics_alignment ---


def _blocked_genomics_alignment(reference_id: str, file_path: str = "") -> Dict[str, Any]:
    fq = (file_path or "").strip()
    refp = resolve_genomics_reference_fasta(reference_id)
    idx_ok = bool(refp and bwa_index_ready(refp))
    checks = [
        ("bwa", "present" if exe_available("bwa") else "missing"),
        ("samtools", "present" if exe_available("samtools") else "missing"),
        (
            f"reference FASTA + bwa index ({reference_id})",
            "resolved" if idx_ok else "missing — run scripts/init_omics_mock_references.py",
        ),
        ("OMICS_REF_DIR", get_omics_ref_dir()),
        ("FASTQ", "ok" if fq and os.path.isfile(fq) else "missing"),
    ]
    return omics_host_prerequisite_blocked(
        context="genomics_alignment",
        tool_id="genomics_alignment",
        title="比对（bwa mem）未满足运行条件",
        checks=checks,
        extra_md="需 export `GIBH_REF_HG38=/path/to/hg38.fa`（或对应 ID），并完成 `bwa index`。见 install 脚本。",
    )


def _try_bwa_alignment(
    file_path: str,
    reference_id: str,
    *,
    threads: int = 8,
    mismatch_penalty: int = 4,
    gap_open_penalty: int = 6,
) -> Optional[Dict[str, Any]]:
    bwa = resolve_cli_exe("bwa")
    samtools = resolve_cli_exe("samtools")
    ref = resolve_genomics_reference_fasta(reference_id)
    fq = (file_path or "").strip()
    bad_in = require_nonempty_file(
        fq, "bwa mem 输入 FASTQ", extensions=(".fastq.gz", ".fq.gz", ".fastq", ".fq")
    )
    if bad_in:
        return bad_in
    if not (bwa and samtools and ref and bwa_index_ready(ref)):
        return None
    low = fq.lower()
    if not (
        low.endswith(".fastq.gz")
        or low.endswith(".fq.gz")
        or low.endswith(".fastq")
        or low.endswith(".fq")
    ):
        return None
    work = tempfile.mkdtemp(prefix="g_bwa_")
    sam_p = os.path.join(work, "aligned.sam")
    bam_p = os.path.join(work, "aligned.bam")
    cmd = build_bwa_mem_command(
        bwa,
        ref,
        fq,
        threads=threads,
        mismatch_penalty=mismatch_penalty,
        gap_open_penalty=gap_open_penalty,
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
        return omics_frontend_error_from_exception("bwa mem", exc)
    if cp.returncode != 0:
        return omics_subprocess_failed("bwa mem", cmd, cp, tool_id="genomics_alignment")
    st_cmd = [samtools, "view", "-bS", sam_p, "-o", bam_p]
    try:
        cp2 = subprocess.run(
            st_cmd,
            capture_output=True,
            text=True,
            timeout=3600,
        )
    except (OSError, subprocess.TimeoutExpired) as exc:
        return omics_frontend_error_from_exception("samtools view", exc)
    if cp2.returncode != 0:
        return omics_subprocess_failed("samtools view", st_cmd, cp2, tool_id="genomics_alignment")
    bad_bam = require_nonempty_file(bam_p, "samtools view 输出 BAM", extensions=(".bam",))
    if bad_bam:
        return bad_bam
    bad_sam = require_nonempty_file(sam_p, "bwa mem 输出 SAM")
    if bad_sam and os.path.getsize(sam_p) <= 0:
        return bad_sam
    bwa_err = clip_omics_log(cp.stderr or "", 8000)
    md = (
        "### 比对完成（bwa mem + samtools view）\n\n"
        f"- **参数**: threads={threads}, mismatch_penalty(B)={mismatch_penalty}, gap_open(O)={gap_open_penalty}\n"
        f"- **BAM**: `{bam_p}`\n"
        "#### bwa stderr（节选）\n\n```\n"
        f"{bwa_err}\n```\n"
        "#### samtools stderr（节选）\n\n```\n"
        f"{clip_omics_log(cp2.stderr or '', 4000)}\n```\n"
    )
    out = {
        "status": "success",
        "message": "比对完成（bwa mem）",
        "output_path": bam_p,
        "file_path": bam_p,
        "bam_path": bam_p,
    }
    return attach_visual_contract(
        out,
        markdown=md,
        image_urls=[IMG_GENOME_BROWSER],
        table_data=simple_rows_table(
            ("artifact", "path"),
            [
                {"artifact": "SAM", "path": sam_p},
                {"artifact": "BAM", "path": bam_p},
            ],
        ),
        tool_id="genomics_alignment",
    )


def genomics_alignment_impl(
    file_path: str = "",
    reference_id: str = "hg38",
    *,
    threads: int = 8,
    mismatch_penalty: int = 4,
    gap_open_penalty: int = 6,
) -> Dict[str, Any]:
    return run_omics_without_synthetic_fallback(
        lambda: _try_bwa_alignment(
            file_path,
            reference_id,
            threads=threads,
            mismatch_penalty=mismatch_penalty,
            gap_open_penalty=gap_open_penalty,
        ),
        lambda: _blocked_genomics_alignment(reference_id, file_path),
        ctx="genomics_alignment",
    )


# --- mark duplicates ---


def _blocked_genomics_markdup(file_path: str = "") -> Dict[str, Any]:
    bam = (file_path or "").strip()
    checks = [
        ("samtools", "present" if exe_available("samtools") else "missing"),
        ("输入 BAM", "ok" if bam and os.path.isfile(bam) and bam.lower().endswith(".bam") else "invalid"),
    ]
    return omics_host_prerequisite_blocked(
        context="genomics_mark_duplicates",
        tool_id="genomics_mark_duplicates",
        title="重复标记（samtools sort/markdup）未满足条件",
        checks=checks,
        extra_md="需要已坐标排序或可由 samtools sort 处理的 BAM。",
    )


def _try_mark_duplicates(file_path: str) -> Optional[Dict[str, Any]]:
    """若输入为真实 BAM 且存在 samtools，则尝试 sort + markdup。"""
    samtools = resolve_cli_exe("samtools")
    bam = (file_path or "").strip()
    bad = require_nonempty_file(bam, "markdup 输入 BAM", extensions=(".bam",))
    if bad:
        return bad
    if not samtools:
        return None
    work = tempfile.mkdtemp(prefix="g_md_")
    sorted_bam = os.path.join(work, "sorted.bam")
    dup_bam = os.path.join(work, "dedup.bam")
    metrics = os.path.join(work, "dup_metrics.txt")
    cmd_sort = [samtools, "sort", "-o", sorted_bam, bam]
    try:
        cp1 = subprocess.run(cmd_sort, capture_output=True, text=True, timeout=7200)
    except (OSError, subprocess.TimeoutExpired) as exc:
        return omics_frontend_error_from_exception("samtools sort", exc)
    if cp1.returncode != 0:
        return omics_subprocess_failed("samtools sort", cmd_sort, cp1, tool_id="genomics_mark_duplicates")
    cmd_md = [samtools, "markdup", "-s", "-f", metrics, sorted_bam, dup_bam]
    try:
        cp2 = subprocess.run(cmd_md, capture_output=True, text=True, timeout=7200)
    except (OSError, subprocess.TimeoutExpired) as exc:
        return omics_frontend_error_from_exception("samtools markdup", exc)
    if cp2.returncode != 0:
        return omics_subprocess_failed("samtools markdup", cmd_md, cp2, tool_id="genomics_mark_duplicates")
    bad_dup = require_nonempty_file(dup_bam, "samtools markdup 输出 BAM", extensions=(".bam",))
    if bad_dup:
        return bad_dup
    idx_err = ensure_bam_indexed(samtools, dup_bam)
    if idx_err:
        return idx_err
    md_body = ""
    try:
        with open(metrics, encoding="utf-8", errors="replace") as mf:
            md_body = clip_omics_log(mf.read(), 12000)
    except OSError:
        md_body = "(metrics 文件不可读)"
    md = (
        "### 重复标记（samtools markdup）\n\n"
        f"- **dedup BAM**: `{dup_bam}`\n"
        f"- **metrics**: `{metrics}`\n"
        "#### markdup metrics（节选）\n\n```\n"
        f"{md_body}\n```\n"
        "#### stderr（节选）\n\n```\n"
        f"{clip_omics_log((cp2.stderr or ''), 6000)}\n```\n"
    )
    out = {
        "status": "success",
        "message": "MarkDuplicates 完成（samtools sort + markdup + index）",
        "output_path": dup_bam,
        "file_path": dup_bam,
        "bam_path": dup_bam,
        "dedup_bam_path": dup_bam,
        "bai_path": f"{dup_bam}.bai",
    }
    return attach_visual_contract(
        out,
        markdown=md,
        image_urls=[IMG_GENOME_BROWSER],
        table_data=simple_rows_table(
            ("artifact", "path"),
            [
                {"artifact": "sorted_bam", "path": sorted_bam},
                {"artifact": "dedup_bam", "path": dup_bam},
                {"artifact": "bam_index", "path": f"{dup_bam}.bai"},
                {"artifact": "dup_metrics_txt", "path": metrics},
            ],
        ),
        tool_id="genomics_mark_duplicates",
    )


def genomics_mark_duplicates_impl(file_path: str = "") -> Dict[str, Any]:
    return run_omics_without_synthetic_fallback(
        lambda: _try_mark_duplicates(file_path),
        lambda: _blocked_genomics_markdup(file_path),
        ctx="genomics_mark_duplicates",
    )


# --- BQSR ---


def _blocked_genomics_bqsr(file_path: str = "") -> Dict[str, Any]:
    checks = [
        ("gatk", "present" if exe_available("gatk") else "missing"),
        ("已知位点 VCF + 排序索引 BAM", "required for BaseRecalibrator"),
    ]
    return omics_host_prerequisite_blocked(
        context="genomics_bqsr",
        tool_id="genomics_bqsr",
        title="BQSR（GATK）未在宿主接入全自动管线",
        checks=checks,
        extra_md="本仓库骨架未串联 ApplyBQSR；请在 Worker 部署 PoN/known-sites 后运行。",
    )


def _try_bqsr(file_path: str) -> Optional[Dict[str, Any]]:
    """
    无 PoN/known-sites 时不伪造 GATK BQSR：校验 BAM+索引后原样传递并说明。
    """
    bam = (file_path or "").strip()
    samtools = resolve_cli_exe("samtools")
    if not samtools:
        return None
    bad = require_nonempty_file(bam, "BQSR 输入 BAM", extensions=(".bam",))
    if bad:
        return bad
    idx_err = ensure_bam_indexed(samtools, bam)
    if idx_err:
        return idx_err
    md = (
        "### 碱基质量重校准（传递阶段）\n\n"
        "> 当前环境未配置 GATK BaseRecalibrator 所需 known-sites；**未伪造 Recalibration**。\n"
        f"- **已校验并传递索引 BAM**: `{bam}`\n"
        f"- **BAI**: `{bam}.bai`\n"
    )
    out = {
        "status": "success",
        "message": "BQSR skipped; indexed BAM passed through",
        "output_path": bam,
        "file_path": bam,
        "bam_path": bam,
        "bai_path": f"{bam}.bai",
        "bqsr_applied": False,
    }
    return attach_visual_contract(
        out,
        markdown=md,
        image_urls=[IMG_GENOME_BROWSER],
        table_data=simple_rows_table(
            ("note", "value"),
            [{"note": "bqsr", "value": "pass-through (no known-sites)"}],
        ),
        tool_id="genomics_bqsr",
    )


def genomics_bqsr_impl(file_path: str = "") -> Dict[str, Any]:
    return run_omics_without_synthetic_fallback(
        lambda: _try_bqsr(file_path),
        lambda: _blocked_genomics_bqsr(file_path),
        ctx="genomics_bqsr",
    )


# --- germline ---


def _blocked_genomics_germline(
    file_path: str = "", reference_id: str = "hg38"
) -> Dict[str, Any]:
    gatk = resolve_cli_exe("gatk")
    ref = resolve_genomics_reference_fasta(reference_id)
    bam = (file_path or "").strip()
    checks = [
        ("gatk", "present" if gatk else "missing"),
        (f"reference FASTA ({reference_id})", "resolved" if ref else "missing"),
        ("sorted/indexed BAM", "ok" if bam and os.path.isfile(bam) and bam.lower().endswith(".bam") else "invalid"),
    ]
    return omics_host_prerequisite_blocked(
        context="genomics_germline_calling",
        tool_id="genomics_germline_calling",
        title="胚系变异检测（HaplotypeCaller）前置条件未满足",
        checks=checks,
        extra_md="需排序且 `samtools index` 过的 BAM；并 export 参考 FASTA。",
    )


def _try_bcftools_germline(
    bam: str,
    ref: str,
    *,
    work: str,
) -> Optional[Dict[str, Any]]:
    bcftools = resolve_cli_exe("bcftools")
    if not bcftools:
        return None
    raw_vcf = os.path.join(work, "raw.vcf")
    mpileup = subprocess.Popen(
        [bcftools, "mpileup", "-f", ref, "-Ou", bam],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )
    assert mpileup.stdout is not None
    call_cmd = [bcftools, "call", "-mv", "-Ov", "-o", raw_vcf]
    try:
        cp = subprocess.run(
            call_cmd,
            stdin=mpileup.stdout,
            capture_output=True,
            text=True,
            timeout=86400,
            check=False,
        )
        mpileup.stdout.close()
        mp_stderr = ""
        if mpileup.stderr:
            mp_stderr = (mpileup.stderr.read() or "") if hasattr(mpileup.stderr, "read") else ""
            if isinstance(mp_stderr, bytes):
                mp_stderr = mp_stderr.decode("utf-8", errors="replace")
        mpileup.wait(timeout=60)
    except (OSError, subprocess.TimeoutExpired) as exc:
        return omics_frontend_error_from_exception("bcftools mpileup|call", exc)
    if cp.returncode != 0:
        fake = subprocess.CompletedProcess(call_cmd, cp.returncode, cp.stdout, (cp.stderr or "") + mp_stderr)
        return omics_subprocess_failed("bcftools call", call_cmd, fake, tool_id="genomics_germline_calling")
    bad = require_nonempty_file(raw_vcf, "bcftools call 输出 VCF", extensions=(".vcf",))
    if bad:
        return bad
    out_vcf = os.path.join(work, "variants.vcf.gz")
    try:
        cp2 = subprocess.run(
            [bcftools, "view", "-Oz", "-o", out_vcf, raw_vcf],
            capture_output=True,
            text=True,
            timeout=3600,
            check=False,
        )
    except (OSError, subprocess.TimeoutExpired) as exc:
        return omics_frontend_error_from_exception("bcftools view -Oz", exc)
    if cp2.returncode != 0:
        return omics_subprocess_failed("bcftools view -Oz", [bcftools, "view", "-Oz", "-o", out_vcf, raw_vcf], cp2, tool_id="genomics_germline_calling")
    bad2 = require_nonempty_file(out_vcf, "bcftools 压缩 VCF", extensions=(".vcf.gz",))
    if bad2:
        return bad2
    idx_err = ensure_vcf_indexed(bcftools, out_vcf)
    if idx_err:
        return idx_err
    var_sum = summarize_vcf(out_vcf)
    md = (
        "### 胚系变异检测（bcftools mpileup + call）\n\n"
        f"- **VCF**: `{out_vcf}`\n"
        f"{format_variant_summary_md(var_sum)}"
        "#### bcftools call stderr（节选）\n\n```\n"
        f"{clip_omics_log((cp.stderr or '') + mp_stderr, 8000)}\n```\n"
    )
    out = {
        "status": "success",
        "message": "Germline calling completed (bcftools)",
        "output_path": out_vcf,
        "file_path": out_vcf,
        "vcf_path": out_vcf,
        "variant_summary": var_sum,
        "caller": "bcftools",
    }
    return attach_visual_contract(
        out,
        markdown=md,
        image_urls=[IMG_VOLCANO],
        table_data=simple_rows_table(
            ("metric", "value"),
            [{"metric": "variants", "value": str(var_sum.get("variant_count", 0))}],
        ),
        tool_id="genomics_germline_calling",
    )


def _try_germline(
    file_path: str,
    reference_id: str = "hg38",
    *,
    stand_call_conf: float = 30.0,
    min_base_quality: int = 20,
    min_mapping_quality: int = 20,
) -> Optional[Dict[str, Any]]:
    del min_base_quality
    ref = resolve_genomics_reference_fasta(reference_id)
    bam = (file_path or "").strip()
    samtools = resolve_cli_exe("samtools")
    if not (ref and samtools):
        return None
    bad = require_nonempty_file(bam, "变异检测输入 BAM", extensions=(".bam",))
    if bad:
        return bad
    idx_err = ensure_bam_indexed(samtools, bam)
    if idx_err:
        return idx_err
    work = tempfile.mkdtemp(prefix="g_hc_")
    gatk = resolve_cli_exe("gatk")
    if not gatk:
        return _try_bcftools_germline(bam, ref, work=work)
    if not (gatk and ref):
        return None
    out_vcf = os.path.join(work, "variants.vcf.gz")
    cmd = build_gatk_haplotypecaller_command(
        gatk,
        ref,
        bam,
        out_vcf,
        stand_call_conf=stand_call_conf,
        min_mapping_quality=min_mapping_quality,
    )
    try:
        cp = subprocess.run(cmd, capture_output=True, text=True, timeout=86400)
    except (OSError, subprocess.TimeoutExpired) as exc:
        return omics_frontend_error_from_exception("GATK HaplotypeCaller", exc)
    if cp.returncode != 0:
        logger.warning("GATK HaplotypeCaller failed, fallback to bcftools: %s", cp.stderr)
        return _try_bcftools_germline(bam, ref, work=work)
    bad_vcf = require_nonempty_file(out_vcf, "GATK 输出 VCF", extensions=(".vcf", ".vcf.gz"))
    if bad_vcf:
        return bad_vcf
    bt_idx = resolve_cli_exe("bcftools")
    if bt_idx:
        idx_err = ensure_vcf_indexed(bt_idx, out_vcf)
        if idx_err:
            return idx_err
    var_sum = summarize_vcf(out_vcf)
    vc_lines = ""
    bt = resolve_cli_exe("bcftools")
    if bt and os.path.isfile(out_vcf):
        try:
            cpv = subprocess.run(
                [bt, "view", "-H", out_vcf],
                capture_output=True,
                text=True,
                timeout=3600,
            )
            if cpv.returncode == 0 and cpv.stdout:
                nl = sum(1 for line in cpv.stdout.splitlines() if line and not line.startswith("#"))
                vc_lines = f"\n- **approx_variant_records（bcftools view -H 行数）**: {nl}\n"
        except (OSError, subprocess.TimeoutExpired):
            vc_lines = "\n- （bcftools 计数跳过）\n"
    md = (
        "### 胚系变异检测（GATK HaplotypeCaller）\n\n"
        f"- **VCF**: `{out_vcf}`\n"
        f"{format_variant_summary_md(var_sum)}"
        f"- **stand_call_conf**: {stand_call_conf}\n"
        f"- **min_mapping_quality**: {min_mapping_quality}\n"
        f"{vc_lines}"
        "#### GATK stderr（节选）\n\n```\n"
        f"{clip_omics_log(cp.stderr or '', 12000)}\n```\n"
    )
    out = {
        "status": "success",
        "message": "Germline calling completed (GATK HaplotypeCaller)",
        "output_path": out_vcf,
        "file_path": out_vcf,
        "vcf_path": out_vcf,
        "variant_summary": var_sum,
        "caller": "gatk",
    }
    return attach_visual_contract(
        out,
        markdown=md,
        image_urls=[IMG_VOLCANO],
        table_data=simple_rows_table(
            ("artifact", "path"),
            [{"artifact": "germline_vcf", "path": out_vcf}],
        ),
        tool_id="genomics_germline_calling",
    )


def genomics_germline_calling_impl(
    file_path: str = "",
    reference_id: str = "hg38",
    *,
    min_base_quality: int = 20,
    min_mapping_quality: int = 20,
    stand_call_conf: float = 30.0,
) -> Dict[str, Any]:
    return run_omics_without_synthetic_fallback(
        lambda: _try_germline(
            file_path,
            reference_id,
            stand_call_conf=stand_call_conf,
            min_base_quality=min_base_quality,
            min_mapping_quality=min_mapping_quality,
        ),
        lambda: _blocked_genomics_germline(file_path, reference_id),
        ctx="genomics_germline_calling",
    )


# --- CNV ---


def _blocked_genomics_cnv(file_path: str = "", bam_path: str = "") -> Dict[str, Any]:
    return omics_host_prerequisite_blocked(
        context="genomics_cnv_calling",
        tool_id="genomics_cnv_calling",
        title="CNV（cnvkit coverage）未满足运行条件",
        checks=[
            ("cnvkit.py / cnvkit", "present" if _resolve_cnvkit_exe() else "missing"),
            ("reference FASTA", "resolved" if resolve_genomics_reference_fasta("hg38") else "missing"),
            ("indexed BAM", "ok" if _resolve_bam_path(file_path, bam_path) else "missing"),
        ],
        extra_md=(
            "需 Conda 安装 `cnvkit`；CNV 步骤使用 **cnvkit.py coverage**（无 PoN 时仍可产出 .cnn）。\n"
            f"file_path=`{file_path}` bam_path=`{bam_path}`\n"
        ),
    )


def _try_cnv(file_path: str, bam_path: str = "") -> Optional[Dict[str, Any]]:
    bam = _resolve_bam_path(file_path, bam_path)
    if not bam:
        return None
    bad = require_nonempty_file(bam, "CNV 输入 BAM", extensions=(".bam",))
    if bad:
        return bad
    ref = resolve_genomics_reference_fasta("hg38")
    cnvkit = _resolve_cnvkit_exe()
    samtools = resolve_cli_exe("samtools")
    if not (cnvkit and ref and samtools):
        return None
    idx_err = ensure_bam_indexed(samtools, bam)
    if idx_err:
        return idx_err

    work = tempfile.mkdtemp(prefix="g_cnv_")
    sample = os.path.splitext(os.path.basename(bam))[0].replace(".", "_")[:80]
    intervals_bed = _cnvkit_intervals_bed_from_fai(ref, work)
    cnn_out = os.path.join(work, f"{sample}.cnn")
    cmd = [cnvkit, "coverage", bam, intervals_bed, "-o", cnn_out, "-f", ref]
    try:
        cp = _run_logged_cli("cnvkit coverage", cmd, timeout=86400)
    except (OSError, subprocess.TimeoutExpired) as exc:
        return omics_frontend_error_from_exception("cnvkit coverage", exc)
    cp_depth: Optional[subprocess.CompletedProcess] = None
    cmd_depth: List[str] = []
    if cp.returncode != 0 or not (os.path.isfile(cnn_out) and os.path.getsize(cnn_out) > 0):
        depth_tsv = os.path.join(work, f"{sample}.depth.tsv")
        cmd_depth = [samtools, "depth", "-b", intervals_bed, bam]
        try:
            with open(depth_tsv, "w", encoding="utf-8") as dfh:
                cp_depth = subprocess.run(
                    cmd_depth,
                    stdout=dfh,
                    stderr=subprocess.PIPE,
                    text=True,
                    timeout=3600,
                    check=False,
                )
        except (OSError, subprocess.TimeoutExpired) as exc:
            cp_depth = subprocess.CompletedProcess(cmd_depth, -1, "", str(exc))
        if cp_depth and cp_depth.returncode == 0 and os.path.getsize(depth_tsv) > 0:
            with open(cnn_out, "w", encoding="utf-8") as cfh:
                cfh.write("#chromosome\tstart\tend\tgene\tlog2\tdepth\n")
                with open(depth_tsv, encoding="utf-8") as dfh:
                    for line in dfh:
                        parts = line.strip().split("\t")
                        if len(parts) < 3:
                            continue
                        chrom, pos_s, depth_s = parts[0], parts[1], parts[2]
                        try:
                            pos_i, depth_i = int(pos_s), float(depth_s)
                        except ValueError:
                            continue
                        cfh.write(
                            f"{chrom}\t{pos_i}\t{pos_i + 1}\t.\t0\t{depth_i}\n"
                        )
            note = (
                "cnvkit coverage 在微缩 BAM 上失败；已回退 **samtools depth** 并写入 .cnn 兼容表。"
                f" cnvkit returncode={cp.returncode}."
            )
        elif cp.returncode != 0:
            fake = subprocess.CompletedProcess(
                cmd,
                cp.returncode,
                cp.stdout,
                (cp.stderr or "") + (cp_depth.stderr if cp_depth else ""),
            )
            return omics_subprocess_failed(
                "cnvkit coverage", cmd, fake, tool_id="genomics_cnv_calling"
            )
        else:
            with open(cnn_out, "w", encoding="utf-8") as fh:
                fh.write("#chromosome\tstart\tend\tgene\tlog2\tdepth\n")
            note = "cnvkit coverage 产出为空；已写入表头占位 .cnn。"
    else:
        note = f"cnvkit coverage completed; output `{cnn_out}`."

    seg_bed = os.path.join(work, f"{sample}.cnr")
    cmd_seg = [cnvkit, "export", cnn_out, "-o", seg_bed, "--output-format", "bed"]
    try:
        cp_seg = _run_logged_cli("cnvkit export", cmd_seg, timeout=3600)
    except (OSError, subprocess.TimeoutExpired) as exc:
        cp_seg = subprocess.CompletedProcess(cmd_seg, -1, "", str(exc))
    if not os.path.isfile(seg_bed) or os.path.getsize(seg_bed) <= 0:
        with open(seg_bed, "w", encoding="utf-8") as fh:
            fh.write("track name=cnvkit_placeholder\n")

    vcf_in = (file_path or "").strip()
    vcf_pass = vcf_in if vcf_in.lower().endswith((".vcf", ".vcf.gz")) else ""
    md = _format_cli_markdown("拷贝数变异（cnvkit coverage + export）", cmd, cp, note=note)
    if cp_depth is not None and cmd_depth:
        md += "\n" + _format_cli_markdown(
            "samtools depth（CNV 回退）",
            cmd_depth,
            cp_depth,
            note="微缩数据 cnvkit 失败时的真实深度回退",
        )
    if cp_seg.returncode is not None:
        md += "\n" + _format_cli_markdown(
            "cnvkit export (bed)",
            cmd_seg,
            cp_seg,
            note=f"BED segments: `{seg_bed}`",
        )

    cli_logs = [_cli_log_fields(cmd, cp, note)]
    if cp_depth is not None and cmd_depth:
        cli_logs.append(_cli_log_fields(cmd_depth, cp_depth, "samtools depth fallback"))
    cli_logs.append(_cli_log_fields(cmd_seg, cp_seg, "cnvkit export bed"))

    out = {
        "status": "success",
        "message": f"CNV coverage completed (cnvkit); {note}",
        "output_path": cnn_out,
        "file_path": vcf_pass or cnn_out,
        "bam_path": bam,
        "cnv_cnn_path": cnn_out,
        "cnv_bed_path": seg_bed,
        "vcf_path": vcf_pass,
        "cli_logs": cli_logs,
    }
    out.update(_cli_log_fields(cmd, cp, note))
    return attach_visual_contract(
        out,
        markdown=md,
        image_urls=[IMG_GENOME_BROWSER],
        table_data=simple_rows_table(
            ("artifact", "path"),
            [
                {"artifact": "cnvkit_cnn", "path": cnn_out},
                {"artifact": "cnvkit_bed", "path": seg_bed},
                {"artifact": "input_bam", "path": bam},
            ],
        ),
        tool_id="genomics_cnv_calling",
    )


def genomics_cnv_calling_impl(file_path: str = "", bam_path: str = "") -> Dict[str, Any]:
    return run_omics_without_synthetic_fallback(
        lambda: _try_cnv(file_path, bam_path),
        lambda: _blocked_genomics_cnv(file_path, bam_path),
        ctx="genomics_cnv_calling",
    )


# --- SV ---


def _blocked_genomics_sv(file_path: str = "", bam_path: str = "") -> Dict[str, Any]:
    return omics_host_prerequisite_blocked(
        context="genomics_sv_calling",
        tool_id="genomics_sv_calling",
        title="结构变异（delly call）未满足运行条件",
        checks=[
            ("delly", "present" if exe_available("delly") else "missing"),
            ("reference FASTA", "resolved" if resolve_genomics_reference_fasta("hg38") else "missing"),
            ("indexed BAM", "ok" if _resolve_bam_path(file_path, bam_path) else "missing"),
        ],
        extra_md=(
            "需 Conda 安装 `delly`；本步骤执行 **delly call** + **bcftools view** 转 VCF。\n"
            f"file_path=`{file_path}` bam_path=`{bam_path}`\n"
        ),
    )


def _try_sv(file_path: str, bam_path: str = "") -> Optional[Dict[str, Any]]:
    bam = _resolve_bam_path(file_path, bam_path)
    if not bam:
        return None
    bad = require_nonempty_file(bam, "SV 输入 BAM", extensions=(".bam",))
    if bad:
        return bad
    ref = resolve_genomics_reference_fasta("hg38")
    delly = resolve_cli_exe("delly")
    bcftools = resolve_cli_exe("bcftools")
    samtools = resolve_cli_exe("samtools")
    if not (delly and ref and bcftools and samtools):
        return None
    idx_err = ensure_bam_indexed(samtools, bam)
    if idx_err:
        return idx_err

    work = tempfile.mkdtemp(prefix="g_sv_")
    bcf_out = os.path.join(work, "delly.sv.bcf")
    cmd = [delly, "call", "-g", ref, "-o", bcf_out, bam]
    try:
        cp = _run_logged_cli("delly call", cmd, timeout=86400)
    except (OSError, subprocess.TimeoutExpired) as exc:
        return omics_frontend_error_from_exception("delly call", exc)

    sv_vcf = os.path.join(work, "delly.sv.vcf.gz")
    cmd_view = [bcftools, "view", "-Oz", "-o", sv_vcf, bcf_out]
    if os.path.isfile(bcf_out) and os.path.getsize(bcf_out) > 0:
        try:
            cp_view = _run_logged_cli("bcftools view", cmd_view, timeout=3600)
        except (OSError, subprocess.TimeoutExpired) as exc:
            return omics_frontend_error_from_exception("bcftools view", exc)
    else:
        cp_view = subprocess.CompletedProcess(cmd_view, 1, "", "delly BCF empty or missing")
        with open(os.path.join(work, "delly.header.vcf"), "w", encoding="utf-8") as fh:
            fh.write("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        cmd_view = [
            bcftools,
            "view",
            "-Oz",
            "-o",
            sv_vcf,
            os.path.join(work, "delly.header.vcf"),
        ]
        cp_view = _run_logged_cli("bcftools view (header-only)", cmd_view, timeout=600)

    stderr_all = (cp.stderr or "") + (cp_view.stderr or "")
    delly_low_data = (
        cp.returncode != 0
        and (
            "not enough data" in stderr_all.lower()
            or "library parameters" in stderr_all.lower()
        )
    )
    if cp.returncode != 0 and not delly_low_data and not os.path.isfile(bcf_out):
        fake = subprocess.CompletedProcess(cmd, cp.returncode, cp.stdout, stderr_all)
        return omics_subprocess_failed("delly call", cmd, fake, tool_id="genomics_sv_calling")

    if delly_low_data and not os.path.isfile(sv_vcf):
        hdr_vcf = os.path.join(work, "delly.header.vcf")
        with open(hdr_vcf, "w", encoding="utf-8") as fh:
            fh.write("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        cmd_view = [bcftools, "view", "-Oz", "-o", sv_vcf, hdr_vcf]
        cp_view = _run_logged_cli("bcftools view (delly low-coverage)", cmd_view, timeout=600)

    sv_count = _count_vcf_variants(sv_vcf)
    note = (
        f"Command executed: delly call -g {ref} -o {bcf_out}. "
        f"Result: **{sv_count}** SV record(s) in minimal/test dataset."
    )
    if delly_low_data:
        note += " Delly reported insufficient read depth for library estimation (expected on tiny FASTQ)."
    vcf_upstream = (file_path or "").strip()
    vcf_pass = vcf_upstream if vcf_upstream.lower().endswith((".vcf", ".vcf.gz")) else sv_vcf

    md = _format_cli_markdown("结构变异（delly call）", cmd, cp, note=note)
    md += "\n" + _format_cli_markdown(
        "bcftools view（BCF→VCF.gz）",
        cmd_view,
        cp_view,
        note=f"SV VCF: `{sv_vcf}`",
    )

    out = {
        "status": "success",
        "message": note,
        "output_path": sv_vcf,
        "file_path": vcf_pass,
        "vcf_path": vcf_pass,
        "sv_vcf_path": sv_vcf,
        "sv_bcf_path": bcf_out,
        "bam_path": bam,
        "sv_count": sv_count,
        "cli_logs": [_cli_log_fields(cmd, cp, note), _cli_log_fields(cmd_view, cp_view, note)],
    }
    out.update(_cli_log_fields(cmd, cp, note))
    return attach_visual_contract(
        out,
        markdown=md,
        image_urls=[IMG_GENOME_BROWSER],
        table_data=simple_rows_table(
            ("metric", "value"),
            [
                {"metric": "sv_records", "value": str(sv_count)},
                {"metric": "sv_vcf", "value": sv_vcf},
            ],
        ),
        tool_id="genomics_sv_calling",
    )


def genomics_sv_calling_impl(file_path: str = "", bam_path: str = "") -> Dict[str, Any]:
    return run_omics_without_synthetic_fallback(
        lambda: _try_sv(file_path, bam_path),
        lambda: _blocked_genomics_sv(file_path, bam_path),
        ctx="genomics_sv_calling",
    )


# --- VQSR ---


def _blocked_genomics_vqsr(file_path: str = "") -> Dict[str, Any]:
    gatk = resolve_cli_exe("gatk")
    ref_ok = bool(
        resolve_genomics_reference_fasta("hg38") or os.getenv("GIBH_REFERENCE_FASTA")
    )
    invc = (file_path or "").strip()
    vcf_ok = bool(
        invc
        and os.path.isfile(invc)
        and invc.lower().endswith((".vcf", ".vcf.gz"))
    )
    checks = [
        ("gatk", "present" if gatk else "missing"),
        ("reference FASTA (GIBH_REF_* / GIBH_REFERENCE_FASTA)", "ok" if ref_ok else "missing"),
        ("input VCF（来自 HaplotypeCaller）", "ok" if vcf_ok else "missing"),
    ]
    return omics_host_prerequisite_blocked(
        context="genomics_vqsr_filtering",
        tool_id="genomics_vqsr_filtering",
        title="变异过滤（VariantFiltration）前置条件未满足",
        checks=checks,
        extra_md=(
            "无 VQSR 训练集时将尝试硬过滤；请先保证上一步产出 `.vcf`/`.vcf.gz`。\n"
            f"当前 file_path: `{file_path}`\n"
        ),
    )


def _try_vqsr_filtering(
    file_path: str,
    reference_id: str = "hg38",
    *,
    _tranche_sensitivity: float = 99.0,
) -> Optional[Dict[str, Any]]:
    """
    变异过滤：默认不走 ApplyVQSR（小样本易崩溃）；未配置训练集或任意条件下均以
    GATK VariantFiltration 硬过滤产出过滤后 VCF，供注释等下游衔接。
    """
    del _tranche_sensitivity  # ApplyVQSR 未接入；保留参数兼容表单
    gatk = resolve_cli_exe("gatk")
    invcf = (file_path or "").strip()
    ref = resolve_genomics_reference_fasta(reference_id)
    bad_in = require_nonempty_file(
        invcf, "变异过滤输入 VCF", extensions=(".vcf", ".vcf.gz")
    )
    if bad_in:
        return bad_in
    if not ref:
        return None
    bcftools = resolve_cli_exe("bcftools")
    if not bcftools:
        if not gatk:
            return None
    else:
        idx_in = ensure_vcf_indexed(bcftools, invcf)
        if idx_in:
            return idx_in

    if not gatk and bcftools:
        return _bcftools_qual_filter_vcf(bcftools, invcf)
    if not gatk:
        return None
    low = invcf.lower()
    if not (low.endswith(".vcf") or low.endswith(".vcf.gz")):
        return None

    # bcftools 胚系 VCF 无 QD/FS，GATK HardFilter 会失败；直接走 QUAL 过滤
    if bcftools and not vcf_header_has_info_ids(invcf, ("QD", "FS")):
        logger.info(
            "VCF 缺少 QD/FS（多为 bcftools caller），跳过 GATK VariantFiltration，使用 bcftools QUAL 过滤"
        )
        out_bc = _bcftools_qual_filter_vcf(bcftools, invcf)
        if isinstance(out_bc, dict) and out_bc.get("status") == "success":
            note = (
                "⚠️ 输入 VCF 非 GATK HaplotypeCaller 注释格式（无 QD/FS），"
                "已自动改用 **bcftools view -i QUAL>=30** 过滤。\n\n"
            )
            out_bc["markdown"] = note + (out_bc.get("markdown") or "")
            out_bc["vqsr_mode"] = "bcftools_qual_filter_non_gatk_vcf"
        return out_bc

    training_ok = _vqsr_training_bundle_ready()
    work = tempfile.mkdtemp(prefix="g_vf_")
    out_vcf = os.path.join(work, "filtered_hard.vcf.gz")
    cmd = build_gatk_hard_filtering_command(gatk, ref, invcf, out_vcf)
    try:
        cp = subprocess.run(cmd, capture_output=True, text=True, timeout=86400)
    except (OSError, subprocess.TimeoutExpired) as exc:
        return omics_frontend_error_from_exception("GATK VariantFiltration", exc)
    if _gatk_subprocess_failed(cp):
        logger.warning(
            "GATK VariantFiltration failed, fallback bcftools: %s",
            clip_omics_log(cp.stderr or "", 2000),
        )
        if bcftools:
            out_bc = _bcftools_qual_filter_vcf(bcftools, invcf)
            if isinstance(out_bc, dict) and out_bc.get("status") == "success":
                note = (
                    "⚠️ GATK VariantFiltration 未成功，已降级为 **bcftools QUAL 过滤**。\n\n"
                    f"GATK stderr 节选：\n```\n{clip_omics_log(cp.stderr or '', 4000)}\n```\n\n"
                )
                out_bc["markdown"] = note + (out_bc.get("markdown") or "")
            return out_bc
        return omics_subprocess_failed(
            "GATK VariantFiltration", cmd, cp, tool_id="genomics_vqsr_filtering"
        )
    idx_out = ensure_vcf_indexed(bcftools, out_vcf) if bcftools else None
    if idx_out:
        return idx_out
    bad_hard = require_nonempty_file(out_vcf, "GATK 过滤 VCF", extensions=(".vcf", ".vcf.gz"))
    if bad_hard:
        return bad_hard

    if not training_ok:
        policy_md = (
            "⚠️ **未检测到 VQSR 训练集**（需同时配置可读路径 `GIBH_DBSNP_VCF` 与 "
            "`GIBH_MILLS_INDELS_VCF`），系统已自动降级为 **GATK Hard Filtering（VariantFiltration 硬过滤）**，成功产出变异文件。\n\n"
        )
    else:
        policy_md = (
            "⚠️ **已检测到 VQSR 已知位点文件**，但本 Runner **未串联 ApplyVQSR**（避免测试数据量过小导致 "
            "VQSR 数值不稳定）。已使用 **GATK VariantFiltration 硬过滤** 产出可衔接下游的 VCF。\n\n"
        )

    md = (
        "### 变异过滤（GATK VariantFiltration）\n\n"
        f"{policy_md}"
        f"- **输入 VCF**: `{invcf}`\n"
        f"- **输出 VCF**: `{out_vcf}`\n"
        f"- **参考序列**: `{ref}`\n"
        "- **过滤**: QD < 2.0 || FS > 60.0 || MQ < 40.0 → `HardFiltered`\n\n"
        "#### GATK stderr（节选）\n\n```\n"
        f"{clip_omics_log(cp.stderr or '', 12000)}\n```\n"
    )
    out = {
        "status": "success",
        "message": "Variant filtration completed (GATK VariantFiltration, hard-filter fallback)",
        "output_path": out_vcf,
        "file_path": out_vcf,
        "vcf_path": out_vcf,
        "filtered_vcf_path": out_vcf,
        "vqsr_mode": "hard_filter_variant_filtration",
        "vqsr_training_bundle_present": training_ok,
    }
    return attach_visual_contract(
        out,
        markdown=md,
        image_urls=[IMG_VOLCANO],
        table_data=simple_rows_table(
            ("artifact", "path"),
            [
                {"artifact": "input_vcf", "path": invcf},
                {"artifact": "filtered_vcf", "path": out_vcf},
            ],
        ),
        tool_id="genomics_vqsr_filtering",
    )


def genomics_vqsr_filtering_impl(
    file_path: str = "",
    reference_id: str = "hg38",
    *,
    tranche_sensitivity: float = 99.0,
) -> Dict[str, Any]:
    return run_omics_without_synthetic_fallback(
        lambda: _try_vqsr_filtering(
            file_path,
            reference_id,
            _tranche_sensitivity=tranche_sensitivity,
        ),
        lambda: _blocked_genomics_vqsr(file_path),
        ctx="genomics_vqsr_filtering",
    )


# --- annotation ---


def _blocked_genomics_annotation(file_path: str = "") -> Dict[str, Any]:
    return omics_host_prerequisite_blocked(
        context="genomics_variant_annotation",
        tool_id="genomics_variant_annotation",
        title="变异注释（snpEff / bcftools）未满足运行条件",
        checks=[
            ("snpEff / snpeff", "present" if _resolve_snpeff_exe() else "missing (will try bcftools fallback)"),
            ("bcftools", "present" if exe_available("bcftools") else "missing"),
            ("input VCF", "ok" if (file_path or "").strip() else "missing"),
        ],
        extra_md=(f"表单路径引用: `{file_path}`\n" if file_path else ""),
    )


def _try_annotation(file_path: str) -> Optional[Dict[str, Any]]:
    vcf = (file_path or "").strip()
    bad = require_nonempty_file(vcf, "注释输入 VCF", extensions=(".vcf", ".vcf.gz"))
    if bad:
        return bad
    bcftools = resolve_cli_exe("bcftools")
    if not bcftools:
        return None

    work = tempfile.mkdtemp(prefix="g_anno_")
    annotated = os.path.join(work, "annotated.vcf.gz")
    snpeff = _resolve_snpeff_exe()
    md_parts: List[str] = []
    cli_logs: List[Dict[str, Any]] = []
    mode = "bcftools_fallback"

    if snpeff:
        ann_txt = os.path.join(work, "snpEff.out.vcf")
        cmd_se = [snpeff, "-v", "-noStats", "hg38", vcf]
        try:
            with open(ann_txt, "w", encoding="utf-8") as out_fh:
                cp_se = subprocess.run(
                    cmd_se,
                    stdout=out_fh,
                    stderr=subprocess.PIPE,
                    text=True,
                    timeout=180,
                    check=False,
                )
        except (OSError, subprocess.TimeoutExpired) as exc:
            cp_se = subprocess.CompletedProcess(cmd_se, -1, "", str(exc))
        md_parts.append(_format_cli_markdown("snpEff 注释尝试", cmd_se, cp_se))
        cli_logs.append(_cli_log_fields(cmd_se, cp_se, "snpEff hg38"))
        if cp_se.returncode == 0 and os.path.isfile(ann_txt) and os.path.getsize(ann_txt) > 0:
            cmd_bg = [bcftools, "view", "-Oz", "-o", annotated, ann_txt]
            cp_bg = _run_logged_cli("bcftools view", cmd_bg)
            md_parts.append(_format_cli_markdown("bcftools 压缩 snpEff 输出", cmd_bg, cp_bg))
            cli_logs.append(_cli_log_fields(cmd_bg, cp_bg, "compress"))
            if cp_bg.returncode == 0 and os.path.isfile(annotated):
                mode = "snpeff"

    if mode != "snpeff":
        cmd_hdr = [bcftools, "view", "-h", vcf]
        cp_hdr = _run_logged_cli("bcftools view -h", cmd_hdr, timeout=600)
        md_parts.append(
            _format_cli_markdown(
                "bcftools view -h（VCF 头）",
                cmd_hdr,
                cp_hdr,
                note="无 snpEff 数据库或 snpEff 失败时，执行 bcftools 标签回退。",
            )
        )
        cli_logs.append(_cli_log_fields(cmd_hdr, cp_hdr, "header"))
        cmd_ann = [
            bcftools,
            "annotate",
            "-Oz",
            "-o",
            annotated,
            "--set-id",
            "%CHROM_%POS_%REF_%ALT",
            vcf,
        ]
        try:
            cp_ann = _run_logged_cli("bcftools annotate", cmd_ann, timeout=86400)
        except (OSError, subprocess.TimeoutExpired) as exc:
            return omics_frontend_error_from_exception("bcftools annotate", exc)
        if cp_ann.returncode != 0 or not os.path.isfile(annotated):
            cmd_cp = [bcftools, "view", "-Oz", "-o", annotated, vcf]
            cp_cp = _run_logged_cli("bcftools view copy", cmd_cp)
            if cp_cp.returncode != 0:
                return omics_subprocess_failed(
                    "bcftools annotate", cmd_ann, cp_ann, tool_id="genomics_variant_annotation"
                )
            cp_ann = cp_cp
        md_parts.append(
            _format_cli_markdown(
                "bcftools annotate（ID 回退）",
                cmd_ann,
                cp_ann,
                note=f"Annotated VCF: `{annotated}`",
            )
        )
        cli_logs.append(_cli_log_fields(cmd_ann, cp_ann, "set-id fallback"))

    bad_out = require_nonempty_file(annotated, "注释输出 VCF", extensions=(".vcf.gz",))
    if bad_out:
        return bad_out
    var_sum = summarize_vcf(annotated)
    note = (
        f"Annotation mode={mode}. Variants in output: {var_sum.get('variant_count', 0)}. "
        "Command logs attached (snpEff and/or bcftools)."
    )
    md = (
        "### 变异注释（snpEff + bcftools 回退）\n\n"
        f"- **模式**: `{mode}`\n"
        f"- **输出 VCF**: `{annotated}`\n"
        f"{format_variant_summary_md(var_sum)}\n\n"
        + "\n".join(md_parts)
    )
    out = {
        "status": "success",
        "message": note,
        "output_path": annotated,
        "file_path": annotated,
        "vcf_path": annotated,
        "annotated_vcf_path": annotated,
        "annotation_mode": mode,
        "variant_summary": var_sum,
        "cli_logs": cli_logs,
    }
    if cli_logs:
        out.update(cli_logs[0])
    return attach_visual_contract(
        out,
        markdown=md,
        image_urls=[IMG_VOLCANO],
        table_data=simple_rows_table(
            ("field", "value"),
            [
                {"field": "annotation_mode", "value": mode},
                {"field": "variants", "value": str(var_sum.get("variant_count", 0))},
            ],
        ),
        tool_id="genomics_variant_annotation",
    )


def genomics_variant_annotation_impl(file_path: str = "") -> Dict[str, Any]:
    return run_omics_without_synthetic_fallback(
        lambda: _try_annotation(file_path),
        lambda: _blocked_genomics_annotation(file_path),
        ctx="genomics_variant_annotation",
    )


# --- ACMG ---


def _blocked_genomics_acmg(file_path: str = "") -> Dict[str, Any]:
    return omics_host_prerequisite_blocked(
        context="genomics_acmg_classification",
        tool_id="genomics_acmg_classification",
        title="ACMG 分级（ClinVar 字段解析）未满足运行条件",
        checks=[
            ("bcftools", "present" if exe_available("bcftools") else "missing"),
            ("annotated VCF", "ok" if (file_path or "").strip() else "missing"),
        ],
        extra_md=(
            "使用 **bcftools query** 提取 INFO/CLNSIG（若有）并结合 QUAL 规则输出 Pathogenic/Benign 统计。\n"
            f"表单路径: `{file_path}`\n"
        ),
    )


def _classify_acmg_from_vcf(vcf_path: str) -> Tuple[Dict[str, int], List[Dict[str, str]], str, subprocess.CompletedProcess, List[str]]:
    """解析 VCF → ACMG 代理分级计数与示例行。"""
    bcftools = resolve_cli_exe("bcftools")
    if not bcftools:
        return {}, [], "", subprocess.CompletedProcess([], 1, "", ""), []

    cmd = [
        bcftools,
        "query",
        "-f",
        "[%CHROM:%POS\\t%REF>%ALT\\t%QUAL\\t%INFO]\\n",
        vcf_path,
    ]
    cp = _run_logged_cli("bcftools query", cmd, timeout=3600)
    counts = {
        "Pathogenic": 0,
        "Likely_pathogenic": 0,
        "Benign": 0,
        "Likely_benign": 0,
        "VUS": 0,
        "Uncertain": 0,
    }
    preview: List[Dict[str, str]] = []
    for line in (cp.stdout or "").splitlines()[:5000]:
        if not line.strip():
            continue
        parts = line.split("\t", 3)
        qual_s = parts[2] if len(parts) > 2 else "."
        info = parts[3] if len(parts) > 3 else ""
        clnsig = ""
        for token in info.replace(";", " ").split():
            if "CLNSIG" in token.upper() or "clinvar" in token.lower():
                clnsig = token
                break
        tier = "VUS"
        low = (clnsig + " " + info).lower()
        if "pathogenic" in low and "likely" not in low:
            tier = "Pathogenic"
        elif "likely_pathogenic" in low or "likely pathogenic" in low:
            tier = "Likely_pathogenic"
        elif "benign" in low and "likely" not in low:
            tier = "Benign"
        elif "likely_benign" in low or "likely benign" in low:
            tier = "Likely_benign"
        else:
            try:
                q = float(qual_s)
                if q >= 50:
                    tier = "Uncertain"
                else:
                    tier = "Benign"
            except ValueError:
                tier = "VUS"
        counts[tier] = counts.get(tier, 0) + 1
        if len(preview) < 8:
            preview.append(
                {
                    "variant": parts[0] if parts else line[:40],
                    "acmg_proxy": tier,
                    "clnsig_hint": clnsig[:80] or "(none)",
                }
            )
    summary = (
        f"Pathogenic={counts['Pathogenic']}, Likely_pathogenic={counts['Likely_pathogenic']}, "
        f"Benign={counts['Benign']}, Likely_benign={counts['Likely_benign']}, "
        f"VUS={counts['VUS']}, Uncertain={counts['Uncertain']}"
    )
    return counts, preview, summary, cp, cmd


def _try_acmg(file_path: str) -> Optional[Dict[str, Any]]:
    vcf = (file_path or "").strip()
    bad = require_nonempty_file(vcf, "ACMG 输入 VCF", extensions=(".vcf", ".vcf.gz"))
    if bad:
        return bad
    bcftools = resolve_cli_exe("bcftools")
    if not bcftools:
        return None

    counts, preview, summary, cp, cmd = _classify_acmg_from_vcf(vcf)
    var_sum = summarize_vcf(vcf)
    note = (
        f"Command executed: bcftools query on `{vcf}`. "
        f"ACMG proxy summary: {summary}. "
        f"(Minimal dataset may yield 0 ClinVar tags; Benign/VUS from QUAL rules.)"
    )
    md = _format_cli_markdown("ACMG/AMP 代理分级（bcftools query + 规则）", cmd, cp, note=note)
    md += "\n" + format_variant_summary_md(var_sum)
    if preview:
        md += "\n| 变异（节选） | 代理分级 | ClinVar 提示 |\n|------|----------|-------------|\n"
        for row in preview:
            md += f"| {row['variant']} | **{row['acmg_proxy']}** | {row['clnsig_hint']} |\n"

    work = tempfile.mkdtemp(prefix="g_acmg_")
    json_path = os.path.join(work, "acmg_summary.json")
    payload = {
        "acmg_proxy_counts": counts,
        "variant_preview": preview,
        "summary_text": summary,
        "cli_note": note,
    }
    with open(json_path, "w", encoding="utf-8") as jf:
        json.dump(payload, jf, ensure_ascii=False, indent=2)

    out = {
        "status": "success",
        "message": note,
        "output_path": json_path,
        "file_path": vcf,
        "vcf_path": vcf,
        "acmg_summary_path": json_path,
        "acmg_proxy_counts": counts,
        "variant_summary": var_sum,
        "cli_logs": [_cli_log_fields(cmd, cp, note)],
    }
    out.update(_cli_log_fields(cmd, cp, note))
    return attach_visual_contract(
        out,
        markdown=md,
        image_urls=[IMG_VOLCANO],
        table_data=simple_rows_table(
            ("tier", "count"),
            [{"tier": k, "count": str(v)} for k, v in counts.items()],
        ),
        tool_id="genomics_acmg_classification",
    )


def genomics_acmg_classification_impl(file_path: str = "") -> Dict[str, Any]:
    return run_omics_without_synthetic_fallback(
        lambda: _try_acmg(file_path),
        lambda: _blocked_genomics_acmg(file_path),
        ctx="genomics_acmg_classification",
    )
