"""
基因组学下游步骤：真实 CLI（subprocess）+ stdout/stderr 写入 Markdown；
缺依赖时返回 **status:error** 与环境诊断表，不再使用仿真生物学占位指标。
"""
from __future__ import annotations

import logging
import os
import shutil
import subprocess
import tempfile
from typing import Any, Dict, List, Optional

from .omics_mock_ui import (
    IMG_DNA_HELIX,
    IMG_GENOME_BROWSER,
    IMG_VOLCANO,
    attach_visual_contract,
    simple_rows_table,
)
from .omics_pipeline_env import (
    clip_omics_log,
    exe_available,
    omics_frontend_error_from_exception,
    omics_host_prerequisite_blocked,
    omics_subprocess_failed,
    resolve_cli_exe,
    resolve_reference_fasta,
    run_omics_without_synthetic_fallback,
    write_temp_mock_artifact,
)

logger = logging.getLogger(__name__)

# 与 scripts/setup_minimal_genomics_test_data.sh 制备路径一致；未 export GIBH_REF_HG38 时作 hg38 兜底。
MINIMAL_HG38_REFERENCE_FALLBACK = "/tmp/omics_test_data/ref/chr21.fa"


def resolve_genomics_reference_fasta(reference_id: str) -> Optional[str]:
    """
    解析参考序列：优先环境变量 GIBH_REF_* / GIBH_REFERENCE_FASTA；
    若仍无解且 reference_id 为 hg38/grch38，则尝试本机一键制备的 chr21 微缩参考。
    """
    ref = resolve_reference_fasta(reference_id)
    if ref:
        return ref
    rid = (reference_id or "").strip().lower().replace(" ", "")
    if rid not in ("hg38", "grch38"):
        return None
    env_or_default = os.environ.get("GIBH_REF_HG38", MINIMAL_HG38_REFERENCE_FALLBACK)
    if env_or_default and os.path.isfile(env_or_default):
        return os.path.abspath(env_or_default)
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
    stderr_excerpt = clip_omics_log(cp.stderr or "", 6000)
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
        extra={"trimmed_fastq": out_fq, "fastp_json": json_rep, "fastp_stderr_excerpt": stderr_excerpt},
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
    checks = [
        ("bwa", "present" if exe_available("bwa") else "missing"),
        ("samtools", "present" if exe_available("samtools") else "missing"),
        (f"reference FASTA ({reference_id})", "resolved" if refp else "missing env GIBH_REF_*"),
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
    if not (bwa and samtools and ref and fq and os.path.isfile(fq)):
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
    if not (samtools and bam and os.path.isfile(bam) and bam.lower().endswith(".bam")):
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
        "message": "MarkDuplicates 完成（samtools）",
        "output_path": dup_bam,
        "file_path": dup_bam,
        "dedup_bam_path": dup_bam,
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


def _try_bqsr(_file_path: str) -> Optional[Dict[str, Any]]:
    return None


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


def _try_germline(
    file_path: str,
    reference_id: str = "hg38",
    *,
    stand_call_conf: float = 30.0,
    min_base_quality: int = 20,
    min_mapping_quality: int = 20,
) -> Optional[Dict[str, Any]]:
    del min_base_quality  # 预留与表单一致；CLI 构建可按需扩展
    gatk = resolve_cli_exe("gatk")
    ref = resolve_genomics_reference_fasta(reference_id)
    bam = (file_path or "").strip()
    if not (gatk and ref and bam and os.path.isfile(bam) and bam.lower().endswith(".bam")):
        return None
    work = tempfile.mkdtemp(prefix="g_hc_")
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
        return omics_subprocess_failed(
            "GATK HaplotypeCaller", cmd, cp, tool_id="genomics_germline_calling"
        )
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


def _blocked_genomics_cnv(file_path: str = "") -> Dict[str, Any]:
    return omics_host_prerequisite_blocked(
        context="genomics_cnv_calling",
        tool_id="genomics_cnv_calling",
        title="CNV 检测未在宿主接入全自动管线",
        checks=[
            ("gatk", "present" if exe_available("gatk") else "missing"),
            ("PoN / intervals / counts table", "required"),
        ],
        extra_md=(f"表单路径引用: `{file_path}`\n" if file_path else ""),
    )


def _try_cnv(_file_path: str) -> Optional[Dict[str, Any]]:
    return None


def genomics_cnv_calling_impl(file_path: str = "") -> Dict[str, Any]:
    return run_omics_without_synthetic_fallback(
        lambda: _try_cnv(file_path),
        lambda: _blocked_genomics_cnv(file_path),
        ctx="genomics_cnv_calling",
    )


# --- SV ---


def _blocked_genomics_sv(file_path: str = "") -> Dict[str, Any]:
    return omics_host_prerequisite_blocked(
        context="genomics_sv_calling",
        tool_id="genomics_sv_calling",
        title="结构变异（Manta/Lumpy 等）未在宿主接入自动化",
        checks=[
            (
                "manta scripts",
                "present"
                if (resolve_cli_exe("configManta.py") or resolve_cli_exe("runWorkflow.py"))
                else "missing",
            ),
            ("完整配置与 BAM", "required"),
        ],
        extra_md=(f"表单路径引用: `{file_path}`\n" if file_path else ""),
    )


def _try_sv(_file_path: str) -> Optional[Dict[str, Any]]:
    return None


def genomics_sv_calling_impl(file_path: str = "") -> Dict[str, Any]:
    return run_omics_without_synthetic_fallback(
        lambda: _try_sv(file_path),
        lambda: _blocked_genomics_sv(file_path),
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
    if not gatk or not ref:
        return None
    if not invcf or not os.path.isfile(invcf):
        return None
    low = invcf.lower()
    if not (low.endswith(".vcf") or low.endswith(".vcf.gz")):
        return None

    training_ok = _vqsr_training_bundle_ready()
    work = tempfile.mkdtemp(prefix="g_vf_")
    out_vcf = os.path.join(work, "filtered_hard.vcf.gz")
    cmd = build_gatk_hard_filtering_command(gatk, ref, invcf, out_vcf)
    try:
        cp = subprocess.run(cmd, capture_output=True, text=True, timeout=86400)
    except (OSError, subprocess.TimeoutExpired) as exc:
        return omics_frontend_error_from_exception("GATK VariantFiltration", exc)
    if cp.returncode != 0:
        return omics_subprocess_failed(
            "GATK VariantFiltration", cmd, cp, tool_id="genomics_vqsr_filtering"
        )

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
        title="变异注释（VEP/SnpEff）未在宿主接入自动化",
        checks=[
            ("vep / snpEff", "至少其一需可用"),
            ("缓存与输入 VCF", "required"),
        ],
        extra_md=(f"表单路径引用: `{file_path}`\n" if file_path else ""),
    )


def _try_annotation(_file_path: str) -> Optional[Dict[str, Any]]:
    return None


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
        title="ACMG 分级未在宿主接入自动化规则引擎",
        checks=[
            ("注释后 VCF / TSV", "required"),
            ("ClinVar / 内部证据库", "deployment-specific"),
        ],
        extra_md=(f"表单路径引用: `{file_path}`\n" if file_path else ""),
    )


def _try_acmg(_file_path: str) -> Optional[Dict[str, Any]]:
    return None


def genomics_acmg_classification_impl(file_path: str = "") -> Dict[str, Any]:
    return run_omics_without_synthetic_fallback(
        lambda: _try_acmg(file_path),
        lambda: _blocked_genomics_acmg(file_path),
        ctx="genomics_acmg_classification",
    )
