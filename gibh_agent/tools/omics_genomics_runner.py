"""
基因组学下游步骤：真实 CLI 骨架（subprocess）+ 缺依赖/失败时的仿真输出。

仅在检测到对应可执行文件、输入文件与（比对类）参考序列可用时执行真实命令；
否则返回结构化 markdown / image_urls / table_data，保证前端可渲染。
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
from .omics_derived_analysis import (
    fastq_bundle_from_path,
    genomics_proxy_acmg,
    genomics_proxy_annotation,
    genomics_proxy_bqsr,
    genomics_proxy_cnv,
    genomics_proxy_germline,
    genomics_proxy_markdup,
    genomics_proxy_sv,
    genomics_proxy_trim,
    genomics_proxy_vqsr,
    genomics_proxy_alignment,
)
from .omics_pipeline_env import (
    degraded_banner_md,
    exe_available,
    resolve_reference_fasta,
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


# --- genomics_read_trimming ---


def _degraded_trimming(file_path: str = "") -> Dict[str, Any]:
    b = fastq_bundle_from_path(file_path)
    if b:
        md, tbl, msg = genomics_proxy_trim(b)
        out = _contract(
            "g_trim",
            msg,
            markdown=md,
            image_urls=[IMG_DNA_HELIX],
            table_data=tbl,
            extra={"omics_analysis_mode": "fastq_stream_proxy", "proxy_metrics": b},
            tool_id="genomics_read_trimming",
        )
        return out
    md = degraded_banner_md("fastp 或可读的配对 FASTQ") + (
        "### 接头修剪与质量裁剪（仿真）\n\n"
        "| 指标 | 仿真值 |\n|------|--------|\n"
        "| 保留比例 | 97.2% |\n"
        "| 修剪后 Q30 | 92.1% |\n"
        "| Adapter 残留率 | 0.08% |\n"
    )
    return _contract(
        "g_trim",
        "Read trimming (simulated)",
        markdown=md,
        image_urls=[IMG_DNA_HELIX],
        table_data=simple_rows_table(
            ("metric", "value"),
            [
                {"metric": "retained_fraction", "value": "0.972"},
                {"metric": "adapter_residual", "value": "0.0008"},
            ],
        ),
        tool_id="genomics_read_trimming",
    )


def _try_fastp_trim(file_path: str, input_dir: str) -> Optional[Dict[str, Any]]:
    fastp = shutil.which("fastp")
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
    ]
    try:
        subprocess.run(cmd, check=True, capture_output=True, text=True, timeout=7200)
        md = (
            "### 接头修剪（fastp）\n\n"
            f"- **输出 FASTQ**: `{out_fq}`\n"
            f"- **JSON 报告**: `{json_rep}`\n"
            "- 完整 MultiQC 聚合建议在 Worker 上执行。\n"
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
            extra={"trimmed_fastq": out_fq, "fastp_json": json_rep},
            tool_id="genomics_read_trimming",
        )
    except (subprocess.CalledProcessError, OSError, subprocess.TimeoutExpired) as exc:
        logger.warning("fastp run failed: %s", exc)
        return None


def genomics_read_trimming_impl(file_path: str = "", input_dir: str = "") -> Dict[str, Any]:
    return run_or_degrade(
        lambda: _try_fastp_trim(file_path, input_dir),
        lambda: _degraded_trimming(file_path or input_dir),
        ctx="genomics_read_trimming",
    )


# --- genomics_alignment ---


def _degraded_alignment(reference_id: str, file_path: str = "") -> Dict[str, Any]:
    b = fastq_bundle_from_path(file_path)
    if b:
        md, tbl, msg = genomics_proxy_alignment(b, reference_id)
        out = _contract(
            "g_aln",
            msg,
            markdown=md,
            image_urls=[IMG_GENOME_BROWSER],
            table_data=tbl,
            extra={"omics_analysis_mode": "fastq_stream_proxy", "proxy_metrics": b},
            tool_id="genomics_alignment",
        )
        out["bam_path"] = out["output_path"]
        return out
    md = degraded_banner_md("BWA / Samtools 与有效参考基因组 FASTA（见 GIBH_REF_*）") + (
        "### 参考基因组比对（仿真）\n\n"
        f"- **参考序列 ID**: `{reference_id}`\n"
        "- **比对率**: **98.7%**\n"
        "- **平均覆盖度**: **35X**\n"
        "- **重复率（预估）**: 12.1%\n\n"
        "| #CHROM | POS | ID | REF | ALT | QUAL | FILTER |\n"
        "|--------|-----|----|----|-----|------|--------|\n"
        "| chr1 | 10497 | . | CC | C | 45 | PASS |\n"
        "| chr17 | 7577548 | rs113488022 | C | T | 892 | PASS |\n"
        "| chr22 | 42522503 | . | G | A | 210 | PASS |\n"
    )
    out = _contract(
        "g_aln",
        f"Alignment simulated ref={reference_id}",
        markdown=md,
        image_urls=[IMG_GENOME_BROWSER],
        table_data=simple_rows_table(
            ("metric", "value"),
            [
                {"metric": "overall_mapped_rate", "value": "98.7%"},
                {"metric": "mean_coverage", "value": "35X"},
                {"metric": "duplicate_rate", "value": "12.1%"},
            ],
        ),
        tool_id="genomics_alignment",
    )
    out["bam_path"] = out["output_path"]
    return out


def _try_bwa_alignment(file_path: str, reference_id: str) -> Optional[Dict[str, Any]]:
    bwa = shutil.which("bwa")
    samtools = shutil.which("samtools")
    ref = resolve_reference_fasta(reference_id)
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
    try:
        with open(sam_p, "wb") as sf:
            subprocess.run(
                [bwa, "mem", "-t", "4", ref, fq],
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
            "### 比对完成（bwa mem + samtools view）\n\n"
            f"- **BAM**: `{bam_p}`\n"
            "- 后续请在 Worker 上完成 `samtools sort/index` 与质控。\n"
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
    except (subprocess.CalledProcessError, OSError, subprocess.TimeoutExpired) as exc:
        logger.warning("bwa mem pipeline failed: %s", exc)
        return None


def genomics_alignment_impl(file_path: str = "", reference_id: str = "hg38") -> Dict[str, Any]:
    return run_or_degrade(
        lambda: _try_bwa_alignment(file_path, reference_id),
        lambda: _degraded_alignment(reference_id, file_path),
        ctx="genomics_alignment",
    )


# --- mark duplicates ---


def _degraded_markdup(file_path: str = "") -> Dict[str, Any]:
    b = fastq_bundle_from_path(file_path)
    if b:
        md, tbl, msg = genomics_proxy_markdup(b)
        out = _contract(
            "g_md",
            msg,
            markdown=md,
            image_urls=[IMG_GENOME_BROWSER],
            table_data=tbl,
            extra={"omics_analysis_mode": "fastq_stream_proxy", "proxy_metrics": b},
            tool_id="genomics_mark_duplicates",
        )
        out["dedup_bam_path"] = out["output_path"]
        return out
    md = degraded_banner_md("Picard MarkDuplicates 或含排序后的真实 BAM") + (
        "### 排序与重复标记（仿真）\n\n"
        "| 统计项 | 仿真值 |\n|--------|--------|\n"
        "| 光学重复 | 7.1% |\n"
        "| PCR 重复 | 5.3% |\n"
    )
    out = _contract(
        "g_md",
        "Sort + MarkDuplicates (simulated)",
        markdown=md,
        image_urls=[IMG_GENOME_BROWSER],
        table_data=simple_rows_table(
            ("category", "fraction"),
            [
                {"category": "optical_dup", "fraction": "0.071"},
                {"category": "pcr_dup", "fraction": "0.053"},
            ],
        ),
        tool_id="genomics_mark_duplicates",
    )
    out["dedup_bam_path"] = out["output_path"]
    return out


def _try_mark_duplicates(file_path: str) -> Optional[Dict[str, Any]]:
    """骨架：若输入为真实 BAM 且存在 samtools，则尝试 sort + markdup。"""
    samtools = shutil.which("samtools")
    bam = (file_path or "").strip()
    if not (samtools and bam and os.path.isfile(bam) and bam.lower().endswith(".bam")):
        return None
    work = tempfile.mkdtemp(prefix="g_md_")
    sorted_bam = os.path.join(work, "sorted.bam")
    dup_bam = os.path.join(work, "dedup.bam")
    metrics = os.path.join(work, "dup_metrics.txt")
    try:
        subprocess.run(
            [samtools, "sort", "-o", sorted_bam, bam],
            check=True,
            timeout=7200,
        )
        subprocess.run(
            [
                samtools,
                "markdup",
                "-s",
                "-f",
                metrics,
                sorted_bam,
                dup_bam,
            ],
            check=True,
            timeout=7200,
        )
        md = (
            "### 重复标记（samtools markdup）\n\n"
            f"- **dedup BAM**: `{dup_bam}`\n"
            f"- **metrics**: `{metrics}`\n"
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
                ],
            ),
            tool_id="genomics_mark_duplicates",
        )
    except (subprocess.CalledProcessError, OSError, subprocess.TimeoutExpired) as exc:
        logger.warning("markdup failed: %s", exc)
        return None


def genomics_mark_duplicates_impl(file_path: str = "") -> Dict[str, Any]:
    return run_or_degrade(
        lambda: _try_mark_duplicates(file_path),
        lambda: _degraded_markdup(file_path),
        ctx="genomics_mark_duplicates",
    )


# --- BQSR ---


def _degraded_bqsr(file_path: str = "") -> Dict[str, Any]:
    b = fastq_bundle_from_path(file_path)
    if b:
        md, tbl, msg = genomics_proxy_bqsr(b)
        return _contract(
            "g_bqsr",
            msg,
            markdown=md,
            image_urls=[IMG_DNA_HELIX],
            table_data=tbl,
            extra={"omics_analysis_mode": "fastq_stream_proxy", "proxy_metrics": b},
            tool_id="genomics_bqsr",
        )
    md = degraded_banner_md("GATK BaseRecalibrator / ApplyBQSR 与已知位点资源") + (
        "### 碱基质量重校准 BQSR（仿真）\n\n"
        "- 已知位点集：dbSNP + Mills + 1000G（占位文案）\n"
        "| Round | max_coord_shift |\n"
        "|-------|-----------------|\n"
        "| 1 | 3.2e-4 |\n"
        "| 2 | 8.1e-5 |\n"
    )
    return _contract(
        "g_bqsr",
        "BQSR (simulated)",
        markdown=md,
        image_urls=[IMG_DNA_HELIX],
        table_data=simple_rows_table(
            ("round", "max_coord_shift"),
            [
                {"round": "1", "max_coord_shift": "3.2e-4"},
                {"round": "2", "max_coord_shift": "8.1e-5"},
            ],
        ),
        tool_id="genomics_bqsr",
    )


def _try_bqsr(_file_path: str) -> Optional[Dict[str, Any]]:
    gatk = shutil.which("gatk")
    if not gatk:
        return None
    if not smoke_subprocess([gatk, "--help"], timeout=180):
        return None
    logger.info(
        "BQSR 完整管线需已知位点 VCF 与排序 BAM；宿主仅验证 `gatk --help`。仿真降级。"
    )
    return None


def genomics_bqsr_impl(file_path: str = "") -> Dict[str, Any]:
    return run_or_degrade(
        lambda: _try_bqsr(file_path),
        lambda: _degraded_bqsr(file_path),
        ctx="genomics_bqsr",
    )


# --- germline ---


def _degraded_germline(file_path: str = "") -> Dict[str, Any]:
    b = fastq_bundle_from_path(file_path)
    if b:
        md, tbl, msg = genomics_proxy_germline(b)
        out = _contract(
            "g_germ",
            msg,
            markdown=md,
            image_urls=[IMG_VOLCANO],
            table_data=tbl,
            extra={"omics_analysis_mode": "fastq_stream_proxy", "proxy_metrics": b},
            tool_id="genomics_germline_calling",
        )
        out["vcf_path"] = out["output_path"]
        return out
    md = degraded_banner_md("GATK HaplotypeCaller / bcftools 与参考序列") + (
        "### 胚系变异检测（仿真）\n\n"
        "| 类型 | 仿真计数 |\n|------|----------|\n"
        "| SNV | 3,842,100 |\n"
        "| Indel | 412,055 |\n\n"
        "| #CHROM | POS | REF | ALT | QUAL | FILTER |\n"
        "|--------|-----|-----|-----|------|--------|\n"
        "| chr1 | 248815 | G | A | 214 | PASS |\n"
        "| chr17 | 41245466 | G | A | 99 | PASS |\n"
    )
    out = _contract(
        "g_germ",
        "Germline calling (simulated)",
        markdown=md,
        image_urls=[IMG_VOLCANO],
        table_data=simple_rows_table(
            ("variant_type", "simulated_count"),
            [
                {"variant_type": "SNV", "simulated_count": "3842100"},
                {"variant_type": "INDEL", "simulated_count": "412055"},
            ],
        ),
        tool_id="genomics_germline_calling",
    )
    out["vcf_path"] = out["output_path"]
    return out


def _try_germline(_file_path: str) -> Optional[Dict[str, Any]]:
    bcftools = shutil.which("bcftools")
    gatk = shutil.which("gatk")
    if bcftools and smoke_subprocess([bcftools, "--help"], timeout=60):
        logger.info("Germline 调用需排序 BAM + 参考序列；宿主仅探测 bcftools。仿真降级。")
    elif gatk and smoke_subprocess([gatk, "--help"], timeout=180):
        logger.info("GATK 可用；完整 HaplotypeCaller 建议在 Worker。仿真降级。")
    return None


def genomics_germline_calling_impl(file_path: str = "") -> Dict[str, Any]:
    return run_or_degrade(
        lambda: _try_germline(file_path),
        lambda: _degraded_germline(file_path),
        ctx="genomics_germline_calling",
    )


# --- CNV ---


def _degraded_cnv(file_path: str = "") -> Dict[str, Any]:
    b = fastq_bundle_from_path(file_path)
    if b:
        md, tbl, msg = genomics_proxy_cnv(b)
        return _contract(
            "g_cnv",
            msg,
            markdown=md,
            image_urls=[IMG_VOLCANO],
            table_data=tbl,
            extra={"omics_analysis_mode": "fastq_stream_proxy", "proxy_metrics": b},
            tool_id="genomics_cnv_calling",
        )
    md = degraded_banner_md("GATK CNV / 深度矩阵管线") + (
        "### 拷贝数变异 CNV（仿真）\n\n"
        "代表性片段：chr17 局部扩增 log2≈0.62。\n"
    )
    return _contract(
        "g_cnv",
        "CNV (simulated)",
        markdown=md,
        image_urls=[IMG_VOLCANO],
        table_data=simple_rows_table(
            ("locus", "log2_ratio"),
            [
                {"locus": "chr17:7.2-7.6Mb", "log2_ratio": "0.62"},
                {"locus": "chrX:42.1-42.3Mb", "log2_ratio": "-0.41"},
            ],
        ),
        tool_id="genomics_cnv_calling",
    )


def _try_cnv(_file_path: str) -> Optional[Dict[str, Any]]:
    if exe_available("gatk"):
        logger.info("CNV 管线需 Panel-of-Normals；宿主不启动重型分析。仿真降级。")
    return None


def genomics_cnv_calling_impl(file_path: str = "") -> Dict[str, Any]:
    return run_or_degrade(
        lambda: _try_cnv(file_path),
        lambda: _degraded_cnv(file_path),
        ctx="genomics_cnv",
    )


# --- SV ---


def _degraded_sv(file_path: str = "") -> Dict[str, Any]:
    b = fastq_bundle_from_path(file_path)
    if b:
        md, tbl, msg = genomics_proxy_sv(b)
        return _contract(
            "g_sv",
            msg,
            markdown=md,
            image_urls=[IMG_GENOME_BROWSER],
            table_data=tbl,
            extra={"omics_analysis_mode": "fastq_stream_proxy", "proxy_metrics": b},
            tool_id="genomics_sv_calling",
        )
    md = degraded_banner_md("Manta / Lumpy 等 SV 管线") + (
        "### 结构变异 SV（仿真）\n\n"
        "| SV 类型 | 仿真命中 |\n|---------|----------|\n"
        "| DEL | 182 |\n"
        "| DUP | 41 |\n"
        "| INV | 6 |\n"
    )
    return _contract(
        "g_sv",
        "SV calling (simulated)",
        markdown=md,
        image_urls=[IMG_GENOME_BROWSER],
        table_data=simple_rows_table(
            ("sv_type", "count"),
            [
                {"sv_type": "DEL", "count": "182"},
                {"sv_type": "DUP", "count": "41"},
            ],
        ),
        tool_id="genomics_sv_calling",
    )


def _try_sv(_file_path: str) -> Optional[Dict[str, Any]]:
    if shutil.which("runWorkflow.py") or shutil.which("configManta.py"):
        logger.info("检测到 Manta 脚本路径不代表完整配置；仿真降级。")
    return None


def genomics_sv_calling_impl(file_path: str = "") -> Dict[str, Any]:
    return run_or_degrade(
        lambda: _try_sv(file_path),
        lambda: _degraded_sv(file_path),
        ctx="genomics_sv_calling",
    )


# --- VQSR ---


def _degraded_vqsr(file_path: str = "") -> Dict[str, Any]:
    b = fastq_bundle_from_path(file_path)
    if b:
        md, tbl, msg = genomics_proxy_vqsr(b)
        return _contract(
            "g_vqsr",
            msg,
            markdown=md,
            image_urls=[IMG_VOLCANO],
            table_data=tbl,
            extra={"omics_analysis_mode": "fastq_stream_proxy", "proxy_metrics": b},
            tool_id="genomics_vqsr_filtering",
        )
    md = degraded_banner_md("GATK VQSR / bcftools filter 与训练资源") + (
        "### VQSR + 位点标准化（仿真）\n\n"
        "- PASS 变异：**3,215,902**\n"
        "- 过滤 LowQual：**214,330**\n"
    )
    return _contract(
        "g_vqsr",
        "VQSR (simulated)",
        markdown=md,
        image_urls=[IMG_VOLCANO],
        table_data=simple_rows_table(
            ("filter", "count"),
            [
                {"filter": "PASS", "count": "3215902"},
                {"filter": "LowQual", "count": "214330"},
            ],
        ),
        tool_id="genomics_vqsr_filtering",
    )


def _try_vqsr(_file_path: str) -> Optional[Dict[str, Any]]:
    bt = shutil.which("bcftools")
    if bt and smoke_subprocess([bt, "--help"], timeout=60):
        logger.info("bcftools 可用；VQSR 过滤需完整 VCF 与训练集；仿真降级。")
    return None


def genomics_vqsr_filtering_impl(file_path: str = "") -> Dict[str, Any]:
    return run_or_degrade(
        lambda: _try_vqsr(file_path),
        lambda: _degraded_vqsr(file_path),
        ctx="genomics_vqsr_filtering",
    )


# --- annotation ---


def _degraded_annotation(file_path: str = "") -> Dict[str, Any]:
    b = fastq_bundle_from_path(file_path)
    if b:
        md, tbl, msg = genomics_proxy_annotation(b)
        out = _contract(
            "g_anno",
            msg,
            markdown=md,
            image_urls=[IMG_VOLCANO],
            table_data=tbl,
            extra={"omics_analysis_mode": "fastq_stream_proxy", "proxy_metrics": b},
            tool_id="genomics_variant_annotation",
        )
        out["annotated_vcf_path"] = out["output_path"]
        return out
    md = degraded_banner_md("VEP / snpEff 与缓存数据库") + (
        "### 变异注释（仿真）\n\n"
        "| Gene | Consequence | AF（gnomAD） |\n"
        "|------|-------------|---------------|\n"
        "| BRCA2 | missense_variant | 0.00012 |\n"
        "| TP53 | synonymous_variant | 0.41 |\n"
        "| EGFR | inframe_insertion | 0.00003 |\n"
    )
    out = _contract(
        "g_anno",
        "Annotation (simulated)",
        markdown=md,
        image_urls=[IMG_VOLCANO],
        table_data=simple_rows_table(
            ("gene", "consequence", "gnomad_af"),
            [
                {"gene": "BRCA2", "consequence": "missense_variant", "gnomad_af": "0.00012"},
                {"gene": "TP53", "consequence": "synonymous_variant", "gnomad_af": "0.41"},
                {"gene": "EGFR", "consequence": "inframe_insertion", "gnomad_af": "0.00003"},
            ],
        ),
        tool_id="genomics_variant_annotation",
    )
    out["annotated_vcf_path"] = out["output_path"]
    return out


def _try_annotation(_file_path: str) -> Optional[Dict[str, Any]]:
    vep = shutil.which("vep")
    if vep:
        logger.info("VEP 可执行文件存在；完整注释需缓存与输入 VCF。仿真降级。")
    return None


def genomics_variant_annotation_impl(file_path: str = "") -> Dict[str, Any]:
    return run_or_degrade(
        lambda: _try_annotation(file_path),
        lambda: _degraded_annotation(file_path),
        ctx="genomics_variant_annotation",
    )


# --- ACMG ---


def _degraded_acmg(file_path: str = "") -> Dict[str, Any]:
    b = fastq_bundle_from_path(file_path)
    if b:
        md, tbl, msg = genomics_proxy_acmg(b)
        return _contract(
            "g_acmg",
            msg,
            markdown=md,
            image_urls=[IMG_DNA_HELIX],
            table_data=tbl,
            extra={"omics_analysis_mode": "fastq_stream_proxy", "proxy_metrics": b},
            tool_id="genomics_acmg_classification",
        )
    md = degraded_banner_md("ACMG 规则引擎 / ClinVar 本地化镜像") + (
        "### ACMG/AMP 致病性分级（仿真）\n\n"
        "| 变异 | 分级 |\n|------|------|\n"
        "| NM_000059.3:c.1234A>G | LP |\n"
        "| NC_000017.10:g.7577121G>A | P |\n"
        "| NM_007294.4:c.5266dup | VUS |\n"
    )
    return _contract(
        "g_acmg",
        "ACMG (simulated)",
        markdown=md,
        image_urls=[IMG_DNA_HELIX],
        table_data=simple_rows_table(
            ("variant", "acmg_tier"),
            [
                {"variant": "BRCA2 c.1234A>G", "acmg_tier": "LP"},
                {"variant": "TP53 g.7577121G>A", "acmg_tier": "P"},
                {"variant": "ATM c.5266dup", "acmg_tier": "VUS"},
            ],
        ),
        tool_id="genomics_acmg_classification",
    )


def _try_acmg(_file_path: str) -> Optional[Dict[str, Any]]:
    return None


def genomics_acmg_classification_impl(file_path: str = "") -> Dict[str, Any]:
    return run_or_degrade(
        lambda: _try_acmg(file_path),
        lambda: _degraded_acmg(file_path),
        ctx="genomics_acmg_classification",
    )
