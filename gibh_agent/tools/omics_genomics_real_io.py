"""
基因组真实管线 I/O：参考序列发现、上游产物校验、VCF/FASTQ 摘要（供专家报告）。
"""
from __future__ import annotations

import glob
import json
import os
import re
import subprocess
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

from .omics_pipeline_env import bwa_index_ready, get_omics_ref_dir, resolve_reference_fasta


def discover_genomics_reference() -> Dict[str, Any]:
    """
    扫描仓库 data/references、/tmp/omics_test_data 等，返回最佳 hg38 FASTA + 索引状态。
    """
    candidates: List[Path] = []
    roots = [
        Path(get_omics_ref_dir()) / "genomics",
        Path("/tmp/omics_test_data/ref"),
        Path(os.environ.get("GIBH_REF_HG38", "")).parent
        if os.environ.get("GIBH_REF_HG38")
        else Path("."),
    ]
    for root in roots:
        if not root.is_dir():
            continue
        for pat in ("hg38.fa", "chr21.fa", "*.fa", "*.fasta"):
            candidates.extend(Path(p) for p in glob.glob(str(root / pat)))

    seen: set[str] = set()
    ranked: List[Tuple[int, str]] = []
    for p in candidates:
        sp = str(p.resolve())
        if sp in seen or not p.is_file() or p.stat().st_size < 1000:
            continue
        seen.add(sp)
        score = 0
        if bwa_index_ready(sp):
            score += 100
        if "hg38" in p.name.lower() or "chr21" in p.name.lower():
            score += 10
        ranked.append((score, sp))
    ranked.sort(key=lambda x: (-x[0], -os.path.getsize(x[1])))

    best = ranked[0][1] if ranked else ""
    env_ref = (os.environ.get("GIBH_REF_HG38") or "").strip()
    resolved = resolve_reference_fasta("hg38")
    fasta = resolved or (env_ref if env_ref and os.path.isfile(env_ref) else best)
    return {
        "fasta": fasta,
        "bwa_index_ready": bwa_index_ready(fasta) if fasta else False,
        "fai": f"{fasta}.fai" if fasta and os.path.isfile(f"{fasta}.fai") else "",
        "bwt": f"{fasta}.bwt" if fasta and os.path.isfile(f"{fasta}.bwt") else "",
        "candidates_scanned": len(ranked),
        "omics_ref_dir": get_omics_ref_dir(),
    }


def artifact_error(
    label: str,
    path: str,
    *,
    reason: str,
    stderr: str = "",
) -> Dict[str, Any]:
    md = (
        f"### {label} — 上游产物校验失败\n\n"
        f"- **路径**: `{path or '(empty)'}`\n"
        f"- **原因**: {reason}\n"
    )
    if stderr.strip():
        md += "\n```\n" + stderr[:8000] + "\n```\n"
    return {
        "status": "error",
        "error_category": "data_issue",
        "message": f"{label}: {reason}",
        "user_message": f"{label} 失败：{reason}",
        "fatal_env_prereq": False,
        "can_skip": False,
        "markdown": md,
    }


def require_nonempty_file(
    path: str,
    label: str,
    *,
    extensions: Optional[Tuple[str, ...]] = None,
) -> Optional[Dict[str, Any]]:
    p = (path or "").strip()
    if not p:
        return artifact_error(label, p, reason="路径为空")
    if not os.path.isfile(p):
        return artifact_error(label, p, reason="文件不存在")
    if os.path.getsize(p) <= 0:
        return artifact_error(label, p, reason="文件大小为 0 字节")
    if extensions:
        low = p.lower()
        if not any(low.endswith(ext) for ext in extensions):
            return artifact_error(
                label,
                p,
                reason=f"扩展名不符合预期 {extensions}",
            )
    return None


def ensure_bam_indexed(samtools_exe: str, bam_path: str) -> Optional[Dict[str, Any]]:
    bam = (bam_path or "").strip()
    err = require_nonempty_file(bam, "BAM（索引前）", extensions=(".bam",))
    if err:
        return err
    bai = f"{bam}.bai"
    if os.path.isfile(bai) and os.path.getsize(bai) > 0:
        return None
    cmd = [samtools_exe, "index", bam]
    try:
        cp = subprocess.run(cmd, capture_output=True, text=True, timeout=7200, check=False)
    except (OSError, subprocess.TimeoutExpired) as exc:
        return artifact_error("samtools index", bam, reason=str(exc))
    if cp.returncode != 0:
        return artifact_error(
            "samtools index",
            bam,
            reason=f"exit {cp.returncode}",
            stderr=cp.stderr or "",
        )
    if not os.path.isfile(bai) or os.path.getsize(bai) <= 0:
        return artifact_error("samtools index", bai, reason="未生成 .bai 或为空")
    return None


def vcf_tabix_index_present(vcf_path: str) -> bool:
    p = (vcf_path or "").strip()
    if not p:
        return False
    if os.path.isfile(f"{p}.tbi") and os.path.getsize(f"{p}.tbi") > 0:
        return True
    if os.path.isfile(f"{p}.csi") and os.path.getsize(f"{p}.csi") > 0:
        return True
    return False


def ensure_vcf_indexed(bcftools_exe: str, vcf_path: str) -> Optional[Dict[str, Any]]:
    """GATK 读取 block-gzip VCF 须 .tbi/.csi；bcftools 产出默认无索引。"""
    vcf = (vcf_path or "").strip()
    err = require_nonempty_file(vcf, "VCF（建索引前）", extensions=(".vcf", ".vcf.gz"))
    if err:
        return err
    if vcf_tabix_index_present(vcf):
        return None
    if not vcf.lower().endswith(".vcf.gz"):
        return None
    cmd = [bcftools_exe, "index", "-t", vcf]
    try:
        cp = subprocess.run(cmd, capture_output=True, text=True, timeout=3600, check=False)
    except (OSError, subprocess.TimeoutExpired) as exc:
        return artifact_error("bcftools index", vcf, reason=str(exc))
    if cp.returncode != 0:
        return artifact_error(
            "bcftools index",
            vcf,
            reason=f"exit {cp.returncode}",
            stderr=cp.stderr or "",
        )
    if not vcf_tabix_index_present(vcf):
        return artifact_error("bcftools index", vcf, reason="未生成 .tbi/.csi")
    return None


def vcf_header_has_info_ids(vcf_path: str, info_ids: Tuple[str, ...]) -> bool:
    """检查 VCF 头是否声明给定 INFO（GATK HardFilter 需 QD/FS 等）。"""
    from .omics_pipeline_env import resolve_cli_exe

    vcf = (vcf_path or "").strip()
    if not vcf or not info_ids:
        return False
    bcftools = resolve_cli_exe("bcftools")
    if not bcftools:
        return False
    try:
        cp = subprocess.run(
            [bcftools, "view", "-h", vcf],
            capture_output=True,
            text=True,
            timeout=600,
            check=False,
        )
    except (OSError, subprocess.TimeoutExpired):
        return False
    if cp.returncode != 0:
        return False
    header = cp.stdout or ""
    for iid in info_ids:
        if f"ID={iid}," not in header:
            return False
    return True


def load_fastp_summary(json_path: str) -> Dict[str, Any]:
    out: Dict[str, Any] = {"json_path": json_path}
    try:
        with open(json_path, encoding="utf-8") as fh:
            data = json.load(fh)
        summ = data.get("summary") or {}
        after = summ.get("after_filtering") or {}
        before = summ.get("before_filtering") or {}
        out["after_filtering"] = after
        out["before_filtering"] = before
        out["total_reads_after"] = after.get("total_reads")
        out["q30_rate_after"] = after.get("q30_rate")
        out["gc_content_after"] = after.get("gc_content")
    except (OSError, json.JSONDecodeError) as exc:
        out["parse_error"] = str(exc)
    return out


def summarize_vcf(vcf_path: str) -> Dict[str, Any]:
    """用 bcftools stats 或简单解析得到变异摘要。"""
    p = (vcf_path or "").strip()
    summary: Dict[str, Any] = {"vcf_path": p, "variant_count": 0}
    if not p or not os.path.isfile(p):
        summary["error"] = "VCF 不存在"
        return summary

    bcftools = None
    for name in ("bcftools",):
        from .omics_pipeline_env import resolve_cli_exe

        bcftools = resolve_cli_exe(name)
        if bcftools:
            break

    if bcftools:
        try:
            cp = subprocess.run(
                [bcftools, "stats", "-s", "-", p],
                capture_output=True,
                text=True,
                timeout=600,
                check=False,
            )
            if cp.returncode == 0 and cp.stdout:
                summary["bcftools_stats_excerpt"] = cp.stdout[:12000]
                m = re.search(r"SN\s+\t0\s+number of records:\s+(\d+)", cp.stdout)
                if m:
                    summary["variant_count"] = int(m.group(1))
                m2 = re.search(r"SN\s+\t0\s+number of SNPs:\s+(\d+)", cp.stdout)
                if m2:
                    summary["snp_count"] = int(m2.group(1))
                m3 = re.search(r"SN\s+\t0\s+number of indels:\s+(\d+)", cp.stdout)
                if m3:
                    summary["indel_count"] = int(m3.group(1))
                titv = re.search(
                    r"Tstv\s+.*?\t(\d+)\s+.*?\t([\d.]+)",
                    cp.stdout,
                    re.DOTALL,
                )
                if titv:
                    summary["tstv_ratio"] = titv.group(2)
                return summary
        except (OSError, subprocess.TimeoutExpired):
            pass

    # 回退：统计非 header 行
    import gzip

    opener = gzip.open if p.lower().endswith(".gz") else open  # type: ignore[assignment]
    try:
        n = 0
        with opener(p, "rt", encoding="utf-8", errors="replace") as fh:  # type: ignore[arg-type]
            for line in fh:
                if line.startswith("#"):
                    continue
                n += 1
                if n > 5_000_000:
                    break
        summary["variant_count"] = n
        summary["count_method"] = "line_parse"
    except OSError as exc:
        summary["error"] = str(exc)
    return summary


def format_variant_summary_md(summary: Dict[str, Any]) -> str:
    if summary.get("error"):
        return f"- **VCF 摘要错误**: {summary['error']}\n"
    lines = [
        f"- **VCF 路径**: `{summary.get('vcf_path', '')}`",
        f"- **变异记录数**: **{summary.get('variant_count', 0)}**",
    ]
    if "snp_count" in summary:
        lines.append(f"- **SNP 数（bcftools stats）**: {summary['snp_count']}")
    if "indel_count" in summary:
        lines.append(f"- **Indel 数**: {summary['indel_count']}")
    if "tstv_ratio" in summary:
        lines.append(f"- **Ti/Tv 比例**: {summary['tstv_ratio']}")
    return "\n".join(lines) + "\n"
