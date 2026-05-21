"""
基因组管线步骤输出：核心指标前置 + Base64 图表 + CLI 日志折叠。
供 omics_genomics_runner / omics_genomics_pipeline_tools 组装 Markdown。
"""
from __future__ import annotations

import base64
import io
import json
import logging
import os
import subprocess
from typing import Any, Dict, List, Optional, Sequence, Tuple

from .omics_pipeline_env import clip_omics_log

logger = logging.getLogger(__name__)

MetricRow = Tuple[str, str]


def metrics_table_md(rows: Sequence[MetricRow], *, title: str = "核心指标摘要") -> str:
    if not rows:
        return ""
    lines = [f"### {title}", "", "| 指标 | 值 |", "| --- | --- |"]
    for k, v in rows:
        lines.append(f"| **{k}** | {v} |")
    lines.append("")
    return "\n".join(lines)


def wrap_cli_logs_markdown(
    stdout: str = "",
    stderr: str = "",
    *,
    summary_label: str = "🔍 点击查看该步骤底层真实运行日志 (CLI Logs)",
    max_chars: int = 12000,
) -> str:
    """将 stdout/stderr 包进 HTML details，避免铺满右栏。"""
    body_parts: List[str] = []
    if stderr and stderr.strip():
        body_parts.append(f"**stderr**\n\n```plaintext\n{clip_omics_log(stderr, max_chars)}\n```")
    if stdout and stdout.strip():
        body_parts.append(f"**stdout**\n\n```plaintext\n{clip_omics_log(stdout, max_chars)}\n```")
    if not body_parts:
        body_parts.append("_（无额外 CLI 文本输出）_")
    inner = "\n\n".join(body_parts)
    return (
        f'\n<details class="omics-cli-logs">\n<summary>{summary_label}</summary>\n\n'
        f"{inner}\n\n</details>\n"
    )


def cli_logs_from_process(
    cp: subprocess.CompletedProcess,
    cmd: Optional[List[str]] = None,
    *,
    note: str = "",
    extra_lines: str = "",
) -> str:
    cmd_s = " ".join(str(c) for c in cmd) if cmd else "(unknown)"
    header = (
        f"- **Command executed**: `{cmd_s}`\n"
        f"- **Exit code**: `{cp.returncode}`\n"
    )
    if note:
        header += f"- **Result note**: {note}\n"
    if extra_lines:
        header += extra_lines.rstrip() + "\n"
    return header + wrap_cli_logs_markdown(cp.stdout or "", cp.stderr or "")


def base64_png_markdown(b64: str, *, alt: str = "chart", width_px: int = 520) -> str:
    if not b64:
        return ""
    return f'\n<div align="center">\n\n![{alt}](data:image/png;base64,{b64})\n\n</div>\n'


def _fig_to_base64(fig) -> str:
    import matplotlib.pyplot as plt

    buf = io.BytesIO()
    fig.savefig(buf, format="png", bbox_inches="tight", dpi=120)
    buf.seek(0)
    b64 = base64.b64encode(buf.read()).decode("utf-8")
    plt.close(fig)
    return b64


def plot_fastp_quality_curve_md(fastp_json_path: str) -> str:
    """从 fastp.json 解析 read quality mean 曲线；失败则返回空。"""
    path = (fastp_json_path or "").strip()
    if not path or not os.path.isfile(path):
        return ""
    try:
        with open(path, encoding="utf-8") as fh:
            data = json.load(fh)
    except (OSError, json.JSONDecodeError) as exc:
        logger.debug("fastp.json 解析失败: %s", exc)
        return ""

    means: Optional[List[float]] = None
    for key in ("read1_before_filtering", "read2_before_filtering", "read1_after_filtering"):
        block = data.get(key) or {}
        qc = block.get("quality_curves") or {}
        raw = qc.get("mean")
        if isinstance(raw, list) and raw:
            means = [float(x) for x in raw]
            break
    if not means:
        return ""

    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        fig, ax = plt.subplots(figsize=(6.5, 3.2))
        positions = list(range(1, len(means) + 1))
        ax.plot(positions, means, color="#2563eb", linewidth=1.8, label="Mean Phred")
        ax.axhline(30, color="#16a34a", linestyle="--", linewidth=1, alpha=0.7, label="Q30")
        ax.axhline(20, color="#ca8a04", linestyle=":", linewidth=1, alpha=0.6, label="Q20")
        ax.set_xlabel("Read position (bp)")
        ax.set_ylabel("Mean quality (Phred)")
        ax.set_title("Per-base sequencing quality (fastp)")
        ax.legend(loc="lower right", fontsize=8)
        ax.grid(True, alpha=0.25)
        ax.set_ylim(bottom=0)
        b64 = _fig_to_base64(fig)
        return base64_png_markdown(b64, alt="fastp_quality_curve")
    except Exception as exc:  # noqa: BLE001
        logger.warning("fastp 质量曲线绘图失败: %s", exc)
        return ""


def plot_streaming_quality_curve_md(fastq_path: str, *, max_reads: int = 1200) -> str:
    """流式抽样 FASTQ 绘制位置平均质量（用于 raw_qc）。"""
    path = (fastq_path or "").strip()
    if not path or not os.path.isfile(path):
        return ""
    import gzip

    lower = path.lower()
    opener = gzip.open if lower.endswith(".gz") else open  # type: ignore[assignment]
    mode = "rt" if lower.endswith(".gz") else "r"
    max_pos = 150
    sums = [0.0] * max_pos
    counts = [0] * max_pos
    n_reads = 0
    try:
        with opener(path, mode, encoding="ascii", errors="replace") as fh:  # type: ignore[arg-type]
            while n_reads < max_reads:
                if not fh.readline():
                    break
                seq = fh.readline()
                fh.readline()
                qual = fh.readline()
                if not qual:
                    break
                qline = qual.rstrip("\n\r")
                for i, qch in enumerate(qline[:max_pos]):
                    q = max(0, ord(qch) - 33)
                    sums[i] += q
                    counts[i] += 1
                n_reads += 1
    except OSError as exc:
        logger.debug("FASTQ 质量曲线抽样失败: %s", exc)
        return ""

    means = [sums[i] / counts[i] if counts[i] else 0.0 for i in range(max_pos)]
    while means and means[-1] == 0:
        means.pop()
    if not means:
        return ""

    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        fig, ax = plt.subplots(figsize=(6.5, 3.2))
        ax.plot(range(1, len(means) + 1), means, color="#7c3aed", linewidth=1.8)
        ax.axhline(30, color="#16a34a", linestyle="--", linewidth=1, alpha=0.7)
        ax.set_xlabel("Read position (bp)")
        ax.set_ylabel("Mean quality (Phred)")
        ax.set_title(f"Estimated per-base quality (sampled {n_reads} reads)")
        ax.grid(True, alpha=0.25)
        ax.set_ylim(bottom=0)
        b64 = _fig_to_base64(fig)
        return base64_png_markdown(b64, alt="fastq_quality_curve")
    except Exception as exc:  # noqa: BLE001
        logger.warning("FASTQ 质量曲线绘图失败: %s", exc)
        return ""


def plot_variant_composition_pie_md(variant_summary: Dict[str, Any]) -> str:
    """变异组成饼图：SNP / Indel / Other（基于 bcftools stats 或计数）。"""
    if not variant_summary:
        return ""
    total = int(variant_summary.get("variant_count") or 0)
    snp = int(variant_summary.get("snp_count") or 0)
    indel = int(variant_summary.get("indel_count") or 0)
    if total <= 0 and snp <= 0 and indel <= 0:
        labels = ["No variants called", "Background"]
        sizes = [1.0, 0.0001]
        colors = ["#94a3b8", "#e2e8f0"]
    else:
        other = max(0, total - snp - indel)
        labels = []
        sizes = []
        colors = []
        if snp > 0:
            labels.append(f"SNP ({snp})")
            sizes.append(float(snp))
            colors.append("#2563eb")
        if indel > 0:
            labels.append(f"Indel ({indel})")
            sizes.append(float(indel))
            colors.append("#dc2626")
        if other > 0:
            labels.append(f"Other ({other})")
            sizes.append(float(other))
            colors.append("#64748b")
        if not sizes:
            labels = [f"Records ({total})"]
            sizes = [float(total or 1)]
            colors = ["#2563eb"]

    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        fig, ax = plt.subplots(figsize=(5.5, 3.5))
        ax.pie(
            sizes,
            labels=labels,
            colors=colors[: len(sizes)],
            autopct="%1.1f%%",
            startangle=90,
            textprops={"fontsize": 9},
        )
        titv = variant_summary.get("tstv_ratio")
        title = "Variant composition"
        if titv:
            title += f" (Ti/Tv ≈ {titv})"
        ax.set_title(title)
        b64 = _fig_to_base64(fig)
        return base64_png_markdown(b64, alt="variant_composition")
    except Exception as exc:  # noqa: BLE001
        logger.warning("变异饼图绘图失败: %s", exc)
        return ""


def assemble_genomics_step_markdown(
    *,
    title: str,
    metrics: Optional[Sequence[MetricRow]] = None,
    chart_md: str = "",
    body_md: str = "",
    cmd: Optional[List[str]] = None,
    cp: Optional[subprocess.CompletedProcess] = None,
    cli_note: str = "",
    cli_extra: str = "",
    footer_note: str = "",
) -> str:
    """
    全模态节点输出契约：指标表 → 图表 → 业务说明 → 折叠 CLI 日志。
    """
    parts: List[str] = [f"## {title}", ""]
    if metrics:
        parts.append(metrics_table_md(metrics))
    if chart_md:
        parts.append(chart_md.strip())
        parts.append("")
    if body_md:
        parts.append(body_md.strip())
        parts.append("")
    if footer_note:
        parts.append(footer_note.strip())
        parts.append("")
    if cp is not None:
        parts.append(
            cli_logs_from_process(cp, cmd, note=cli_note, extra_lines=cli_extra).strip()
        )
    return "\n".join(parts).strip() + "\n"
