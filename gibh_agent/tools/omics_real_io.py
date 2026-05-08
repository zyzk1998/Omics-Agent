"""
组学原始文件真实 I/O 与轻量级统计（纯 Python / 标准库），供 pipeline 首步质控使用。

设计原则：不依赖容器内 FastQC/fastp；保证数据真实流经读取逻辑并产出可核验指标。
"""
from __future__ import annotations

import gzip
import logging
import os
import xml.etree.ElementTree as ET
from typing import Any, Dict, Optional

logger = logging.getLogger(__name__)

_DEFAULT_FASTQ_STREAM_CAP = 100_000


def _local_tag(elem_tag: str) -> str:
    if "}" in elem_tag:
        return elem_tag.rsplit("}", 1)[-1]
    return elem_tag


def compute_fastq_stats(file_path: str, max_reads: Optional[int] = None) -> Dict[str, Any]:
    """
    流式读取 FASTQ（支持 .gz），统计 reads 数、碱基总数、GC%、平均读长。
    """
    path = os.path.abspath(file_path)
    if not os.path.isfile(path):
        raise FileNotFoundError(path)

    lower = path.lower()
    use_gzip = lower.endswith(".gz")
    open_fn = gzip.open if use_gzip else open  # type: ignore[assignment]
    mode = "rt" if use_gzip else "r"

    n_reads = 0
    total_bases = 0
    gc_bases = 0
    len_sum = 0

    with open_fn(path, mode, encoding="ascii", errors="replace") as fh:  # type: ignore[misc]
        while True:
            header = fh.readline()
            if not header:
                break
            seq = fh.readline()
            plus = fh.readline()
            qual = fh.readline()
            if not seq or not qual:
                raise ValueError(f"FASTQ 不完整或截断（文件: {path}）")
            seq = seq.strip()
            if not seq:
                continue
            n_reads += 1
            slen = len(seq)
            len_sum += slen
            total_bases += slen
            su = seq.upper()
            gc_bases += sum(1 for b in su if b in ("G", "C"))

            if max_reads is not None and n_reads >= max_reads:
                break

    gc_frac = (gc_bases / total_bases) if total_bases else 0.0
    mean_len = (len_sum / n_reads) if n_reads else 0.0

    return {
        "n_reads": n_reads,
        "total_bases": int(total_bases),
        "gc_fraction": round(gc_frac, 6),
        "gc_percent": round(gc_frac * 100.0, 4),
        "mean_read_length": round(mean_len, 4),
        "input_path": path,
        "compression": "gzip" if use_gzip else "none",
    }


def compute_fastq_stream_bundle(
    file_path: str,
    max_reads: Optional[int] = None,
) -> Dict[str, Any]:
    """
    单次流式扫描 FASTQ：碱基统计 + 质量值统计（Phred），可选最多读取 max_reads 条以控制耗时。
    用于无 BWA/GATK 时的下游「轻量代理」指标；非完整文件时应查看 truncated。
    """
    path = os.path.abspath(file_path)
    if not os.path.isfile(path):
        raise FileNotFoundError(path)
    lower = path.lower()
    use_gzip = lower.endswith(".gz")
    open_fn = gzip.open if use_gzip else open  # type: ignore[assignment]
    mode = "rt" if use_gzip else "r"
    cap = max_reads if max_reads is not None else _DEFAULT_FASTQ_STREAM_CAP

    n_reads = 0
    total_bases = 0
    gc_bases = 0
    len_sum = 0
    qual_sum = 0
    qual_bases = 0
    q30_bases = 0
    n_bases = 0
    truncated_early = False

    with open_fn(path, mode, encoding="ascii", errors="replace") as fh:  # type: ignore[misc]
        while True:
            header = fh.readline()
            if not header:
                break
            seq = fh.readline()
            plus = fh.readline()
            qual = fh.readline()
            if not seq or not qual:
                raise ValueError(f"FASTQ 不完整或截断（文件: {path}）")
            seq = seq.strip()
            if not seq:
                continue
            n_reads += 1
            slen = len(seq)
            len_sum += slen
            total_bases += slen
            su = seq.upper()
            gc_bases += sum(1 for b in su if b in ("G", "C"))
            n_bases += sum(1 for b in su if b == "N")

            qline = qual.rstrip("\n\r")
            for qch in qline[:slen]:
                q = ord(qch) - 33
                if q < 0:
                    q = 0
                qual_sum += q
                qual_bases += 1
                if q >= 30:
                    q30_bases += 1

            if cap is not None and n_reads >= cap:
                truncated_early = True
                break

    gc_frac = (gc_bases / total_bases) if total_bases else 0.0
    mean_len = (len_sum / n_reads) if n_reads else 0.0
    mean_qual = (qual_sum / qual_bases) if qual_bases else 0.0
    q30_frac = (q30_bases / qual_bases) if qual_bases else 0.0
    n_frac = (n_bases / total_bases) if total_bases else 0.0

    truncated = bool(cap is not None and truncated_early)

    return {
        "n_reads": n_reads,
        "total_bases": int(total_bases),
        "gc_fraction": round(gc_frac, 6),
        "gc_percent": round(gc_frac * 100.0, 4),
        "mean_read_length": round(mean_len, 4),
        "input_path": path,
        "compression": "gzip" if use_gzip else "none",
        "qual_mean": round(mean_qual, 4),
        "q30_fraction": round(q30_frac, 6),
        "n_fraction": round(n_frac, 6),
        "stream_cap": cap,
        "truncated": truncated,
    }


def compute_mzml_basic_stats(file_path: str) -> Dict[str, Any]:
    """
    使用 iterparse 扫描 mzML：谱图数、色谱数、文件大小（不加载整文件到内存）。
    """
    path = os.path.abspath(file_path)
    if not os.path.isfile(path):
        raise FileNotFoundError(path)
    lower = path.lower()
    if not lower.endswith(".mzml"):
        raise ValueError(f"当前实现仅支持 .mzML 文件: {path}")

    size_mb = os.path.getsize(path) / (1024.0 * 1024.0)
    spectrum_count = 0
    chromatogram_count = 0

    # iterparse：大文件友好
    for event, elem in ET.iterparse(path, events=("end",)):
        tag = _local_tag(elem.tag)
        if tag == "spectrum":
            spectrum_count += 1
        elif tag == "chromatogram":
            chromatogram_count += 1
        elem.clear()

    return {
        "file_size_mb": round(size_mb, 6),
        "spectrum_count": spectrum_count,
        "chromatogram_count": chromatogram_count,
        "input_path": path,
    }


def compute_mzml_derived_stats(file_path: str) -> Dict[str, Any]:
    """
    单次 iterparse 汇总 centroid 谱：峰数（defaultArrayLength）、TIC、保留时间窗、m/z 扫描窗等。
    用于无搜库引擎时的下游定量/差异代理指标（基于文件实测字段）。
    """
    path = os.path.abspath(file_path)
    if not os.path.isfile(path):
        raise FileNotFoundError(path)
    if not path.lower().endswith(".mzml"):
        raise ValueError(f"当前实现仅支持 .mzML 文件: {path}")

    total_peaks = 0
    spec_n = 0
    tics: list[float] = []
    bp_mzs: list[float] = []
    rt_vals: list[float] = []
    low_mzs: list[float] = []
    hi_mzs: list[float] = []

    for _event, elem in ET.iterparse(path, events=("end",)):
        if _local_tag(elem.tag) != "spectrum":
            elem.clear()
            continue
        spec_n += 1
        try:
            total_peaks += int(elem.get("defaultArrayLength") or 0)
        except ValueError:
            pass
        tic_v = bp_v = rt_v = lowmz_v = himz_v = None
        for cv in elem.iter():
            if _local_tag(cv.tag) != "cvParam":
                continue
            acc = cv.get("accession")
            val = cv.get("value")
            name = (cv.get("name") or "").lower()
            if not val:
                continue
            try:
                fv = float(val)
            except ValueError:
                continue
            if acc == "MS:1000285":
                tic_v = fv
            elif acc == "MS:1000504":
                bp_v = fv
            elif acc == "MS:1000016" and "time" in name:
                rt_v = fv
            elif acc == "MS:1000528":
                lowmz_v = fv
            elif acc == "MS:1000527":
                himz_v = fv
        if tic_v is not None:
            tics.append(tic_v)
        if bp_v is not None:
            bp_mzs.append(bp_v)
        if rt_v is not None:
            rt_vals.append(rt_v)
        if lowmz_v is not None:
            low_mzs.append(lowmz_v)
        if himz_v is not None:
            hi_mzs.append(himz_v)
        elem.clear()

    def _mean(xs: list[float]) -> float:
        return sum(xs) / len(xs) if xs else 0.0

    sum_tic = sum(tics)
    size_mb = os.path.getsize(path) / (1024.0 * 1024.0)
    return {
        "input_path": path,
        "spectrum_count": spec_n,
        "total_centroid_peaks": total_peaks,
        "mean_peaks_per_spectrum": round(total_peaks / spec_n, 4) if spec_n else 0.0,
        "tic_sum": round(sum_tic, 2),
        "tic_mean": round(_mean(tics), 4),
        "tic_min": round(min(tics), 4) if tics else 0.0,
        "tic_max": round(max(tics), 4) if tics else 0.0,
        "base_peak_mz_mean": round(_mean(bp_mzs), 6),
        "rt_span_sec": round(max(rt_vals) - min(rt_vals), 4) if len(rt_vals) > 1 else 0.0,
        "rt_min": round(min(rt_vals), 4) if rt_vals else 0.0,
        "rt_max": round(max(rt_vals), 4) if rt_vals else 0.0,
        "mz_window_low_min": round(min(low_mzs), 6) if low_mzs else 0.0,
        "mz_window_high_max": round(max(hi_mzs), 6) if hi_mzs else 0.0,
        "file_size_mb": round(size_mb, 6),
    }


def format_fastq_qc_markdown(stats: Dict[str, Any]) -> str:
    return (
        "## 原始测序数据质控（真实统计）\n\n"
        "以下指标由服务进程**直接读取 FASTQ/FASTQ.GZ** 计算，非占位 Mock。\n\n"
        "| 指标 | 值 |\n"
        "| --- | --- |\n"
        f"| Reads | **{stats['n_reads']}** |\n"
        f"| Total bases | {stats['total_bases']} |\n"
        f"| Mean read length | {stats['mean_read_length']} |\n"
        f"| GC content | **{stats['gc_percent']}%** |\n"
        f"| Input | `{stats.get('input_path', '')}` |\n"
    )


def format_mzml_qc_markdown(stats: Dict[str, Any]) -> str:
    return (
        "## 原始质谱文件检视（真实统计）\n\n"
        "以下指标由 **XML iterparse** 扫描 mzML 得到（谱图/色谱计数），非占位 Mock。\n\n"
        "| 指标 | 值 |\n"
        "| --- | --- |\n"
        f"| Spectrum count | **{stats['spectrum_count']}** |\n"
        f"| Chromatogram count | {stats['chromatogram_count']} |\n"
        f"| File size (MB) | {stats['file_size_mb']} |\n"
        f"| Input | `{stats.get('input_path', '')}` |\n"
    )
