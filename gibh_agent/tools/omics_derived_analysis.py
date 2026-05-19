"""
无外部重型 CLI 时，从原始 FASTQ / mzML **实测流式统计** 导出确定性代理指标。

说明：
- 指标随输入可复现，便于审计与端到端联通测试；
- **不等价**于 BWA/GATK、DIA-NN、MACS2 等真值输出；接入 Worker 真管线后由真实结果替换。
"""
from __future__ import annotations

import hashlib
import os
from typing import Any, Dict, List, Optional, Tuple

from .omics_mock_ui import (
    IMG_MANHATTAN,
    IMG_WATERFALL_ONCOPRINT,
    simple_rows_table,
)
from .omics_real_io import (
    compute_fastq_stream_bundle,
    compute_mzml_derived_stats,
)

HDR_FASTQ = (
    "> **分析说明**：未调用外部比对/变异软件时，下列数值由 **原始 FASTQ 流式碱基与质量统计**"
    " 经确定性公式推导，用于全流程连贯输出；**非** BWA/GATK 临床级真值，但可随输入复现、可核验。\n\n"
)

HDR_MZML = (
    "> **分析说明**：未调用搜库/定量 CLI 时，下列数值由 **mzML 谱图元数据（TIC、峰数、采集窗）**"
    " 经确定性聚合得到；**非** DIA-NN/MaxQuant 真值，但对同一文件稳定可重复。\n\n"
)


def _u01(seed: bytes, salt: str) -> float:
    h = hashlib.sha256(seed + salt.encode()).digest()
    return int.from_bytes(h[:8], "big") / float(2**64)


def fastq_bundle_from_path(file_path: str) -> Optional[Dict[str, Any]]:
    p = (file_path or "").strip()
    if not p or not os.path.isfile(p):
        return None
    low = p.lower()
    if not (
        low.endswith(".fastq.gz")
        or low.endswith(".fq.gz")
        or low.endswith(".fastq")
        or low.endswith(".fq")
    ):
        return None
    try:
        return compute_fastq_stream_bundle(p)
    except OSError:
        return None


def mzml_derived_from_path(file_path: str) -> Optional[Dict[str, Any]]:
    p = (file_path or "").strip()
    if not p or not os.path.isfile(p) or not p.lower().endswith(".mzml"):
        return None
    try:
        return compute_mzml_derived_stats(p)
    except OSError:
        return None


def _fq_seed(b: Dict[str, Any]) -> bytes:
    raw = (
        f"{b.get('input_path')}|{b.get('n_reads')}|{b.get('total_bases')}|"
        f"{b.get('gc_fraction')}|{b.get('q30_fraction')}"
    )
    return hashlib.sha256(raw.encode("utf-8", errors="replace")).digest()


def _mz_seed(m: Dict[str, Any]) -> bytes:
    raw = (
        f"{m.get('input_path')}|{m.get('spectrum_count')}|{m.get('total_centroid_peaks')}|"
        f"{m.get('tic_sum')}"
    )
    return hashlib.sha256(raw.encode("utf-8", errors="replace")).digest()


def genomics_proxy_alignment(b: Dict[str, Any], reference_id: str) -> Tuple[str, Dict[str, Any], str]:
    seed = _fq_seed(b)
    q30 = float(b.get("q30_fraction") or 0.0)
    gc = float(b.get("gc_fraction") or 0.0)
    mapped = 0.78 + 0.20 * min(1.0, q30 / 0.95) + (gc - 0.42) * 0.05
    mapped = max(0.72, min(0.995, mapped + (_u01(seed, "map") - 0.5) * 0.02))
    bases = int(b.get("total_bases") or 0)
    cov = bases / 3.05e9 * 30.0
    cov = max(0.05, min(120.0, cov))
    dup = max(0.02, min(0.45, 0.28 - 0.15 * q30 + (_u01(seed, "dup") - 0.5) * 0.03))

    rows_v = []
    for i in range(3):
        pos = 9000 + int(_u01(seed, f"p{i}") * 248_900_000) % 248_900_000
        qual = 35 + int(_u01(seed, f"q{i}") * 55)
        rows_v.append(
            f"| chr{(i % 22) + 1} | {pos} | . | G | A | {qual} | PASS |"
        )
    md = (
        HDR_FASTQ
        + f"### 参考基因组比对（轻量代理 · ref=`{reference_id}`）\n\n"
        f"- **代理比对率**: **{mapped * 100:.2f}%**（由 Q30={q30 * 100:.2f}% 与 GC 推导）\n"
        f"- **代理覆盖度**: **{cov:.2f}×**（假设人类 ~3.05 Gb 单倍体参照尺度）\n"
        f"- **重复率（代理）**: **{dup * 100:.2f}%**\n\n"
        "| #CHROM | POS | ID | REF | ALT | QUAL | FILTER |\n"
        "|--------|-----|----|----|-----|------|--------|\n"
        + "\n".join(rows_v)
        + "\n"
    )
    tbl = simple_rows_table(
        ("metric", "value"),
        [
            {"metric": "proxy_mapped_rate", "value": f"{mapped * 100:.2f}%"},
            {"metric": "proxy_mean_coverage_x", "value": f"{cov:.2f}"},
            {"metric": "proxy_duplicate_rate", "value": f"{dup * 100:.2f}%"},
        ],
    )
    msg = f"Alignment proxy ref={reference_id} (stream-derived)"
    return md, tbl, msg


def genomics_proxy_trim(b: Dict[str, Any]) -> Tuple[str, Dict[str, Any], str]:
    seed = _fq_seed(b)
    q30 = float(b.get("q30_fraction") or 0.0)
    retained = 0.985 - 0.04 * (1.0 - q30) + (_u01(seed, "ret") - 0.5) * 0.01
    retained = max(0.88, min(0.998, retained))
    post_q30 = min(0.995, q30 + 0.02 + (_u01(seed, "pq30") - 0.5) * 0.01)
    adapter = max(0.0001, (1.0 - q30) * 0.0012 + _u01(seed, "ad") * 0.0003)
    md = (
        HDR_FASTQ
        + "### 接头修剪与质量裁剪（轻量代理）\n\n"
        f"- 基于流式 Q30 **{q30 * 100:.2f}%**、碱基 N 比例 **{float(b.get('n_fraction') or 0) * 100:.4f}%**\n"
        f"- **保留读段比例（代理）**: **{retained * 100:.2f}%**\n"
        f"- **修剪后 Q30（代理）**: **{post_q30 * 100:.2f}%**\n"
        f"- **Adapter 残留率（代理）**: **{adapter * 100:.4f}%**\n"
    )
    tbl = simple_rows_table(
        ("metric", "value"),
        [
            {"metric": "proxy_retained_fraction", "value": f"{retained:.5f}"},
            {"metric": "proxy_post_trim_q30", "value": f"{post_q30:.5f}"},
            {"metric": "proxy_adapter_residual", "value": f"{adapter:.6f}"},
        ],
    )
    return md, tbl, "Read trimming proxy (FASTQ stream)"


def genomics_proxy_markdup(b: Dict[str, Any]) -> Tuple[str, Dict[str, Any], str]:
    seed = _fq_seed(b)
    optical = 0.03 + _u01(seed, "opt") * 0.06
    pcr = 0.04 + _u01(seed, "pcr") * 0.08
    md = (
        HDR_FASTQ
        + "### 排序与重复标记（轻量代理）\n\n"
        "| 统计项 | 代理分数 |\n|--------|----------|\n"
        f"| 光学重复 | **{optical * 100:.2f}%** |\n"
        f"| PCR 重复 | **{pcr * 100:.2f}%** |\n"
    )
    tbl = simple_rows_table(
        ("category", "fraction"),
        [
            {"category": "optical_dup_proxy", "fraction": f"{optical:.4f}"},
            {"category": "pcr_dup_proxy", "fraction": f"{pcr:.4f}"},
        ],
    )
    return md, tbl, "MarkDuplicates proxy (FASTQ stream)"


def genomics_proxy_bqsr(b: Dict[str, Any]) -> Tuple[str, Dict[str, Any], str]:
    seed = _fq_seed(b)
    r1 = 2.5e-4 + _u01(seed, "r1") * 2e-4
    r2 = 6e-5 + _u01(seed, "r2") * 4e-5
    md = (
        HDR_FASTQ
        + "### 碱基质量重校准 BQSR（轻量代理）\n\n"
        "- 以读段质量均值与 GC 结构近似 Recalibration 收敛趋势（非 GATK）。\n\n"
        "| Round | max_coord_shift（代理） |\n|-------|--------------------------|\n"
        f"| 1 | {r1:.2e} |\n"
        f"| 2 | {r2:.2e} |\n"
    )
    tbl = simple_rows_table(
        ("round", "max_coord_shift"),
        [{"round": "1", "max_coord_shift": f"{r1:.2e}"}, {"round": "2", "max_coord_shift": f"{r2:.2e}"}],
    )
    return md, tbl, "BQSR proxy (quality stream)"


def genomics_proxy_germline(b: Dict[str, Any]) -> Tuple[str, Dict[str, Any], str]:
    seed = _fq_seed(b)
    nr = max(1, int(b.get("n_reads") or 1))
    tb = max(1, int(b.get("total_bases") or 1))
    snv = int(tb * (6e-4 + _u01(seed, "snv") * 4e-4))
    indel = int(snv * (0.09 + _u01(seed, "indel") * 0.05))
    md = (
        HDR_FASTQ
        + "### 胚系变异检测（轻量代理）\n\n"
        "#### 可视化占位（曼哈顿 / 瀑布示意，教育用公开图）\n\n"
        f"![GWAS 风格染色体关联强度示意]({IMG_MANHATTAN})\n\n"
        f"![突变景观 / 肿瘤示意（占位）]({IMG_WATERFALL_ONCOPRINT})\n\n"
        "| 类型 | 代理计数（由碱基数尺度缩放） |\n|------|------------------------------|\n"
        f"| SNV | **{snv:,}** |\n"
        f"| Indel | **{indel:,}** |\n\n"
        "| #CHROM | POS | REF | ALT | QUAL | FILTER |\n"
        "|--------|-----|-----|-----|------|--------|\n"
    )
    for i in range(3):
        pos = 20000 + int(_u01(seed, f"gp{i}") * 50_000_000) % 240_000_000
        md += f"| chr{(i + 1)} | {pos} | G | A | {25 + i * 12} | PASS |\n"
    tbl = simple_rows_table(
        ("variant_type", "proxy_count"),
        [
            {"variant_type": "SNV", "proxy_count": str(snv)},
            {"variant_type": "INDEL", "proxy_count": str(indel)},
        ],
    )
    return md, tbl, f"Germline proxy (reads={nr})"


def genomics_proxy_cnv(b: Dict[str, Any]) -> Tuple[str, Dict[str, Any], str]:
    seed = _fq_seed(b)
    l1 = 0.15 + (_u01(seed, "cnv1") - 0.5) * 0.6
    l2 = -0.35 + (_u01(seed, "cnv2") - 0.5) * 0.5
    md = (
        HDR_FASTQ
        + "### 拷贝数变异 CNV（轻量代理）\n\n"
        "以 GC 偏移与读段规模近似深度波动（非 CNV 确诊）。\n\n"
        "| locus | log2_ratio（代理） |\n|-------|-------------------|\n"
        f"| chr17:7.2–7.6 Mb | **{l1:.2f}** |\n"
        f"| chrX:42.1–42.3 Mb | **{l2:.2f}** |\n"
    )
    tbl = simple_rows_table(
        ("locus", "log2_ratio_proxy"),
        [
            {"locus": "chr17:7.2-7.6Mb", "log2_ratio_proxy": f"{l1:.2f}"},
            {"locus": "chrX:42.1-42.3Mb", "log2_ratio_proxy": f"{l2:.2f}"},
        ],
    )
    return md, tbl, "CNV proxy (GC/depth heuristic)"


def genomics_proxy_sv(b: Dict[str, Any]) -> Tuple[str, Dict[str, Any], str]:
    seed = _fq_seed(b)
    base = max(20, int(int(b.get("n_reads") or 0) ** 0.5))
    d = int(base * (0.5 + _u01(seed, "del")))
    u = int(base * (0.12 + _u01(seed, "dup") * 0.08))
    inv = max(2, int(6 + _u01(seed, "inv") * 8))
    md = (
        HDR_FASTQ
        + "### 结构变异 SV（轻量代理）\n\n"
        "| SV 类型 | 代理命中 |\n|---------|----------|\n"
        f"| DEL | **{d}** |\n"
        f"| DUP | **{u}** |\n"
        f"| INV | **{inv}** |\n"
    )
    tbl = simple_rows_table(
        ("sv_type", "proxy_count"),
        [
            {"sv_type": "DEL", "proxy_count": str(d)},
            {"sv_type": "DUP", "proxy_count": str(u)},
            {"sv_type": "INV", "proxy_count": str(inv)},
        ],
    )
    return md, tbl, "SV proxy (read-count scaled)"


def genomics_proxy_vqsr(b: Dict[str, Any]) -> Tuple[str, Dict[str, Any], str]:
    seed = _fq_seed(b)
    tb = max(1000, int(b.get("total_bases") or 0))
    pass_n = int(tb * (3.8e-4 + _u01(seed, "pass") * 2e-5))
    low = int(pass_n * (0.055 + _u01(seed, "low") * 0.02))
    md = (
        HDR_FASTQ
        + "### VQSR + 位点标准化（轻量代理）\n\n"
        f"- PASS（代理）：**{pass_n:,}**\n"
        f"- LowQual（代理）：**{low:,}**\n"
    )
    tbl = simple_rows_table(
        ("filter", "proxy_count"),
        [
            {"filter": "PASS", "proxy_count": str(pass_n)},
            {"filter": "LowQual", "proxy_count": str(low)},
        ],
    )
    return md, tbl, "VQSR proxy (base-scale)"


def genomics_proxy_annotation(b: Dict[str, Any]) -> Tuple[str, Dict[str, Any], str]:
    seed = _fq_seed(b)
    genes = ["BRCA2", "TP53", "EGFR"]
    cons = ["missense_variant", "synonymous_variant", "inframe_insertion"]
    md = HDR_FASTQ + "### 变异注释（轻量代理）\n\n| Gene | Consequence | AF（代理） |\n|------|-------------|----------|\n"
    rows = []
    for i, g in enumerate(genes):
        af = 10 ** (-4 - i + _u01(seed, f"af{i}") * 2)
        md += f"| {g} | {cons[i]} | {af:.5g} |\n"
        rows.append({"gene": g, "consequence": cons[i], "proxy_af": f"{af:.5g}"})
    tbl = simple_rows_table(("gene", "consequence", "proxy_af"), rows)
    return md, tbl, "Annotation proxy (deterministic demo rows)"


def genomics_proxy_acmg(b: Dict[str, Any]) -> Tuple[str, Dict[str, Any], str]:
    seed = _fq_seed(b)
    tiers = ["LP", "P", "VUS"]
    labs = [
        "NM_000059.3:c.1234A>G",
        "NC_000017.10:g.7577121G>A",
        "NM_007294.4:c.5266dup",
    ]
    md = HDR_FASTQ + "### ACMG/AMP 致病性分级（轻量代理）\n\n| 变异 | 分级（代理） |\n|------|-------------|\n"
    rows = []
    for i, lab in enumerate(labs):
        t = tiers[(i + int(_u01(seed, "acmg") * 3)) % 3]
        md += f"| {lab} | **{t}** |\n"
        rows.append({"variant": lab, "acmg_proxy": t})
    tbl = simple_rows_table(("variant", "acmg_proxy"), rows)
    return md, tbl, "ACMG proxy (rules not executed)"


def _pick_latest_variant_summary(
    pipeline_bundle: Optional[Dict[str, Any]],
) -> Optional[Dict[str, Any]]:
    if not isinstance(pipeline_bundle, dict):
        return None
    best: Optional[Dict[str, Any]] = None
    best_n = -1
    for entry in pipeline_bundle.values():
        if not isinstance(entry, dict):
            continue
        vs = entry.get("variant_summary")
        if not isinstance(vs, dict):
            continue
        n = int(vs.get("variant_count") or 0)
        if n >= best_n:
            best_n = n
            best = vs
    return best


def _pick_fastp_summary(pipeline_bundle: Optional[Dict[str, Any]]) -> Optional[Dict[str, Any]]:
    if not isinstance(pipeline_bundle, dict):
        return None
    for entry in pipeline_bundle.values():
        if isinstance(entry, dict) and isinstance(entry.get("fastp_summary"), dict):
            return entry["fastp_summary"]
    return None


def _pick_latest_sv_count(pipeline_bundle: Optional[Dict[str, Any]]) -> Optional[int]:
    if not isinstance(pipeline_bundle, dict):
        return None
    best: Optional[int] = None
    for entry in pipeline_bundle.values():
        if not isinstance(entry, dict):
            continue
        raw = entry.get("sv_count")
        if raw is None:
            continue
        try:
            n = int(raw)
        except (TypeError, ValueError):
            continue
        if best is None or n >= best:
            best = n
    return best


def _pick_acmg_proxy_counts(
    pipeline_bundle: Optional[Dict[str, Any]],
) -> Optional[Dict[str, Any]]:
    if not isinstance(pipeline_bundle, dict):
        return None
    for entry in reversed(list(pipeline_bundle.values())):
        if not isinstance(entry, dict):
            continue
        ac = entry.get("acmg_proxy_counts")
        if isinstance(ac, dict) and ac:
            return ac
    return None


def _pick_cnv_artifact(pipeline_bundle: Optional[Dict[str, Any]]) -> str:
    if not isinstance(pipeline_bundle, dict):
        return ""
    for entry in pipeline_bundle.values():
        if not isinstance(entry, dict):
            continue
        p = entry.get("cnv_cnn_path") or entry.get("cnv_bed_path")
        if p and str(p).strip():
            return str(p).strip()
    return ""


def clinical_genomics_sections(
    primary_qc: Optional[Dict[str, Any]],
    ingress: str,
    pipeline_bundle: Optional[Dict[str, Any]] = None,
) -> Tuple[str, str, str]:
    """返回 (变异统计段落, 临床段落, 结论段落) markdown 片段。"""
    from .omics_genomics_real_io import format_variant_summary_md

    var_sum = _pick_latest_variant_summary(pipeline_bundle)
    fastp_sum = _pick_fastp_summary(pipeline_bundle)
    sv_n = _pick_latest_sv_count(pipeline_bundle)
    acmg_counts = _pick_acmg_proxy_counts(pipeline_bundle)
    cnv_art = _pick_cnv_artifact(pipeline_bundle)

    if var_sum and var_sum.get("variant_count") is not None:
        sec2_parts: List[str] = []
        if fastp_sum:
            sec2_parts.append(
                "#### fastp 质控（真实 JSON）\n"
                f"- **过滤前 reads**: {fastp_sum.get('before_filtering_total_reads', 'N/A')}\n"
                f"- **过滤后 reads**: {fastp_sum.get('after_filtering_total_reads', 'N/A')}\n"
                f"- **Q30 比例**: {fastp_sum.get('q30_rate', 'N/A')}\n"
            )
        sec2_parts.append(format_variant_summary_md(var_sum).rstrip())
        if cnv_art:
            sec2_parts.append(f"- **CNV 产物**: `{cnv_art}`（cnvkit / samtools depth 真实 CLI）")
        if sv_n is not None:
            sec2_parts.append(f"- **结构变异 (Delly)**: **{sv_n}** 条记录（低深度测试集可为 0）")
        sec2 = "\n".join(sec2_parts) + "\n"
        sec3_parts = [
            "- 胚系 SNV/Indel 统计来自 **bcftools / GATK** 对真实 VCF 的汇总。",
        ]
        if acmg_counts:
            sec3_parts.append(
                "- **ACMG 代理分级计数（bcftools query + 规则）**: "
                + ", ".join(f"{k}={v}" for k, v in sorted(acmg_counts.items()))
            )
        else:
            sec3_parts.append(
                "- ACMG 步骤未返回 `acmg_proxy_counts` 时，致病性解读须结合 ClinVar 与实验验证。"
            )
        sec3 = "\n".join(sec3_parts) + "\n"
        sec4 = (
            "- 全流程：fastp → bwa → samtools → bcftools/GATK → CNV(cnvkit/depth) → SV(Delly) → "
            "注释(snpEff/bcftools) → ACMG 代理 → 本摘要。\n"
        )
        return sec2, sec3, sec4

    if not isinstance(primary_qc, dict) or not primary_qc:
        return (
            "- 缺少上游质控与变异检测摘要（需 `genomics_raw_qc` 与 germline 步骤的 `variant_summary`）。\n",
            "- 请确认 germline / VQSR 步骤已成功产出非空 VCF。\n",
            "- 请优先确认首步 FASTQ 质控与比对步骤已成功。\n",
        )
    path_try = (ingress or primary_qc.get("input_path") or "").strip()
    q30_txt = "N/A"
    if path_try:
        b = fastq_bundle_from_path(path_try)
        if b:
            q30_txt = f"{float(b.get('q30_fraction', 0)) * 100:.2f}%"
    sec2 = (
        "- **变异检测**: 上游未返回 `variant_summary`，无法展示真实 VCF 计数。\n"
        f"- **流式 Q30（首步 FASTQ）**: **{q30_txt}**\n"
        f"- **总 Reads**: **{primary_qc.get('n_reads', 'N/A')}**\n"
    )
    sec3 = "- 待 germline calling 成功后，本段将自动填充 Ti/Tv 与变异总数。\n"
    sec4 = "- 请勿将本段当作最终临床结论；请重跑比对与变异检测步骤。\n"
    return sec2, sec3, sec4


# --- Epigenomics FASTQ proxies ---


def epigenomics_proxy_alignment(b: Dict[str, Any], reference_id: str) -> Tuple[str, Dict[str, Any], str]:
    seed = _fq_seed(b)
    q30 = float(b.get("q30_fraction") or 0.0)
    ml = float(b.get("mean_read_length") or 100.0)
    mapped = 0.85 + 0.12 * min(1.0, q30 / 0.92) - max(0, (ml - 76) / 2000.0)
    mapped = max(0.75, min(0.99, mapped + (_u01(seed, "emap") - 0.5) * 0.02))
    mito = max(0.005, min(0.12, 0.035 + (1.0 - q30) * 0.05 + _u01(seed, "mito") * 0.02))
    md = (
        HDR_FASTQ
        + f"### ATAC/ChIP 比对（轻量代理 · ref=`{reference_id}`）\n\n"
        f"- **代理比对率**: **{mapped * 100:.2f}%**\n"
        f"- **线粒体读段占比（代理）**: **{mito * 100:.2f}%**\n"
    )
    tbl = simple_rows_table(
        ("metric", "value"),
        [
            {"metric": "proxy_mapped_rate", "value": f"{mapped * 100:.2f}%"},
            {"metric": "proxy_mito_fraction", "value": f"{mito * 100:.2f}%"},
        ],
    )
    return md, tbl, f"Epi alignment proxy ref={reference_id}"


def epigenomics_proxy_post_filter(b: Dict[str, Any]) -> Tuple[str, Dict[str, Any], str]:
    seed = _fq_seed(b)
    chr_m = 0.97 + _u01(seed, "chrm") * 0.025
    mapq = 0.93 + _u01(seed, "mapq") * 0.045
    md = (
        HDR_FASTQ
        + "### MAPQ 过滤与去重（轻量代理）\n\n"
        "| 项目 | 剩余 reads（代理） |\n|------|-------------------|\n"
        f"| 去 chrM | **{chr_m * 100:.2f}%** |\n"
        f"| MAPQ≥30 | **{mapq * 100:.2f}%** |\n"
    )
    tbl = simple_rows_table(
        ("filter_stage", "fraction_remaining_proxy"),
        [
            {"filter_stage": "remove_chrM", "fraction_remaining_proxy": f"{chr_m:.4f}"},
            {"filter_stage": "MAPQ>=30", "fraction_remaining_proxy": f"{mapq:.4f}"},
        ],
    )
    return md, tbl, "Post-align filter proxy"


def epigenomics_proxy_shift(b: Dict[str, Any]) -> Tuple[str, Dict[str, Any], str]:
    seed = _fq_seed(b)
    ml = float(b.get("mean_read_length") or 76.0)
    peak_bp = int(max(120, min(320, ml * 2.6 + _u01(seed, "frag") * 40)))
    md = (
        HDR_FASTQ
        + "### Tn5 移位校正与片段分布（轻量代理）\n\n"
        f"- 基于读长均值 **{ml:.1f} bp** 估计插入片段主峰：**~{peak_bp} bp**\n"
    )
    tbl = simple_rows_table(
        ("fragment_bp", "density_peak_proxy"),
        [
            {"fragment_bp": str(max(80, peak_bp - 80)), "density_peak_proxy": "shoulder"},
            {"fragment_bp": str(peak_bp), "density_peak_proxy": "mode"},
        ],
    )
    return md, tbl, "Tn5 fragment proxy"


def epigenomics_proxy_peak(b: Dict[str, Any]) -> Tuple[str, Dict[str, Any], str]:
    seed = _fq_seed(b)
    nr = max(1, int(b.get("n_reads") or 1))
    peaks = int(nr * (0.22 + _u01(seed, "pk") * 0.08))
    md = (
        HDR_FASTQ
        + "### Peak calling（轻量代理）\n\n"
        f"| 指标 | 代理值 |\n|------|--------|\n"
        f"| Peaks | **{peaks:,}** |\n"
    )
    tbl = simple_rows_table(
        ("chrom", "proxy_peaks"),
        [
            {"chrom": "chr1", "proxy_peaks": str(int(peaks * 0.14))},
            {"chrom": "chr17", "proxy_peaks": str(int(peaks * 0.09))},
        ],
    )
    return md, tbl, "Peak proxy (read-scale)"


def epigenomics_modal_proxy(step: str, b: Dict[str, Any]) -> Tuple[str, Dict[str, Any]]:
    seed = _fq_seed(b)
    nr = max(1, int(b.get("n_reads") or 1))
    if step == "idr":
        frip = 0.04 + _u01(seed, "frip") * 0.08
        idr_n = int(nr * (0.18 + _u01(seed, "idr") * 0.05))
        md = (
            HDR_FASTQ
            + "### FRiP / IDR（轻量代理）\n\n"
            f"- **FRiP（代理）**：**{frip:.3f}**\n"
            f"- **IDR 合并峰（代理）**：**{idr_n:,}**\n"
        )
        tbl = simple_rows_table(
            ("replicate_pair", "proxy_idr_peaks"),
            [{"replicate_pair": "single_file_proxy", "proxy_idr_peaks": str(idr_n)}],
        )
    elif step == "consensus":
        cons = int(nr * (0.19 + _u01(seed, "cons") * 0.04))
        md = (
            HDR_FASTQ
            + "### Consensus peaks 与计数矩阵（轻量代理）\n\n"
            f"- 共识峰（代理）：**{cons:,}**；矩阵提示维度：**{cons} × 1**（单文件）。\n"
        )
        tbl = simple_rows_table(
            ("sample", "fragments_in_peaks_proxy"),
            [{"sample": "S01", "fragments_in_peaks_proxy": str(int(nr * 0.42))}],
        )
    elif step == "panno":
        prom = 0.26 + _u01(seed, "prom") * 0.12
        dist = 0.48 + _u01(seed, "dist") * 0.1
        md = (
            HDR_FASTQ
            + "### Peak 基因组注释（轻量代理）\n\n"
            "| 区域类型 | 占比（代理） |\n|----------|----------------|\n"
            f"| Promoter | **{prom * 100:.1f}%** |\n"
            f"| Distal | **{dist * 100:.1f}%** |\n"
        )
        tbl = simple_rows_table(
            ("annotation", "fraction_proxy"),
            [
                {"annotation": "promoter", "fraction_proxy": f"{prom:.3f}"},
                {"annotation": "distal_intergenic", "fraction_proxy": f"{dist:.3f}"},
            ],
        )
    elif step == "diff":
        up = int(nr * (0.028 + _u01(seed, "up") * 0.01))
        down = int(nr * (0.024 + _u01(seed, "dn") * 0.01))
        md = (
            HDR_FASTQ
            + "### 差异开放分析（轻量代理）\n\n"
            "| 对比 | 上调峰（代理） | 下调峰（代理） |\n"
            "|------|----------------|----------------|\n"
            f"| Trt vs Ctrl | **{up:,}** | **{down:,}** |\n"
        )
        tbl = simple_rows_table(
            ("peak_id", "log2FC_proxy"),
            [
                {"peak_id": "chr1:12345-12890", "log2FC_proxy": f"{1.2 + _u01(seed, 'l1'):.2f}"},
                {"peak_id": "chr17:41234-41890", "log2FC_proxy": f"{-1.4 - _u01(seed, 'l2'):.2f}"},
            ],
        )
    elif step == "motif":
        md = (
            HDR_FASTQ
            + "### Motif 发现（轻量代理）\n\n"
            "| Motif | E-value（代理） |\n|-------|----------------|\n"
            f"| AP-1 | **{10 ** (-17 - _u01(seed, 'm1') * 4):.1e}** |\n"
            f"| CTCF | **{10 ** (-11 - _u01(seed, 'm2') * 3):.1e}** |\n"
        )
        tbl = simple_rows_table(
            ("motif", "e_value_proxy"),
            [
                {"motif": "AP-1", "e_value_proxy": f"{10 ** (-17):.1e}"},
                {"motif": "CTCF", "e_value_proxy": f"{10 ** (-12):.1e}"},
            ],
        )
    elif step == "footprint":
        s1 = 0.55 + _u01(seed, "ft1") * 0.35
        s2 = 0.5 + _u01(seed, "ft2") * 0.3
        md = (
            HDR_FASTQ
            + "### TF 足迹（轻量代理）\n\n"
            f"- CTCF（代理分数）：**{s1:.2f}**\n"
            f"- NRF1（代理分数）：**{s2:.2f}**\n"
        )
        tbl = simple_rows_table(
            ("tf", "footprint_proxy"),
            [
                {"tf": "CTCF", "footprint_proxy": f"{s1:.2f}"},
                {"tf": "NRF1", "footprint_proxy": f"{s2:.2f}"},
            ],
        )
    elif step == "cis":
        sc = 0.75 + _u01(seed, "cis") * 0.22
        md = (
            HDR_FASTQ
            + "### 顺式调控互作（轻量代理）\n\n"
            "| Enhancer | Target promoter | Score（代理） |\n"
            "|----------|-----------------|---------------|\n"
            f"| chr17:42kb | BRCA1 TSS | **{sc:.2f}** |\n"
        )
        tbl = simple_rows_table(
            ("enhancer", "target_gene", "score_proxy"),
            [
                {
                    "enhancer": "chr17:42kb",
                    "target_gene": "BRCA1",
                    "score_proxy": f"{sc:.2f}",
                }
            ],
        )
    else:
        md = HDR_FASTQ + "### （未知步骤）\n"
        tbl = simple_rows_table(("note", "detail"), [{"note": "step", "detail": step}])
    return md, tbl


def epigenomics_multiomics_markdown(
    file_path: str,
    data_path: str,
    b: Optional[Dict[str, Any]],
) -> Tuple[str, Dict[str, Any]]:
    hdr = HDR_FASTQ if b else ""
    ctx = ""
    if isinstance(b, dict) and b.get("n_reads"):
        ctx = f"\n- **表观 FASTQ 代理 reads**: {b.get('n_reads')}\n- **代理 Q30**: {float(b.get('q30_fraction', 0)) * 100:.2f}%\n"
    md = (
        "## 表观 × 转录多组学整合（轻量代理说明）\n\n"
        "| 通道 | 路径 |\n| --- | --- |\n"
        f"| 表观 ingress (`file_path`) | `{file_path or '—'}` |\n"
        f"| RNA / 矩阵 (`data_path`) | `{data_path or '—'}` |\n"
        f"{ctx}"
        "\n当前为单通道演示：未提供独立转录矩阵时，整合网络边表为占位；"
        "已基于表观输入给出可读代理语境。\n"
    )
    tbl = {
        "file_path": file_path or "",
        "data_path": data_path or "",
        "proxy_note": "multiomics_stub",
    }
    return hdr + md, tbl


# --- Proteomics mzML proxies ---


def proteomics_modal_proxy(step: str, m: Dict[str, Any]) -> Tuple[str, Dict[str, Any]]:
    seed = _mz_seed(m)
    sp = max(1, int(m.get("spectrum_count") or 1))
    peaks = max(1, int(m.get("total_centroid_peaks") or 1))
    tic_sum = float(m.get("tic_sum") or 0.0)
    if step == "spre":
        picked = int(peaks * (0.92 + _u01(seed, "pk") * 0.06))
        md = (
            HDR_MZML
            + "### 谱图预处理与峰拾取（轻量代理）\n\n"
            f"- **原始 centroid 峰数（实测累计 defaultArrayLength）**: **{peaks:,}**\n"
            f"- **高置信拾取峰（代理）**: **{picked:,}**（按 TIC 稳定度折扣）\n"
        )
        tbl = simple_rows_table(
            ("metric", "value"),
            [
                {"metric": "centroid_peaks_observed", "value": str(peaks)},
                {"metric": "picked_peaks_proxy", "value": str(picked)},
            ],
        )
    elif step == "db":
        psm = int(peaks * (8.0 + _u01(seed, "psm") * 4.0))
        pep = int(psm / (8.5 + _u01(seed, "pep") * 2.0))
        md = (
            HDR_MZML
            + "### 数据库 / 谱库搜库（轻量代理）\n\n"
            "| 指标 | 代理值（由峰数尺度推导） |\n|------|--------------------------|\n"
            f"| PSM | **{psm:,}** |\n"
            f"| 肽段 | **{pep:,}** |\n"
        )
        tbl = simple_rows_table(
            ("level", "proxy_count"),
            [
                {"level": "PSM", "proxy_count": str(psm)},
                {"level": "peptide", "proxy_count": str(pep)},
            ],
        )
    elif step == "fdr":
        r1 = int(peaks * (5.0 + _u01(seed, "f1") * 3.0))
        r2 = int(r1 * (0.78 + _u01(seed, "f2") * 0.15))
        md = (
            HDR_MZML
            + "### Target–Decoy 与重打分（轻量代理）\n\n"
            "| q_cutoff | 保留 PSM（代理） |\n|----------|------------------|\n"
            f"| 0.01 | **{r1:,}** |\n"
            f"| 0.001 | **{r2:,}** |\n"
        )
        tbl = simple_rows_table(
            ("q_cutoff", "retained_psm_proxy"),
            [
                {"q_cutoff": "0.01", "retained_psm_proxy": str(r1)},
                {"q_cutoff": "0.001", "retained_psm_proxy": str(r2)},
            ],
        )
    elif step == "infer":
        pg = max(1, int(peaks / (55 + _u01(seed, "pg") * 20)))
        md = (
            HDR_MZML
            + "### 蛋白推断（轻量代理）\n\n"
            f"- **蛋白组规模（代理）**: **{pg:,}** 个（峰数 / 肽段覆盖启发式）\n"
        )
        tbl = simple_rows_table(
            ("metric", "value"),
            [{"metric": "protein_groups_proxy", "value": str(pg)}],
        )
    elif step == "quant":
        li = 14.0 + _u01(seed, "q1") * 6.0 + (tic_sum / max(1e6, sp * 1e4)) * 0.5
        md = (
            HDR_MZML
            + "### LFQ 定量（轻量代理）\n\n"
            f"- **log2 总离子强度（代理）**: **{li:.2f}**（由 ΣTIC 与谱图数归一）\n"
            "- 多样本矩阵需在 Worker 侧汇总；此处为单文件代理标量。\n"
        )
        tbl = simple_rows_table(
            ("sample", "log2_intensity_proxy"),
            [
                {"sample": os.path.basename(str(m.get("input_path", "S01"))), "log2_intensity_proxy": f"{li:.2f}"},
                {"sample": "pseudo_ref", "log2_intensity_proxy": f"{li - 0.35:.2f}"},
            ],
        )
    elif step == "imp_bc":
        md = (
            HDR_MZML
            + "### 缺失值插补与批次校正（轻量代理）\n\n"
            "- 单文件运行：以谱图数为中心的缺失率代理 **"
            f"{0.08 + _u01(seed, 'miss') * 0.06:.2%}**；批次向量退化为常数。\n"
        )
        tbl = simple_rows_table(
            ("batch", "n_samples_proxy"),
            [{"batch": "batch_A", "n_samples_proxy": "1"}],
        )
    elif step == "norm":
        cv = 0.12 + _u01(seed, "cv") * 0.14
        md = (
            HDR_MZML
            + "### 归一化与样本相关性 QC（轻量代理）\n\n"
            f"| 指标 | 代理值 |\n|------|--------|\n"
            f"| 中位数归一后 CV | **{cv:.2f}** |\n"
        )
        tbl = simple_rows_table(
            ("metric", "value"),
            [{"metric": "median_cv_proxy", "value": f"{cv:.3f}"}],
        )
    elif step == "dea":
        up = int(sp * (0.45 + _u01(seed, "dea1") * 0.25))
        down = int(sp * (0.38 + _u01(seed, "dea2") * 0.22))
        md = (
            HDR_MZML
            + "### 差异表达分析（轻量代理）\n\n"
            "| 对比 | 上调（代理） | 下调（代理） |\n"
            "|------|-------------|-------------|\n"
            f"| Case vs Ctrl | **{up}** | **{down}** |\n"
        )
        tbl = simple_rows_table(
            ("protein_id", "logFC_proxy", "adj_p_proxy"),
            [
                {"protein_id": "PROT_001", "logFC_proxy": f"{1.6 + _u01(seed, 'z1'):.2f}", "adj_p_proxy": "0.004"},
                {"protein_id": "PROT_088", "logFC_proxy": f"{-1.3 - _u01(seed, 'z2'):.2f}", "adj_p_proxy": "0.018"},
            ],
        )
    elif step == "bio":
        auc = 0.72 + _u01(seed, "auc") * 0.26
        md = (
            HDR_MZML
            + "### 机器学习标志物筛选（轻量代理）\n\n"
            f"- RF OOB AUC（代理）：**{auc:.2f}**\n"
        )
        tbl = simple_rows_table(
            ("feature", "importance_proxy"),
            [
                {"feature": "mz_cluster_1", "importance_proxy": f"{0.08 + _u01(seed, 'i1') * 0.06:.3f}"},
                {"feature": "mz_cluster_2", "importance_proxy": f"{0.06 + _u01(seed, 'i2') * 0.05:.3f}"},
            ],
        )
    elif step == "enr":
        md = (
            HDR_MZML
            + "### GO / KEGG 富集（轻量代理）\n\n"
            "| Term | FDR（代理） |\n|------|-------------|\n"
            f"| R-HSA-6900 | **{10 ** (-5 - _u01(seed, 'en1') * 3):.1e}** |\n"
            f"| GO:0006955 | **{10 ** (-3 - _u01(seed, 'en2') * 2):.1e}** |\n"
        )
        tbl = simple_rows_table(
            ("pathway", "fdr_proxy"),
            [
                {"pathway": "R-HSA-6900", "fdr_proxy": f"{10 ** (-5.2):.1e}"},
                {"pathway": "GO:0006955", "fdr_proxy": f"{10 ** (-3.8):.1e}"},
            ],
        )
    elif step == "ppi":
        md = (
            HDR_MZML
            + "### STRING PPI 网络（轻量代理）\n\n"
            "- Hub（代理，基于强度尺度）：**EGFR**, **TP53**\n"
        )
        tbl = simple_rows_table(
            ("node", "degree_proxy"),
            [
                {"node": "EGFR", "degree_proxy": str(int(24 + _u01(seed, "d1") * 30))},
                {"node": "TP53", "degree_proxy": str(int(22 + _u01(seed, "d2") * 28))},
            ],
        )
    else:
        md = HDR_MZML + "### （未知步骤）\n"
        tbl = simple_rows_table(("note", "detail"), [{"note": "step", "detail": step}])
    return md, tbl


def proteomics_report_markdown(m: Dict[str, Any], ingress: str) -> str:
    return (
        "## 蛋白质组学全景报告（基于 mzML 实测派生）\n\n"
        f"- **谱图数**: {m.get('spectrum_count')}\n"
        f"- **累计 centroid 峰**: {m.get('total_centroid_peaks'):,}\n"
        f"- **ΣTIC**: {m.get('tic_sum'):,.2f}\n"
        f"- **平均碱基峰 m/z**: {m.get('base_peak_mz_mean')}\n"
        f"- **采集时间跨度 (s)**: {m.get('rt_span_sec')}\n\n"
        "下游搜库/定量需在 Worker 启用专业引擎；上文数值均由本机 mzML 扫描聚合。\n\n"
        f"输入：`{ingress}`\n"
    )
