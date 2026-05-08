"""
组学 Mock 步骤的「可视化契约」：与右栏 / 工作台渲染对齐。

成功返回除 status/message 外，建议始终包含：
- markdown：人类可读摘要（可含 Markdown 表格）
- image_urls：公开示例图 URL + 可选 matplotlib data URI（教育/占位）
- table_data：结构化预览（columns/rows 或领域自定义键）

重型算子仍须在 Worker/TaaS；此处仅为 Mock 层的 UI 占位，满足前端富文本展示。
"""
from __future__ import annotations

from typing import Any, Dict, List, Optional, Sequence

# Wikimedia Commons — 教育/示意图占位（HTTP 200 稳定优先）
IMG_DNA_HELIX = (
    "https://upload.wikimedia.org/wikipedia/commons/thumb/e/e6/"
    "DNA_double_helix_horizontal.png/640px-DNA_double_helix_horizontal.png"
)
IMG_MS_SCHEMATIC = (
    "https://upload.wikimedia.org/wikipedia/commons/thumb/4/4f/"
    "Mass_spectrometry_schematic.svg/500px-Mass_spectrometry_schematic.svg.png"
)
IMG_CHIP_WORKFLOW = (
    "https://upload.wikimedia.org/wikipedia/commons/thumb/8/8c/"
    "ChIP-seq_workflow.svg/500px-ChIP-seq_workflow.svg.png"
)
IMG_VOLCANO = (
    "https://upload.wikimedia.org/wikipedia/commons/thumb/2/2f/"
    "Volcano_plot.svg/480px-Volcano_plot.svg.png"
)
IMG_PCA = (
    "https://upload.wikimedia.org/wikipedia/commons/thumb/0/0d/"
    "Principal_component_analysis.svg/480px-Principal_component_analysis.svg.png"
)
# Phred / 测序质量示意（质控步骤）
IMG_PHRED_QUALITY = (
    "https://upload.wikimedia.org/wikipedia/commons/thumb/6/64/"
    "Quality_scores_plot.png/640px-Quality_scores_plot.png"
)
# IGV / 基因组浏览器类教育示意图（比对与覆盖）
IMG_GENOME_BROWSER = (
    "https://upload.wikimedia.org/wikipedia/commons/thumb/3/3a/"
    "Chromosome_Structure.svg/640px-Chromosome_Structure.svg.png"
)
# GWAS / 曼哈顿图（变异关联示意）
IMG_MANHATTAN = (
    "https://upload.wikimedia.org/wikipedia/commons/thumb/7/72/"
    "Manhattan_plot_from_a_GWAS_of_kidney_stone_disease_%28white_background%29.jpg/"
    "800px-Manhattan_plot_from_a_GWAS_of_kidney_stone_disease_%28white_background%29.jpg"
)
# 瀑布图 / oncoprint 风格（体细胞突变瀑布示意）
IMG_WATERFALL_ONCOPRINT = (
    "https://upload.wikimedia.org/wikipedia/commons/thumb/6/65/"
    "Tumor_schematic.svg/640px-Tumor_schematic.svg.png"
)
# 差异表达：聚类热图（公开 PCA 结构图作样本聚类教育占位；矩阵热图由 synthetic heatmap 兜底）
IMG_HEATMAP_CLUSTERED = IMG_PCA
# 通路富集：谱图/强度示意（详细气泡形态由 synthetic bubble 生成）
IMG_PATHWAY_BUBBLE = IMG_MS_SCHEMATIC
# Peak 基因组分布：经典饼图示意（Wikimedia 稳定缩略图）
IMG_PEAK_GENOME_PIE = (
    "https://upload.wikimedia.org/wikipedia/commons/thumb/7/72/"
    "Pie_chart.svg/500px-Pie_chart.svg.png"
)
# Motif：序列 Logo（经典 DNA logo 教育图）
IMG_SEQUENCE_LOGO = (
    "https://upload.wikimedia.org/wikipedia/commons/thumb/e/e4/"
    "DNA_sequence_logo.png/640px-DNA_sequence_logo.png"
)
# CNV：染色体尺度变异示意（与 synthetic cnv_heatmap 搭配）
IMG_CNV_PROFILE_HEATMAP = IMG_MANHATTAN

PLACEHOLDER_IMAGES_GENOMICS = {
    "qc": IMG_DNA_HELIX,
    "align": IMG_GENOME_BROWSER,
    "variant": IMG_VOLCANO,
    "report": IMG_DNA_HELIX,
}
PLACEHOLDER_IMAGES_PROTEOMICS = {
    "qc": IMG_MS_SCHEMATIC,
    "quant": IMG_PCA,
    "dea": IMG_VOLCANO,
}
PLACEHOLDER_IMAGES_EPIGENOMICS = {
    "qc": IMG_DNA_HELIX,
    "peak": IMG_CHIP_WORKFLOW,
    "diff": IMG_VOLCANO,
}


def _preset_for_tool(tool_id: str) -> Dict[str, Any]:
    """按 tool_id 返回 {urls: [], synthetic: kind}；未知工具使用通用预设。"""
    m = OMICS_TOOL_VISUAL_PRESETS.get(tool_id)
    if m:
        return m
    return {"urls": [IMG_DNA_HELIX, IMG_VOLCANO], "synthetic": "generic"}


def enrich_omics_visual_urls(tool_id: Optional[str], existing: Optional[Sequence[str]] = None) -> List[str]:
    """
    合并：调用方已有 URL + 工具预设公开图 + matplotlib 合成图（兜底绝不返回空列表）。
    """
    from .omics_visual_synthetic import synthetic_png_data_uri

    tid = (tool_id or "").strip()
    raw = [x for x in (existing or ()) if isinstance(x, str) and x.strip()]
    preset = _preset_for_tool(tid)
    urls = list(preset.get("urls") or [])
    syn_kind = str(preset.get("synthetic") or "generic")

    out: List[str] = []
    for u in raw + urls:
        if u not in out:
            out.append(u)
    seed = abs(hash(tid)) % (2**31) if tid else 42
    syn = synthetic_png_data_uri(syn_kind, seed=seed)
    if syn and syn not in out:
        out.append(syn)
    if len(out) < 2:
        out.extend([IMG_DNA_HELIX, IMG_VOLCANO])
    # 去重保序，限制长度避免 SSE 过大
    dedup: List[str] = []
    for u in out:
        if u not in dedup:
            dedup.append(u)
    return dedup[:10]


# tool_id -> 公开示意图（2 张）+ 合成图类型（与 omics_visual_synthetic.plot_kind 对齐）
OMICS_TOOL_VISUAL_PRESETS: Dict[str, Dict[str, Any]] = {
    # --- Genomics ---
    "genomics_raw_qc": {"urls": [IMG_PHRED_QUALITY, IMG_DNA_HELIX], "synthetic": "coverage"},
    "genomics_read_trimming": {"urls": [IMG_PHRED_QUALITY, IMG_GENOME_BROWSER], "synthetic": "coverage"},
    "genomics_alignment": {"urls": [IMG_GENOME_BROWSER, IMG_DNA_HELIX], "synthetic": "coverage"},
    "genomics_mark_duplicates": {"urls": [IMG_GENOME_BROWSER, IMG_VOLCANO], "synthetic": "generic"},
    "genomics_bqsr": {"urls": [IMG_PHRED_QUALITY, IMG_GENOME_BROWSER], "synthetic": "coverage"},
    "genomics_germline_calling": {"urls": [IMG_MANHATTAN, IMG_WATERFALL_ONCOPRINT], "synthetic": "manhattan"},
    "genomics_cnv_calling": {
        "urls": [IMG_CNV_PROFILE_HEATMAP, IMG_GENOME_BROWSER],
        "synthetic": "cnv_heatmap",
    },
    "genomics_sv_calling": {"urls": [IMG_MANHATTAN, IMG_DNA_HELIX], "synthetic": "manhattan"},
    "genomics_vqsr_filtering": {"urls": [IMG_VOLCANO, IMG_MANHATTAN], "synthetic": "generic"},
    "genomics_variant_annotation": {"urls": [IMG_VOLCANO, IMG_MANHATTAN], "synthetic": "waterfall"},
    "genomics_acmg_classification": {"urls": [IMG_WATERFALL_ONCOPRINT, IMG_VOLCANO], "synthetic": "waterfall"},
    "genomics_clinical_reporting": {"urls": [IMG_DNA_HELIX, IMG_MANHATTAN], "synthetic": "generic"},
    # --- Proteomics ---
    "proteomics_raw_qc_conversion": {"urls": [IMG_MS_SCHEMATIC, IMG_PCA], "synthetic": "spectra"},
    "proteomics_spectrum_preprocessing": {"urls": [IMG_MS_SCHEMATIC, IMG_VOLCANO], "synthetic": "spectra"},
    "proteomics_database_search": {"urls": [IMG_MS_SCHEMATIC, IMG_VOLCANO], "synthetic": "spectra"},
    "proteomics_fdr_rescoring": {"urls": [IMG_VOLCANO, IMG_MS_SCHEMATIC], "synthetic": "generic"},
    "proteomics_protein_inference": {"urls": [IMG_PCA, IMG_VOLCANO], "synthetic": "waterfall"},
    "proteomics_quantification": {"urls": [IMG_PCA, IMG_MS_SCHEMATIC], "synthetic": "spectra"},
    "proteomics_imputation_batch_correction": {"urls": [IMG_PCA, IMG_VOLCANO], "synthetic": "generic"},
    "proteomics_normalization_qc": {"urls": [IMG_PCA, IMG_VOLCANO], "synthetic": "generic"},
    "proteomics_differential_analysis": {
        "urls": [IMG_VOLCANO, IMG_HEATMAP_CLUSTERED],
        "synthetic": "heatmap",
    },
    "proteomics_biomarker_discovery": {"urls": [IMG_VOLCANO, IMG_PCA], "synthetic": "generic"},
    "proteomics_functional_enrichment": {
        "urls": [IMG_PATHWAY_BUBBLE, IMG_VOLCANO],
        "synthetic": "bubble",
    },
    "proteomics_ppi_network_analysis": {"urls": [IMG_VOLCANO, IMG_PCA], "synthetic": "generic"},
    "proteomics_clinical_reporting": {"urls": [IMG_MS_SCHEMATIC, IMG_VOLCANO], "synthetic": "spectra"},
    # --- Epigenomics ---
    "epigenomics_raw_qc_trimming": {"urls": [IMG_PHRED_QUALITY, IMG_CHIP_WORKFLOW], "synthetic": "coverage"},
    "epigenomics_alignment": {"urls": [IMG_GENOME_BROWSER, IMG_CHIP_WORKFLOW], "synthetic": "coverage"},
    "epigenomics_post_align_filtering": {"urls": [IMG_GENOME_BROWSER, IMG_VOLCANO], "synthetic": "coverage"},
    "epigenomics_shift_fragment_analysis": {"urls": [IMG_CHIP_WORKFLOW, IMG_GENOME_BROWSER], "synthetic": "peak"},
    "epigenomics_peak_calling": {
        "urls": [IMG_PEAK_GENOME_PIE, IMG_CHIP_WORKFLOW],
        "synthetic": "pie",
    },
    "epigenomics_reproducibility_idr": {"urls": [IMG_CHIP_WORKFLOW, IMG_VOLCANO], "synthetic": "peak"},
    "epigenomics_consensus_peak_counting": {"urls": [IMG_VOLCANO, IMG_CHIP_WORKFLOW], "synthetic": "generic"},
    "epigenomics_peak_annotation": {"urls": [IMG_CHIP_WORKFLOW, IMG_GENOME_BROWSER], "synthetic": "generic"},
    "epigenomics_diff_accessibility": {"urls": [IMG_VOLCANO, IMG_PCA], "synthetic": "waterfall"},
    "epigenomics_motif_discovery": {
        "urls": [IMG_SEQUENCE_LOGO, IMG_DNA_HELIX],
        "synthetic": "logo",
    },
    "epigenomics_tf_footprinting": {"urls": [IMG_CHIP_WORKFLOW, IMG_VOLCANO], "synthetic": "peak"},
    "epigenomics_cis_regulatory_interactions": {"urls": [IMG_CHIP_WORKFLOW, IMG_GENOME_BROWSER], "synthetic": "generic"},
    "epigenomics_multiomics_integration": {"urls": [IMG_VOLCANO, IMG_PCA], "synthetic": "generic"},
}


def attach_visual_contract(
    result: Dict[str, Any],
    *,
    markdown: str,
    image_urls: Optional[Sequence[str]] = None,
    table_data: Optional[Dict[str, Any]] = None,
    tool_id: Optional[str] = None,
) -> Dict[str, Any]:
    """将可视化字段写入工具返回 dict；若传入 tool_id 或 result['tool_id']，自动 enrich image_urls。"""
    result["markdown"] = markdown
    tid = tool_id or (result.get("tool_id") if isinstance(result.get("tool_id"), str) else None)
    result["image_urls"] = enrich_omics_visual_urls(tid, image_urls)
    if tid and "tool_id" not in result:
        result["tool_id"] = tid
    if table_data is not None:
        result["table_data"] = table_data
    inner = result.setdefault("data", {})
    if isinstance(inner, dict):
        inner["markdown"] = markdown
        inner["image_urls"] = result["image_urls"]
        if tid:
            inner["tool_id"] = tid
        if table_data is not None:
            inner["table_data"] = table_data
    return result


def simple_rows_table(
    columns: Sequence[str], rows: Sequence[Dict[str, Any]]
) -> Dict[str, Any]:
    return {"columns": list(columns), "rows": list(rows)}
