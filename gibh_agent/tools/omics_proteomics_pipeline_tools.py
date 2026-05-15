"""
蛋白组学全流程原子工具

与 ProteomicsWorkflow 的 tool_id 一一对应。
首步 `proteomics_raw_qc_conversion`：对 mzML 做真实谱图/色谱计数与体积统计。
"""
from __future__ import annotations

import json
import logging
import os
import re
import tempfile
from typing import Any, Dict, List, Optional

from ..core.tool_registry import registry
from ..core.utils import safe_tool_execution
from .omics_derived_analysis import (
    mzml_derived_from_path,
    proteomics_report_markdown,
)
from .omics_pipeline_env import modal_tool_with_degradation, write_temp_mock_artifact
from .omics_mock_ui import (
    IMG_HEATMAP_CLUSTERED,
    IMG_MS_SCHEMATIC,
    IMG_PATHWAY_BUBBLE,
    IMG_PCA,
    IMG_VOLCANO,
    attach_visual_contract,
    simple_rows_table,
)
from .omics_real_io import compute_mzml_basic_stats, format_mzml_qc_markdown

logger = logging.getLogger(__name__)


def build_proteomics_database_search_cli(
    mzml_path: str,
    *,
    fragment_tol_da: float,
    missed_cleavages: int,
    peptide_fdr: float,
) -> List[str]:
    """
    DDA/DIA 搜库管线占位 CLI（MaxQuant / DIA-NN / MSFragger 语义对齐，供 Worker 与单测校验）。
    宿主环境未必安装二进制；仅保证参数字符串与业界工具字段一致。
    """
    return [
        "diann",
        "--f",
        mzml_path,
        "--frag-mass-tolerance",
        str(fragment_tol_da),
        "--missed-cleavages",
        str(missed_cleavages),
        "--peptide-fdr",
        str(peptide_fdr),
    ]


def _proteomics_derived_or_modal(
    prefix: str,
    msg: str,
    step_key: str,
    file_path: str,
    *,
    markdown_sim: str,
    image_urls: List[str],
    table_data: Dict[str, Any],
    tool_human_label: str,
    candidate_exes: Optional[List[str]] = None,
    tool_id: str = "",
) -> Dict[str, Any]:
    del step_key  # 保留签名兼容；下游禁止 mzML 统计代理冒充搜库真值
    return modal_tool_with_degradation(
        prefix,
        msg,
        markdown_sim=markdown_sim,
        image_urls=image_urls,
        table_data=table_data,
        tool_human_label=tool_human_label,
        candidate_exes=candidate_exes,
        real_runner=None,
        tool_id=tool_id or None,
    )


def _mock(
    prefix: str,
    msg: str,
    *,
    md: Optional[str] = None,
    image_urls: Optional[List[str]] = None,
    table_data: Optional[Dict[str, Any]] = None,
) -> Dict[str, Any]:
    fd, path = tempfile.mkstemp(prefix=f"{prefix}_", suffix=".mock.txt")
    os.close(fd)
    with open(path, "w", encoding="utf-8") as fh:
        fh.write(msg + "\n")
    out: Dict[str, Any] = {
        "status": "success",
        "message": msg,
        "output_path": path,
        "file_path": path,
    }
    md_body = md or f"### 蛋白质组步骤（Mock）\n\n{msg}\n"
    imgs = image_urls if image_urls is not None else [IMG_MS_SCHEMATIC]
    tbl = (
        table_data
        if table_data is not None
        else simple_rows_table(
            ("stage", "detail"),
            ({"stage": prefix, "detail": msg[:200]},),
        )
    )
    return attach_visual_contract(
        out, markdown=md_body, image_urls=imgs, table_data=tbl
    )


@registry.register(
    name="proteomics_raw_qc_conversion",
    description="mzML 原始谱图检视（真实 spectrum/chromatogram 计数）；RAW 转换需仪器厂商工具链",
    category="Proteomics",
    output_type="file_path",
)
@safe_tool_execution
def proteomics_raw_qc_conversion(
    file_path: str = "", mz_tolerance_ppm: int = 10
) -> Dict[str, Any]:
    src = (file_path or "").strip()
    if not src:
        return {"status": "error", "message": "必须提供 file_path（mzML）"}
    path = os.path.abspath(src)
    if not os.path.isfile(path):
        return {"status": "error", "message": f"输入不是可读文件: {path}"}
    low = path.lower()
    if low.endswith(".raw"):
        return {
            "status": "error",
            "message": "当前算子直接支持 mzML；请将 RAW 经 msconvert 转为 mzML 后再跑全流程。",
        }
    if not low.endswith(".mzml"):
        return {"status": "error", "message": f"首步质控当前仅支持 .mzML，收到: {path}"}

    stats = compute_mzml_basic_stats(path)
    md = format_mzml_qc_markdown(stats)

    safe_base = re.sub(r"[^\w.\-]+", "_", os.path.basename(path))[:120]
    qdir = os.path.abspath(os.getenv("RESULTS_DIR", os.path.join(os.getcwd(), "results")))
    qdir = os.path.join(qdir, "omics_qc_reports")
    os.makedirs(qdir, exist_ok=True)
    json_path = os.path.join(qdir, f"proteomics_raw_qc_{safe_base}.json")
    try:
        with open(json_path, "w", encoding="utf-8") as jf:
            json.dump({"tool_id": "proteomics_raw_qc_conversion", "qc_metrics": stats}, jf, ensure_ascii=False, indent=2)
    except OSError as exc:
        logger.warning("QC JSON 写入失败: %s", exc)
        json_path = ""

    out = {
        "status": "success",
        "message": (
            f"mzML 检视完成：{stats['spectrum_count']} spectra，"
            f"{stats['chromatogram_count']} chromatograms"
        ),
        "markdown": md,
        "qc_metrics": stats,
        "summary": (
            f"{stats['spectrum_count']} spectra, "
            f"{stats['file_size_mb']} MB"
        ),
        "output_path": json_path or path,
        "file_path": path,
    }
    return attach_visual_contract(
        out,
        markdown=md,
        image_urls=[IMG_MS_SCHEMATIC],
        table_data=simple_rows_table(
            ("metric", "value"),
            [
                {"metric": "spectrum_count", "value": str(stats["spectrum_count"])},
                {"metric": "chromatogram_count", "value": str(stats["chromatogram_count"])},
                {"metric": "file_size_mb", "value": str(stats["file_size_mb"])},
            ],
        ),
        tool_id="proteomics_raw_qc_conversion",
    )


@registry.register(
    name="proteomics_spectrum_preprocessing",
    description="基线扣除、平滑、同位素解卷积与峰拾取（Mock）",
    category="Proteomics",
    output_type="file_path",
)
@safe_tool_execution
def proteomics_spectrum_preprocessing(
    file_path: str = "", snr_threshold: float = 3.0
) -> Dict[str, Any]:
    return _proteomics_derived_or_modal(
        "p_spre",
        "Spectrum preprocessing / peak picking",
        "spre",
        file_path,
        markdown_sim=(
            "### 谱图预处理与峰拾取（仿真）\n\n"
            "- 平滑 / 基线校正（占位流程）\n"
            "- **拾取峰数（仿真）**：18,420\n"
        ),
        image_urls=[IMG_MS_SCHEMATIC],
        table_data=simple_rows_table(
            ("metric", "value"),
            [{"metric": "peaks_detected", "value": "18420"}],
        ),
        tool_human_label="msconvert / OpenMS 等预处理 CLI",
        candidate_exes=["msconvert", "FileConverter"],
        tool_id="proteomics_spectrum_preprocessing",
    )


@registry.register(
    name="proteomics_database_search",
    description="DDA/DIA 搜库（MaxQuant/DIA-NN：fragment_tol、missed_cleavages、peptide_fdr）",
    category="Proteomics",
    output_type="file_path",
)
@safe_tool_execution
def proteomics_database_search(
    file_path: str = "",
    fragment_tol_da: float = 0.05,
    missed_cleavages: int = 2,
    peptide_fdr: float = 0.01,
) -> Dict[str, Any]:
    return _proteomics_derived_or_modal(
        "p_db",
        "Database / spectral library search",
        "db",
        file_path,
        markdown_sim=(
            "### 数据库 / 谱库搜库（仿真）\n\n"
            "| 指标 | 仿真值 |\n|------|--------|\n"
            "| PSM | 1,253,400 |\n| 肽段 | 153,200 |\n"
        ),
        image_urls=[IMG_MS_SCHEMATIC],
        table_data=simple_rows_table(
            ("level", "count"),
            [
                {"level": "PSM", "count": "1253400"},
                {"level": "peptide", "count": "153200"},
            ],
        ),
        tool_human_label="DIA-NN / Comet / MSFragger 等搜库引擎",
        candidate_exes=["diann", "comet", "msfragger"],
        tool_id="proteomics_database_search",
    )


@registry.register(
    name="proteomics_fdr_rescoring",
    description="Target-Decoy + Percolator 等 FDR 与重打分（Mock）",
    category="Proteomics",
    output_type="file_path",
)
@safe_tool_execution
def proteomics_fdr_rescoring(
    file_path: str = "", target_fdr: float = 0.01
) -> Dict[str, Any]:
    return _proteomics_derived_or_modal(
        "p_fdr",
        "FDR rescoring",
        "fdr",
        file_path,
        markdown_sim=(
            "### Target–Decoy 与 Percolator 重打分（仿真）\n\n"
            "- 全局 q-value：≤ **1%**（占位）\n"
        ),
        image_urls=[IMG_VOLCANO],
        table_data=simple_rows_table(
            ("q_cutoff", "retained_psm"),
            [
                {"q_cutoff": "0.01", "retained_psm": "1180200"},
                {"q_cutoff": "0.001", "retained_psm": "992400"},
            ],
        ),
        tool_human_label="Percolator / PeptideProphet",
        candidate_exes=["percolator"],
        tool_id="proteomics_fdr_rescoring",
    )


@registry.register(
    name="proteomics_protein_inference",
    description="最大简约蛋白推断与定性分组（Mock）",
    category="Proteomics",
    output_type="file_path",
)
@safe_tool_execution
def proteomics_protein_inference(
    file_path: str = "", razor_min_peptides: int = 2
) -> Dict[str, Any]:
    return _proteomics_derived_or_modal(
        "p_inf",
        "Protein inference",
        "infer",
        file_path,
        markdown_sim=(
            "### 蛋白推断（仿真）\n\n"
            "- 推断蛋白组：**8,932** 个\n"
        ),
        image_urls=[IMG_MS_SCHEMATIC],
        table_data=simple_rows_table(
            ("metric", "value"),
            [{"metric": "protein_groups", "value": "8932"}],
        ),
        tool_human_label="ProteinProphet / MaxQuant 蛋白分组（GUI 或容器）",
        candidate_exes=["ProteinProphet", "maxquant"],
        tool_id="proteomics_protein_inference",
    )


@registry.register(
    name="proteomics_quantification",
    description="LFQ/TMT 等丰度定量（Mock）",
    category="Proteomics",
    output_type="file_path",
)
@safe_tool_execution
def proteomics_quantification(
    file_path: str = "", lfq_ratio_type: str = "Median"
) -> Dict[str, Any]:
    return _proteomics_derived_or_modal(
        "p_quant",
        "Quantification",
        "quant",
        file_path,
        markdown_sim=(
            "### LFQ 定量（仿真）\n\n"
            "下图 PCA 为公开示意图；真实管线对接矩阵后替换。\n"
        ),
        image_urls=[IMG_PCA],
        table_data=simple_rows_table(
            ("sample", "log2_intensity"),
            [
                {"sample": "A1", "log2_intensity": "18.42"},
                {"sample": "B1", "log2_intensity": "18.07"},
            ],
        ),
        tool_human_label="FlashLFQ / MaxQuant LFQ 模块",
        candidate_exes=["FlashLFQ", "maxquant"],
        tool_id="proteomics_quantification",
    )


@registry.register(
    name="proteomics_imputation_batch_correction",
    description="缺失值插补与 ComBat 等批次校正（Mock）",
    category="Proteomics",
    output_type="file_path",
)
@safe_tool_execution
def proteomics_imputation_batch_correction(
    file_path: str = "", knn_neighbors: int = 5
) -> Dict[str, Any]:
    return _proteomics_derived_or_modal(
        "p_imp_bc",
        "Imputation & batch correction",
        "imp_bc",
        file_path,
        markdown_sim="### 缺失值插补与批次校正（仿真）\n\nKNN 插补 + ComBat 占位完成。\n",
        image_urls=[IMG_PCA],
        table_data=simple_rows_table(
            ("batch", "n_samples"),
            [
                {"batch": "2025-01", "n_samples": "12"},
                {"batch": "2025-02", "n_samples": "10"},
            ],
        ),
        tool_human_label="R-commons / ComBatSeq（Rscript）",
        candidate_exes=["Rscript"],
        tool_id="proteomics_imputation_batch_correction",
    )


@registry.register(
    name="proteomics_normalization_qc",
    description="中位数/分位数标准化与 PCA、相关性 QC（Mock）",
    category="Proteomics",
    output_type="file_path",
)
@safe_tool_execution
def proteomics_normalization_qc(
    file_path: str = "", norm_method: str = "quantile"
) -> Dict[str, Any]:
    return _proteomics_derived_or_modal(
        "p_norm_qc",
        "Normalization & QC",
        "norm",
        file_path,
        markdown_sim=(
            "### 归一化与样本相关性 QC（仿真）\n\n"
            "| 指标 | 值 |\n|------|-----|\n"
            "| 中位数归一后 CV | 0.18 |\n"
        ),
        image_urls=[IMG_PCA],
        table_data=simple_rows_table(
            ("metric", "value"),
            [{"metric": "median_cv", "value": "0.18"}],
        ),
        tool_human_label="R / Python 归一化脚本环境",
        candidate_exes=["Rscript", "python3"],
        tool_id="proteomics_normalization_qc",
    )


@registry.register(
    name="proteomics_differential_analysis",
    description="Limma 等差异表达（Mock）",
    category="Proteomics",
    output_type="file_path",
)
@safe_tool_execution
def proteomics_differential_analysis(
    file_path: str = "",
    group_column: str = "",
    log2fc_cutoff: float = 1.0,
) -> Dict[str, Any]:
    out = _proteomics_derived_or_modal(
        "p_dea",
        "Differential analysis",
        "dea",
        file_path,
        markdown_sim=(
            "### 差异表达分析（仿真 / Limma）\n\n"
            "| 对比 | 上调 | 下调 |\n"
            "|------|------|------|\n"
            "| Case vs Ctrl | 412 | 355 |\n"
        ),
        image_urls=[IMG_VOLCANO, IMG_HEATMAP_CLUSTERED],
        table_data=simple_rows_table(
            ("protein_id", "logFC", "adj_p"),
            [
                {"protein_id": "P12345", "logFC": "1.82", "adj_p": "0.003"},
                {"protein_id": "Q98765", "logFC": "-1.41", "adj_p": "0.011"},
            ],
        ),
        tool_human_label="Limma / DEP 等差异分析环境（R）",
        candidate_exes=["Rscript"],
        tool_id="proteomics_differential_analysis",
    )
    if group_column:
        out["group_column"] = group_column
    return out


@registry.register(
    name="proteomics_biomarker_discovery",
    description="随机森林/SVM/LASSO 标志物筛选（Mock）",
    category="Proteomics",
    output_type="file_path",
)
@safe_tool_execution
def proteomics_biomarker_discovery(
    file_path: str = "", n_top_features: int = 50
) -> Dict[str, Any]:
    return _proteomics_derived_or_modal(
        "p_bio",
        "Biomarker ML",
        "bio",
        file_path,
        markdown_sim=(
            "### 机器学习标志物筛选（仿真）\n\n"
            "- RF OOB AUC：**0.91**\n"
        ),
        image_urls=[IMG_VOLCANO],
        table_data=simple_rows_table(
            ("feature", "importance"),
            [
                {"feature": "PROT_001", "importance": "0.142"},
                {"feature": "PROT_088", "importance": "0.098"},
            ],
        ),
        tool_human_label="Python/R 机器学习栈（sklearn / caret）",
        candidate_exes=["python3", "Rscript"],
        tool_id="proteomics_biomarker_discovery",
    )


@registry.register(
    name="proteomics_functional_enrichment",
    description="GO/KEGG/Reactome/GSEA（Mock）",
    category="Proteomics",
    output_type="file_path",
)
@safe_tool_execution
def proteomics_functional_enrichment(
    file_path: str = "", enrich_padj: float = 0.05
) -> Dict[str, Any]:
    return _proteomics_derived_or_modal(
        "p_enr",
        "Functional enrichment",
        "enr",
        file_path,
        markdown_sim=(
            "### GO / KEGG 富集（仿真）\n\n"
            "| Term | FDR |\n|------|-----|\n"
            "| R-HSA-6900 | 2.1e-5 |\n"
        ),
        image_urls=[IMG_PATHWAY_BUBBLE, IMG_VOLCANO],
        table_data=simple_rows_table(
            ("pathway", "fdr"),
            [
                {"pathway": "R-HSA-6900", "fdr": "2.1e-5"},
                {"pathway": "GO:0006955", "fdr": "4.3e-4"},
            ],
        ),
        tool_human_label="clusterProfiler / gProfiler / STRING API 客户端",
        candidate_exes=["Rscript"],
        tool_id="proteomics_functional_enrichment",
    )


@registry.register(
    name="proteomics_ppi_network_analysis",
    description="STRING PPI 与 Hub 识别（Mock）",
    category="Proteomics",
    output_type="file_path",
)
@safe_tool_execution
def proteomics_ppi_network_analysis(
    file_path: str = "", string_score_cutoff: int = 400
) -> Dict[str, Any]:
    return _proteomics_derived_or_modal(
        "p_ppi",
        "PPI network analysis",
        "ppi",
        file_path,
        markdown_sim=(
            "### STRING PPI 网络（仿真）\n\n"
            "- Hub 蛋白：**EGFR**, **TP53**（示例）\n"
        ),
        image_urls=[IMG_MS_SCHEMATIC],
        table_data=simple_rows_table(
            ("node", "degree"),
            [
                {"node": "EGFR", "degree": "42"},
                {"node": "TP53", "degree": "38"},
            ],
        ),
        tool_human_label="Cytoscape / stringApp 或本地网络脚本",
        candidate_exes=["cytoscape", "Rscript"],
        tool_id="proteomics_ppi_network_analysis",
    )


@registry.register(
    name="proteomics_clinical_reporting",
    description="多维全景标准化报告（Mock；PDF）",
    category="Proteomics",
    output_type="file_path",
)
@safe_tool_execution
def proteomics_clinical_reporting(
    file_path: str = "", report_depth: str = "standard"
) -> Dict[str, Any]:
    out = _mock("p_report", "Clinical / panorama reporting placeholder")
    out["report_path"] = out["output_path"]
    out["pdf_url"] = None
    out["html_url"] = None
    depth_line = f"- **报告深度**: `{report_depth}`\n\n"
    m = mzml_derived_from_path(file_path)
    if m:
        out["markdown"] = depth_line + proteomics_report_markdown(m, file_path or "")
        out["omics_analysis_mode"] = "mzml_derived_proxy"
        out["derived_mzml_metrics"] = m
    else:
        out["markdown"] = depth_line + (
            "## 蛋白质组学全景报告（Mock）\n\n"
            "- **差异蛋白**：占位 Top-N 列表（真实管线对接 Limma/矩阵）。\n"
            "- **通路**：GO/KEGG 富集占位。\n\n"
            "输入：`{inp}`\n"
        ).format(inp=file_path or "（未指定）")
    out["image_urls"] = [IMG_MS_SCHEMATIC]
    out["json_url"] = None
    out["table_data"] = {
        "proteins_preview": [],
        "pathway_hits": [],
        "ingress_file_path": file_path or "",
    }
    out["data"] = {
        "report_stage": "proteomics_clinical_reporting",
        "ingress_file_path": file_path or "",
        "pdf_url": out["pdf_url"],
        "html_url": out["html_url"],
        "table_data": out["table_data"],
    }
    return attach_visual_contract(
        out,
        markdown=out["markdown"],
        image_urls=out["image_urls"],
        table_data=out["table_data"],
        tool_id="proteomics_clinical_reporting",
    )
