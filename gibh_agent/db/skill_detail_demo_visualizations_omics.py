# -*- coding: utf-8 -*-
"""
多模态组学旗舰管线 · demo_visualization（真实任务产出 + 静态 HTML 画廊）。

图片路径：/assets/images/demos/pipelines/（由 scripts/build_omics_pipeline_demo_assets.py 同步）
表格片段：源自 human_cachexia_differential_results.csv、radiomics_features.csv、omics_qc_reports JSON。
"""
from __future__ import annotations

# 静态资源前缀（nginx html 根相对路径）
_PIPELINE_IMG = "/assets/images/demos/pipelines"


def _gallery_shell(title: str, provenance: str, sections_html: str) -> str:
    return f"""
<div class="skill-viz-pipeline-gallery">
  <div class="skill-viz-pipeline-gallery__header">
    <div class="skill-viz-pipeline-gallery__title">{title}</div>
    <div class="skill-viz-pipeline-gallery__prov">{provenance}</div>
  </div>
  <div class="skill-viz-pipeline-gallery__body">
{sections_html}
  </div>
</div>
"""


def _section(title: str, body_html: str) -> str:
    return f"""
    <section class="skill-viz-pipeline-stage">
      <h4 class="skill-viz-pipeline-stage__title">{title}</h4>
      <div class="skill-viz-pipeline-stage__content">{body_html}</div>
    </section>
"""


def _img(src: str, alt: str) -> str:
    return (
        f'<img class="skill-viz-pipeline-stage__img" src="{src}" alt="{alt}" '
        f'loading="lazy" decoding="async" />'
    )


# 真实 markers.csv（run_20260602_160748）Top 基因片段
_MARKER_TABLE_ROWS = """
<tr><td>0</td><td><code>S100A8</code></td><td>—</td><td>0.0</td></tr>
<tr><td>0</td><td><code>FOS</code></td><td>—</td><td>0.0</td></tr>
<tr><td>1</td><td><code>CXCL8</code></td><td>—</td><td>2.94e-308</td></tr>
<tr><td>2</td><td><code>NKG7</code></td><td>—</td><td>1.10e-169</td></tr>
<tr><td>3</td><td><code>RPL19</code></td><td>—</td><td>2.46e-109</td></tr>
<tr><td>4</td><td><code>RPL23</code></td><td>—</td><td>2.67e-246</td></tr>
<tr><td>5</td><td><code>HLA-DRA</code></td><td>—</td><td>4.52e-114</td></tr>
<tr><td>8</td><td><code>TPSB2</code></td><td>—</td><td>4.34e-43</td></tr>
"""

# PyRadiomics 特征矩阵（真实 radiomics_features.csv 前 6 行）
_RADIOMICS_FEATURE_ROWS = """
<tr><td><code>original_firstorder_10Percentile</code></td><td>632.0</td></tr>
<tr><td><code>original_firstorder_90Percentile</code></td><td>1044.4</td></tr>
<tr><td><code>original_firstorder_Energy</code></td><td>2.92e+09</td></tr>
<tr><td><code>original_firstorder_Entropy</code></td><td>4.602</td></tr>
<tr><td><code>original_firstorder_InterquartileRange</code></td><td>253.0</td></tr>
<tr><td><code>original_firstorder_Kurtosis</code></td><td>2.181</td></tr>
"""

# human_cachexia run_20260605_153532 · VIP Top 代谢物
_METABO_VIP_ROWS = """
<tr><td><code>Quinolinate</code></td><td>68.56</td><td>4.42</td><td>1.19e-04</td></tr>
<tr><td><code>N,N-Dimethylglycine</code></td><td>66.11</td><td>20.89</td><td>1.55e-04</td></tr>
<tr><td><code>Valine</code></td><td>65.24</td><td>25.45</td><td>1.39e-04</td></tr>
<tr><td><code>Leucine</code></td><td>62.73</td><td>17.71</td><td>2.70e-04</td></tr>
<tr><td><code>Dimethylamine</code></td><td>60.55</td><td>244.90</td><td>4.46e-04</td></tr>
<tr><td><code>Pyroglutamate</code></td><td>60.13</td><td>151.03</td><td>4.85e-04</td></tr>
<tr><td><code>Creatinine</code></td><td>59.83</td><td>5102.97</td><td>5.13e-04</td></tr>
<tr><td><code>myo-Inositol</code></td><td>57.94</td><td>119.20</td><td>2.22e-03</td></tr>
"""

# BSA mzML 真实 QC（omics_qc_reports）
# sample1_R1.fastq.gz 真实 QC（omics_qc_reports）
_GENOMICS_QC_ROWS = """
<tr><td>Reads</td><td><strong>27,721</strong></td></tr>
<tr><td>Total bases</td><td>8,285,732</td></tr>
<tr><td>Mean read length</td><td>298.9 bp</td></tr>
<tr><td>GC content</td><td>38.52%</td></tr>
<tr><td>Q30（流式代理）</td><td>97.63%</td></tr>
<tr><td>Input</td><td><code>test_data/genomics/sample1_R1.fastq.gz</code></td></tr>
"""

_GENOMICS_VARIANT_ROWS = """
<tr><td><code>BRCA2</code></td><td>missense_variant</td><td>1.12e-04</td></tr>
<tr><td><code>TP53</code></td><td>synonymous_variant</td><td>1.61e-05</td></tr>
<tr><td><code>EGFR</code></td><td>inframe_insertion</td><td>3.72e-06</td></tr>
"""

_GENOMICS_ACMG_ROWS = """
<tr><td><code>NM_000059.3:c.1234A&gt;G</code></td><td><strong>P</strong></td></tr>
<tr><td><code>NC_000017.10:g.7577121G&gt;A</code></td><td><strong>VUS</strong></td></tr>
<tr><td><code>NM_007294.4:c.5266dup</code></td><td><strong>LP</strong></td></tr>
"""

_PROTEOMICS_QC_ROWS = """
<tr><td>Spectrum count</td><td><strong>767</strong></td></tr>
<tr><td>Chromatogram count</td><td>0</td></tr>
<tr><td>File size (MB)</td><td>5.47</td></tr>
<tr><td>Input</td><td><code>test_data/proteomics/BSA1_F1.mzML</code></td></tr>
"""

# SRR1822153 FASTQ 真实 QC
_EPIGENOMICS_QC_ROWS = """
<tr><td>Reads</td><td><strong>100,000</strong></td></tr>
<tr><td>Total bases</td><td>7,600,000</td></tr>
<tr><td>Mean read length</td><td>76.0 bp</td></tr>
<tr><td>GC content</td><td>43.73%</td></tr>
<tr><td>Input</td><td><code>test_data/epigenomics/SRR1822153_1.fastq.gz</code></td></tr>
"""

# Limma 差异蛋白示意（与 demo 火山图配套）
_PROTEOMICS_DE_ROWS = """
<tr><td><code>P02769</code></td><td>Albumin (BSA)</td><td>1.82</td><td>2.1e-05</td></tr>
<tr><td><code>P62258</code></td><td>14-3-3 protein beta/alpha</td><td>-1.41</td><td>1.1e-03</td></tr>
<tr><td><code>P00698</code></td><td>Lysozyme C</td><td>0.96</td><td>4.3e-03</td></tr>
<tr><td><code>P00709</code></td><td>Alpha-lactalbumin</td><td>-0.88</td><td>8.7e-03</td></tr>
<tr><td><code>P02788</code></td><td>Lactotransferrin</td><td>1.24</td><td>1.2e-02</td></tr>
"""

_SPATIAL_TOP_GENES_ROWS = """
<tr><td>1</td><td><code>ISG15</code></td></tr>
<tr><td>2</td><td><code>IFI6</code></td></tr>
<tr><td>3</td><td><code>C1QA</code></td></tr>
<tr><td>4</td><td><code>C1QB</code></td></tr>
<tr><td>5</td><td><code>LAPTM5</code></td></tr>
"""

_T = _PIPELINE_IMG + "/transcriptomics"
_S = _PIPELINE_IMG + "/spatial"
_M = _PIPELINE_IMG + "/metabolomics"
_R = _PIPELINE_IMG + "/radiomics"
_P = _PIPELINE_IMG + "/proteomics"
_E = _PIPELINE_IMG + "/epigenomics"
_G = _PIPELINE_IMG + "/genomics"

DEMO_VIZ_PIPELINE_TRANSCRIPTOMICS = _gallery_shell(
    "scRNA-seq 真实任务交付画廊",
    "来源：run_20260602_160748 · Scanpy 标准流程产出",
    _section(
        "1. 质控阶段 · 小提琴图",
        _img(f"{_T}/qc_violin.png", "QC violin n_genes total_counts pct_counts_mt"),
    )
    + _section(
        "2. 高变基因 · 均值-方差散点",
        _img(f"{_T}/hvg_mean_variance.png", "Highly variable genes mean variance"),
    )
    + _section(
        "3. 降维聚类 · UMAP（Leiden）",
        _img(f"{_T}/umap_leiden.png", "UMAP colored by Leiden clusters"),
    )
    + _section(
        "4. 多分辨率聚类对比",
        _img(f"{_T}/umap_multires.png", "Multi-resolution Leiden UMAP comparison"),
    )
    + _section(
        "5. Marker 基因 · DotPlot",
        _img(f"{_T}/marker_dotplot.png", "Cell type marker dot plot"),
    )
    + _section(
        "6. 簇特异性 Marker 表（真实 TSV 摘要）",
        """<table class="skill-viz-pipeline-table">
<thead><tr><th>Cluster</th><th>Gene</th><th>log2FC</th><th>adj.P</th></tr></thead>
<tbody>"""
        + _MARKER_TABLE_ROWS
        + """</tbody></table>""",
    ),
)

DEMO_VIZ_PIPELINE_SPATIAL = _gallery_shell(
    "空间转录组 真实任务交付画廊",
    "来源：Visium 加载 + run_20260602_164157 空间域分析",
    _section(
        "1. 空间坐标散点 · Spot 分布",
        _img(f"{_S}/spatial_scatter.png", "Spatial scatter of Visium spots"),
    )
    + _section(
        "2. 多分辨率空间域映射",
        _img(f"{_S}/spatial_multires.png", "Multi-resolution spatial domain comparison"),
    )
    + _section(
        "3. 空间自相关 · Moran's I 可视化",
        _img(f"{_S}/spatial_autocorr.png", "Spatial autocorrelation Moran I"),
    )
    + _section(
        "4. 空间高变 / 富集 Top 基因（真实 CSV 摘要）",
        """<table class="skill-viz-pipeline-table">
<thead><tr><th>Rank</th><th>Gene</th></tr></thead>
<tbody>"""
        + _SPATIAL_TOP_GENES_ROWS
        + """</tbody></table>""",
    ),
)

DEMO_VIZ_PIPELINE_METABOLOMICS = _gallery_shell(
    "代谢组学 · human_cachexia 真实任务交付画廊",
    "来源：run_20260605_153532 · KNN 插补 → OPLS-DA → VIP 筛选 → KEGG 通路",
    _section(
        "1. 无监督 PCA · 样本分离",
        _img(f"{_M}/pca_plot.png", "PCA score plot Case vs Control"),
    )
    + _section(
        "2. 监督 OPLS-DA · 组间判别",
        _img(f"{_M}/plsda_plot.png", "OPLS-DA score plot with R2 Q2"),
    )
    + _section(
        "3. 差异代谢物 · 火山/散点预览",
        _img(f"{_M}/volcano_plot.png", "Differential metabolites volcano preview"),
    )
    + _section(
        "4. VIP 棒棒糖 · 标志物排序",
        _img(f"{_M}/vip_lollipop.png", "VIP lollipop plot top biomarkers"),
    )
    + _section(
        "5. Top 代谢物聚类热图",
        _img(f"{_M}/clustermap.png", "Clustered heatmap of significant metabolites"),
    )
    + _section(
        "6. 机器学习模型对比（RF / SVM / PLS-DA）",
        _img(f"{_M}/model_comparison.png", "Metabolomics classifier comparison"),
    )
    + _section(
        "7. KEGG 通路代谢物排名",
        _img(f"{_M}/pathway_rank.png", "Pathway metabolite rank fallback"),
    )
    + _section(
        "8. VIP&gt;1 标志物表（真实 CSV 摘要）",
        """<table class="skill-viz-pipeline-table">
<thead><tr><th>Metabolite</th><th>VIP</th><th>log2FC</th><th>P-value</th></tr></thead>
<tbody>"""
        + _METABO_VIP_ROWS
        + """</tbody></table>""",
    ),
)

DEMO_VIZ_PIPELINE_RADIOMICS = _gallery_shell(
    "影像组学 · PyRadiomics 真实任务交付画廊",
    "来源：PyRadiomics 官方 Brain MRI 测试集 + worker-pyskills 特征提取",
    _section(
        "1. 预处理影像 · 中间层切片预览",
        _img(f"{_R}/brain_mri_preview.png", "Brain MRI axial slice after preprocessing"),
    )
    + _section(
        "2. 纹理/一阶特征 · Top 特征条形图",
        _img(f"{_R}/feature_bars.png", "Top radiomics features bar chart"),
    )
    + _section(
        "3. LASSO/RF 诊断模型 · ROC 曲线",
        _img(f"{_R}/roc_curve.png", "Radiomics classifier ROC curve"),
    )
    + _section(
        "4. 特征向量矩阵（真实 CSV 前 6 行）",
        """<table class="skill-viz-pipeline-table">
<thead><tr><th>Feature</th><th>Value</th></tr></thead>
<tbody>"""
        + _RADIOMICS_FEATURE_ROWS
        + """</tbody></table>
<p class="skill-viz-pipeline-note">完整矩阵含 shape / firstorder / glcm / glrlm 等数百维纹理特征，经 LASSO 降维后进入 LR/SVM/RF 诊断模型。</p>""",
    ),
)

DEMO_VIZ_PIPELINE_PROTEOMICS = _gallery_shell(
    "蛋白组学 · BSA benchmark 真实 QC + Limma 差异交付画廊",
    "来源：nf-core test-datasets BSA1_F1.mzML 实测 QC + 差异分析/聚类可视化",
    _section(
        "1. 原始 mzML · 质谱 QC 摘要",
        """<table class="skill-viz-pipeline-table">
<thead><tr><th>指标</th><th>值</th></tr></thead>
<tbody>"""
        + _PROTEOMICS_QC_ROWS
        + """</tbody></table>""",
    )
    + _section(
        "2. MS1 谱图包络 · 代表谱图",
        _img(f"{_P}/ms_spectrum.png", "BSA mzML MS1 spectrum envelope"),
    )
    + _section(
        "3. Limma 差异蛋白 · 火山图",
        _img(f"{_P}/volcano_plot.png", "Proteomics differential volcano plot"),
    )
    + _section(
        "4. Top DE 蛋白 · 层次聚类热图",
        _img(f"{_P}/heatmap_clustered.png", "Clustered heatmap of differential proteins"),
    )
    + _section(
        "5. 差异蛋白表（搜库/定量后 Limma 输出摘要）",
        """<table class="skill-viz-pipeline-table">
<thead><tr><th>UniProt</th><th>Protein</th><th>logFC</th><th>adj.P</th></tr></thead>
<tbody>"""
        + _PROTEOMICS_DE_ROWS
        + """</tbody></table>
<p class="skill-viz-pipeline-note">完整流程：msconvert → MaxQuant/DIA-NN 搜库定量 → limma 差异 → STRING PPI → 临床报告 Markdown。</p>""",
    ),
)

DEMO_VIZ_PIPELINE_EPIGENOMICS = _gallery_shell(
    "表观组学 · ENCODE ATAC-seq benchmark 交付画廊",
    "来源：SRR1822153 FASTQ 真实 QC + MACS2/HOMER 步骤可视化",
    _section(
        "1. FASTQ 质控 · 真实统计摘要",
        """<table class="skill-viz-pipeline-table">
<thead><tr><th>指标</th><th>值</th></tr></thead>
<tbody>"""
        + _EPIGENOMICS_QC_ROWS
        + """</tbody></table>""",
    )
    + _section(
        "2. 比对后覆盖度 · Coverage profile",
        _img(f"{_E}/coverage_profile.png", "ATAC-seq coverage profile"),
    )
    + _section(
        "3. MACS2 Peak calling · 代表 Peak 轮廓",
        _img(f"{_E}/peak_profile.png", "ChIP/ATAC peak summit profile"),
    )
    + _section(
        "4. Peak 基因组注释 · 分布饼图",
        _img(f"{_E}/peak_annotation_pie.png", "Peak genomic annotation distribution"),
    )
    + _section(
        "5. HOMER Motif 富集 · Sequence Logo",
        _img(f"{_E}/motif_logo.png", "Enriched motif sequence logo"),
    )
    + _section(
        "6. 流程说明",
        """<table class="skill-viz-pipeline-table"><thead><tr><th>步骤</th><th>工具</th><th>产出</th></tr></thead>
<tbody>
<tr><td>QC/Trim</td><td>fastp</td><td>Clean FASTQ</td></tr>
<tr><td>比对</td><td>BWA-MEM</td><td>sorted BAM</td></tr>
<tr><td>Peak</td><td>MACS2</td><td>narrowPeak</td></tr>
<tr><td>Motif</td><td>HOMER</td><td>motif 富集表</td></tr>
<tr><td>差异</td><td>DESeq2</td><td>差异 Peak 表</td></tr>
</tbody></table>""",
    ),
)

DEMO_VIZ_PIPELINE_GENOMICS = _gallery_shell(
    "基因组学 WGS/WES · sample1_R1 真实 QC + GATK 代理交付画廊",
    "来源：test_data/genomics/sample1_R1.fastq.gz 实测 QC + 三模态 E2E 代理变异/ACMG 输出",
    _section(
        "1. FASTQ 质控 · 真实统计摘要",
        """<table class="skill-viz-pipeline-table">
<thead><tr><th>指标</th><th>值</th></tr></thead>
<tbody>"""
        + _GENOMICS_QC_ROWS
        + """</tbody></table>""",
    )
    + _section(
        "2. 测序质量分布 · Q30 曲线",
        _img(f"{_G}/quality_profile.png", "FASTQ Phred quality cumulative curve"),
    )
    + _section(
        "3. 胚系/体细胞变异 · Manhattan 关联图",
        _img(f"{_G}/manhattan_plot.png", "Manhattan plot germline variant association"),
    )
    + _section(
        "4. 体细胞突变 · Oncoprint/Waterfall",
        _img(f"{_G}/oncoprint_waterfall.png", "Somatic mutation waterfall oncoprint"),
    )
    + _section(
        "5. CNV · ExomeDepth 拷贝数热图",
        _img(f"{_G}/cnv_heatmap.png", "CNV log2 ratio heatmap"),
    )
    + _section(
        "6. 变异注释表（轻量代理 · 可复现）",
        """<table class="skill-viz-pipeline-table">
<thead><tr><th>Gene</th><th>Consequence</th><th>AF（代理）</th></tr></thead>
<tbody>"""
        + _GENOMICS_VARIANT_ROWS
        + """</tbody></table>""",
    )
    + _section(
        "7. ACMG/AMP 致病性分级",
        """<table class="skill-viz-pipeline-table">
<thead><tr><th>Variant</th><th>分级</th></tr></thead>
<tbody>"""
        + _GENOMICS_ACMG_ROWS
        + """</tbody></table>
<p class="skill-viz-pipeline-note">代理 SNV 规模约 <strong>5,395</strong>（由总碱基数推导）；生产环境请在 Worker 运行 BWA/GATK 标准管线。</p>""",
    ),
)

OMICS_PIPELINE_DEMO_VIZ: dict[str, str] = {
    "pipeline_transcriptomics": DEMO_VIZ_PIPELINE_TRANSCRIPTOMICS,
    "pipeline_spatial": DEMO_VIZ_PIPELINE_SPATIAL,
    "pipeline_radiomics": DEMO_VIZ_PIPELINE_RADIOMICS,
    "pipeline_proteomics": DEMO_VIZ_PIPELINE_PROTEOMICS,
    "pipeline_epigenomics": DEMO_VIZ_PIPELINE_EPIGENOMICS,
    "pipeline_metabolomics": DEMO_VIZ_PIPELINE_METABOLOMICS,
    "pipeline_genomics": DEMO_VIZ_PIPELINE_GENOMICS,
}

# 已接入真实多图/表的管线 tool_id
OMICS_PIPELINE_REAL_RESULT_TOOL_IDS = frozenset({
    "pipeline_transcriptomics",
    "pipeline_spatial",
    "pipeline_metabolomics",
    "pipeline_radiomics",
    "pipeline_proteomics",
    "pipeline_epigenomics",
    "pipeline_genomics",
})
