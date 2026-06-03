# -*- coding: utf-8 -*-
"""
多模态组学旗舰管线 · demo_visualization（真实任务产出 + 静态 HTML 画廊）。

图片路径：/assets/images/demos/pipelines/（由 scripts/build_omics_pipeline_demo_assets.py 同步）
表格片段：源自 run_20260602_160748/markers.csv、radiomics_features.csv、spatial_pathway CSV。
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

# PyRadiomics 特征矩阵（uploads/.../radiomics_features.csv 前 5 行）
_RADIOMICS_FEATURE_ROWS = """
<tr><td><code>original_firstorder_10Percentile</code></td><td>632.0</td></tr>
<tr><td><code>original_firstorder_90Percentile</code></td><td>1044.4</td></tr>
<tr><td><code>original_firstorder_Energy</code></td><td>2.92e+09</td></tr>
<tr><td><code>original_firstorder_Entropy</code></td><td>4.602</td></tr>
<tr><td><code>original_firstorder_InterquartileRange</code></td><td>253.0</td></tr>
"""

# spatial_pathway_enrichment_top_genes_fallback.csv
_SPATIAL_TOP_GENES_ROWS = """
<tr><td>1</td><td><code>ISG15</code></td></tr>
<tr><td>2</td><td><code>IFI6</code></td></tr>
<tr><td>3</td><td><code>C1QA</code></td></tr>
<tr><td>4</td><td><code>C1QB</code></td></tr>
<tr><td>5</td><td><code>LAPTM5</code></td></tr>
"""

_T = _PIPELINE_IMG + "/transcriptomics"
_S = _PIPELINE_IMG + "/spatial"

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

DEMO_VIZ_PIPELINE_RADIOMICS = _gallery_shell(
    "影像组学 · PyRadiomics 特征提取预览",
    "来源：worker-pyskills 产出 radiomics_features.csv（真实单病例特征向量）",
    _section(
        "1. 一阶统计特征矩阵（前 5 行）",
        """<table class="skill-viz-pipeline-table">
<thead><tr><th>Feature</th><th>Value</th></tr></thead>
<tbody>"""
        + _RADIOMICS_FEATURE_ROWS
        + """</tbody></table>
<p class="skill-viz-pipeline-note">完整矩阵含 shape / firstorder / glcm 等数百维纹理特征，经 LASSO 降维后进入 LR/SVM/RF 诊断模型。</p>""",
    ),
)

# 未完全落地管线：保留专业示意 HTML（无虚假真实图路径）
def _placeholder_gallery(title: str, prov: str, stages: list[tuple[str, str]]) -> str:
    body = "".join(_section(t, c) for t, c in stages)
    return _gallery_shell(title, prov, body)


DEMO_VIZ_PIPELINE_PROTEOMICS = _placeholder_gallery(
    "蛋白组学全流程 · 标准交付结构",
    "编排 DAG 已发布；真实搜库图随 MaxQuant/DIANN 任务挂载",
    [
        ("规划步骤", """<table class="skill-viz-pipeline-table"><thead><tr><th>步骤</th><th>工具</th><th>产出</th></tr></thead>
<tbody>
<tr><td>RAW QC</td><td>msconvert</td><td>质谱质量报告</td></tr>
<tr><td>搜库</td><td>MaxQuant / DIANN</td><td>peptide/protein 矩阵</td></tr>
<tr><td>差异</td><td>limma</td><td>火山图 PNG</td></tr>
<tr><td>PPI</td><td>STRING</td><td>互作网络 HTML</td></tr>
</tbody></table>"""),
        ("说明", "<p class=\"skill-viz-pipeline-note\">环境就绪后，本栏将替换为真实火山图、热图与蛋白表（与 data/results 同步）。</p>"),
    ],
)

DEMO_VIZ_PIPELINE_EPIGENOMICS = _placeholder_gallery(
    "表观组学全流程 · 标准交付结构",
    "ATAC/ChIP 13 步 DAG；Peak/Motif 图待任务写入",
    [
        ("规划步骤", """<table class="skill-viz-pipeline-table"><thead><tr><th>步骤</th><th>工具</th><th>产出</th></tr></thead>
<tbody>
<tr><td>比对</td><td>BWA</td><td>BAM</td></tr>
<tr><td>Peak</td><td>MACS2</td><td>narrowPeak</td></tr>
<tr><td>Motif</td><td>HOMER</td><td>motif 富集表</td></tr>
<tr><td>差异</td><td>DESeq2</td><td>差异 Peak 表</td></tr>
</tbody></table>"""),
        ("说明", "<p class=\"skill-viz-pipeline-note\">组蛋白修饰（H3K27ac 等）分支与多组学整合节点已在工作流元数据中配置。</p>"),
    ],
)

DEMO_VIZ_PIPELINE_METABOLOMICS = _placeholder_gallery(
    "代谢组学 · 标志物发现交付结构",
    "OPLS-DA / VIP 主链路可执行；图表将绑定 run 目录",
    [
        ("规划步骤", """<table class="skill-viz-pipeline-table"><thead><tr><th>步骤</th><th>方法</th><th>产出</th></tr></thead>
<tbody>
<tr><td>插补</td><td>KNN</td><td>完整矩阵</td></tr>
<tr><td>模型</td><td>OPLS-DA</td><td>R2/Q2 得分图</td></tr>
<tr><td>筛选</td><td>VIP&gt;1</td><td>标志物表</td></tr>
<tr><td>通路</td><td>KEGG</td><td>富集 TSV</td></tr>
</tbody></table>"""),
        ("说明", "<p class=\"skill-viz-pipeline-note\">human_cachexia 等 benchmark 跑通后，VIP 棒棒糖与火山图将自动同步至本画廊。</p>"),
    ],
)

DEMO_VIZ_PIPELINE_GENOMICS = _placeholder_gallery(
    "基因组学 WGS/WES · 标准交付结构",
    "GATK Best Practices 12 步；VCF/报告随 FASTQ 任务生成",
    [
        ("规划步骤", """<table class="skill-viz-pipeline-table"><thead><tr><th>步骤</th><th>工具</th><th>产出</th></tr></thead>
<tbody>
<tr><td>QC</td><td>FastQC</td><td>MultiQC HTML</td></tr>
<tr><td>变异</td><td>HaplotypeCaller</td><td>gVCF/VCF</td></tr>
<tr><td>CNV</td><td>ExomeDepth</td><td>CNV 表</td></tr>
<tr><td>报告</td><td>ACMG</td><td>临床 Markdown</td></tr>
</tbody></table>"""),
        ("说明", "<p class=\"skill-viz-pipeline-note\">真实胚系/体细胞变异表与 TMB 摘要将引用 omics_genomics_report_ui 产出，禁止手写虚假坐标。</p>"),
    ],
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
    "pipeline_radiomics",
})
