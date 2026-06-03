# -*- coding: utf-8 -*-
"""特色科研流程 · STED-EC / 时空动力学 demo_visualization（真实教程产出画廊）。"""
from __future__ import annotations

_STED = "/assets/demo/spatiotemporal_tutorial"


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


_DRIVER_GENES_ROWS = """
<tr><td><code>MKI67</code></td><td>day_3</td><td>0.0012</td></tr>
<tr><td><code>TOP2A</code></td><td>day_3</td><td>0.0028</td></tr>
<tr><td><code>CDK1</code></td><td>day_5</td><td>0.0041</td></tr>
<tr><td><code>STMN1</code></td><td>day_5</td><td>0.0067</td></tr>
<tr><td><code>HMGB2</code></td><td>day_7</td><td>0.0093</td></tr>
"""

DEMO_VIZ_STED_EC_TRAJECTORY = _gallery_shell(
    "STED-EC 基础轨迹推断 · 真实交付画廊",
    "来源：spatiotemporal_tutorial · Moscot OTT 四步 DAG 产出",
    _section(
        "1. 数据校验 · 各时间点细胞数",
        _img(f"{_STED}/sted_ec_validation_cells_per_time.png", "Cells per time point validation"),
    )
    + _section(
        "2. 时序标准化 · 清洗后分布",
        _img(f"{_STED}/sted_ec_format_cells_per_time.png", "Formatted time series cell counts"),
    )
    + _section(
        "3. 时空 UMAP · 按时间点着色",
        _img(f"{_STED}/umap_time.png", "UMAP colored by time"),
    )
    + _section(
        "4. 细胞类型 UMAP 与演化动态",
        _img(f"{_STED}/umap_celltype.png", "UMAP colored by cell type")
        + _img(f"{_STED}/cell_type_dynamics.png", "Cell type dynamics stacked bar"),
    ),
)

DEMO_VIZ_SPATIOTEMPORAL_DYNAMICS = _gallery_shell(
    "单细胞时空动力学完全体 · 真实交付画廊",
    "来源：spatiotemporal_tutorial · 六步工具 DAG + 专家解读",
    _section(
        "1. 数据校验 · 各时间点细胞数",
        _img(f"{_STED}/sted_ec_validation_cells_per_time.png", "Cells per time point validation"),
    )
    + _section(
        "2. 最优传输轨迹 · 时空 UMAP",
        _img(f"{_STED}/umap_time.png", "UMAP colored by time"),
    )
    + _section(
        "3. 细胞类型演化",
        _img(f"{_STED}/umap_celltype.png", "UMAP colored by cell type")
        + _img(f"{_STED}/cell_type_dynamics.png", "Cell type dynamics"),
    )
    + _section(
        "4. 驱动基因候选（示意 TSV 摘要）",
        """<table class="skill-viz-pipeline-table">
<thead><tr><th>Gene</th><th>Group</th><th>adj.P</th></tr></thead>
<tbody>"""
        + _DRIVER_GENES_ROWS
        + """</tbody></table>""",
    )
    + _section(
        "5. 通路富集 · 条形图",
        _img(f"{_STED}/sted_ec_pathway_enrichment_bar.png", "Pathway enrichment bar plot"),
    )
    + _section(
        "6. 通路富集 · 气泡图",
        _img(f"{_STED}/sted_ec_pathway_enrichment_bubble.png", "Pathway enrichment bubble plot"),
    ),
)

RESEARCH_FLOW_DEMO_VIZ = {
    "sted_ec_trajectory": DEMO_VIZ_STED_EC_TRAJECTORY,
    "spatiotemporal_dynamics_sc": DEMO_VIZ_SPATIOTEMPORAL_DYNAMICS,
    "spatiotemporal_dynamics": DEMO_VIZ_SPATIOTEMPORAL_DYNAMICS,
}
