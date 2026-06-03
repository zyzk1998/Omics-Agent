# -*- coding: utf-8 -*-
"""
特色科研流程 · 2 条 STED-EC 旗舰技能 detailed_spec（对标组学管线 task 展示深度）。

前端 Featured 卡片 id → tool_id：sted_ec_trajectory / spatiotemporal_dynamics_sc
"""
from __future__ import annotations

from typing import Any, Dict, List, Optional

from gibh_agent.db.skill_detail_demo_visualizations_sted import RESEARCH_FLOW_DEMO_VIZ
from gibh_agent.db.skill_pipeline_specs import build_prompt_engineering_guide

SKILL_NAME_TO_RESEARCH_FLOW_TOOL_ID: Dict[str, str] = {
    "单细胞时空轨迹推断 (STED-EC)": "sted_ec_trajectory",
    "单细胞时空动力学分析": "spatiotemporal_dynamics_sc",
}

_MERMAID_CLASS_DEF = (
    "    classDef core fill:#e1f5fe,stroke:#0288d1,stroke-width:2px;\n"
    "    classDef algo fill:#f3e5f5,stroke:#8e24aa,stroke-width:2px;\n"
    "    classDef branchNode fill:#fff3e0,stroke:#f57c00,stroke-width:2px;"
)


def _workflow_topology_section(title: str, mermaid_body: str) -> str:
    return f"{title}\n\n```mermaid\n{mermaid_body.strip()}\n```"


def _mermaid_diagram(body: str) -> str:
    return "graph LR\n" + _MERMAID_CLASS_DEF + "\n" + body.strip()


_MERMAID_STED_EC = _mermaid_diagram("""
    H5AD["AnnData h5ad"] --> Validate["sted_ec_data_validation<br/>数据与元数据校验"]
    Validate --> Format["sted_ec_time_series_formatting<br/>时间序列标准化"]
    Format --> Moscot["sted_ec_moscot_trajectory<br/>Moscot 最优传输"]
    Moscot --> Plot["sted_ec_plot_trajectory<br/>时空 UMAP 可视化"]

    Plot --> Outputs{"交付物"}
    Outputs --> UMAPTime["时空 UMAP PNG"]
    Outputs --> UMAPType["细胞类型 UMAP"]
    Outputs --> Dynamics["演化堆叠柱状图"]
    Outputs --> ZipPkg["sted_ec_results.zip"]

    class Validate,Format,Plot core;
    class Moscot,UMAPTime,UMAPType,Dynamics,ZipPkg algo;
    class Outputs branchNode;
""")

_MERMAID_SPATIOTEMPORAL = _mermaid_diagram("""
    H5AD["AnnData h5ad"] --> Validate["数据与元数据校验"]
    Validate --> Format["时间序列标准化"]
    Format --> Moscot["Moscot OTT<br/>转移概率矩阵"]
    Moscot --> Plot["轨迹可视化<br/>UMAP 与动态图"]

    Plot --> Driver{"下游挖掘"}
    Driver --> DriverGenes["sted_ec_driver_gene_extraction<br/>Wilcoxon 驱动基因"]
    Driver --> SkipEnrich["跳过富集"]

    DriverGenes --> Enrich["sted_ec_pathway_enrichment<br/>KEGG / GO BP"]
    SkipEnrich --> Report["Orchestrator 专家解读"]

    Enrich --> VizOut{"可视化输出"}
    VizOut --> BarPlot["富集条形图"]
    VizOut --> BubblePlot["富集气泡图"]
    Enrich --> Report
    BarPlot --> Report
    BubblePlot --> Report

    class Validate,Format,Plot,DriverGenes core;
    class Moscot,Enrich,BarPlot,BubblePlot,Report algo;
    class Driver,VizOut branchNode;
""")


def _spec(
    tool_id: str,
    description_long: str,
    usage_hint: str,
    inputs: List[Dict[str, Any]],
    outputs: List[Dict[str, Any]],
    parameters_table: List[Dict[str, Any]],
    query_examples: List[Any],
    workflow_highlights: Optional[List[str]] = None,
) -> Dict[str, Any]:
    out: Dict[str, Any] = {
        "tool_id": tool_id,
        "description_long": description_long,
        "usage_hint": usage_hint,
        "inputs": inputs,
        "outputs": outputs,
        "parameters_table": parameters_table,
        "query_examples": query_examples,
        "demo_visualization": RESEARCH_FLOW_DEMO_VIZ[tool_id],
    }
    if workflow_highlights:
        out["workflow_highlights"] = workflow_highlights
    return out


_SPEC_STED_EC = _spec(
    tool_id="sted_ec_trajectory",
    description_long=(
        "基于最优传输理论（Optimal Transport）的单细胞时空轨迹推断引擎（STED-EC 通道 A · 四步 DAG）。"
        "整合 Moscot OTT 与 Scanpy 生态：自动校验 h5ad 时间列与细胞类型探针（缺失时触发聚类 + Marker 自动标注），"
        "清洗并标准化时间标签，计算相邻时间点细胞转移概率矩阵，输出时空 UMAP、细胞类型演化堆叠图与 ZIP 打包下载。"
        "适合技能广场演示与快速轨迹推断；完整驱动基因 / 通路 / 专家报告请使用「单细胞时空动力学分析」完全体。\n\n"
        + _workflow_topology_section("### 🧬 STED-EC 基础轨迹工作流拓扑图", _MERMAID_STED_EC)
    ),
    usage_hint=(
        "上传含 time_key（如 day）的 h5ad；细胞类型列可选（cell_type/celltype 等自动探针）。"
        "编排器路由 sted_ec_trajectory → STED_EC 四步卡片，勿手写 JSON 工作流块。"
    ),
    inputs=[
        {
            "name": "时空表达矩阵",
            "type": ".h5ad",
            "required": True,
            "description": "AnnData，obs 须含时间序列列（如 day）；可选 cell_type 列。",
        },
        {
            "name": "元数据列名",
            "type": "obs 列",
            "required": False,
            "description": "time_key、cell_type_key；缺省自动嗅探 day 与常见细胞类型列。",
        },
    ],
    outputs=[
        {
            "name": "validation_qc",
            "type": "PNG",
            "description": "各时间点细胞数校验图 sted_ec_validation_cells_per_time.png。",
        },
        {
            "name": "trajectory_umap",
            "type": "PNG + h5ad",
            "description": "时空 UMAP、细胞类型 UMAP、演化堆叠图；写出 sted_ec_after_umap.h5ad。",
        },
        {
            "name": "sted_ec_results_zip",
            "type": "ZIP",
            "description": "sted_ec_report_images 打包，含全部可视化 PNG。",
        },
    ],
    parameters_table=[
        {"name": "file_path", "type": "path", "required": True, "description": "h5ad 路径（左侧资产注入）。"},
        {"name": "time_key", "type": "string", "required": False, "description": "时间列名，默认 day。"},
        {"name": "cell_type_key", "type": "string", "required": False, "description": "细胞类型列；空则自动探针/标注。"},
        {"name": "n_pcs", "type": "int", "required": False, "description": "Moscot PCA 维数（默认 40）。"},
        {"name": "epsilon", "type": "float", "required": False, "description": "OT 正则化系数（默认 0.001）。"},
    ],
    query_examples=build_prompt_engineering_guide(
        "请帮我执行单细胞时空轨迹推断分析。基于最优传输理论探索细胞在时间序列上的发育与分化轨迹，采用默认参数即可。",
        "请修改默认工作流：time_key 改为 timepoint；Moscot epsilon 设为 0.005；n_pcs=50；"
        "跳过自动细胞类型标注，强制使用 obs 中的 cell_type 列。",
        "基于刚才的 UMAP 图，请解读 day_3 到 day_7 之间主要扩增的 cell type，"
        "并列出转移矩阵中概率最高的 3 条 cell-cell 迁移路径。",
    ),
    workflow_highlights=["Moscot OTT", "四步 DAG", "时空 UMAP", "细胞类型演化", "ZIP 交付"],
)

_SPEC_SPATIOTEMPORAL = _spec(
    tool_id="spatiotemporal_dynamics_sc",
    description_long=(
        "单细胞时空动力学分析完全体（SPATIOTEMPORAL_DYNAMICS 通道 B · 六步工具 DAG + 系统专家解读）。"
        "在 STED-EC 四步轨迹推断基础上，延伸 Wilcoxon 驱动基因挖掘与 gseapy Enrichr 通路富集（KEGG / GO BP），"
        "输出 driver_genes.csv、富集条形图/气泡图，并由 Orchestrator 统一生成与 RNA/代谢组同级的 AnalysisSummary 专家 Markdown。"
        "支持 kidney/lung 等演示 h5ad 与 mini_sted_ec 快速验收。\n\n"
        + _workflow_topology_section("### 🧬 时空动力学完全体工作流拓扑图", _MERMAID_SPATIOTEMPORAL)
    ),
    usage_hint=(
        "自然语言触发「时空动力学」或点击本卡片；六步 DAG 至通路富集，专家报告由编排器 Reporting 自动生成。"
        "可与左侧 sted-ec.h5ad / mini_sted_ec.h5ad 演示数据联动。"
    ),
    inputs=[
        {
            "name": "时空表达矩阵",
            "type": ".h5ad",
            "required": True,
            "description": "含时间序列与可选细胞类型注释的 AnnData。",
        },
        {
            "name": "分组/时间点",
            "type": "obs 列",
            "required": True,
            "description": "time_key 必填；cell_type_key 可自动探针。",
        },
    ],
    outputs=[
        {
            "name": "trajectory_bundle",
            "type": "PNG + ZIP",
            "description": "时空 UMAP、细胞类型演化、sted_ec_results.zip。",
        },
        {
            "name": "driver_genes_csv",
            "type": "TSV",
            "description": "按时间点/细胞类型 Wilcoxon Top 驱动基因表。",
        },
        {
            "name": "pathway_enrichment",
            "type": "PNG + CSV",
            "description": "KEGG/GO BP 富集条形图、气泡图与 enrichment_results.csv。",
        },
        {
            "name": "expert_report",
            "type": "Markdown",
            "description": "Orchestrator 生成的专家解读与数据诊断报告。",
        },
    ],
    parameters_table=[
        {"name": "file_path", "type": "path", "required": True, "description": "h5ad 路径。"},
        {"name": "time_key", "type": "string", "required": False, "description": "时间列名。"},
        {"name": "top_n", "type": "int", "required": False, "description": "驱动基因 Top N（默认 100）。"},
        {"name": "organism", "type": "string", "required": False, "description": "通路富集物种（默认 Human）。"},
        {"name": "groupby_mode", "type": "enum", "required": False, "description": "auto | time | cell_type。"},
    ],
    query_examples=build_prompt_engineering_guide(
        "帮我执行单细胞时空动力学分析。基于最优传输理论探索细胞在时间序列上的发育与分化轨迹，"
        "执行数据校验、时序标准化、轨迹推断、动态可视化、驱动基因挖掘及通路富集分析全流程。",
        "请加大 Moscot max_iterations=30；驱动基因 top_n=200；通路富集仅跑 KEGG；"
        "跳过专家报告前的冗余 QC 图重复输出。",
        "请解读 Top 5 驱动基因对应的 KEGG 通路，并结合 cell_type 演化图说明 day_5 免疫微环境重塑的关键节点。",
    ),
    workflow_highlights=[
        "Moscot OTT",
        "驱动基因",
        "KEGG/GO 富集",
        "六步完全体",
        "专家解读",
    ],
)

RESEARCH_FLOW_SPECS_BY_TOOL_ID: Dict[str, Dict[str, Any]] = {
    "sted_ec_trajectory": _SPEC_STED_EC,
    "spatiotemporal_dynamics_sc": _SPEC_SPATIOTEMPORAL,
    "spatiotemporal_dynamics": _SPEC_SPATIOTEMPORAL,
}
