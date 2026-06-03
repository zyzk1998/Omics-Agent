# -*- coding: utf-8 -*-
"""
多模态组学 · 7 条旗舰管线 detailed_spec（大师级说明书 + 结构化 Prompt 工程指南）。

通过技能广场 **name** → pipeline tool_id 解析；与 [Skill_Route] / [Omics_Route] 执行路由解耦。
结构见 docs/技能详情 Demo 页面规范.md §2.1（query_examples 支持对象数组）。
"""
from __future__ import annotations

from typing import Any, Dict, List, Optional

from gibh_agent.db.skill_detail_demo_visualizations_omics import OMICS_PIPELINE_DEMO_VIZ


def _workflow_topology_section(title: str, mermaid_body: str) -> str:
    """将 Mermaid 工作流蓝图嵌入 description_long Markdown（节点 ID 纯英文，标签双引号）。"""
    return f"{title}\n\n```mermaid\n{mermaid_body.strip()}\n```"


def _mermaid_diagram(body: str) -> str:
    """classDef 置于节点定义之前，避免 Mermaid 10 解析顺序问题。"""
    return "graph LR\n" + _MERMAID_CLASS_DEF + "\n" + body.strip()


_MERMAID_CLASS_DEF = (
    "    classDef core fill:#e1f5fe,stroke:#0288d1,stroke-width:2px;\n"
    "    classDef algo fill:#f3e5f5,stroke:#8e24aa,stroke-width:2px;\n"
    "    classDef branchNode fill:#fff3e0,stroke:#f57c00,stroke-width:2px;"
)

_MERMAID_TRANSCRIPTOMICS = _mermaid_diagram("""
    Raw["Raw FASTQ"] --> CR["Cell Ranger<br/>比对定量"]
    CR --> QC["Scanpy QC<br/>质控过滤"]
    QC --> Scrublet["Scrublet<br/>双细胞去除"]

    Scrublet --> Batch{"批次效应校正"}
    Batch -->|推荐| Harmony["Harmony"]
    Batch -->|图神经网络| BBKNN["BBKNN"]
    Batch -->|深度生成模型| scVI["scVI"]

    Harmony --> Annotate["CellTypist + Marker<br/>智能细胞注释"]
    BBKNN --> Annotate
    scVI --> Annotate

    Annotate --> Downstream{"高级下游分析"}
    Downstream --> DEG["差异基因分析 (DEG)"]
    Downstream --> CellChat["CellChat 细胞通讯"]
    Downstream --> Trajectory["Monocle3 / PAGA<br/>发育轨迹推断"]

    class CR,QC,Scrublet,Annotate core;
    class Harmony,BBKNN,scVI,DEG,CellChat,Trajectory algo;
    class Batch,Downstream branchNode;
""")

_MERMAID_SPATIAL = _mermaid_diagram("""
    Visium["Visium / Stereo-seq"] --> Validate["Squidpy<br/>空间数据校验"]
    Validate --> QC["Scanpy QC<br/>质控过滤"]
    QC --> SVG["Moran I<br/>空间高变基因"]
    SVG --> SGraph["空间邻域图<br/>Squidpy"]
    SGraph --> Cluster["Leiden<br/>空间域聚类"]
    Cluster --> MultiRes["多分辨率<br/>物理映射对比"]

    MultiRes --> Deconv{"空间解卷积"}
    Deconv -->|贝叶斯推断| C2L["Cell2location"]
    Deconv -->|稳健回归| RCTD["RCTD"]
    Deconv -->|跳过| Annotate["空间域<br/>生物学注释"]

    C2L --> Annotate
    RCTD --> Annotate

    Annotate --> Downstream{"空间下游分析"}
    Downstream --> SpaDEG["空间差异表达"]
    Downstream --> Niche["微环境识别"]
    Downstream --> Enrich["通路富集"]

    class Validate,QC,SVG,SGraph,Cluster,MultiRes,Annotate core;
    class C2L,RCTD,SpaDEG,Niche,Enrich algo;
    class Deconv,Downstream branchNode;
""")

_MERMAID_PROTEOMICS = _mermaid_diagram("""
    RawMS["RAW / mzML"] --> MSQC["谱图质控<br/>MSFragger / DIA-NN"]
    MSQC --> Mode{"采集模式"}
    Mode -->|DDA| DDA["DDA 搜库路径"]
    Mode -->|DIA| DIA["DIA 谱图库路径"]

    DDA --> EngineDDA{"搜库引擎"}
    EngineDDA -->|Labelfree| MQ["MaxQuant"]
    EngineDDA -->|开放搜库| FP["FragPipe"]

    DIA --> DIANN["DIANN<br/>谱图定量"]

    MQ --> Quant{"定量策略"}
    FP --> Quant
    DIANN --> Quant

    Quant -->|LFQ| LFQ["LFQ 蛋白矩阵"]
    Quant -->|TMT| TMT["TMT 通道定量"]
    Quant -->|DIA| DIAQuant["DIA 前体离子定量"]

    LFQ --> Impute["KNN / MinProb<br/>缺失值插补"]
    TMT --> Impute
    DIAQuant --> Impute

    Impute --> Diff["差异蛋白分析"]
    Diff --> Volcano["火山图 / 热图"]
    Volcano --> PPI["STRING PPI<br/>互作网络"]
    Volcano --> EnrichP["GO / KEGG 富集"]

    class MSQC,Impute,Diff,Volcano core;
    class MQ,FP,DIANN,PPI,EnrichP algo;
    class Mode,EngineDDA,Quant branchNode;
""")

_MERMAID_EPIGENOMICS = _mermaid_diagram("""
    FASTQ["FASTQ 读段"] --> FastQC["FastQC / fastp<br/>质控修剪"]
    FastQC --> Align["BWA-MEM<br/>比对"]
    Align --> Filter["Duplicate 过滤<br/>Tn5 校正"]

    Filter --> Assay{"实验类型"}
    Assay -->|染色质开放| ATAC["ATAC-seq"]
    Assay -->|组蛋白修饰| ChIP["ChIP-seq"]

    ATAC --> PeakA["MACS2<br/>Peak Calling"]
    ChIP --> PeakC["MACS2 + Input<br/>Peak Calling"]

    PeakA --> IDR["IDR<br/>重现性峰集"]
    PeakC --> IDR

    IDR --> Annot["HOMER / ChIPseeker<br/>Peak 注释"]
    Annot --> Motif["HOMER TF Motif<br/>富集分析"]

    Motif --> DiffPeak["差异开放区域"]
    DiffPeak --> Footprint["足迹分析<br/>顺式调控"]
    DiffPeak --> Integrate["表观转录整合"]

    class FastQC,Align,Filter,IDR,Annot core;
    class PeakA,PeakC,Motif,Footprint,Integrate algo;
    class Assay branchNode;
""")

_MERMAID_METABOLOMICS = _mermaid_diagram("""
    Matrix["代谢物矩阵"] --> Audit["零值审计<br/>稀疏性校验"]
    Audit --> ImputeM["KNN / QRILC<br/>缺失值插补"]
    ImputeM --> Norm["Log2 + Pareto<br/>数据归一化"]

    Norm --> Model{"多元统计建模"}
    Model -->|无监督| PCA["PCA 探索"]
    Model -->|有监督| OPLSDA["OPLS-DA 判别"]

    PCA --> Compare["模型效能对比<br/>R2 / Q2"]
    OPLSDA --> Compare

    Compare --> VIP["VIP 评分筛选"]
    VIP --> BioMarker["标志物表<br/>VIP 大于 1 且 P 小于 0.05"]
    BioMarker --> VolcanoM["火山图 / 热图"]
    VolcanoM --> Pathway["KEGG / HMDB<br/>代谢通路富集"]

    class Audit,ImputeM,Norm,VIP,BioMarker core;
    class PCA,OPLSDA,VolcanoM,Pathway algo;
    class Model branchNode;
""")

_MERMAID_RADIOMICS = _mermaid_diagram("""
    DICOM["CT / MRI + ROI Mask"] --> ROIVal["数据与 ROI<br/>校验"]
    ROIVal --> Preproc["图像预处理<br/>重采样 / 归一化"]
    Preproc --> PyRad["PyRadiomics<br/>特征提取"]

    PyRad --> VarFilter["方差过滤<br/>低方差剔除"]
    VarFilter --> DimRed{"特征降维"}
    DimRed -->|线性投影| PCA_R["PCA 降维"]
    DimRed -->|稀疏筛选| LASSO["LASSO 特征选择"]

    PCA_R --> ModelBuild["训练集划分<br/>分层采样"]
    LASSO --> ModelBuild

    ModelBuild --> Algo{"分类器选择"}
    Algo -->|逻辑回归| LR["Logistic Regression"]
    Algo -->|支持向量机| SVM["SVM"]
    Algo -->|集成学习| RF["Random Forest"]

    LR --> Eval["ROC / AUC 对比"]
    SVM --> Eval
    RF --> Eval

    Eval --> RadScore["Rad-Score<br/>风险概率"]
    Eval --> FeatImp["特征权重<br/>相关性热图"]

    class ROIVal,Preproc,PyRad,VarFilter,ModelBuild core;
    class PCA_R,LASSO,LR,SVM,RF,RadScore,FeatImp algo;
    class DimRed,Algo branchNode;
""")

_MERMAID_GENOMICS = _mermaid_diagram("""
    WES["FASTQ / BAM"] --> QC_G["FastQC / fastp<br/>质控修剪"]
    QC_G --> AlignG["BWA-MEM<br/>比对"]
    AlignG --> Dedup["MarkDuplicates<br/>去重"]
    Dedup --> BQSR["BQSR<br/>碱基质量重校准"]

    BQSR --> Variant{"变异检测"}
    Variant -->|胚系| HC["HaplotypeCaller<br/>SNP / Indel"]
    Variant -->|体细胞| Mutect["Mutect2"]

    HC --> VQSR["VQSR / Filter<br/>质量过滤"]
    Mutect --> VQSR

    BQSR --> CNVPath["CNV / SV<br/>结构变异"]

    VQSR --> AnnotG["VEP / ANNOVAR<br/>功能注释"]
    CNVPath --> AnnotG

    AnnotG --> ACMG["ACMG 致病性<br/>分类"]
    AnnotG --> TMB["TMB 汇总"]
    ACMG --> Report["临床可读报告"]
    TMB --> Report

    class QC_G,AlignG,Dedup,BQSR,VQSR,AnnotG core;
    class HC,Mutect,CNVPath,ACMG,TMB algo;
    class Variant branchNode;
""")

# 广场卡片 Skill.name → 静态 tool_id
SKILL_NAME_TO_PIPELINE_TOOL_ID: Dict[str, str] = {
    "转录组学标准全流程": "pipeline_transcriptomics",
    "空间域识别与聚类": "pipeline_spatial",
    "差异标志物发现": "pipeline_metabolomics",
    "多算法诊断模型构建": "pipeline_radiomics",
    "基因组学全流程": "pipeline_genomics",
    "蛋白组学全流程": "pipeline_proteomics",
    "表观组学全流程": "pipeline_epigenomics",
}


def build_prompt_engineering_guide(
    zero_shot: str,
    dynamic_routing: str,
    post_analysis: str,
) -> List[Dict[str, str]]:
    """高阶 Prompt 教程：三阶调用范式（与前端 skill-detail-prompt-guide 对齐）。"""
    return [
        {
            "tier": "zero_shot",
            "label": "一键化标准执行（Standard Auto-pipeline）",
            "hint": "零参数自动化启动标准 DAG；数据从左侧资产栏注入或由编排器绑定 file_path / h5ad。",
            "prompt": zero_shot.strip(),
        },
        {
            "tier": "dynamic_routing",
            "label": "动态规划与参数精调（Dynamic Routing & Param Tuning）",
            "hint": "可跳过步骤、替换算法模块（如 Harmony 替代 BBKNN），在工作流卡片或对话中覆盖 QC/聚类/批次校正参数。",
            "prompt": dynamic_routing.strip(),
        },
        {
            "tier": "post_analysis",
            "label": "结果追问与深度挖掘（Post-analysis Deep Dive）",
            "hint": "基于已产出的 UMAP/聚类/Marker 表进行亚群下钻、共表达评估与补充可视化。",
            "prompt": post_analysis.strip(),
        },
    ]


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
        "demo_visualization": OMICS_PIPELINE_DEMO_VIZ[tool_id],
    }
    if workflow_highlights:
        out["workflow_highlights"] = workflow_highlights
    return out


# ---------------------------------------------------------------------------
# 1. 转录组学标准全流程（标杆）
# ---------------------------------------------------------------------------
_SCRNA_ZERO = (
    "请按标准单细胞转录组（scRNA-seq）全流程分析我左侧资产中的表达矩阵，"
    "采用默认质控、双细胞去除、降维聚类与细胞注释参数执行即可。"
)
_SCRNA_DYNAMIC = (
    "请修改默认工作流：跳过数据比对环节直接从 h5ad 矩阵开始；在质控步骤，请将线粒体基因保留阈值收紧至 < 5%，"
    "并过滤掉表达基因数少于 200 的细胞；在降维聚类前，强制使用 Harmony 算法对 'batch' 字段进行批次效应矫正，"
    "Leiden 聚类分辨率设为 0.8。"
)
_SCRNA_POST = (
    "基于刚才生成的 UMAP 图，请帮我单独提取出 'CD8+ T cells' 亚群，重新进行细分聚类，"
    "寻找耗竭 T 细胞 (Exhausted T) 的特征 Marker。顺便帮我评估一下 PD-1 (PDCD1) 和 CTLA4 "
    "在这些亚群里的共表达分布，输出小提琴图。"
)

_SPEC_TRANSCRIPTOMICS = _spec(
    tool_id="pipeline_transcriptomics",
    description_long=(
        "端到端的工业级单细胞转录组（scRNA-seq）分析引擎。无缝整合 Cell Ranger 与 Scanpy/Seurat 底层生态。"
        "支持从原始 FASTQ 测序数据的自动比对定量，到严苛的质控过滤、双细胞去除（Scrublet）、"
        "多算法批次效应校正（Harmony/BBKNN/scVI）。内置智能大模型驱动的细胞类型自动注释"
        "（CellTypist 融合 Marker 知识库），并一键延伸至差异表达分析（DEG）、细胞通讯（CellChat）"
        "及发育轨迹推断（Monocle3/PAGA）。\n\n"
        + _workflow_topology_section("### 🧬 标准分析工作流拓扑图", _MERMAID_TRANSCRIPTOMICS)
    ),
    usage_hint=(
        "旗舰管线：由 SOP Planner 动态编排步骤卡片。可在对话中「魔改」DAG（跳过比对、替换 Harmony、调整 resolution），"
        "勿在助手回复中手写 JSON 工作流块。"
    ),
    inputs=[
        {
            "name": "表达数据",
            "type": "10x 目录 / .h5ad / .rds / .mtx+.tsv",
            "required": True,
            "description": "10x Genomics 标准目录、AnnData h5ad、Seurat rds 或稀疏矩阵三件套。",
        },
        {
            "name": "元数据",
            "type": "obs 列 / 样本表",
            "required": False,
            "description": "batch、sample、condition 等批次与分组字段，供 Harmony/差异分析使用。",
        },
    ],
    outputs=[
        {
            "name": "qc_report",
            "type": "交互式 HTML",
            "description": "线粒体/核糖体/基因数分布、Scrublet 双细胞评分与过滤前后细胞数。",
        },
        {
            "name": "embedding_plots",
            "type": "PNG / HTML",
            "description": "多分辨率 UMAP/t-SNE、批次校正前后对比、Leiden 聚类着色。",
        },
        {
            "name": "marker_and_deg",
            "type": "TSV + DotPlot",
            "description": "簇特异性 Marker、Wilcoxon/t-test DEG 表、气泡图与火山图。",
        },
        {
            "name": "cellchat_enrichment",
            "type": "网络图 + 表格",
            "description": "CellChat 细胞通讯拓扑、GO/KEGG 富集结果与专家解读 Markdown。",
        },
    ],
    parameters_table=[
        {"name": "file_path", "type": "path", "required": False, "description": "FASTQ 压缩包、10x 目录或 h5ad 路径。"},
        {"name": "min_genes", "type": "int", "required": False, "description": "每细胞最少检测基因数（默认 200，可对话收紧）。"},
        {"name": "max_pct_mt", "type": "float", "required": False, "description": "线粒体基因比例上限（如 5% 即 0.05）。"},
        {"name": "batch_key", "type": "string", "required": False, "description": "obs 中批次列名，供 Harmony/BBKNN/scVI。"},
        {"name": "integration_method", "type": "enum", "required": False, "description": "harmony | bbknn | scvi — 可对话指定替换默认。"},
        {"name": "resolution", "type": "float", "required": False, "description": "Leiden 聚类分辨率（如 0.8）。"},
        {"name": "doublet_method", "type": "string", "required": False, "description": "Scrublet 双细胞检测与过滤阈值。"},
    ],
    query_examples=build_prompt_engineering_guide(_SCRNA_ZERO, _SCRNA_DYNAMIC, _SCRNA_POST),
    workflow_highlights=[
        "Cell Ranger 定量 → Scanpy QC",
        "Scrublet 双细胞",
        "Harmony / BBKNN / scVI",
        "CellTypist + Marker",
        "DEG · CellChat · Monocle3/PAGA",
    ],
)

# ---------------------------------------------------------------------------
# 2. 空间组学
# ---------------------------------------------------------------------------
_SPEC_SPATIAL = _spec(
    tool_id="pipeline_spatial",
    description_long=(
        "空间转录组工业级分析引擎，原生支持 10x Visium、Stereo-seq 及带 H&E 底图的空间 AnnData。"
        "整合 Squidpy/Scanpy 空间图构建、Moran's I 空间自相关高变基因（SVG）筛选、"
        "多分辨率空间域 Leiden 聚类与组织切片物理映射。支持 Cell2location/RCTD 空间解卷积，"
        "将单细胞参考图谱映射至 spot 水平，识别肿瘤核心、浸润边缘与基质微环境，"
        "并延伸空间差异表达与通路富集。\n\n"
        + _workflow_topology_section("### 🗺️ 空间域识别工作流拓扑图", _MERMAID_SPATIAL)
    ),
    usage_hint="上传含 uns['spatial'] 或坐标矩阵的 h5ad；可要求单细胞参考 + 空间切片联合解卷积。H&E 缺失时自动降级为坐标散点图。",
    inputs=[
        {"name": "空间表达矩阵", "type": ".h5ad / Visium 目录", "required": True, "description": "含空间坐标与可选 H&E 图像的 AnnData。"},
        {"name": "单细胞参考（可选）", "type": ".h5ad", "required": False, "description": "用于 Cell2location/RCTD 解卷积的参考图谱。"},
    ],
    outputs=[
        {"name": "spatial_domain_maps", "type": "PNG / HTML", "description": "多分辨率空间域 1×3 对比、UMAP 与 H&E 双屏联动。"},
        {"name": "svg_table", "type": "TSV", "description": "Moran's I 排序的空间高变基因列表。"},
        {"name": "deconvolution_matrix", "type": "TSV + 热图", "description": "spot × cell type 丰度矩阵与组织微环境注释。"},
        {"name": "enrichment_report", "type": "Markdown", "description": "空间域通路富集与生物学解读。"},
    ],
    parameters_table=[
        {"name": "file_path", "type": "path", "required": True, "description": "Visium 输出目录或空间 h5ad。"},
        {"name": "n_neighbors", "type": "int", "required": False, "description": "空间邻域图邻居数。"},
        {"name": "resolution", "type": "float", "required": False, "description": "空间域 Leiden 分辨率（0.3/0.5/0.8 多分辨率对比）。"},
        {"name": "deconv_method", "type": "enum", "required": False, "description": "cell2location | rctd | 无（仅聚类）。"},
    ],
    query_examples=build_prompt_engineering_guide(
        "请对我左侧资产中的 Visium 空间转录组数据执行标准空间域识别与聚类管线，采用默认参数即可。",
        "请跳过原始比对，直接从空间 h5ad 开始；Moran's I 筛选 Top 300 SVG；空间域聚类分辨率用 0.5 和 0.8 做对比；"
        "对 'batch' 切片批次做 Harmony 整合后再聚类。",
        "请帮我把单细胞数据映射到这张空间切片上做空间解卷积，参考图谱用同一项目的 scRNA h5ad；"
        "输出肿瘤核心区 vs 浸润边缘的 cell type 比例条形图，并对 SPP1+ 基质区做空间热图。",
    ),
    workflow_highlights=["Visium/Stereo-seq", "Moran's I SVG", "空间域聚类", "Cell2location/RCTD", "微环境识别"],
)

# ---------------------------------------------------------------------------
# 3–7. 其余旗舰管线
# ---------------------------------------------------------------------------
_SPEC_PROTEOMICS = _spec(
    tool_id="pipeline_proteomics",
    description_long=(
        "工业级蛋白质组学 DDA/DIA 分析引擎，兼容 MaxQuant、DIANN、FragPipe 搜库定量生态。"
        "覆盖 RAW/mzML 质控、谱图预处理、数据库搜库、FDR 重打分、蛋白推断与 LFQ/TMT 定量、"
        "缺失值 KNN/MinProb 插补、批次校正与差异蛋白分析。输出火山图、热图、"
        "GO/KEGG 通路富集及 STRING PPI 蛋白互作网络，支撑生物标志物发现与机制假说。\n\n"
        + _workflow_topology_section("### 🧪 蛋白组定量工作流拓扑图", _MERMAID_PROTEOMICS)
    ),
    usage_hint="首行含 [Omics_Route: proteomics]；默认 Human UniProt + Trypsin + LFQ，可在对话中切换 DIA/TMT。",
    inputs=[
        {"name": "质谱原始数据", "type": "RAW / mzML", "required": True, "description": "DDA 或 DIA 下机文件，可多样本批量。"},
        {"name": "实验设计", "type": "分组表", "required": True, "description": "样本 ID、condition、batch 列。"},
    ],
    outputs=[
        {"name": "protein_matrix", "type": "TSV", "description": "蛋白 × 样本定量矩阵（插补后）。"},
        {"name": "volcano_heatmap", "type": "PNG", "description": "差异蛋白火山图与聚类热图。"},
        {"name": "ppi_network", "type": "HTML/TSV", "description": "STRING PPI 拓扑与 hub 蛋白列表。"},
        {"name": "enrichment_tables", "type": "TSV", "description": "GO/KEGG 富集结果。"},
    ],
    parameters_table=[
        {"name": "file_path", "type": "path", "required": True, "description": "RAW/mzML 或搜库结果目录。"},
        {"name": "search_engine", "type": "enum", "required": False, "description": "maxquant | diann | fragpipe。"},
        {"name": "quant_method", "type": "enum", "required": False, "description": "LFQ | TMT | DIA。"},
        {"name": "imputation", "type": "string", "required": False, "description": "缺失值插补策略（KNN 等）。"},
        {"name": "fdr_threshold", "type": "float", "required": False, "description": "蛋白/肽段 FDR 阈值。"},
    ],
    query_examples=build_prompt_engineering_guide(
        "请启动蛋白质组学标准全流程：对左侧 RAW/mzML 资产完成搜库、LFQ 定量与差异蛋白分析。",
        "请使用 DIANN 替代 MaxQuant 搜库；缺失值用 KNN 插补；差异分析对比 Tumor vs Normal，FDR<0.01；跳过 PPI 步骤。",
        "基于火山图，请列出 Top 20 上调膜蛋白并绘制 STRING PPI 子网络，解读与免疫检查点相关的 hub。",
    ),
    workflow_highlights=["DIA/DDA", "MaxQuant/DIANN", "缺失值插补", "差异火山图", "PPI 网络"],
)

_SPEC_EPIGENOMICS = _spec(
    tool_id="pipeline_epigenomics",
    description_long=(
        "表观遗传组学 ATAC-seq / ChIP-seq 工业引擎，遵循 ENCODE 与 GATK 系预处理最佳实践。"
        "覆盖 FASTQ 质控、BWA 比对、线粒体/PCR  duplicate 过滤、Tn5 移位校正、MACS2 Peak Calling、"
        "IDR 重现性峰集、Peak 注释（HOMER/ChIPseeker）、差异开放区域、HOMER TF motif 富集、"
        "足迹分析与顺式调控推断；组蛋白修饰（H3K27ac/H3K4me3）与 ATAC 可分支编排。\n\n"
        + _workflow_topology_section("### 🧬 表观组 Peak 工作流拓扑图", _MERMAID_EPIGENOMICS)
    ),
    usage_hint="首行含 [Omics_Route: epigenomics]；默认 hg38；末步可与 RNA 矩阵做多组学整合（data_path 占位）。",
    inputs=[
        {"name": "测序数据", "type": "FASTQ / BAM", "required": True, "description": "ATAC 或 ChIP-seq 双端/单端读段。"},
        {"name": "对照样本", "type": "IgG / Input", "required": False, "description": "ChIP-seq 背景对照。"},
    ],
    outputs=[
        {"name": "peak_bed", "type": "BED", "description": "MACS2 narrowPeak / broadPeak 与 IDR 共有峰。"},
        {"name": "motif_enrichment", "type": "HTML/TSV", "description": "TF motif 富集与靶基因预测。"},
        {"name": "diff_open_regions", "type": "TSV + 热图", "description": "条件间差异 Peak 与注释。"},
        {"name": "integration_report", "type": "Markdown", "description": "表观–转录联合解读（可选）。"},
    ],
    parameters_table=[
        {"name": "file_path", "type": "path", "required": True, "description": "FASTQ 或 BAM 路径。"},
        {"name": "assay_type", "type": "enum", "required": False, "description": "ATAC | H3K27ac | H3K4me3 | 其他 ChIP。"},
        {"name": "macs_qvalue", "type": "float", "required": False, "description": "MACS2 q-value 阈值。"},
        {"name": "genome_build", "type": "string", "required": False, "description": "参考基因组（默认 hg38）。"},
    ],
    query_examples=build_prompt_engineering_guide(
        "请对我左侧资产中的 ATAC-seq FASTQ 执行表观遗传组学标准全流程（质控至 Peak/Motif）。",
        "请按 H3K27ac ChIP-seq 分支跑 Peak；MACS2 qvalue=0.01；跳过多组学整合；比对后仅保留 MAPQ>30。",
        "基于差异开放区域，请做 NFKB1 motif 富集并列出 Top 靶基因启动子，结合通路解读炎症应答。",
    ),
    workflow_highlights=["ATAC/ChIP", "MACS2 Peak", "IDR", "TF Motif", "组蛋白修饰分支"],
)

_SPEC_METABOLOMICS = _spec(
    tool_id="pipeline_metabolomics",
    description_long=(
        "非靶向代谢组学生物标志物发现引擎。支持高维代谢物矩阵的零值审计、KNN/QRILC 缺失插补、"
        "Log2 转换与 Pareto/Auto-scaling 归一化。内置 PCA 与 OPLS-DA 有监督模型、"
        "置换检验与 VIP 评分筛选（VIP>1 & P<0.05），输出火山图、VIP 棒棒糖图、Clustermap，"
        "并对接 KEGG/HMDB 代谢通路富集，服务临床诊断与机制研究。\n\n"
        + _workflow_topology_section("### 🔬 代谢标志物发现工作流拓扑图", _MERMAID_METABOLOMICS)
    ),
    usage_hint="须明确 Disease vs Control 等两组以上标签；PLS-DA 前系统自动校验分组有效性。",
    inputs=[
        {"name": "代谢物矩阵", "type": "CSV / TSV", "required": True, "description": "样本 × 代谢物丰度表。"},
        {"name": "样本分组", "type": "Metadata", "required": True, "description": "condition、batch 等列。"},
    ],
    outputs=[
        {"name": "model_comparison", "type": "PNG", "description": "PCA vs OPLS-DA 1×3 模型效能对比与 R2/Q2。"},
        {"name": "biomarker_table", "type": "TSV", "description": "VIP、P-value、FC 联合筛选标志物。"},
        {"name": "pathway_enrichment", "type": "TSV", "description": "KEGG 代谢通路富集。"},
    ],
    parameters_table=[
        {"name": "file_path", "type": "path", "required": True, "description": "代谢物矩阵路径。"},
        {"name": "group_column", "type": "string", "required": True, "description": "分组列名。"},
        {"name": "vip_threshold", "type": "float", "required": False, "description": "VIP  cutoff（默认 1.0）。"},
        {"name": "p_threshold", "type": "float", "required": False, "description": "显著性阈值。"},
    ],
    query_examples=build_prompt_engineering_guide(
        "请对左侧非靶向代谢组矩阵执行差异标志物发现全流程，分组为 Tumor vs Control。",
        "请用 OPLS-DA 替代 PLS-DA；VIP>1.5 且 P<0.01；归一化改为 Pareto scaling；跳过通路富集。",
        "请解读 Top 5 VIP 代谢物的 KEGG 通路，并绘制其在两组间的箱线图与相关性热图。",
    ),
    workflow_highlights=["OPLS-DA", "VIP 筛选", "缺失插补", "代谢通路富集"],
)

_SPEC_RADIOMICS = _spec(
    tool_id="pipeline_radiomics",
    description_long=(
        "医学影像组学（Radiomics）诊断建模引擎。基于 PyRadiomics 从 CT/MRI ROI mask 提取"
        "一阶统计、形状与纹理特征（GLCM/GLRLM 等），经方差过滤与 LASSO 降维筛选核心特征，"
        "并行训练 Logistic Regression、SVM、Random Forest 分类器，输出多算法 ROC/AUC 对比、"
        "校准曲线、Rad-Score 与特征权重/相关性热图，满足临床可解释性与严格 train/test 隔离。\n\n"
        + _workflow_topology_section("### 🏥 影像组学诊断建模工作流拓扑图", _MERMAID_RADIOMICS)
    ),
    usage_hint="上传影像组学特征 CSV 并指明 label 列；禁止在对话中泄露测试集标签到训练步骤。",
    inputs=[
        {"name": "特征矩阵", "type": "CSV", "required": True, "description": "PyRadiomics 或自定义 radiomics 特征表。"},
        {"name": "临床标签", "type": "列名", "required": True, "description": "二分类诊断 label（如 Malignant/Benign）。"},
    ],
    outputs=[
        {"name": "lasso_path", "type": "PNG", "description": "LASSO 系数路径与入选特征。"},
        {"name": "roc_comparison", "type": "PNG", "description": "LR/SVM/RF 多算法 ROC 与 AUC 标注。"},
        {"name": "feature_importance", "type": "PNG + TSV", "description": "权重棒棒糖图与相关热图。"},
    ],
    parameters_table=[
        {"name": "file_path", "type": "path", "required": True, "description": "特征 CSV 路径。"},
        {"name": "label_column", "type": "string", "required": True, "description": "分类标签列。"},
        {"name": "test_size", "type": "float", "required": False, "description": "测试集比例（默认 0.2）。"},
        {"name": "random_state", "type": "int", "required": False, "description": "随机种子，保证可重复。"},
    ],
    query_examples=build_prompt_engineering_guide(
        "请基于左侧 PyRadiomics 特征矩阵构建多算法诊断模型，分类标签列名为 label。",
        "请加大 LASSO 惩罚（更小 alpha）、仅用 Random Forest 与 SVM、test_size=0.3、输出 SHAP 风格特征排名。",
        "基于 ROC 图，请给出 Rad-Score 切点与灵敏度/特异度，并解释 Top 3 纹理特征的临床含义。",
    ),
    workflow_highlights=["PyRadiomics", "LASSO", "LR/SVM/RF", "ROC/AUC"],
)

_SPEC_GENOMICS = _spec(
    tool_id="pipeline_genomics",
    description_long=(
        "基因组学 WGS/WES 胚系/肿瘤变异分析引擎，对齐 GATK Best Practices。"
        "主线覆盖 FastQC 质控、fastp 修剪、BWA-MEM 比对、MarkDuplicates、BQSR、"
        "HaplotypeCaller 胚系 SNP/Indel、CNV（ExomeDepth 等）、结构变异、VQSR、"
        "VEP/ANNOVAR 注释、ACMG 致病性分类与 TMB 汇总，输出临床可读报告与 IGV 友好 VCF。\n\n"
        + _workflow_topology_section("### 🧬 基因组变异分析工作流拓扑图", _MERMAID_GENOMICS)
    ),
    usage_hint="首行含 [Omics_Route: genomics]；默认 GRCh38；助手勿编造变异坐标，由执行器产出真实 VCF。",
    inputs=[
        {"name": "测序数据", "type": "FASTQ / BAM / VCF", "required": True, "description": "WGS/WES 原始或中间文件。"},
        {"name": "参考基因组", "type": "FASTA 索引", "required": False, "description": "默认 hg38，可由环境变量注入。"},
    ],
    outputs=[
        {"name": "vcf_annotated", "type": "VCF", "description": "VQSR 后注释变异集。"},
        {"name": "cnv_sv_summary", "type": "TSV", "description": "CNV/SV 计数与区域摘要。"},
        {"name": "tmb_acmg_report", "type": "Markdown/HTML", "description": "TMB、ACMG 分类与临床报告。"},
    ],
    parameters_table=[
        {"name": "file_path", "type": "path", "required": True, "description": "FASTQ/BAM/VCF。"},
        {"name": "input_dir", "type": "path", "required": False, "description": "多样本 FASTQ 目录。"},
        {"name": "min_read_length", "type": "int", "required": False, "description": "fastp 最短读长过滤。"},
        {"name": "report_style", "type": "string", "required": False, "description": "临床报告版式。"},
    ],
    query_examples=build_prompt_engineering_guide(
        "请对左侧 WES FASTQ 执行基因组学标准全流程，产出胚系变异注释与临床可读报告。",
        "已有 BAM，请跳过 trim/align，从 BQSR 开始；胚系变异用 GATK HaplotypeCaller；补充 CNV 与 TMB 汇总。",
        "请解读报告中的 BRAF V600E 与 TP53 变异 ACMG 分类，并列出同源重组修复相关基因突变。",
    ),
    workflow_highlights=["WGS/WES", "GATK", "SNP/Indel", "CNV/SV", "TMB", "ACMG"],
)

OMICS_PIPELINE_SPECS_BY_TOOL_ID: Dict[str, Dict[str, Any]] = {
    "pipeline_transcriptomics": _SPEC_TRANSCRIPTOMICS,
    "pipeline_spatial": _SPEC_SPATIAL,
    "pipeline_proteomics": _SPEC_PROTEOMICS,
    "pipeline_epigenomics": _SPEC_EPIGENOMICS,
    "pipeline_metabolomics": _SPEC_METABOLOMICS,
    "pipeline_radiomics": _SPEC_RADIOMICS,
    "pipeline_genomics": _SPEC_GENOMICS,
}
