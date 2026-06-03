# -*- coding: utf-8 -*-
"""
无 [Skill_Route] 占位技能 · 完整 detailed_spec（生物医药 + 化学）。

通过技能广场 **name** → tool_id 映射解析；不占用执行器路由名。
"""
from __future__ import annotations

from typing import Any, Dict, List, Optional

from gibh_agent.db.skill_detail_demo_visualizations_unrouted import UNROUTED_VIZ

# 广场卡片 Skill.name → 静态 tool_id（与 UNROUTED_VIZ / SPECS 键一致）
SKILL_NAME_TO_UNROUTED_TOOL_ID: Dict[str, str] = {
    "Evo2": "evo2",
    "ESM3": "esm3",
    "AlphaFold2": "alphafold2",
    "AlphaFold2-Multimer": "alphafold2_multimer",
    "DiffSBDD": "diffsbdd",
    "OligoFormer": "oligoformer",
    "MSA search": "msa_search",
    "ESMFold": "esmfold",
    "ProtGPT2": "protgpt2",
    "RFdiffusion": "rfdiffusion",
    "BioGPT": "biogpt",
    "BioGraph": "biograph",
    "蛋白质资料提取工具": "protein_info_extractor",
    "蛋白质结构渲染工具": "protein_structure_renderer",
    "RNA二级结构可视化工具": "rna_secondary_structure_viz",
    "抗体人源化": "antibody_humanization",
    "分析细菌生长曲线": "bacterial_growth_curve",
    "LigandMPNN": "ligandmpnn",
    "ESM-Variants": "esm_variants",
    "抗体序列生成": "antibody_sequence_generation",
    "AlphaFold数据库查询": "alphafold_db_query",
    "UCSC 基因组浏览器查询": "ucsc_genome_browser_query",
    "查询GSEA支持的数据库工具": "gsea_supported_databases_query",
    "分子可视化工具": "chem_molecule_viewer",
    "LAMMPS": "lammps_md",
    "CP2K": "cp2k_qc",
}


def _io(kind: str) -> tuple[List[Dict[str, str]], List[Dict[str, str]]]:
    if kind == "db_query":
        return (
            [{"name": "检索词 / ID", "type": "string", "required": True, "description": "基因名、UniProt、坐标或库名。"}],
            [{"name": "markdown", "type": "Markdown + 表格", "description": "命中记录与字段摘要。"}],
        )
    if kind == "viz":
        return (
            [
                {"name": "序列 / 结构 / 表格", "type": "文本或文件", "required": True, "description": "待可视化的输入数据。"},
            ],
            [{"name": "image / html", "type": "PNG / SVG / HTML", "description": "图表或交互式预览。"}],
        )
    if kind == "chem_sim":
        return (
            [{"name": "输入 deck / 拓扑", "type": "文件", "required": True, "description": "模拟输入文件与力场/基组说明。"}],
            [{"name": "log / trajectory", "type": "文件 + 摘要", "description": "轨迹、能量、RDF 等分析结果。"}],
        )
    if kind == "text_ai":
        return (
            [{"name": "生物医学文本", "type": "自然语言", "required": True, "description": "摘要、文献段落或对话上下文。"}],
            [{"name": "markdown", "type": "Markdown", "description": "NER / 关系 / 摘要等结构化输出。"}],
        )
    # predict / analysis default
    return (
        [
            {"name": "序列 / 结构 / 参数", "type": "文本或文件", "required": True, "description": "模型输入与可选超参数。"},
        ],
        [{"name": "report", "type": "Markdown + 文件", "description": "预测结果、评分表与可视化附件。"}],
    )


def _params(tool_id: str, extra: Optional[List[Dict[str, Any]]] = None) -> List[Dict[str, Any]]:
    base = [
        {"name": "user_request", "type": "string", "required": True, "description": "研究目标、物种背景与输出要求。"},
        {"name": "context", "type": "string", "required": False, "description": "序列、结构文件路径或补充材料。"},
    ]
    if extra:
        return base + extra
    return base


def _spec(
    tool_id: str,
    description_long: str,
    usage_hint: str,
    kind: str,
    query_example: str,
    parameters_table: Optional[List[Dict[str, Any]]] = None,
) -> Dict[str, Any]:
    inputs, outputs = _io(kind)
    return {
        "tool_id": tool_id,
        "description_long": description_long,
        "usage_hint": usage_hint,
        "inputs": inputs,
        "outputs": outputs,
        "parameters_table": parameters_table or _params(tool_id),
        "query_examples": [query_example],
        "demo_visualization": UNROUTED_VIZ[tool_id],
    }


_UNROUTED_META: List[Dict[str, Any]] = [
    {
        "tool_id": "evo2",
        "kind": "predict",
        "description_long": "Evo2 生物学基础模型：整合长程基因组上下文，对单核苷酸变化保持敏感性，可用于变异效应与调控元件分析（接入规划中）。",
        "usage_hint": "本项为能力占位说明；正式接入后请提供基因组坐标、参考序列与待评估变异列表。",
        "query_example": "您好。我希望使用 Evo2 评估 chr17 TP53 外显子区 c.524G>A 在长程基因组上下文中的似然变化，并输出碱基分辨率热图摘要。",
    },
    {
        "tool_id": "esm3",
        "kind": "predict",
        "description_long": "ESM3 多模态蛋白基础模型：联合序列、结构与功能提示进行蛋白质设计与进化模拟（接入规划中）。",
        "usage_hint": "请说明设计约束（长度、折叠类型或功能关键词）与是否需要结构条件输入。",
        "query_example": "您好。请使用 ESM3 在激酶折叠骨架约束下生成约 120 aa 的新序列，并报告与模板结构的 RMSD 估计。",
    },
    {
        "tool_id": "alphafold2",
        "kind": "predict",
        "description_long": "AlphaFold2：由氨基酸序列预测单链蛋白质三维结构，输出 PDB/mmCIF 及 per-residue pLDDT 置信度（接入规划中）。",
        "usage_hint": "请提供 FASTA 或单条蛋白序列；长序列或多结构域建议注明分割策略。",
        "query_example": "您好。请对 EGFR 激酶域序列执行 AlphaFold2 结构预测，并汇总 pLDDT 分区与低置信区段。",
    },
    {
        "tool_id": "alphafold2_multimer",
        "kind": "predict",
        "description_long": "AlphaFold2-Multimer：预测蛋白复合物三维组装，提供 ipTM / pTM 等复合物置信度指标（接入规划中）。",
        "usage_hint": "请给出各链序列及化学计量比；界面分析需配合后续对接或实验验证。",
        "query_example": "您好。请预测 EGFR–ERBB2 异源二聚体复合物结构，并报告 ipTM 与界面 pLDDT。",
    },
    {
        "tool_id": "diffsbdd",
        "kind": "predict",
        "description_long": "DiffSBDD：基于扩散模型的结构驱动药物设计，在蛋白结合口袋约束下生成匹配配体候选（接入规划中）。",
        "usage_hint": "需提供受体 PDB 与口袋定义；生成配体须后续对接与 ADMET 评估。",
        "query_example": "您好。请针对激酶 ATP 结合口袋生成 Top-10 候选配体，并给出对接分数排序表。",
    },
    {
        "tool_id": "oligoformer",
        "kind": "predict",
        "description_long": "OligoFormer：面向特定 mRNA 靶标的 siRNA 自动设计与效率推荐（接入规划中）。",
        "usage_hint": "请提供靶 mRNA 序列或 RefSeq ID 及物种信息。",
        "query_example": "您好。请为 TP53 mRNA 3′UTR 设计 3 条 siRNA 候选，并比较双链 ΔG 与 off-target 风险。",
    },
    {
        "tool_id": "msa_search",
        "kind": "analysis",
        "description_long": "MSA search：将查询蛋白序列与数据库比对，生成多序列比对（MSA）供进化与保守性分析（接入规划中）。",
        "usage_hint": "请粘贴 FASTA 或提供序列文件；大型 MSA 可能耗时较长。",
        "query_example": "您好。请对 412 aa 查询序列执行 MSA 检索，并返回命中数、保守块与 CLUSTAL 格式节选。",
    },
    {
        "tool_id": "esmfold",
        "kind": "predict",
        "description_long": "ESMFold：快速单序列蛋白结构预测，适合大规模筛选与结构初探（接入规划中）。",
        "usage_hint": "输入单条氨基酸序列即可；精度略低于 AlphaFold2 但速度更快。",
        "query_example": "您好。请对单域抗体 VH 序列执行 ESMFold 预测并报告 pLDDT 均值与 PDB 输出路径。",
    },
    {
        "tool_id": "protgpt2",
        "kind": "predict",
        "description_long": "ProtGPT2：蛋白质语言模型，可用于从头序列生成与序列续写（接入规划中）。",
        "usage_hint": "提供种子序列片段与期望生成长度；注意生成序列需实验或结构预测验证。",
        "query_example": "您好。请以给定 20 aa 提示续写 100 aa 蛋白序列，并报告 perplexity 与低复杂度区段。",
    },
    {
        "tool_id": "rfdiffusion",
        "kind": "predict",
        "description_long": "RFdiffusion：蛋白质结合剂骨架生成模型，在靶标结构约束下设计可折叠 binder（接入规划中）。",
        "usage_hint": "需提供靶蛋白结构与结合位点；输出骨架需后续序列设计与验证。",
        "query_example": "您好。请针对 PD-L1 β-sheet 界面设计 3 个 binder 骨架并比较界面残基数与可折叠性评分。",
    },
    {
        "tool_id": "biogpt",
        "kind": "text_ai",
        "description_long": "BioGPT：面向生物医学文本的预训练语言模型，支持 NER、关系抽取、摘要与对话（接入规划中）。",
        "usage_hint": "请粘贴英文生物医学段落并说明任务类型（实体识别 / 关系 / 摘要）。",
        "query_example": "您好。请对下列摘要执行药物–靶点关系抽取，并以表格列出实体对与置信度。\n\n…aspirin inhibits COX-1…",
    },
    {
        "tool_id": "biograph",
        "kind": "viz",
        "description_long": "BioGraph：基因表达与生物信息学专用可视化，可生成热图、PCA、火山图等（接入规划中）。",
        "usage_hint": "请上传表达矩阵或差异分析结果表，并指定图表类型。",
        "query_example": "您好。请基于 48×12 基因表达矩阵绘制样本聚类热图与 PCA 双标图。",
    },
    {
        "tool_id": "protein_info_extractor",
        "kind": "db_query",
        "description_long": "蛋白质资料提取：根据 UniProt/蛋白质 ID 汇总功能注释、亚细胞定位、组织表达与疾病关联（接入规划中）。",
        "usage_hint": "请提供 UniProt 登录号或基因名；多 ID 可批量提交。",
        "query_example": "您好。请提取 P04637 (TP53) 的功能注释、亚细胞定位、组织表达与疾病关联摘要表。",
    },
    {
        "tool_id": "protein_structure_renderer",
        "kind": "viz",
        "description_long": "蛋白质结构渲染：自 PDB/mmCIF 导入坐标，生成 Cartoon / Surface / 配体相互作用等视图（接入规划中）。",
        "usage_hint": "请上传 PDB 或提供 PDB ID；可指定链 ID 与渲染模式。",
        "query_example": "您好。请渲染 PDB 1M17 的 Cartoon 与配体 sticks 视图，并导出 PNG。",
    },
    {
        "tool_id": "rna_secondary_structure_viz",
        "kind": "viz",
        "description_long": "RNA 二级结构可视化：由碱基序列生成点括号表示与弧线图 / 径向图（接入规划中）。",
        "usage_hint": "请提供 RNA 序列（A/U/G/C）；可与 RNAfold 结果联动。",
        "query_example": "您好。请为序列 GGGGAUAGGUUCAACCUCCUU 绘制二级结构弧线图并标注 MFE。",
    },
    {
        "tool_id": "antibody_humanization",
        "kind": "predict",
        "description_long": "抗体人源化：基于 Sapiens 天然抗体库或 CDR 移植深度学习流程，降低免疫原性（接入规划中）。",
        "usage_hint": "请提供鼠源或嵌合抗体 VH/VL 序列；输出需结合人源性评分复核。",
        "query_example": "您好。请对给定鼠源抗体 VH/VL 执行 CDR 移植人源化，并列出回复突变建议。",
    },
    {
        "tool_id": "bacterial_growth_curve",
        "kind": "analysis",
        "description_long": "细菌生长曲线分析：基于 OD600 时间序列拟合 Logistic / Gompertz 等模型，估计 μmax、延迟期与倍增时间（接入规划中）。",
        "usage_hint": "请上传含 time 与 OD600 列的 CSV；注明菌株与培养条件。",
        "query_example": "您好。请拟合附件 OD600 生长曲线，报告 μmax、延迟期 λ 与倍增时间 td 及 95% CI。",
    },
    {
        "tool_id": "ligandmpnn",
        "kind": "predict",
        "description_long": "LigandMPNN：在配体、核酸等非蛋白环境约束下设计蛋白序列（接入规划中）。",
        "usage_hint": "需提供蛋白–配体复合物 PDB 与固定配体原子；输出序列需结构预测验证。",
        "query_example": "您好。请在固定配体环境下重设计结合位点 15 个氨基酸，并报告 ΔΔG 变化排序。",
    },
    {
        "tool_id": "esm_variants",
        "kind": "predict",
        "description_long": "ESM-Variants：交互式评估蛋白序列氨基酸突变的似然变化，辅助致病性初筛（接入规划中）。",
        "usage_hint": "请提供野生型序列与待扫描位点；Δ log-likelihood 需结合其他证据解读。",
        "query_example": "您好。请扫描 TP53 蛋白第 156、220、273 位点常见突变并输出 Δ log-likelihood 热图摘要。",
    },
    {
        "tool_id": "antibody_sequence_generation",
        "kind": "predict",
        "description_long": "抗体序列生成：对抗体序列进行定向突变并实时监测疏水性、免疫原性等指标（接入规划中）。",
        "usage_hint": "请提供 WT VH/VL 与突变策略（CDR 或 FR 区）。",
        "query_example": "您好。请在 CDR-H3 区域生成 2 个变体并比较疏水性与免疫原性评分相对 WT 的变化。",
    },
    {
        "tool_id": "alphafold_db_query",
        "kind": "db_query",
        "description_long": "AlphaFold 数据库查询：按 UniProt 检索预测结构条目、置信度统计与 PDB/mmCIF 下载方式（接入规划中）。",
        "usage_hint": "请提供 UniProt 登录号；可批量查询并比较 pLDDT 分布。",
        "query_example": "您好。请查询 UniProt 登录号 P04637（p53）在 AlphaFold 结构数据库中的三维预测模型条目、置信度分区及下载方式。",
    },
    {
        "tool_id": "ucsc_genome_browser_query",
        "kind": "db_query",
        "description_long": "UCSC 基因组浏览器查询：按坐标或基因名检索 hg38/hg19 等装配下的注释轨道（接入规划中）。",
        "usage_hint": "请说明基因组装配版本、染色体坐标或基因符号。",
        "query_example": "您好。请在 hg38 上检索 chr17 TP53 基因座附近的 RefSeq、ClinVar 与 ENCODE 轨道摘要。",
    },
    {
        "tool_id": "gsea_supported_databases_query",
        "kind": "db_query",
        "description_long": "查询 GSEA / GSEApy 支持的基因集数据库名称、物种覆盖与典型用途（接入规划中）。",
        "usage_hint": "可指定物种或分析类型（ORA / GSEA）以筛选推荐库。",
        "query_example": "您好。请列出 GSEApy 常用的人类基因集数据库（MSigDB Hallmark、KEGG、GO BP）及基因集规模。",
    },
    {
        "tool_id": "chem_molecule_viewer",
        "kind": "viz",
        "description_long": "分子可视化工具：基于 RDKit.js 从 SMILES、SDF、Mol 等格式读取结构并生成交互式 2D 图像（接入规划中）。",
        "usage_hint": "请提供 SMILES 或上传分子文件；可指定导出 PNG 分辨率。",
        "query_example": "您好。请读取阿司匹林 SMILES 并生成交互式 2D 结构图，导出 400×400 PNG。",
    },
    {
        "tool_id": "lammps_md",
        "kind": "chem_sim",
        "description_long": "LAMMPS：经典分子动力学引擎，模拟液体、固体或气体粒子系综的平衡与输运性质（接入规划中）。",
        "usage_hint": "需准备 data 文件、力场与 input 脚本；大规模模拟建议 HPC 提交。",
        "query_example": "您好。请对水盒子体系执行 NVT→NPT 2 ns 模拟，并报告密度收敛与 O–H RDF 峰值。",
        "parameters_table": _params(
            "lammps_md",
            [{"name": "file_path", "type": "string", "required": True, "description": "LAMMPS input / data 文件路径。"}],
        ),
    },
    {
        "tool_id": "cp2k_qc",
        "kind": "chem_sim",
        "description_long": "CP2K：量子化学与固体物理模拟包，支持 DFT、几何优化与 ab initio 分子动力学（接入规划中）。",
        "usage_hint": "需提供 CP2K 输入文件与基组/泛函说明；计算量随体系大小急剧增长。",
        "query_example": "您好。请对给定小分子执行 B3LYP/6-31G* 单点能计算，并输出总能量与最大力收敛报告。",
        "parameters_table": _params(
            "cp2k_qc",
            [{"name": "file_path", "type": "string", "required": True, "description": "CP2K 输入 deck 路径。"}],
        ),
    },
]

SKILL_DETAILED_SPECS_UNROUTED: Dict[str, Dict[str, Any]] = {
    m["tool_id"]: _spec(
        tool_id=m["tool_id"],
        description_long=m["description_long"],
        usage_hint=m["usage_hint"],
        kind=m["kind"],
        query_example=m["query_example"],
        parameters_table=m.get("parameters_table"),
    )
    for m in _UNROUTED_META
}

assert set(SKILL_NAME_TO_UNROUTED_TOOL_ID.values()) == set(SKILL_DETAILED_SPECS_UNROUTED.keys())
