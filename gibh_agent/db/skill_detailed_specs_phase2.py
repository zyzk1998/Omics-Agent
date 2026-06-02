# -*- coding: utf-8 -*-
"""技能详情 Phase 2：第二批 10 项 detailed_spec。"""
from __future__ import annotations

from gibh_agent.db.skill_detail_demo_visualizations import (
    DEMO_VIZ_BEPIPRED,
    DEMO_VIZ_BLAST_HITS,
    DEMO_VIZ_COVER_LETTER,
    DEMO_VIZ_MOL_CONVERT,
    DEMO_VIZ_PUBMED,
    DEMO_VIZ_RESUME,
    DEMO_VIZ_SEQ_CONVERT,
    DEMO_VIZ_SMILES_CID,
    DEMO_VIZ_SUBSTRUCTURE,
    DEMO_VIZ_UNIPROT,
)

SKILL_DETAILED_SPECS_PHASE2: dict = {
    "pubmed_query": {
        "tool_id": "pubmed_query",
        "description_long": (
            "在 **PubMed** 检索生物医学文献：支持关键词、MeSH 风格短语与组合检索，"
            "返回标题、作者、期刊、PMID 与摘要片段的 Markdown 表格，便于快速精读与引用整理。"
        ),
        "usage_hint": "需 launch-skills 微服务联网；未改写检索词时使用 Demo JSON。",
        "inputs": [{"name": "检索式", "type": "string", "required": True, "description": "PubMed 检索关键词或短语。"}],
        "outputs": [
            {"name": "markdown", "type": "Markdown + 表格", "description": "文献命中列表与摘要摘录。"},
        ],
        "parameters_table": [
            {"name": "query", "type": "string", "required": True, "description": "检索关键词。"},
            {"name": "limit", "type": "integer", "required": False, "description": "返回条数，默认 10。"},
        ],
        "query_examples": [
            (
                "您好。我需要在 PubMed 检索相关文献。\n\n"
                "```json\n"
                '{"query": "single cell RNA-seq immunotherapy", "limit": 8}\n'
                "```"
            ),
        ],
        "demo_visualization": DEMO_VIZ_PUBMED,
    },
    "uniprot_query": {
        "tool_id": "uniprot_query",
        "description_long": (
            "查询 **UniProt** 蛋白条目：获取官方蛋白名、基因名、功能注释、亚细胞定位、"
            "结构域与相关疾病关联，适合靶点调研与序列注释前的信息核对。"
        ),
        "usage_hint": "支持基因符号（如 TP53）或 UniProt accession；需联网。",
        "inputs": [{"name": "检索词", "type": "string", "required": True, "description": "蛋白名、基因名或 accession。"}],
        "outputs": [{"name": "markdown", "type": "Markdown + 表格", "description": "蛋白注释摘要与关键字段。"}],
        "parameters_table": [
            {"name": "query", "type": "string", "required": True, "description": "检索关键词。"},
            {"name": "limit", "type": "integer", "required": False, "description": "返回条数，默认 5。"},
        ],
        "query_examples": [
            (
                "您好。我需要查询 UniProt 蛋白注释。\n\n"
                "```json\n"
                '{"query": "TP53", "limit": 5}\n'
                "```"
            ),
        ],
        "demo_visualization": DEMO_VIZ_UNIPROT,
    },
    "bepipred3_prediction": {
        "tool_id": "bepipred3_prediction",
        "description_long": (
            "基于蛋白语言模型的 **B 细胞表位预测（BepiPred-3.0）**："
            "对输入蛋白序列逐残基打分，识别潜在线性表位区段，"
            "并可在工作台渲染表位分数柱状图与阈值筛选结果。"
        ),
        "usage_hint": "输入 FASTA 或纯序列文本；远程计算可能需数分钟。",
        "inputs": [{"name": "蛋白序列", "type": "FASTA / 文本", "required": True, "description": "单条或多条蛋白序列。"}],
        "outputs": [
            {"name": "markdown", "type": "Markdown", "description": "预测摘要与高置信区段。"},
            {"name": "chart_html", "type": "HTML", "description": "表位分数可视化（若启用）。"},
        ],
        "parameters_table": [
            {"name": "sequence_text", "type": "string", "required": False, "description": "序列文本或 FASTA。"},
            {"name": "sequence_or_path", "type": "string", "required": False, "description": "FASTA 文件路径（与 sequence_text 二选一）。"},
        ],
        "query_examples": [
            (
                "您好。请对下列蛋白序列执行 BepiPred-3.0 表位预测。\n\n"
                "```json\n"
                '{"sequence_text": ">7lj4_B\\nMKTIIALSYIFCLVFA\\n", "threshold": "top20"}\n'
                "```"
            ),
        ],
        "demo_visualization": DEMO_VIZ_BEPIPRED,
    },
    "seq_format_converter": {
        "tool_id": "seq_format_converter",
        "description_long": (
            "核酸序列 **格式互转**：支持 FASTA ↔ GenBank 等常见格式，"
            "自动补全 record id、描述行与碱基分组，便于下游提交 NCBI 或流水线输入。"
        ),
        "usage_hint": "粘贴 FASTA 或使用 `sequence_or_path` 指向本地文件。",
        "inputs": [{"name": "序列", "type": "FASTA / GenBank 文本", "required": True, "description": "待转换的序列内容。"}],
        "outputs": [{"name": "markdown", "type": "Markdown + code", "description": "转换后序列与格式说明。"}],
        "parameters_table": [
            {"name": "sequence_text", "type": "string", "required": False, "description": "序列文本。"},
            {"name": "sequence_or_path", "type": "string", "required": False, "description": "序列文件路径。"},
            {"name": "input_format", "type": "string", "required": True, "description": "fasta / genbank。"},
            {"name": "output_format", "type": "string", "required": True, "description": "目标格式。"},
            {"name": "record_id", "type": "string", "required": False, "description": "GenBank 记录 ID。"},
        ],
        "query_examples": [
            (
                "您好。我需要将 FASTA 序列转换为 GenBank 格式。\n\n"
                "```json\n"
                '{"sequence_text": ">demo_gene\\nATGCGTACGTTAGCTAGCTAGCTAGCTAG\\n", '
                '"input_format": "fasta", "output_format": "genbank", "record_id": "demo_gene"}\n'
                "```"
            ),
        ],
        "demo_visualization": DEMO_VIZ_SEQ_CONVERT,
    },
    "rdkit_mol_format_convert": {
        "tool_id": "rdkit_mol_format_convert",
        "description_long": (
            "使用 RDKit 进行**小分子结构格式转换**：SMILES、InChI、Mol 块等格式互转，"
            "适合化学信息学预处理与数据库提交前的格式统一。"
        ),
        "usage_hint": "在 `input_text` 中粘贴 SMILES 或 Mol 块；指定 input/output 格式。",
        "inputs": [{"name": "结构文本", "type": "SMILES / Mol", "required": True, "description": "待转换分子表示。"}],
        "outputs": [{"name": "markdown", "type": "Markdown", "description": "转换结果与验证状态。"}],
        "parameters_table": [
            {"name": "input_text", "type": "string", "required": True, "description": "输入结构字符串。"},
            {"name": "input_format", "type": "string", "required": True, "description": "smiles / inchi / mol 等。"},
            {"name": "output_format", "type": "string", "required": True, "description": "目标格式。"},
        ],
        "query_examples": [
            (
                "您好。我需要将阿司匹林 SMILES 转为 InChI。\n\n"
                "```json\n"
                '{"input_text": "CC(=O)Oc1ccccc1C(=O)O", "input_format": "smiles", "output_format": "inchi"}\n'
                "```"
            ),
        ],
        "demo_visualization": DEMO_VIZ_MOL_CONVERT,
    },
    "nucleotide_sequence_blast": {
        "tool_id": "nucleotide_sequence_blast",
        "description_long": (
            "对核酸序列执行 **BLASTn** 检索（支持远程 NCBI）："
            "返回命中 accession、描述、E-value 与比对区间，用于同源性与物种来源快速判断。"
        ),
        "usage_hint": "远程 blastn 通常需 1–3 分钟；Demo 使用 core_nt + blastn-short 加速。",
        "inputs": [{"name": "核酸序列", "type": "FASTA / 文本", "required": True, "description": "查询序列。"}],
        "outputs": [{"name": "markdown", "type": "Markdown + 表格", "description": "BLAST 命中摘要。"}],
        "parameters_table": [
            {"name": "sequence_text", "type": "string", "required": False, "description": "序列文本。"},
            {"name": "sequence_or_path", "type": "string", "required": False, "description": "FASTA 路径。"},
            {"name": "database", "type": "string", "required": False, "description": "默认 core_nt。"},
            {"name": "max_target_seqs", "type": "integer", "required": False, "description": "命中条数上限。"},
        ],
        "query_examples": [
            (
                "您好。我需要对核酸序列做 blastn 检索。\n\n"
                "```json\n"
                '{"sequence_text": "ATGCGTACGTTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAG", '
                '"use_remote": true, "database": "core_nt", "blast_task": "blastn-short", "max_target_seqs": 5}\n'
                "```"
            ),
        ],
        "demo_visualization": DEMO_VIZ_BLAST_HITS,
    },
    "tailored_resume": {
        "tool_id": "tailored_resume",
        "description_long": (
            "根据 **职位描述（JD）** 与候选人背景生成 ATS 友好的 **Markdown 简历**，"
            "并给出优势亮点与待补差距建议，适合求职投递前的快速定制。"
        ),
        "usage_hint": "JD 放 `user_request`；履历要点放 `context`；纯 Prompt 技能，秒级返回。",
        "inputs": [
            {"name": "职位描述", "type": "自然语言", "required": True, "description": "目标岗位与技能要求。"},
            {"name": "候选人背景", "type": "文本", "required": True, "description": "经历、项目与技能清单。"},
        ],
        "outputs": [{"name": "markdown", "type": "Markdown", "description": "完整简历 + 求职建议。"}],
        "parameters_table": [
            {"name": "user_request", "type": "string", "required": True, "description": "JD 摘要。"},
            {"name": "context", "type": "string", "required": True, "description": "候选人履历。"},
        ],
        "query_examples": [
            (
                "您好。请根据职位描述生成定制化 Markdown 简历。\n\n"
                "```json\n"
                '{"user_request": "Senior Data Analyst；SQL/Python/可视化/A/B 测试；医疗行业优先", '
                '"context": "5 年 RetailCo 数据分析；Tableau/Power BI；HealthPlus 实习"}\n'
                "```"
            ),
        ],
        "demo_visualization": DEMO_VIZ_RESUME,
    },
    "journal_cover_letter": {
        "tool_id": "journal_cover_letter",
        "description_long": (
            "根据论文题目、核心亮点与目标期刊风格，生成 **Cover Letter Markdown 草稿**，"
            "包含编辑称呼、研究意义、与期刊 scope 的契合点及推荐审稿人占位（若用户提供）。"
        ),
        "usage_hint": "期刊与文章类型写 `user_request`；Title/Highlights 写 `context`。",
        "inputs": [
            {"name": "投稿信息", "type": "自然语言", "required": True, "description": "目标期刊、文章类型、语言。"},
            {"name": "论文亮点", "type": "文本", "required": True, "description": "Title、Highlights、创新点。"},
        ],
        "outputs": [{"name": "markdown", "type": "Markdown", "description": "Cover Letter 正文。"}],
        "parameters_table": [
            {"name": "user_request", "type": "string", "required": True, "description": "期刊与体裁。"},
            {"name": "context", "type": "string", "required": True, "description": "题目与亮点。"},
        ],
        "query_examples": [
            (
                "您好。请起草期刊 Cover Letter。\n\n"
                "```json\n"
                '{"user_request": "目标期刊：Nature Communications；Article type: Research Article；英文", '
                '"context": "Title: scRNA atlas under PD-1 blockade\\nHighlights: trajectory; validation; biomarkers"}\n'
                "```"
            ),
        ],
        "demo_visualization": DEMO_VIZ_COVER_LETTER,
    },
    "rdkit_substructure_search": {
        "tool_id": "rdkit_substructure_search",
        "description_long": (
            "使用 RDKit 在候选分子中执行**子结构匹配**："
            "给定 SMARTS/SMILES 子结构，判断目标分子是否包含该片段，"
            "支持单分子或 SMILES 列表批量筛选。"
        ),
        "usage_hint": "苯环 Demo：`c1ccccc1`；需 launch-skills 本地 RDKit。",
        "inputs": [
            {"name": "子结构", "type": "SMILES", "required": True, "description": "查询子结构。"},
            {"name": "目标分子", "type": "SMILES / 列表", "required": True, "description": "待筛选化合物。"},
        ],
        "outputs": [{"name": "markdown", "type": "Markdown + 表格", "description": "匹配结果与命中分子。"}],
        "parameters_table": [
            {"name": "substructure_smiles", "type": "string", "required": True, "description": "子结构 SMILES。"},
            {"name": "target_smiles", "type": "string", "required": False, "description": "单个目标 SMILES。"},
            {"name": "smiles_list", "type": "array", "required": False, "description": "批量 SMILES。"},
        ],
        "query_examples": [
            (
                "您好。我需要在候选分子中做苯环子结构筛选。\n\n"
                "```json\n"
                '{"substructure_smiles": "c1ccccc1", "target_smiles": "c1ccccc1CCO"}\n'
                "```"
            ),
        ],
        "demo_visualization": DEMO_VIZ_SUBSTRUCTURE,
    },
    "smiles_to_cid": {
        "tool_id": "smiles_to_cid",
        "description_long": (
            "将 **SMILES** 解析为 **PubChem Compound ID (CID)**："
            "支持多 CID 歧义时返回候选列表，便于后续 PubChem/ChEMBL 交叉检索。"
        ),
        "usage_hint": "Demo 使用阿司匹林 SMILES；需联网访问 PubChem。",
        "inputs": [{"name": "SMILES", "type": "string", "required": True, "description": "小分子 SMILES。"}],
        "outputs": [{"name": "markdown", "type": "Markdown + 表格", "description": "CID 与化合物名称。"}],
        "parameters_table": [
            {"name": "smiles", "type": "string", "required": True, "description": "输入 SMILES。"},
            {"name": "max_cids", "type": "integer", "required": False, "description": "最多返回 CID 数。"},
        ],
        "query_examples": [
            (
                "您好。我需要将阿司匹林 SMILES 转为 PubChem CID。\n\n"
                "```json\n"
                '{"smiles": "CC(=O)Oc1ccccc1C(=O)O", "max_cids": 10}\n'
                "```"
            ),
        ],
        "demo_visualization": DEMO_VIZ_SMILES_CID,
    },
}
