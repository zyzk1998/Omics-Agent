# -*- coding: utf-8 -*-
"""
首发 10 项核心技能 prompt_template（与 seed_skills / patch_launch_skills_batch 同源）。

首行 [Skill_Route: <tool_id>] 须与 gibh_agent/skills/skill_*.py 中 skill_id 完全一致。
每模板须含 **广场一键体验 JSON**（与 gibh_agent/skills/launch_skill_demos.py 一致），供 SkillAgent 与广场直接发送。
"""
from __future__ import annotations

LAUNCH_CORE_PROMPTS_BY_SKILL_NAME: dict[str, str] = {
    "ChEMBL药物检索": """[Skill_Route: chembl_drug_search]
您好。我需要在 ChEMBL 中检索药物记录。

**广场一键体验（助手侧：用户未提供其它检索词、未上传附件时，必须将下列 JSON 原样写入工具调用参数，禁止留空 query）**
```json
{"query": "aspirin", "search_type": "drug", "limit": 10}
```

**可复制一句直接发送**
请按上述 JSON 检索 ChEMBL 中包含 aspirin 的药物，最多 10 条。

（助手侧：必填 `query`；用户给出 CHEMBL ID 时可用 `search_type=molecule`。）
""",
    "子结构搜索化合物": """[Skill_Route: rdkit_substructure_search]
您好。我需要在候选分子中做子结构（苯环）筛选。

**广场一键体验（助手侧：未改需求时须写入下列 JSON）**
```json
{"substructure_smiles": "c1ccccc1", "target_smiles": "c1ccccc1CCO"}
```

**可复制一句直接发送**
请用苯环子结构 `c1ccccc1` 筛选分子 `c1ccccc1CCO` 是否匹配。

（助手侧：`substructure_smiles` 与 `target_smiles` 或 `smiles_list` / `file_path` 二选一以上。）
""",
    "核酸序列比对": """[Skill_Route: nucleotide_sequence_blast]
您好。我需要对核酸序列做 blastn 检索。

**广场一键体验（助手侧：未上传 FASTA 时须写入下列 JSON）**
```json
{"sequence_text": "ATGCGTACGTTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAG", "use_remote": true, "database": "core_nt", "blast_task": "blastn-short", "max_target_seqs": 5}
```

**可复制一句直接发送**
请对以上核酸序列执行远程 blastn，返回 5 条命中。

**说明**：远程 NCBI 检索通常需 **1–3 分钟**（已自动使用 `core_nt` + `blastn-short` 加速；全库 `nt` 可能需 5–15 分钟）。

（助手侧：`sequence_text` 或 `sequence_or_path` 二选一。）
""",
    "蛋白质序列比对": """[Skill_Route: protein_sequence_blast]
您好。我需要对蛋白序列做 blastp 检索。

**广场一键体验（助手侧：未上传 FASTA 时须写入下列 JSON）**
```json
{"sequence_text": "MKTIIALSYIFCLVFA", "use_remote": true, "database": "swissprot", "blast_task": "blastp-short", "max_target_seqs": 5}
```

**可复制一句直接发送**
请对序列 MKTIIALSYIFCLVFA 执行远程 blastp。

（助手侧：`sequence_text` 或 `sequence_or_path` 二选一。）
""",
    "通过SMILES获取CID": """[Skill_Route: smiles_to_cid]
您好。我需要将阿司匹林 SMILES 转为 PubChem CID。

**广场一键体验（助手侧：未改 SMILES 时须写入下列 JSON）**
```json
{"smiles": "CC(=O)Oc1ccccc1C(=O)O", "max_cids": 10}
```

**可复制一句直接发送**
请将 SMILES `CC(=O)Oc1ccccc1C(=O)O` 解析为 PubChem CID。

（助手侧：必填 `smiles`。）
""",
    "MHC关联表位检索": """[Skill_Route: mhc_epitope_search]
您好。我需要在 IEDB 检索 HLA-A*02:01 相关 MHC 表位实验。

**广场一键体验（助手侧：用户未提供等位基因或肽段时，必须写入下列 JSON，禁止两项皆空）**
```json
{"mhc_allele": "HLA-A*02:01", "limit": 10}
```

**可复制一句直接发送**
请检索 mhc_allele 为 HLA-A*02:01 的 IEDB 记录，limit=10。

（助手侧：`mhc_allele` 与 `peptide_sequence` 至少一项非空。）
""",
    "ChIPAtlas实验获取工具": """[Skill_Route: chipatlas_experiment_search]
您好。我需要获取 ChIP-Atlas 实验 SRX018625 的元数据。

**广场一键体验（助手侧：未改 accession 时须写入下列 JSON）**
```json
{"expid": "SRX018625", "limit": 10}
```

**可复制一句直接发送**
请获取 expid=SRX018625 的 ChIP-Atlas 实验详情。

（助手侧：`expid` 或 `query` 至少一项；关键词示例 query: H3K4me3 HeLa。）
""",
    "3D分子结构渲染工具": """[Skill_Route: rdkit_3d_mol_render]
您好。我需要由阿司匹林 SMILES 生成 3D MOL 块。

**广场一键体验（助手侧：未改结构时须写入下列 JSON）**
```json
{"smiles": "CC(=O)Oc1ccccc1C(=O)O", "output_format": "mol_block"}
```

**可复制一句直接发送**
请对 `CC(=O)Oc1ccccc1C(=O)O` 生成 3D mol_block。

（助手侧：必填 `smiles`。）
""",
    "检索相似小分子": """[Skill_Route: chembl_similar_molecules]
您好。我需要在 ChEMBL 中检索与阿司匹林结构相似的小分子。

**广场一键体验（助手侧：未改结构时须写入下列 JSON）**
```json
{"smiles": "CC(=O)Oc1ccccc1C(=O)O", "similarity_threshold": 70, "limit": 10}
```

**可复制一句直接发送**
请按 SMILES 阿司匹林在 ChEMBL 做 70% 相似性检索，limit=10。

（助手侧：`smiles` 或 `chembl_id` 二选一。）
""",
    "分子格式转换工具": """[Skill_Route: rdkit_mol_format_convert]
您好。我需要将阿司匹林 SMILES 转为 InChI。

**广场一键体验（助手侧：未改输入时须写入下列 JSON）**
```json
{"input_text": "CC(=O)Oc1ccccc1C(=O)O", "input_format": "smiles", "output_format": "inchi"}
```

**可复制一句直接发送**
请将上述 SMILES 转为 InChI 表示。

（助手侧：必填 `input_text` 及格式字段。）
""",
    "PubMed数据库查询": """[Skill_Route: pubmed_query]
您好。我需要在 PubMed 检索相关文献。

**广场一键体验（助手侧：须将下列 JSON 原样写入工具调用）**
```json
{"query": "single cell RNA-seq immunotherapy", "limit": 8}
```

（助手侧：必填 `query`；可选 `limit` 默认 10。）
""",
    "UniProt数据库查询": """[Skill_Route: uniprot_query]
您好。我需要查询 UniProt 蛋白注释。

**广场一键体验（助手侧：须将下列 JSON 原样写入工具调用）**
```json
{"query": "TP53", "limit": 5}
```

（助手侧：必填 `query`（蛋白名/基因名/accession）。）
""",
    "序列格式转换工具": """[Skill_Route: seq_format_converter]
您好。我需要将 FASTA 序列转换为 GenBank 格式。

**广场一键体验（助手侧：须将下列 JSON 原样写入工具调用）**
```json
{"sequence_text": ">demo_gene\\nATGCGTACGTTAGCTAGCTAGCTAGCTAG\\n", "input_format": "fasta", "output_format": "genbank", "record_id": "demo_gene"}
```

（助手侧：`sequence_text` 或 `sequence_or_path`；`input_format` / `output_format` 为 fasta 或 genbank。）
""",
    "分子量计算工具": """[Skill_Route: calc_molecular_weight]
您好。我需要计算小分子分子量与理化描述符。

**广场一键体验（助手侧：须将下列 JSON 原样写入工具调用）**
```json
{"smiles": "CC(=O)Oc1ccccc1C(=O)O"}
```

（助手侧：必填 `smiles`。）
""",
    "ClinVar 数据库查询": """[Skill_Route: clinvar_query]
您好。我需要在 ClinVar 检索变异临床意义。

**广场一键体验（助手侧：须将下列 JSON 原样写入工具调用）**
```json
{"query": "BRCA1[gene]", "limit": 8}
```

（助手侧：必填 `query`（基因名、rsID 或 HGVS）；可选 `limit`。）
""",
    "dbSNP 数据库查询": """[Skill_Route: dbsnp_query]
您好。我需要在 dbSNP 检索 SNP 位点信息。

**广场一键体验（助手侧：须将下列 JSON 原样写入工具调用）**
```json
{"query": "rs699", "limit": 5}
```

（助手侧：必填 `query`（rsID 或基因名）。）
""",
    "GWAS catalog 数据库查询": """[Skill_Route: gwas_catalog_query]
您好。我需要在 GWAS Catalog 检索性状关联或研究。

**广场一键体验（助手侧：须将下列 JSON 原样写入工具调用）**
```json
{"query": "type 2 diabetes", "limit": 8}
```

（助手侧：必填 `query`（rsID 或疾病/性状名）。）
""",
    "Reactome 数据库查询": """[Skill_Route: reactome_query]
您好。我需要检索 Reactome 人类通路。

**广场一键体验（助手侧：须将下列 JSON 原样写入工具调用）**
```json
{"query": "apoptosis", "limit": 8}
```

（助手侧：必填 `query`（通路关键词）。）
""",
    "InterPro 数据库查询": """[Skill_Route: interpro_query]
您好。我需要在 InterPro 检索蛋白结构域。

**广场一键体验（助手侧：须将下列 JSON 原样写入工具调用）**
```json
{"query": "kinase", "limit": 8}
```

（助手侧：必填 `query`（结构域/IPR/关键词）。）
""",
    "GEO 数据库查询": """[Skill_Route: geo_query]
您好。我需要检索 GEO 表达数据集元数据（不下载矩阵）。

**广场一键体验（助手侧：须将下列 JSON 原样写入工具调用）**
```json
{"query": "breast cancer", "limit": 8}
```

（助手侧：必填 `query`（关键词或 GSE accession）。）
""",
    "基因蛋白信息查询器": """[Skill_Route: gene_protein_info_query]
您好。我需要查询人类基因 Ensembl 注释。

**广场一键体验（助手侧：须将下列 JSON 原样写入工具调用）**
```json
{"query": "TP53"}
```

（助手侧：必填 `query`（人类基因符号）。）
""",
    "获取mRNA序列工具": """[Skill_Route: mrna_sequence_fetch]
您好。我需要获取人类基因 mRNA（cDNA）序列。

**广场一键体验（助手侧：须将下列 JSON 原样写入工具调用）**
```json
{"sequence_or_path": "TP53"}
```

（助手侧：必填 `sequence_or_path`（基因符号）或 `query`。）
""",
    "药物标识符交叉检索": """[Skill_Route: drug_id_crossref]
您好。我需要药物标识符多库交叉引用。

**广场一键体验（助手侧：须将下列 JSON 原样写入工具调用）**
```json
{"query": "aspirin"}
```

（助手侧：`query`（药名）或 `smiles` 至少一项。）
""",
    "代谢物名称检索": """[Skill_Route: refmet_metabolite_search]
您好。我需要在 RefMet 检索代谢物标准化名称。

**广场一键体验（助手侧：须将下列 JSON 原样写入工具调用）**
```json
{"query": "cholesterol", "match_threshold": 0.8}
```

（助手侧：必填 `query`（代谢物名）。）
""",
    "RNAcentral 按编号查询": """[Skill_Route: rnacentral_query]
您好。我需要查询 RNAcentral RNA 条目。

**广场一键体验（助手侧：须将下列 JSON 原样写入工具调用）**
```json
{"query": "URS000065FE5A"}
```

（助手侧：必填 `query`（URS 编号或关键词）。）
""",
    "Ensembl GO术语后代查询": """[Skill_Route: ensembl_go_descendants]
您好。我需要查询 GO 术语的子术语与后代。

**广场一键体验（助手侧：须将下列 JSON 原样写入工具调用）**
```json
{"query": "GO:0006915", "limit": 15}
```

（助手侧：必填 `query`（GO ID）。）
""",
    "OpenTargets_按ChEMBL ID获取父分子和子分子": """[Skill_Route: opentargets_chembl_hierarchy]
您好。我需要 Open Targets 药物分子父子层级。

**广场一键体验（助手侧：须将下列 JSON 原样写入工具调用）**
```json
{"query": "CHEMBL25"}
```

（助手侧：必填 `query`（ChEMBL ID）。）
""",
    "FDA药品标签字段检索": """[Skill_Route: fda_drug_label_search]
您好。我需要检索 FDA 药品标签关键字段。

**广场一键体验（助手侧：须将下列 JSON 原样写入工具调用）**
```json
{"query": "aspirin", "limit": 3}
```

（助手侧：必填 `query`（品牌名/通用名）。）
""",
}

LAUNCH_CORE_SKILL_SEEDS: list[dict[str, str]] = [
    {
        "name": "ChEMBL药物检索",
        "main_category": "化学",
        "sub_category": "信息检索",
        "description": "根据药物名称或 ChEMBL ID 检索 ChEMBL 药物与分子信息（REST API）。",
    },
    {
        "name": "子结构搜索化合物",
        "main_category": "化学",
        "sub_category": "数据分析",
        "description": "使用 RDKit 在目标分子或 SMILES 列表中执行子结构匹配筛选。",
    },
    {
        "name": "核酸序列比对",
        "main_category": "生物医药",
        "sub_category": "数据分析",
        "description": "使用 NCBI BLAST blastn 对核酸序列进行同源性搜索（支持远程 nt 库）。",
    },
    {
        "name": "蛋白质序列比对",
        "main_category": "生物医药",
        "sub_category": "数据分析",
        "description": "使用 NCBI BLAST blastp 对蛋白序列进行同源性搜索（支持远程 nr 库）。",
    },
    {
        "name": "通过SMILES获取CID",
        "main_category": "化学",
        "sub_category": "信息检索",
        "description": "通过 SMILES 查询 PubChem 化合物 CID（PUG REST）。",
    },
    {
        "name": "MHC关联表位检索",
        "main_category": "生物医药",
        "sub_category": "数据分析",
        "description": "基于 IEDB Query API 检索 MHC 相关表位实验记录。",
    },
    {
        "name": "ChIPAtlas实验获取工具",
        "main_category": "生物医药",
        "sub_category": "信息检索",
        "description": "检索 ChIP-Atlas 实验元数据或按 SRX/GSM 获取详情。",
    },
    {
        "name": "检索相似小分子",
        "main_category": "化学",
        "sub_category": "数据分析",
        "description": "通过 ChEMBL 相似性 API 检索结构类似小分子（需 SMILES 或 ChEMBL ID）。",
    },
    {
        "name": "序列格式转换工具",
        "main_category": "生物医药",
        "sub_category": "数据处理",
        "description": "BioPython SeqIO：FASTA / GenBank 等序列格式互转。",
    },
    {
        "name": "ClinVar 数据库查询",
        "main_category": "生物医药",
        "sub_category": "信息检索",
        "description": "ClinVar E-utilities 变异临床意义检索，表格 + Markdown。",
    },
    {
        "name": "dbSNP 数据库查询",
        "main_category": "生物医药",
        "sub_category": "信息检索",
        "description": "dbSNP E-utilities SNP 位点与基因关联检索。",
    },
    {
        "name": "GWAS catalog 数据库查询",
        "main_category": "生物医药",
        "sub_category": "信息检索",
        "description": "GWAS Catalog REST 性状—SNP 关联与研究检索。",
    },
    {
        "name": "Reactome 数据库查询",
        "main_category": "生物医药",
        "sub_category": "信息检索",
        "description": "Reactome ContentService 人类通路检索。",
    },
    {
        "name": "InterPro 数据库查询",
        "main_category": "生物医药",
        "sub_category": "信息检索",
        "description": "InterPro REST 蛋白结构域与家族检索。",
    },
    {
        "name": "GEO 数据库查询",
        "main_category": "生物医药",
        "sub_category": "信息检索",
        "description": "NCBI GEO DataSets 元数据检索（不下载表达矩阵）。",
    },
    {
        "name": "基因蛋白信息查询器",
        "main_category": "生物医药",
        "sub_category": "信息检索",
        "description": "Ensembl REST 人类基因注释与转录本概览。",
    },
    {
        "name": "获取mRNA序列工具",
        "main_category": "生物医药",
        "sub_category": "信息检索",
        "description": "Ensembl REST 按基因符号获取代表性 cDNA 序列（FASTA）。",
    },
    {
        "name": "药物标识符交叉检索",
        "main_category": "化学",
        "sub_category": "信息检索",
        "description": "PubChem / ChEMBL / Registry 药物标识符交叉引用宽表。",
    },
    {
        "name": "代谢物名称检索",
        "main_category": "化学",
        "sub_category": "信息检索",
        "description": "RefMet 代谢物标准化名称与分类检索。",
    },
    {
        "name": "RNAcentral 按编号查询",
        "main_category": "生物医药",
        "sub_category": "信息检索",
        "description": "RNAcentral URS 条目与交叉 ID 检索。",
    },
    {
        "name": "Ensembl GO术语后代查询",
        "main_category": "生物医药",
        "sub_category": "信息检索",
        "description": "GO 术语子术语与 QuickGO 后代 ID 列表。",
    },
    {
        "name": "OpenTargets_按ChEMBL ID获取父分子和子分子",
        "main_category": "化学",
        "sub_category": "信息检索",
        "description": "Open Targets GraphQL 药物分子父子层级。",
    },
    {
        "name": "FDA药品标签字段检索",
        "main_category": "化学",
        "sub_category": "信息检索",
        "description": "openFDA 药品标签适应症/警告/用法字段检索。",
    },
]
