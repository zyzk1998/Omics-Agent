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
    "通过SMILES获取CID": """[Skill_Route: pubchem_smiles_to_cid]
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
]
