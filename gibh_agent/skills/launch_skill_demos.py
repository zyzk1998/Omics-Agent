# -*- coding: utf-8 -*-
"""
首发核心技能 — 广场「一键体验」默认参数（化学/BLAST 10 项 + Prompt 软技能 10 项）。

- ``LAUNCH_ISOLATED_TOOL_IDS``：仅这些 tool_id 经 HTTP 委托 ``launch-skills``（BLAST/RDKit/ChEMBL）。
- ``LAUNCH_SKILL_DEMO_ARGS``：含 Prompt 软技能（含阵列一 Batch1 五件套）；演示默认值在 api-server 本地执行，**不**进 launch-skills 白名单。

SkillAgent 填参留空时回退；各 skill execute 在必填项为空时二次兜底。
"""
from __future__ import annotations

from typing import Any, Dict, FrozenSet

# 仅隔离栈（services/launch-skills/main.py _LAUNCH_SKILL_CLASSES）内执行的 tool_id
LAUNCH_ISOLATED_TOOL_IDS: FrozenSet[str] = frozenset(
    {
        "chembl_drug_search",
        "chembl_similar_molecules",
        "pubchem_smiles_to_cid",
        "mhc_epitope_search",
        "chipatlas_experiment_search",
        "rdkit_substructure_search",
        "rdkit_3d_mol_render",
        "rdkit_mol_format_convert",
        "nucleotide_sequence_blast",
        "protein_sequence_blast",
        # 阵列二 Batch1 · 轻量 DB / RDKit / BioPython
        "pubmed_query",
        "uniprot_query",
        "smiles_to_cid",
        "seq_format_converter",
        "calc_molecular_weight",
        # 阵列二 Batch2a · P0 DB 检索
        "clinvar_query",
        "dbsnp_query",
        "gwas_catalog_query",
        "reactome_query",
        "interpro_query",
        "geo_query",
        "gene_protein_info_query",
        "mrna_sequence_fetch",
        # 阵列二 Batch2b · P2 ID / 标签 / ontology
        "drug_id_crossref",
        "refmet_metabolite_search",
        "rnacentral_query",
        "ensembl_go_descendants",
        "opentargets_chembl_hierarchy",
        "fda_drug_label_search",
    }
)

# tool_id → 演示用 kwargs（须与 @registry / BaseSkill.execute 形参一致）
LAUNCH_SKILL_DEMO_ARGS: Dict[str, Dict[str, Any]] = {
    "chembl_drug_search": {
        "query": "aspirin",
        "search_type": "drug",
        "limit": 10,
    },
    "rdkit_substructure_search": {
        "substructure_smiles": "c1ccccc1",
        "target_smiles": "c1ccccc1CCO",
    },
    "nucleotide_sequence_blast": {
        "sequence_text": "ATGCGTACGTTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAG",
        "use_remote": True,
        "database": "core_nt",
        "blast_task": "blastn-short",
        "max_target_seqs": 5,
    },
    "protein_sequence_blast": {
        "sequence_text": "MKTIIALSYIFCLVFA",
        "use_remote": True,
        "database": "swissprot",
        "blast_task": "blastp-short",
        "max_target_seqs": 5,
    },
    "pubchem_smiles_to_cid": {
        "smiles": "CC(=O)Oc1ccccc1C(=O)O",
        "max_cids": 10,
    },
    "mhc_epitope_search": {
        "mhc_allele": "HLA-A*02:01",
        "limit": 10,
    },
    "chipatlas_experiment_search": {
        "expid": "SRX018625",
        "limit": 10,
    },
    "rdkit_3d_mol_render": {
        "smiles": "CC(=O)Oc1ccccc1C(=O)O",
        "output_format": "mol_block",
    },
    "chembl_similar_molecules": {
        "smiles": "CC(=O)Oc1ccccc1C(=O)O",
        "similarity_threshold": 70,
        "limit": 10,
    },
    "rdkit_mol_format_convert": {
        "input_text": "CC(=O)Oc1ccccc1C(=O)O",
        "input_format": "smiles",
        "output_format": "inchi",
    },
    "ppt_outline": {
        "user_request": (
            "单细胞时空动力学分析：研究背景、最优传输轨迹推断、"
            "驱动基因与通路富集、核心结论与展望"
        ),
        "context": "",
    },
    "mindmap_gen": {
        "user_request": (
            "单细胞时空动力学：数据校验→时序标准化→轨迹推断→"
            "驱动基因→通路富集→专家解读"
        ),
        "context": "",
    },
    "weekly_report_writer": {
        "user_request": (
            "周期 2026-05-19 至 2026-05-25；按 SKILL 工作流起草周报，"
            "含个人状态同步与团队对外同步"
        ),
        "context": "上一份报告：待用户补充路径；重点关注项目进展与 🟥/🟨 待办继承",
    },
    "tailored_resume": {
        "user_request": "Senior Data Analyst；SQL/Python/可视化/A/B 测试；医疗行业优先",
        "context": "5 年 RetailCo 数据分析；Tableau/Power BI；HealthPlus 实习；Business Analytics 学士",
    },
    "academic_poster_generator": {
        "user_request": "单细胞轨迹推断论文 → 会议海报故事板（四段式 + 至少 3 张配图方案）",
        "context": "摘要：最优传输轨迹推断；结果含 UMAP 与驱动基因（见用户附件）",
    },
    "academic_abstract_refiner": {
        "user_request": "将下列英文草稿精炼为中英双语 SCI 风格单段摘要（无小标题）",
        "context": (
            "Background: Single-cell RNA sequencing enables high-resolution profiling, but batch effects "
            "limit cross-cohort integration. Methods: We developed OT-Integrate, an optimal-transport "
            "framework with graph smoothing for atlas-level alignment. Results: On 12 public datasets, "
            "OT-Integrate improved cluster purity by 15% (ARI +0.12) and preserved rare cell states. "
            "Conclusions: The method supports reproducible cross-study comparison for translational "
            "single-cell analysis. (Demo 演示原文，可直接体验润色效果)"
        ),
    },
    "email_manager": {
        "user_request": "write email to client about project delay — professional, milestones + revised date",
        "context": "",
    },
    "pdf_extractor": {
        "user_request": "提取 PDF 章节大纲、表格摘要与关键结论",
        "context": "",
        "file_path": "",
    },
    "deep_research": {
        "user_request": "间歇性禁食对代谢健康的影响：益处、风险与适用人群",
        "context": "侧重近五年综述与临床指南；标注证据等级",
    },
    "blueprint_drafter": {
        "user_request": (
            "绘制转录组差异分析流程蓝图：原始 FASTQ → 质控(QC) → 比对定量 → "
            "基因表达矩阵 → 差异表达分析(DEG) → 火山图/热图可视化"
        ),
        "context": "扁平工程图；中文节点标签；突出质控与统计步骤",
    },
    "diff_expr_interpreter": {
        "user_request": "解读下列 Demo DEG 结果并给出验证建议",
        "context": (
            "对比：处理组 vs 对照组；阈值 padj<0.05, |log2FC|>1\n"
            "TP53\t2.1\t1.2e-8\nCDKN1A\t1.8\t3.4e-6\nMYC\t-1.5\t0.002\n"
            "IL6\t1.2\t0.01\nSTAT3\t0.9\t0.03"
        ),
    },
    "go_kegg_narrative": {
        "user_request": "撰写富集结果叙事，突出免疫与代谢主题",
        "context": (
            "GO:0006955 immune response\tpadj=2.1e-5\tGeneRatio=12/180\n"
            "KEGG:hsa04668 TNF signaling\tpadj=8.3e-4\tGeneRatio=9/180\n"
            "GO:0006091 energy metabolism\tpadj=0.012\tGeneRatio=7/180"
        ),
    },
    "single_cell_checklist": {
        "user_request": (
            "人肺腺癌免疫治疗前后 scRNA-seq；比较响应者与非响应者；计划 10x 3' v3"
        ),
        "context": "预期每组 6 例生物学重复；关注 T 细胞耗竭与髓系抑制细胞",
    },
    "journal_cover_letter": {
        "user_request": (
            "目标期刊：Nature Communications；Article type: Research Article；语言：英文"
        ),
        "context": (
            "Title: Single-cell atlas reveals T cell exhaustion dynamics under PD-1 blockade\n"
            "Highlights: (1) New exhaustion trajectory (2) Cross-cohort validation "
            "(3) Proposed biomarker panel"
        ),
    },
    "review_rebuttal_outline": {
        "user_request": "稿件：scRNA-seq 免疫治疗研究；拟回复 Major revision",
        "context": (
            "Reviewer 1: (1) 缺少独立队列验证 (2) 双胞阈值未说明 "
            "(3) 部分标志基因未做 IHC\n"
            "Reviewer 2: 统计方法需澄清伪重复"
        ),
    },
    "acmg_variant_interpretation": {
        "user_request": "BRCA1 胚系变异 ACMG 风格草稿；GRCh38",
        "context": (
            "c.5266dupC (p.Gln1756Profs*74); 家族乳腺癌史; "
            "gnomAD 频率未提供; 预测：蛋白截断"
        ),
    },
    "pipeline_selection_memo": {
        "user_request": "人肿瘤 bulk RNA-seq；20 对肿瘤/癌旁；差异+通路",
        "context": "32 核 128G；偏好开源；无 WGS",
    },
    "omics_metadata_review": {
        "user_request": "RNA-seq 肿瘤/癌旁 metadata 审查",
        "context": (
            "SampleID\tGroup\tBatch\tRIN\n"
            "S1\tTumor\tB1\t8.2\nS2\tNormal\tB1\t\nS3\tTumor\tB2\t7.5"
        ),
    },
    "qpcr_primer_design_guide": {
        "user_request": "人源 IL6 mRNA 相对定量；SYBR Green；组织 RNA",
        "context": "内参候选 GAPDH、ACTB；无现成引物",
    },
    "oral_presentation_outline": {
        "user_request": "组会口头报告 12 分钟；听众：实验室研究生",
        "context": (
            "题目：PD-1 阻断下 T 细胞耗竭；结构：背景→方法→结果→验证→展望"
        ),
    },
    "clinical_trial_protocol_skeleton": {
        "user_request": "II 期单臂免疫联合化疗；晚期 NSCLC；主要终点 ORR",
        "context": "EGFR 野生型；次要 PFS、安全性",
    },
    "drug_repositioning_memo": {
        "user_request": "二甲双胍用于 IPF 抗纤维化重定位假说",
        "context": "AMPK 激活；动物模型部分有效；TGF-β/Smad 推测",
    },
    "protein_function_hypothesis": {
        "user_request": "TP53 功能假说；肿瘤抑制",
        "context": "DNA 结合域；核定位；细胞周期检查点",
    },
    "experiment_failure_postmortem": {
        "user_request": "Western 目标条带弱且无特异性",
        "context": "一抗 1:1000 过夜；冻融组织；预期 55 kDa；多条非特异",
    },
    "literature_matrix_notes": {
        "user_request": "免疫治疗预测生物标志物；3 篇代表文献",
        "context": "A: scRNA n=24; B: bulk meta; C: 空间组试点",
    },
    "multi_omics_storyline": {
        "user_request": "肺癌免疫耐药；基因组+转录组+蛋白组",
        "context": "WGS TP53 突变；RNA 耗竭 T 细胞；蛋白 STAT3-p 升高",
    },
    "ethics_consent_checklist": {
        "user_request": "前瞻性生物样本库；健康人抽血；可识别",
        "context": "存储 5 年；可能测序；无干预",
    },
    "stats_method_advisor": {
        "user_request": "两组独立；连续表达量；非正态",
        "context": "n=18 vs 20；可能有批次；需多重校正",
    },
    "pubmed_query": {
        "query": "single cell RNA-seq immunotherapy",
        "limit": 8,
    },
    "uniprot_query": {
        "query": "TP53",
        "limit": 5,
    },
    "smiles_to_cid": {
        "smiles": "CC(=O)Oc1ccccc1C(=O)O",
        "max_cids": 10,
    },
    "seq_format_converter": {
        "sequence_text": ">demo_gene\nATGCGTACGTTAGCTAGCTAGCTAGCTAG\n",
        "input_format": "fasta",
        "output_format": "genbank",
        "record_id": "demo_gene",
    },
    "calc_molecular_weight": {
        "smiles": "CC(=O)Oc1ccccc1C(=O)O",
    },
    "clinvar_query": {
        "query": "BRCA1[gene]",
        "limit": 8,
    },
    "dbsnp_query": {
        "query": "rs699",
        "limit": 5,
    },
    "gwas_catalog_query": {
        "query": "type 2 diabetes",
        "limit": 8,
    },
    "reactome_query": {
        "query": "apoptosis",
        "limit": 8,
    },
    "interpro_query": {
        "query": "kinase",
        "limit": 8,
    },
    "geo_query": {
        "query": "breast cancer",
        "limit": 8,
    },
    "gene_protein_info_query": {
        "query": "TP53",
    },
    "mrna_sequence_fetch": {
        "sequence_or_path": "TP53",
    },
    "drug_id_crossref": {
        "query": "aspirin",
    },
    "refmet_metabolite_search": {
        "query": "cholesterol",
        "match_threshold": 0.8,
    },
    "rnacentral_query": {
        "query": "URS000065FE5A",
    },
    "ensembl_go_descendants": {
        "query": "GO:0006915",
        "limit": 15,
    },
    "opentargets_chembl_hierarchy": {
        "query": "CHEMBL25",
    },
    "fda_drug_label_search": {
        "query": "aspirin",
        "limit": 3,
    },
}


def _is_blank(val: Any) -> bool:
    if val is None:
        return True
    if isinstance(val, str):
        return not val.strip()
    return False


def apply_launch_demo_defaults(skill_id: str, params: Dict[str, Any]) -> Dict[str, Any]:
    """将演示默认值填入仍为空的字段（不覆盖用户已提供的非空值）。"""
    demo = LAUNCH_SKILL_DEMO_ARGS.get((skill_id or "").strip()) or {}
    out = dict(params)
    for key, default in demo.items():
        if _is_blank(out.get(key)):
            out[key] = default
    return out


def get_demo_arg(skill_id: str, param: str, default: Any = "") -> Any:
    demo = LAUNCH_SKILL_DEMO_ARGS.get((skill_id or "").strip()) or {}
    return demo.get(param, default)
