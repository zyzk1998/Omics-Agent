# -*- coding: utf-8 -*-
"""生物医药大类 · 剩余技能 demo_visualization（Phase 1.3 批量具象化）。"""
from __future__ import annotations

from typing import Dict


def _report(title: str, body: str, footnote: str = "") -> str:
    foot = (
        f'<div style="color:#64748b;font-size:11px;margin-top:10px;">{footnote}</div>'
        if footnote
        else ""
    )
    return f"""
<div class="skill-viz-report" style="font-family:system-ui,-apple-system,sans-serif;font-size:12px;line-height:1.55;padding:14px;background:#fff;border:1px solid #e2e8f0;border-radius:8px;">
  <div style="font-weight:700;font-size:14px;color:#1e40af;margin-bottom:10px;">{title}</div>
  {body}
  {foot}
</div>"""


def _table(headers: list, rows: list) -> str:
    th = "".join(
        f'<th style="padding:6px 8px;text-align:left;border:1px solid #e2e8f0;">{h}</th>'
        for h in headers
    )
    trs = []
    for row in rows:
        tds = "".join(
            f'<td style="padding:6px 8px;border:1px solid #e2e8f0;">{c}</td>' for c in row
        )
        trs.append(f"<tr>{tds}</tr>")
    return (
        '<table style="width:100%;border-collapse:collapse;font-size:11px;margin-bottom:8px;">'
        f'<thead><tr style="background:#f1f5f9;">{th}</tr></thead>'
        f"<tbody>{''.join(trs)}</tbody></table>"
    )


BIOPHARMA_REMAINDER_VIZ: Dict[str, str] = {
    "acmg_variant_interpretation": _report(
        "ACMG 风格变异解读 · TP53 c.524G&gt;A (p.Arg175His)",
        _table(
            ["证据项", "强度", "说明"],
            [
                ["PS1", "强", "同氨基酸改变已知致病变异"],
                ["PM2", "中", "gnomAD 人群频率极低"],
                ["PP3", "支持", "多种 in silico 预测有害"],
                ["BP4", " benign支持", "ClinVar 部分条目冲突需复核"],
            ],
        )
        + '<p style="margin:8px 0 0;color:#334155;"><strong>综合分类：</strong> <span style="color:#dc2626;">Likely Pathogenic</span>（需功能实验确认）</p>',
        "* 示意草案；正式报告须逐条引用 PMID / 数据库 accession。",
    ),
    "cell_cycle_phase_duration_estimation": _report(
        "EdU/BrdU 双脉冲 · 细胞周期时相估计",
        _table(
            ["时相", "时长 (h)", "95% CI", "备注"],
            [
                ["G1", "5.8", "5.2–6.4", "拟合收敛"],
                ["S", "7.6", "7.0–8.2", "EdU 掺入峰对齐"],
                ["G2/M", "3.9", "3.4–4.5", "—"],
                ["死亡率", "0.018", "—", "每 h 分数"],
            ],
        ),
    ),
    "chipatlas_experiment_search": _report(
        "ChIPAtlas 实验检索 · H3K27ac / K562",
        _table(
            ["Experiment ID", "Antibody", "Cell", "Peaks"],
            [
                ["SRX1234567", "H3K27ac", "K562", "48,210"],
                ["SRX1234588", "H3K27ac", "K562", "51,902"],
                ["SRX1234601", "H3K4me3", "K562", "22,145"],
            ],
        ),
    ),
    "clinical_trial_protocol_skeleton": _report(
        "临床试验方案骨架 · Phase II 单臂",
        "<ol style='margin:0;padding-left:1.2em;color:#334155;'>"
        "<li><strong>标题与注册号占位</strong> — CONSORT 2025 扩展项</li>"
        "<li><strong>主要终点</strong> — ORR @ 12 周 (RECIST 1.1)</li>"
        "<li><strong>次要终点</strong> — PFS, DOR, 安全性 (CTCAE v5.0)</li>"
        "<li><strong>入排标准</strong> — ECOG 0–1；可测量病灶 ≥1；既往线数 ≤2</li>"
        "<li><strong>统计计划</strong> — Simon 两阶段；α=0.05, β=0.20；目标 ORR 35%</li>"
        "<li><strong>安全性监测</strong> — DSMB 章程附件；SAE 24h 上报</li>"
        "<li><strong>数据管理</strong> — eCRF / 源数据核查计划</li></ol>"
        + _table(
            ["里程碑", "时间窗"],
            [["首例入组", "M0"], ["中期分析", "M6"], ["数据库锁库", "M18"]],
        ),
    ),
    "clinvar_query": _report(
        "ClinVar 检索 · rs1042522 (TP53)",
        _table(
            ["Variation ID", "Clinical significance", "Review status", "Condition"],
            [
                ["VCV000012345", "Pathogenic / Likely pathogenic", "reviewed by expert panel", "Li-Fraumeni syndrome"],
                ["VCV000098765", "Uncertain significance", "criteria provided", "Hereditary cancer predisposition"],
            ],
        ),
    ),
    "dbsnp_query": _report(
        "dbSNP 检索 · rs1042522",
        _table(
            ["rsID", "Gene", "Alleles", "MAF (EUR)", "Functional class"],
            [
                ["rs1042522", "TP53", "C/T", "0.28", "missense_variant"],
                ["rs17878362", "TP53", "del/—", "0.02", "inframe_deletion"],
            ],
        ),
    ),
    "dna_sequence_generate": _report(
        "DNA 续写生成 · 提示 120 nt → 续写 100 nt",
        "<div style='font-family:ui-monospace,monospace;font-size:10px;background:#f8fafc;padding:10px;border-radius:6px;'>"
        "<div style='color:#64748b;margin-bottom:4px;'>Prompt 末 40 nt</div>"
        "<div style='word-break:break-all;color:#475569;'>…GACGCAGAAGACGGTGATTTCTGCATTTCCATCTGAG</div>"
        "<div style='color:#64748b;margin:8px 0 4px;'>续写片段 (temperature=0.8)</div>"
        "<div style='word-break:break-all;color:#1e40af;'>GTACCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAG</div>"
        "</div>"
        + _table(
            ["指标", "值"],
            [["平均 step 概率", "0.41"], ["GC%", "52.1"], ["终止密码子", "未出现"]],
        ),
    ),
    "drug_repositioning_memo": _report(
        "药物重定位假说 · Metformin → 抗肿瘤免疫",
        "<p style='margin:0 0 8px;color:#334155;'><strong>假说：</strong>AMPK 激活 → mTOR 抑制 → CD8+ T 细胞代谢重编程增强 PD-1 阻断响应。</p>"
        + _table(
            ["证据层", "要点", "置信"],
            [
                ["机制", "AMPK/mTOR/HIF1α 轴", "中"],
                ["临床前", "MC38 同系模型联合 anti-PD-1", "中"],
                ["临床", "回顾性队列 n=412", "低-中"],
            ],
        ),
    ),
    "ensembl_go_descendants": _report(
        "Ensembl GO 后代术语 · GO:0006955 (immune response)",
        _table(
            ["GO ID", "Term", "Depth"],
            [
                ["GO:0002250", "adaptive immune response", "2"],
                ["GO:0002449", "lymphocyte mediated immunity", "3"],
                ["GO:0002456", "T cell mediated immunity", "4"],
            ],
        ),
    ),
    "ethics_consent_checklist": _report(
        "伦理与知情同意要点清单",
        "<ul style='margin:0 0 10px;padding-left:1.2em;color:#334155;'>"
        "<li>☑ 研究方案已通过 IRB/伦理委员会批准（批件号占位）</li>"
        "<li>☑ 知情同意书含可撤回条款与数据二次使用说明</li>"
        "<li>☐ 脆弱人群额外保护（未成年人/孕妇/认知障碍，如适用）</li>"
        "<li>☑ 遗传样本出境/共享合规说明（人类遗传资源条例）</li>"
        "<li>☑ 数据脱敏与访问权限分级</li></ul>"
        + _table(
            ["检查项", "状态", "负责人"],
            [
                ["ICF 版本与日期", "v2.1 / 2026-01", "PI"],
                ["风险获益说明", "已审", "伦理秘书"],
                ["替代治疗说明", "已纳入", "PI"],
            ],
        ),
    ),
    "experiment_failure_postmortem": _report(
        "实验失败复盘 · Western Blot 无条带",
        _table(
            ["维度", "发现", "行动项"],
            [
                ["假设", "一抗种属不匹配", "更换兔源一抗"],
                ["操作", "转膜电流过高", "250 mA × 90 min → 100 V 限压"],
                ["对照", "内参 GAPDH 正常", "排除上样/转膜全局失败"],
            ],
        ),
    ),
    "gene_protein_info_query": _report(
        "基因蛋白注释 · TP53 (Homo sapiens)",
        _table(
            ["字段", "值"],
            [
                ["Official symbol", "TP53"],
                ["Chromosome", "17p13.1"],
                ["UniProt", "P04637"],
                ["Function", "Cell cycle arrest, apoptosis, DNA repair"],
                ["Subcellular", "Nucleus, cytosol"],
            ],
        ),
    ),
    "geo_query": _report(
        "GEO 检索 · single-cell RNA-seq PD-1 blockade",
        _table(
            ["Accession", "Title", "Samples", "Platform"],
            [
                ["GSE123456", "scRNA atlas under PD-1 blockade", "12", "GPL24676"],
                ["GSE123457", "T cell exhaustion dynamics", "8", "GPL24676"],
            ],
        ),
    ),
    "gseapy_analysis": _report(
        "GSEA 富集 · MSigDB Hallmark",
        _table(
            ["Term", "NES", "FDR q-val", "Leading edge"],
            [
                ["HALLMARK_INTERFERON_GAMMA_RESPONSE", "2.14", "0.001", "42/200"],
                ["HALLMARK_G2M_CHECKPOINT", "1.87", "0.008", "38/200"],
                ["HALLMARK_OXIDATIVE_PHOSPHORYLATION", "-1.92", "0.012", "35/200"],
            ],
        ),
    ),
    "gwas_catalog_query": _report(
        "GWAS Catalog · Type 2 diabetes",
        _table(
            ["rsID", "Gene", "P-value", "OR", "Trait"],
            [
                ["rs7903146", "TCF7L2", "5×10⁻⁴⁸", "1.37", "Type 2 diabetes"],
                ["rs1801282", "PPARG", "2×10⁻¹⁵", "0.86", "Type 2 diabetes"],
            ],
        ),
    ),
    "interpro_query": _report(
        "InterPro · P04637 (p53)",
        _table(
            ["InterPro ID", "Entry name", "Type"],
            [
                ["IPR002117", "p53 tumour suppressor family", "Family"],
                ["IPR008967", "p53 DNA-binding domain", "Domain"],
                ["IPR010991", "p53 tetramerisation domain", "Domain"],
            ],
        ),
    ),
    "itc_binding_thermodynamic_analysis": _report(
        "ITC 一位点结合拟合 · 蛋白–配体",
        _table(
            ["参数", "值", "单位", "误差"],
            [
                ["K<sub>d</sub>", "2.4", "μM", "±0.3"],
                ["ΔH", "-8.6", "kcal/mol", "±0.4"],
                ["ΔS", "-12.1", "cal/mol/K", "±1.2"],
                ["ΔG", "-7.2", "kcal/mol", "—"],
                ["n", "0.98", "—", "±0.05"],
            ],
        ),
        "拟合模型：One set of sites；χ² = 1.2×10³",
    ),
    "journal_cover_letter": _report(
        "Cover Letter 草稿 · Nature Communications",
        "<div style='font-family:Georgia,serif;font-size:12px;line-height:1.65;color:#1e293b;'>"
        "<p style='margin:0 0 8px;'>Dear Editor,</p>"
        "<p style='margin:0 0 8px;'>We submit <em>Single-cell atlas reveals T cell exhaustion dynamics under PD-1 blockade</em> "
        "for consideration as a Research Article in <strong>Nature Communications</strong>.</p>"
        "<p style='margin:0 0 8px;'><strong>Highlights:</strong> (1) cross-cohort validation (2) trajectory inference (3) biomarker panel</p>"
        "<p style='margin:0;'>This work is original, not under consideration elsewhere, and all authors approve submission.</p>"
        "</div>"
        + _table(
            ["附件", "说明"],
            [["Main manuscript", "PDF + 源文件"], ["Supplementary", "12 figures, 3 tables"], ["Author list", "ORCID 已链接"]],
        ),
    ),
    "literature_matrix_notes": _report(
        "文献矩阵精读 · PD-1 / T cell exhaustion",
        _table(
            ["PMID", "Design", "Key finding", "Limitation"],
            [
                ["38123456", "scRNA-seq, n=24", "Tex 轨迹分叉", "单中心"],
                ["37890123", "Bulk + IHC", "TOX+ 亚群", "无功能验证"],
            ],
        ),
    ),
    "mrna_sequence_fetch": _report(
        "mRNA 序列检索 · EGFR (NCBI)",
        "<div style='font-family:ui-monospace,monospace;font-size:10px;background:#f8fafc;padding:10px;border-radius:6px;'>"
        "<div style='color:#64748b;'>&gt;NM_005228.5 Homo sapiens EGFR mRNA, complete cds</div>"
        "<div style='word-break:break-all;color:#334155;margin-top:4px;'>ATGCGACCCTCCGGGACGGCCGGGGCAGCGCTCCTGGCGCTGCTGGCTGCGCTCTGCCCGGCGAGTCGGGCTCT…</div>"
        "<div style='color:#64748b;margin-top:6px;font-family:system-ui;'>长度 5,988 nt · CDS 4,401 nt · 3′ UTR 1,587 nt</div>"
        "</div>"
        + _table(
            ["字段", "值"],
            [
                ["Gene", "EGFR"],
                ["RefSeq", "NM_005228.5"],
                ["Protein", "NP_005219.2 (1210 aa)"],
                ["Chromosome", "7p11.2"],
            ],
        ),
    ),
    "multi_omics_storyline": _report(
        "多组学整合故事线 · 肿瘤免疫",
        "<div style='display:flex;gap:6px;flex-wrap:wrap;margin-bottom:8px;'>"
        "<span style='padding:4px 10px;background:#dbeafe;border-radius:4px;font-size:11px;'>转录组 DEG</span>"
        "<span style='color:#94a3b8;'>→</span>"
        "<span style='padding:4px 10px;background:#d1fae5;border-radius:4px;font-size:11px;'>蛋白验证</span>"
        "<span style='color:#94a3b8;'>→</span>"
        "<span style='padding:4px 10px;background:#fef3c7;border-radius:4px;font-size:11px;'>代谢通路</span>"
        "<span style='color:#94a3b8;'>→</span>"
        "<span style='padding:4px 10px;background:#fce7f3;border-radius:4px;font-size:11px;'>临床关联</span>"
        "</div>"
        "<p style='margin:0;color:#475569;'>核心叙事：耗竭 Tex 亚群与糖酵解通量升高耦合，提示联合代谢干预假说。</p>",
    ),
    "nucleotide_sequence_blast": _report(
        "核酸 BLAST · 命中摘要",
        _table(
            ["Accession", "Description", "Identity", "E-value", "Coverage"],
            [
                ["NM_001126112.2", "Homo sapiens TP53", "99.2%", "0.0", "100%"],
                ["XM_024452123.1", "Pan troglodytes TP53", "98.1%", "1.2e-180", "98%"],
                ["NR_046018.2", "TP53 transcript variant", "97.5%", "3.4e-165", "95%"],
            ],
        ),
    ),
    "omics_metadata_review": _report(
        "组学 Metadata 规范审查",
        _table(
            ["字段", "状态", "建议"],
            [
                ["sample_id", "✓ 唯一", "—"],
                ["batch", "⚠ 缺失", "补充测序批次/flowcell"],
                ["condition", "✓", "对照/处理标签清晰"],
                ["organism", "✓ NCBI Taxon 9606", "—"],
            ],
        ),
    ),
    "oral_presentation_outline": _report(
        "口头报告讲稿 · 15 min",
        _table(
            ["时段", "Slide", "要点"],
            [
                ["0–2 min", "1–2", "背景与科学问题"],
                ["2–8 min", "3–6", "方法 + 关键结果图"],
                ["8–13 min", "7–8", "机制讨论与局限"],
                ["13–15 min", "9", "总结与致谢 / Q&A"],
            ],
        ),
    ),
    "pipeline_selection_memo": _report(
        "Pipeline 选型备忘录 · 单细胞 RNA-seq",
        _table(
            ["环节", "推荐", "备选", "理由"],
            [
                ["比对", "Cell Ranger", "STARsolo", "与 10x 官方兼容"],
                ["QC", "scDblFinder", "DoubletDetection", "社区基准较新"],
                ["整合", "Harmony", "scVI", "批次效应中等规模"],
            ],
        ),
    ),
    "protease_kinetics_analysis": _report(
        "蛋白酶动力学 · 荧光底物拟合",
        _table(
            ["参数", "估计值", "单位", "R²"],
            [
                ["K<sub>m</sub>", "18.3", "μM", "0.991"],
                ["V<sub>max</sub>", "245", "RFU/min", "—"],
                ["k<sub>cat</sub>/K<sub>m</sub>", "1.2×10⁵", "M⁻¹s⁻¹", "—"],
            ],
        ),
    ),
    "protein_function_hypothesis": _report(
        "蛋白功能假说推演 · 未特征蛋白 Q8NXXX",
        "<p style='margin:0 0 8px;color:#334155;'><strong>假说：</strong>含 DUFxxxx 结构域，可能参与 RNA 结合与应激颗粒组装。</p>"
        + _table(
            ["证据", "来源", "权重"],
            [
                ["InterPro DUF", "结构域注释", "高"],
                ["共表达模块", "GTEx 脑组织", "中"],
                ["PPI 预测", "G3BP1, TIA1", "中"],
            ],
        ),
    ),
    "protein_sequence_blast": _report(
        "蛋白 BLAST · 命中摘要",
        _table(
            ["Accession", "Organism", "Identity", "E-value", "Length"],
            [
                ["NP_000537.3", "Homo sapiens p53", "100%", "0.0", "393 aa"],
                ["XP_023812345.1", "Macaca mulatta", "96%", "2e-250", "393 aa"],
            ],
        ),
    ),
    "protein_sequence_conservation_analysis": _report(
        "多序列保守性 · TP53 同源比对",
        _table(
            ["物种", "Identity", "保守位点 (%)", "系统发育距离"],
            [
                ["Human", "100%", "—", "0"],
                ["Mouse", "78.5%", "62%", "0.21"],
                ["Zebrafish", "45.2%", "38%", "0.89"],
            ],
        )
        + "<div style='font-size:11px;color:#64748b;margin-top:6px;'>附：NJ 系统发育树与 Shannon 熵轨迹图（工作台）。</div>",
    ),
    "qpcr_primer_design_guide": _report(
        "qPCR 引物设计说明 · ACTB 内参",
        _table(
            ["引物", "序列 (5′→3′)", "Tm", "产物 (bp)"],
            [
                ["ACTB-F", "CTCTTCCAGCCTTCCTTCCT", "60.1°C", "154"],
                ["ACTB-R", "AGCACTGTGTTGGCGTACAG", "59.8°C", "—"],
            ],
        )
        + "<p style='margin:8px 0 0;color:#334155;'>建议：扩增效率 95–105%；溶解曲线单峰；无引物二聚体。</p>",
    ),
    "reactome_query": _report(
        "Reactome 通路 · Apoptosis",
        _table(
            ["Pathway ID", "Name", "Species", "Events"],
            [
                ["R-HSA-109581", "Apoptosis", "Homo sapiens", "128"],
                ["R-HSA-75153", "Apoptotic execution phase", "Homo sapiens", "42"],
            ],
        ),
    ),
    "review_rebuttal_outline": _report(
        "Rebuttal 要点 · Reviewer #2",
        _table(
            ["意见", "回应策略", "新增数据"],
            [
                ["样本量不足", "补充独立队列 n=32", "Figure R1"],
                ["批次效应", "Harmony 整合 + 敏感性分析", "Supp Fig S3"],
                ["机制深度", "CRISPR 验证关键基因", "Figure R2"],
            ],
        ),
    ),
    "rna_secondary_structure_analysis": _report(
        "RNA 二级结构 · 碱基配对概率",
        "<div style='font-family:ui-monospace;font-size:10px;background:#f8fafc;padding:8px;border-radius:6px;margin-bottom:8px;'>"
        "<div>结构: <code>(((((.....)))))....</code></div>"
        "<div>MFE: −12.8 kcal/mol · 配对概率均值: 0.78</div></div>"
        + _table(
            ["茎区", "位置", "平均 p(i,j)"],
            [
                ["Stem 1", "1–5 ⟷ 16–20", "0.92"],
                ["Stem 2", "6–8 ⟷ 12–14", "0.71"],
            ],
        ),
    ),
    "rnacentral_query": _report(
        "RNAcentral · URS0000123456",
        _table(
            ["字段", "值"],
            [
                ["URS ID", "URS0000123456"],
                ["RNA type", "miRNA"],
                ["Description", "Homo sapiens miR-21-5p"],
                ["Xrefs", "miRBase:MIMAT0000076"],
            ],
        ),
    ),
    "seq_format_converter": _report(
        "序列格式转换 · FASTA → GenBank",
        "<div style='font-family:ui-monospace,monospace;font-size:10px;background:#f8fafc;padding:10px;border-radius:6px;'>"
        "<div style='color:#1e40af;'>LOCUS       demo_gene 42 bp ds-DNA linear UNK 01-JAN-2026</div>"
        "<div>DEFINITION  Synthetic demo sequence for format conversion.</div>"
        "<div>ACCESSION   demo_gene</div>"
        "<div>VERSION     demo_gene.1</div>"
        "<div>KEYWORDS    .</div>"
        "<div>SOURCE      synthetic DNA</div>"
        "<div>ORIGIN</div>"
        "<div>        1 atgcgtacgt tagctagcta gctagctagc tag</div>"
        "<div>//</div></div>"
        + _table(
            ["转换", "输入", "输出"],
            [["FASTA → GenBank", "1 条序列", "LOCUS + ORIGIN 块"], ["校验", "ATGC 字母表", "通过"]],
        ),
    ),
    "stats_method_advisor": _report(
        "统计方法选择建议 · 两组成连续变量",
        _table(
            ["条件", "推荐", "备选"],
            [
                ["正态 + 方差齐", "Independent t-test", "Welch t-test"],
                ["非正态", "Mann-Whitney U", "Permutation test"],
                ["配对样本", "Paired t-test", "Wilcoxon signed-rank"],
            ],
        )
        + "<p style='margin:8px 0 0;color:#334155;'>请先报告 Shapiro-Wilk 与 Levene 检验结果再定稿。</p>",
    ),
    "go_kegg_narrative": _report(
        "GO/KEGG 富集叙事 · 免疫激活模块",
        _table(
            ["Term", "Database", "padj", "GeneRatio", "解读要点"],
            [
                ["GO:0006955 immune response", "GO", "2.1e-5", "12/180", "上游 IFN 信号增强"],
                ["KEGG:hsa04668 TNF signaling", "KEGG", "8.3e-4", "9/180", "与耗竭 Tex 表型一致"],
                ["GO:0002250 adaptive immune response", "GO", "1.2e-3", "8/180", "T 细胞受体信号富集"],
            ],
        )
        + "<p style='margin:8px 0 0;color:#475569;'>叙事建议：将 TNF/IFN 通路与临床响应亚组关联，避免仅罗列 term。</p>",
    ),
    "uniprot_query": _report(
        "UniProt 检索 · EGFR (P00533)",
        _table(
            ["字段", "值"],
            [
                ["Entry", "P00533"],
                ["Protein names", "Epidermal growth factor receptor"],
                ["Gene", "EGFR"],
                ["Organism", "Homo sapiens"],
                ["Function", "Receptor tyrosine kinase; EGF binding"],
                ["Subcellular location", "Cell membrane, Endosome"],
            ],
        )
        + _table(
            ["结构域", "位置"],
            [["EGF receptor L domain", "1–350"], ["Protein kinase domain", "712–979"]],
        ),
    ),
    "single_cell_checklist": _report(
        "单细胞实验设计检查清单",
        "<ul style='margin:0 0 10px;padding-left:1.2em;color:#334155;'>"
        "<li>☑ 样本量 ≥ 3 生物学重复 / 条件</li>"
        "<li>☑ 批次效应评估与整合策略预注册</li>"
        "<li>☑ 双胞检测工具选定 (scDblFinder / DoubletDetection)</li>"
        "<li>☑ 参考面板 / 线粒体比例 QC 阈值</li>"
        "<li>☐ 细胞周期回归策略（是否保留 cycling 信号）</li>"
        "<li>☑ 下游 DE/轨迹/通讯分析路径一致</li></ul>"
        + _table(
            ["QC 指标", "建议阈值"],
            [["nFeature_RNA", "200–6000"], ["percent.mt", "< 15%"], ["doublet rate", "< 8%"]],
        ),
    ),
}
