# -*- coding: utf-8 -*-
"""无 [Skill_Route] 占位技能 · demo_visualization 具象化（生物医药 + 化学）。"""
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


UNROUTED_VIZ: Dict[str, str] = {
    "evo2": _report(
        "Evo2 · 长程基因组上下文与变异效应（示意）",
        _table(
            ["窗口", "序列长度", "SNV 位点", "Δlog-likelihood", "解读"],
            [
                ["chr17:7674000–7676000", "2 kb", "c.524G&gt;A", "+2.8", "局部上下文敏感"],
                ["TP53 外显子区", "180 bp", "同义突变", "+0.1", "中性倾向"],
                ["启动子近端", "512 bp", "Indel −3 bp", "+4.1", "需实验验证"],
            ],
        )
        + '<div style="font-family:ui-monospace;font-size:10px;background:#f8fafc;padding:8px;border-radius:6px;color:#475569;">'
        "…GACGCAGAAGACGGTGATTTCTGCATTTCCATCTGAG…</div>",
        "* 占位预览；正式接入后将展示真实模型打分与碱基概率热图。",
    ),
    "esm3": _report(
        "ESM3 · 多模态蛋白生成（示意）",
        _table(
            ["模态", "输入", "输出摘要"],
            [
                ["序列", "骨架约束 120 aa", "IFVQLVES… 设计序列"],
                ["结构", "Cα 坐标片段", "RMSD 1.8 Å vs 模板"],
                ["功能提示", "激酶活性", "ATP 结合位点保守性 0.82"],
            ],
        )
        + '<p style="margin:8px 0 0;color:#334155;">联合序列–结构–功能条件采样 · 新型蛋白骨架探索</p>',
    ),
    "alphafold2": _report(
        "AlphaFold2 · 单链结构预测 · EGFR 激酶域",
        _table(
            ["指标", "值", "说明"],
            [
                ["pLDDT 均值", "88.4", "高置信"],
                ["pTM", "0.91", "拓扑可信"],
                ["可建模残基", "298 / 310", "96%"],
                ["PDB 输出", "AF-EGFR-kinase.pdb", "含 B-factor=pLDDT"],
            ],
        )
        + '<div style="height:70px;background:linear-gradient(90deg,#dc2626 0%,#fbbf24 35%,#22c55e 70%,#2563eb 100%);border-radius:6px;margin-top:8px;opacity:0.85;"></div>'
        + '<p style="margin:6px 0 0;font-size:11px;color:#64748b;">pLDDT 色带：低 → 高置信（示意）</p>',
    ),
    "alphafold2_multimer": _report(
        "AlphaFold2-Multimer · 复合物预测",
        _table(
            ["链", "UniProt", "ipTM", "pLDDT", "界面 pLDDT"],
            [
                ["A", "P00533 (EGFR)", "0.78", "86.2", "72.4"],
                ["B", "P04626 (ERBB2)", "0.78", "84.9", "71.8"],
            ],
        )
        + '<p style="margin:8px 0 0;color:#334155;"><strong>复合物 ipTM 0.78</strong> · 二聚界面置信度良好（示意）</p>',
    ),
    "diffsbdd": _report(
        "DiffSBDD · 结构驱动配体生成",
        _table(
            ["Rank", "SMILES (节选)", "对接分数", "口袋匹配"],
            [
                ["1", "Cc1nc(N)nc2c1ccc(O)c2", "−8.4 kcal/mol", "0.91"],
                ["2", "COc1cc2ncnc(N)c2cc1F", "−7.9 kcal/mol", "0.87"],
                ["3", "Clc1ccc2[nH]cnc2c1", "−7.2 kcal/mol", "0.84"],
            ],
        )
        + '<p style="margin:8px 0 0;color:#64748b;">基于口袋 3D 网格的扩散生成 · 需受体 PDB + 结合位点</p>',
    ),
    "oligoformer": _report(
        "OligoFormer · siRNA 设计推荐",
        _table(
            ["排名", "正义链 (19nt)", "ΔG 双链", "GC%", "Off-target 风险"],
            [
                ["1", "GCAUGAGUCUAGCUAGCUA", "−28.4 kcal/mol", "47%", "低"],
                ["2", "UAGCUAGCUAGCUAGCAUG", "−26.1 kcal/mol", "42%", "中"],
                ["3", "CUGAGUCUAGCUAGCAUGA", "−25.3 kcal/mol", "53%", "低"],
            ],
        )
        + '<p style="margin:8px 0 0;color:#334155;">靶 mRNA：<code>NM_000546.5</code> (TP53) 3′UTR 区域（示意）</p>',
    ),
    "msa_search": _report(
        "MSA search · 多序列比对摘要",
        _table(
            ["统计量", "值"],
            [
                ["查询长度", "412 aa"],
                ["命中序列数", "248"],
                ["有效比对列", "385"],
                ["平均同一性", "34.2%"],
                ["保守块", "12 个"],
            ],
        )
        + '<div style="font-family:ui-monospace;font-size:9px;background:#f8fafc;padding:8px;border-radius:6px;margin-top:6px;color:#475569;">'
        "CLUSTAL 格式 MSA 节选 · 输出 FASTA / Stockholm</div>",
    ),
    "esmfold": _report(
        "ESMFold · 快速折叠 · 单域抗体 VH",
        _table(
            ["指标", "值"],
            [
                ["序列长度", "118 aa"],
                ["推理耗时", "4.2 s (GPU 示意)"],
                ["pLDDT 均值", "82.6"],
                ["输出", "ESMFold-VH.pdb"],
            ],
        ),
    ),
    "protgpt2": _report(
        "ProtGPT2 · 从头蛋白序列生成",
        '<div style="font-family:ui-monospace;font-size:10px;background:#f8fafc;padding:10px;border-radius:6px;margin-bottom:8px;">'
        "<div style='color:#64748b;margin-bottom:4px;'>Prompt 末 20 aa</div>"
        "<div style='color:#475569;'>…MKTAYIAKQRQISFVKSHFSRQDILDLWIYHTQGYFPD</div>"
        "<div style='color:#64748b;margin:8px 0 4px;'>生成 100 aa (T=0.8)</div>"
        "<div style='word-break:break-all;color:#1e40af;'>"
        "MKTAYIAKQRQISFVKSHFSRQDILDLWIYHTQGYFPDQKPGKAPKLLIYWA"
        "SFLYSGVPSRFSGSRSGTDFTLTISSLQPEDFATYYCQQYSTVPWTFGQGTKVEIK"
        "</div></div>"
        + _table(
            ["Perplexity", "结构倾向", "低复杂度区"],
            [["2.1", "α/β 混合", "无"]],
        ),
    ),
    "rfdiffusion": _report(
        "RFdiffusion · 结合剂骨架设计",
        _table(
            ["设计 #", "靶标", "界面残基", "RMSD", "可折叠性"],
            [
                ["1", "PD-L1 β-sheet", "18", "1.2 Å", "0.89"],
                ["2", "PD-L1 β-sheet", "22", "1.5 Å", "0.85"],
                ["3", "PD-L1 β-sheet", "16", "1.8 Å", "0.81"],
            ],
        )
        + '<p style="margin:8px 0 0;color:#64748b;">输出骨架 PDB → 后续 ProteinMPNN / AlphaFold 验证流程</p>',
    ),
    "biogpt": _report(
        "BioGPT · 生物医学文本理解",
        _table(
            ["任务", "输入片段", "输出"],
            [
                ["NER", "…aspirin inhibits COX-1…", "Drug: aspirin · Gene: COX-1"],
                ["关系", "metformin treats T2DM", "treats(metformin, T2DM)"],
                ["摘要", "Long abstract…", "3-sentence summary"],
            ],
        ),
    ),
    "biograph": _report(
        "BioGraph · 基因表达可视化",
        '<div style="display:grid;grid-template-columns:1fr 1fr;gap:8px;margin-bottom:8px;">'
        '<div style="height:80px;background:linear-gradient(180deg,#1e40af,#dc2626,#fbbf24);border-radius:6px;display:flex;align-items:center;justify-content:center;color:#fff;font-size:11px;">热图 (48×12 genes)</div>'
        '<div style="height:80px;background:#eff6ff;border-radius:6px;display:flex;align-items:center;justify-content:center;color:#1e40af;font-size:11px;">PCA · PC1 42% PC2 18%</div>'
        "</div>"
        + _table(
            ["图表", "用途"],
            [["Cluster heatmap", "样本 × 基因"], ["Volcano overlay", "DEG 标注"]],
        ),
    ),
    "protein_info_extractor": _report(
        "蛋白质资料提取 · P04637 (TP53)",
        _table(
            ["字段", "内容"],
            [
                ["蛋白名", "Cellular tumor antigen p53"],
                ["亚细胞定位", "Nucleus · Cytoplasm · Mitochondrion"],
                ["组织表达", "Ubiquitous · 高：spleen, colon"],
                ["疾病关联", "Li-Fraumeni · 多种肿瘤"],
                ["功能", "转录因子 · DNA 损伤应答"],
            ],
        ),
    ),
    "protein_structure_renderer": _report(
        "蛋白质 3D 结构渲染 · PDB 1M17",
        _table(
            ["渲染模式", "配色", "输出"],
            [
                ["Cartoon", "二级结构", "1M17_cartoon.png"],
                ["Surface", "疏水性", "1M17_surface.png"],
                ["Sticks + H-bond", "配体 ERK 抑制剂", "1M17_ligand.png"],
            ],
        )
        + '<div style="height:80px;background:#eff6ff;border:1px dashed #93c5fd;border-radius:8px;display:flex;align-items:center;justify-content:center;color:#1e40af;font-size:11px;margin-top:8px;">WebGL / PyMOL 风格预览占位</div>',
    ),
    "rna_secondary_structure_viz": _report(
        "RNA 二级结构可视化 · 点括号图",
        '<div style="font-family:ui-monospace;font-size:11px;background:#f8fafc;padding:10px;border-radius:6px;margin-bottom:8px;">'
        "<div>序列: GGGGAUAGGUUCAACCUCCUU</div>"
        "<div style='color:#1e40af;margin-top:4px;'>结构: (((((.....)))))....</div>"
        "<div style='color:#64748b;margin-top:4px;'>MFE: −12.8 kcal/mol · 配对碱基: 10 / 22</div></div>"
        + _table(
            ["茎区", "位置", "类型"],
            [
                ["Stem 1", "1–5 ⟷ 16–20", "GC 丰富"],
                ["环区", "6–15", "内环 9 nt"],
            ],
        )
        + '<div style="height:70px;background:#fff;border:1px solid #e2e8f0;border-radius:6px;display:flex;align-items:center;justify-content:center;color:#64748b;font-size:11px;">'
        "弧线图 / 径向布局 · SVG / PNG · 碱基配色按配对概率</div>",
    ),
    "antibody_humanization": _report(
        "抗体人源化 · CDR 移植与回复突变",
        _table(
            ["链", "原序列同一性", "人源化后", "CDR 保留"],
            [
                ["VH", "78%", "94%", "H1/H2/H3 移植"],
                ["VL", "81%", "96%", "L1/L2/L3 移植"],
            ],
        )
        + _table(
            ["回复突变位点", "原因"],
            [["VH-F68L", "界面 packing"], ["VL-Y36F", "构象稳定"]],
        ),
    ),
    "bacterial_growth_curve": _report(
        "细菌生长曲线拟合 · OD600",
        _table(
            ["参数", "估计值", "95% CI"],
            [
                ["最大比生长速率 μmax", "0.42 h⁻¹", "0.38–0.46"],
                ["延迟期 λ", "1.8 h", "1.5–2.1"],
                ["倍增时间 td", "1.65 h", "1.51–1.82"],
                ["最大密度 ODmax", "2.14", "2.08–2.20"],
            ],
        )
        + '<div style="height:60px;background:linear-gradient(180deg,#fff 0%,#dbeafe 100%);border-radius:6px;margin-top:8px;border:1px solid #e2e8f0;"></div>',
    ),
    "ligandmpnn": _report(
        "LigandMPNN · 配体环境序列设计",
        _table(
            ["位点", "原 aa", "设计 aa", "ΔΔG (示意)"],
            [
                ["A:45", "Leu", "Phe", "−0.8"],
                ["A:82", "Val", "Ile", "−0.3"],
                ["A:103", "Ser", "Thr", "−0.5"],
            ],
        )
        + '<p style="margin:8px 0 0;color:#334155;">固定配体 + 骨架 · 显式考虑非蛋白环境</p>',
    ),
    "esm_variants": _report(
        "ESM-Variants · 氨基酸突变效应",
        _table(
            ["位点", "突变", "Δ log-likelihood", "预测影响"],
            [
                ["156", "R→H", "+3.2", "有害倾向"],
                ["220", "V→M", "+0.4", "中性"],
                ["273", "R→C", "+5.1", "强有害"],
            ],
        )
        + '<div style="height:40px;background:linear-gradient(90deg,#22c55e,#fbbf24,#dc2626);border-radius:4px;margin-top:8px;"></div>',
    ),
    "antibody_sequence_generation": _report(
        "抗体序列生成 · 突变扫描",
        _table(
            ["变体", "CDR-H3 变化", "疏水性", "免疫原性评分"],
            [
                ["WT", "—", "0.42", "0.31"],
                ["M1", "YYD→YYG", "0.38", "0.28"],
                ["M2", "YYD→YFD", "0.45", "0.35"],
            ],
        ),
    ),
    "alphafold_db_query": _report(
        "AlphaFold DB · P04637 (p53)",
        _table(
            ["字段", "值"],
            [
                ["UniProt", "P04637"],
                ["AF 模型 ID", "AF-P04637-F1"],
                ["序列长度", "393 aa"],
                ["pLDDT 均值", "84.2"],
                ["下载", "PDB / mmCIF / PAED"],
            ],
        )
        + '<p style="margin:8px 0 0;color:#334155;">置信度分区：高 (pLDDT&gt;90) 58% · 低 &lt;50 仅 4%</p>',
    ),
    "ucsc_genome_browser_query": _report(
        "UCSC 基因组浏览器 · chr17 TP53 区域",
        _table(
            ["Track", "装配", "坐标", "摘要"],
            [
                ["RefSeq Genes", "hg38", "chr17:7668402-7687550", "TP53 外显子结构"],
                ["ClinVar", "hg38", "chr17:7675088", "Pathogenic 变异簇"],
                ["ENCODE H3K27ac", "hg38", "chr17:7676000", "启动子活跃"],
            ],
        ),
    ),
    "gsea_supported_databases_query": _report(
        "GSEA 支持的数据库列表",
        _table(
            ["库名", "物种", "基因集数", "用途"],
            [
                ["MSigDB Hallmark", "Human", "50", "通路 hallmark"],
                ["KEGG", "Human/Mouse", "186+", "代谢与信号"],
                ["GO Biological Process", "Multi", "7500+", "功能注释"],
                ["Reactome", "Human", "1500+", "反应通路"],
            ],
        ),
    ),
    "chem_molecule_viewer": _report(
        "分子可视化 · RDKit.js 2D 结构",
        _table(
            ["输入格式", "渲染", "交互"],
            [
                ["SMILES", "2D 键线式", "缩放 / 导出 PNG"],
                ["SDF / Mol", "立体化学保留", "原子 hover 提示"],
                ["InChI", "标准化后绘制", "—"],
            ],
        )
        + '<div style="height:90px;background:#fff;border:1px solid #e2e8f0;border-radius:8px;display:flex;align-items:center;justify-content:center;margin-top:8px;">'
        '<svg width="120" height="60" xmlns="http://www.w3.org/2000/svg">'
        '<polygon points="60,8 88,24 88,52 60,68 32,52 32,24" fill="none" stroke="#2563eb" stroke-width="1.5"/>'
        '<text x="60" y="42" text-anchor="middle" font-size="9" fill="#64748b">RDKit.js</text></svg></div>',
    ),
    "lammps_md": _report(
        "LAMMPS · 分子动力学模拟摘要",
        _table(
            ["阶段", "系综", "步数", "关键观测"],
            [
                ["NVT 平衡", "298 K", "50 ps", "T 收敛 ±2 K"],
                ["NPT 生产", "298 K, 1 bar", "2 ns", "ρ = 0.998 g/cm³"],
                ["RDF 分析", "O–H", "—", "峰值 2.8 Å"],
            ],
        )
        + '<p style="margin:8px 0 0;color:#64748b;">输入：data 文件 + 力场 · 输出：traj.dcd / log.lammps</p>',
    ),
    "cp2k_qc": _report(
        "CP2K · 量子化学 / AIMD 摘要",
        _table(
            ["计算类型", "方法", "基组", "结果"],
            [
                ["单点能", "DFT B3LYP", "6-31G*", "E = −232.451 Ha"],
                ["几何优化", "DFT PBE", "TZV2P", "max force &lt; 1e-4"],
                ["AIMD", "Born-Oppenheimer", "300 K", "5 ps 轨迹"],
            ],
        )
        + '<p style="margin:8px 0 0;color:#334155;">固体 / 液体 / 气相体系 · 需 CP2K 输入 deck 与 HPC 资源</p>',
    ),
}
