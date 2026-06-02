# -*- coding: utf-8 -*-
"""技能详情抽屉 · demo_visualization 静态 HTML/路径（Phase 1 首批 10 项）。"""
from __future__ import annotations

# 轻量级内联 HTML，由前端 trusted 渲染；禁止用户 UGC 直接写入此字段。
DEMO_VIZ_BLUEPRINT_RNASEQ = """
<div class="skill-viz-blueprint" style="font-family:system-ui,-apple-system,sans-serif;padding:20px;background:#f8fafc;border-radius:10px;">
  <div style="font-size:13px;font-weight:600;color:#64748b;margin-bottom:14px;">RNA-seq 差异分析流程蓝图（示意）</div>
  <div style="display:flex;align-items:center;justify-content:center;flex-wrap:wrap;gap:10px;">
    <div style="padding:10px 14px;background:#dbeafe;border:2px solid #3b82f6;border-radius:8px;font-size:12px;color:#1e40af;">原始 FASTQ</div>
    <span style="color:#94a3b8;">→</span>
    <div style="padding:10px 14px;background:#e0e7ff;border:2px solid #6366f1;border-radius:8px;font-size:12px;color:#3730a3;">质控 QC</div>
    <span style="color:#94a3b8;">→</span>
    <div style="padding:10px 14px;background:#d1fae5;border:2px solid #10b981;border-radius:8px;font-size:12px;color:#065f46;">表达矩阵</div>
    <span style="color:#94a3b8;">→</span>
    <div style="padding:10px 14px;background:#fef3c7;border:2px solid #f59e0b;border-radius:8px;font-size:12px;color:#92400e;">差异分析 DEG</div>
    <span style="color:#94a3b8;">→</span>
    <div style="padding:10px 14px;background:#fce7f3;border:2px solid #ec4899;border-radius:8px;font-size:12px;color:#9d174d;">火山图 / 热图</div>
  </div>
</div>
"""

DEMO_VIZ_PPT_OUTLINE = """
<table style="width:100%;border-collapse:collapse;font-family:system-ui;font-size:12px;">
  <thead><tr style="background:#f1f5f9;"><th style="padding:8px;text-align:left;border:1px solid #e2e8f0;">页</th>
  <th style="padding:8px;text-align:left;border:1px solid #e2e8f0;">标题</th><th style="padding:8px;text-align:left;border:1px solid #e2e8f0;">要点</th></tr></thead>
  <tbody>
  <tr><td style="padding:8px;border:1px solid #e2e8f0;">1</td><td style="padding:8px;border:1px solid #e2e8f0;">研究背景</td><td style="padding:8px;border:1px solid #e2e8f0;">单细胞时空动力学 · 科学问题</td></tr>
  <tr><td style="padding:8px;border:1px solid #e2e8f0;">2</td><td style="padding:8px;border:1px solid #e2e8f0;">方法</td><td style="padding:8px;border:1px solid #e2e8f0;">最优传输轨迹推断 · 质控</td></tr>
  <tr><td style="padding:8px;border:1px solid #e2e8f0;">3</td><td style="padding:8px;border:1px solid #e2e8f0;">核心发现</td><td style="padding:8px;border:1px solid #e2e8f0;">驱动基因 · 通路富集</td></tr>
  <tr><td style="padding:8px;border:1px solid #e2e8f0;">4</td><td style="padding:8px;border:1px solid #e2e8f0;">结论与展望</td><td style="padding:8px;border:1px solid #e2e8f0;">验证计划 · 后续工作</td></tr>
  </tbody>
</table>
"""

DEMO_VIZ_PUBMED = """
<table style="width:100%;border-collapse:collapse;font-family:system-ui;font-size:11px;">
  <thead><tr style="background:#f1f5f9;"><th style="padding:6px 8px;text-align:left;border:1px solid #e2e8f0;">PMID</th>
  <th style="padding:6px 8px;text-align:left;border:1px solid #e2e8f0;">Title</th><th style="padding:6px 8px;text-align:left;border:1px solid #e2e8f0;">Journal</th></tr></thead>
  <tbody>
  <tr><td style="padding:6px 8px;border-bottom:1px solid #e2e8f0;">38123456</td><td style="padding:6px 8px;border-bottom:1px solid #e2e8f0;">Single-cell atlas of T cell exhaustion…</td><td style="padding:6px 8px;border-bottom:1px solid #e2e8f0;">Nat Med</td></tr>
  <tr><td style="padding:6px 8px;border-bottom:1px solid #e2e8f0;">37890123</td><td style="padding:6px 8px;border-bottom:1px solid #e2e8f0;">scRNA-seq in immunotherapy response</td><td style="padding:6px 8px;border-bottom:1px solid #e2e8f0;">Cell</td></tr>
  </tbody>
</table>
"""

DEMO_VIZ_UNIPROT = """
<table style="width:100%;border-collapse:collapse;font-family:system-ui;font-size:11px;">
  <thead><tr style="background:#f1f5f9;"><th style="padding:6px 8px;text-align:left;border:1px solid #e2e8f0;">Entry</th>
  <th style="padding:6px 8px;text-align:left;border:1px solid #e2e8f0;">Protein</th><th style="padding:6px 8px;text-align:left;border:1px solid #e2e8f0;">Function</th></tr></thead>
  <tbody>
  <tr><td style="padding:6px 8px;border-bottom:1px solid #e2e8f0;">P04637</td><td style="padding:6px 8px;border-bottom:1px solid #e2e8f0;">Cellular tumor antigen p53</td><td style="padding:6px 8px;border-bottom:1px solid #e2e8f0;">Tumor suppressor</td></tr>
  </tbody>
</table>
"""

DEMO_VIZ_BEPIPRED = """
<div style="font-family:system-ui;font-size:11px;padding:8px;background:#fff;">
  <div style="font-weight:600;margin-bottom:8px;color:#334155;">BepiPred-3.0 Epitope Score（示意）</div>
  <div style="display:flex;align-items:flex-end;gap:2px;height:80px;border-bottom:1px solid #cbd5e1;padding-bottom:4px;">
    <div style="flex:1;background:linear-gradient(#6366f1,#a855f7);height:45%;border-radius:2px 2px 0 0;"></div>
    <div style="flex:1;background:linear-gradient(#6366f1,#eab308);height:72%;border-radius:2px 2px 0 0;"></div>
    <div style="flex:1;background:linear-gradient(#6366f1,#eab308);height:88%;border-radius:2px 2px 0 0;"></div>
    <div style="flex:1;background:linear-gradient(#6366f1,#a855f7);height:55%;border-radius:2px 2px 0 0;"></div>
    <div style="flex:1;background:linear-gradient(#6366f1,#6366f1);height:35%;border-radius:2px 2px 0 0;"></div>
  </div>
  <div style="margin-top:6px;color:#64748b;">— 阈值线（Top 20%）</div>
</div>
"""

DEMO_VIZ_SEQ_CONVERT = """
<div style="font-family:ui-monospace,monospace;font-size:10px;line-height:1.5;background:#f8fafc;padding:10px;border-radius:6px;">
  <div style="color:#64748b;margin-bottom:4px;">FASTA → GenBank</div>
  <div style="color:#1e40af;">LOCUS       demo_gene 42 bp ds-DNA linear UNK</div>
  <div>ORIGIN</div>
  <div>        1 atgcgtacgt tagctagcta gctagctagc tag</div>
  <div>//</div>
</div>
"""

DEMO_VIZ_MOL_CONVERT = """
<div style="font-family:ui-monospace,monospace;font-size:11px;padding:10px;background:#f8fafc;border-radius:6px;">
  <div style="color:#64748b;margin-bottom:4px;">SMILES → InChI</div>
  <div style="word-break:break-all;color:#1e293b;">InChI=1S/C9H8O4/c1-6(10)13-8-5-3-2-4-7(8)9(11)12/h2-5H,1H3,(H,11,12)</div>
</div>
"""

DEMO_VIZ_BLAST_HITS = """
<table style="width:100%;border-collapse:collapse;font-family:system-ui;font-size:11px;">
  <thead><tr style="background:#1e293b;color:#f8fafc;">
  <th style="padding:6px 8px;">Accession</th><th style="padding:6px 8px;">Description</th><th style="padding:6px 8px;">E-value</th></tr></thead>
  <tbody style="color:#334155;">
  <tr><td style="padding:5px 8px;border-bottom:1px solid #e2e8f0;">NM_001234</td><td style="padding:5px 8px;border-bottom:1px solid #e2e8f0;">Homo sapiens …</td><td style="padding:5px 8px;border-bottom:1px solid #e2e8f0;">2.3e-45</td></tr>
  <tr><td style="padding:5px 8px;border-bottom:1px solid #e2e8f0;">XM_567890</td><td style="padding:5px 8px;border-bottom:1px solid #e2e8f0;">Mus musculus …</td><td style="padding:5px 8px;border-bottom:1px solid #e2e8f0;">1.1e-12</td></tr>
  </tbody>
</table>
"""

DEMO_VIZ_RESUME = """
<div style="font-family:system-ui;font-size:12px;padding:12px;background:#fff;border:1px solid #e2e8f0;border-radius:8px;">
  <div style="font-weight:700;font-size:14px;">Jane Doe</div>
  <div style="color:#64748b;margin:4px 0 8px;">Senior Data Analyst · Healthcare</div>
  <div style="font-weight:600;margin-bottom:4px;">Summary</div>
  <div style="color:#475569;line-height:1.5;">5+ years SQL/Python · A/B testing · Tableau</div>
</div>
"""

DEMO_VIZ_COVER_LETTER = """
<div style="font-family:Georgia,serif;font-size:12px;line-height:1.6;padding:14px;background:#fff;border:1px solid #e2e8f0;border-radius:8px;color:#1e293b;">
  <p style="margin:0 0 8px;">Dear Editor,</p>
  <p style="margin:0 0 8px;">We submit our manuscript entitled <em>scRNA atlas under PD-1 blockade</em> for consideration as a Research Article in <strong>Nature Communications</strong>.</p>
  <p style="margin:0;color:#64748b;font-family:system-ui;font-size:11px;">（示意草稿 · 实际内容由技能生成）</p>
</div>
"""

DEMO_VIZ_SUBSTRUCTURE = """
<div style="font-family:system-ui;font-size:12px;padding:12px;background:#fff;border:1px solid #e2e8f0;border-radius:8px;">
  <div style="margin-bottom:8px;"><span style="padding:4px 8px;background:#dbeafe;border-radius:4px;">子结构: 苯环 c1ccccc1</span></div>
  <div><span style="padding:4px 8px;background:#d1fae5;border-radius:4px;">✓ c1ccccc1CCO 匹配</span></div>
</div>
"""

DEMO_VIZ_SMILES_CID = """
<table style="width:100%;border-collapse:collapse;font-family:system-ui;font-size:11px;">
  <thead><tr style="background:#f1f5f9;"><th style="padding:6px 8px;text-align:left;">CID</th><th style="padding:6px 8px;text-align:left;">Name</th></tr></thead>
  <tbody><tr><td style="padding:6px 8px;border-bottom:1px solid #e2e8f0;">2244</td><td style="padding:6px 8px;border-bottom:1px solid #e2e8f0;">Aspirin</td></tr></tbody>
</table>
"""

DEMO_VIZ_MINDMAP = """
<div style="font-family:system-ui;font-size:12px;line-height:1.6;padding:12px;background:#fff;border:1px solid #e2e8f0;border-radius:8px;">
  <div style="text-align:center;font-weight:600;color:#2563eb;margin-bottom:10px;">单细胞时空动力学</div>
  <div style="display:grid;grid-template-columns:1fr 1fr;gap:8px;">
    <div style="padding:8px;background:#eff6ff;border-radius:6px;">数据校验</div>
    <div style="padding:8px;background:#eff6ff;border-radius:6px;">时序标准化</div>
    <div style="padding:8px;background:#f0fdf4;border-radius:6px;">轨迹推断</div>
    <div style="padding:8px;background:#f0fdf4;border-radius:6px;">驱动基因</div>
    <div style="padding:8px;background:#fefce8;border-radius:6px;">通路富集</div>
    <div style="padding:8px;background:#fefce8;border-radius:6px;">专家解读</div>
  </div>
</div>
"""

DEMO_VIZ_DEG_TABLE = """
<table style="width:100%;border-collapse:collapse;font-family:ui-monospace,monospace;font-size:11px;">
  <thead><tr style="background:#1e293b;color:#f8fafc;">
  <th style="padding:6px 8px;">gene</th><th style="padding:6px 8px;">log2FC</th><th style="padding:6px 8px;">padj</th></tr></thead>
  <tbody style="color:#334155;">
  <tr style="background:#fef2f2;"><td style="padding:5px 8px;border-bottom:1px solid #e2e8f0;">TP53</td><td style="padding:5px 8px;border-bottom:1px solid #e2e8f0;">2.1</td><td style="padding:5px 8px;border-bottom:1px solid #e2e8f0;">1.2e-8</td></tr>
  <tr style="background:#fef2f2;"><td style="padding:5px 8px;border-bottom:1px solid #e2e8f0;">CDKN1A</td><td style="padding:5px 8px;border-bottom:1px solid #e2e8f0;">1.8</td><td style="padding:5px 8px;border-bottom:1px solid #e2e8f0;">3.4e-6</td></tr>
  <tr style="background:#eff6ff;"><td style="padding:5px 8px;border-bottom:1px solid #e2e8f0;">MYC</td><td style="padding:5px 8px;border-bottom:1px solid #e2e8f0;">-1.5</td><td style="padding:5px 8px;border-bottom:1px solid #e2e8f0;">0.002</td></tr>
  </tbody>
</table>
"""

DEMO_VIZ_WEEKLY_REPORT = """
<div style="font-family:system-ui;font-size:12px;padding:14px;background:#fff;border:1px solid #e2e8f0;border-radius:8px;">
  <div style="font-weight:600;margin-bottom:8px;">📋 周报结构预览</div>
  <div style="color:#64748b;margin-bottom:6px;">个人状态同步</div>
  <ul style="margin:0 0 12px 1.2em;padding:0;color:#334155;"><li>✅ 项目 A 里程碑完成</li><li>🟥 实验 B 需复盘</li></ul>
  <div style="color:#64748b;margin-bottom:6px;">团队对外同步</div>
  <ul style="margin:0 0 0 1.2em;padding:0;color:#334155;"><li>里程碑 / 风险 / 协作需求</li></ul>
</div>
"""

DEMO_VIZ_DEEP_RESEARCH = """
<div style="font-family:system-ui;font-size:12px;padding:14px;background:#fff;border:1px solid #e2e8f0;border-radius:8px;">
  <div style="font-weight:600;color:#1e40af;margin-bottom:8px;">Executive Summary</div>
  <p style="margin:0 0 10px;color:#475569;line-height:1.5;">间歇性禁食与代谢健康：益处、风险与适用人群…</p>
  <div style="display:flex;gap:8px;flex-wrap:wrap;">
    <span style="padding:4px 10px;background:#dbeafe;border-radius:999px;font-size:11px;">共识</span>
    <span style="padding:4px 10px;background:#fef3c7;border-radius:999px;font-size:11px;">争议</span>
    <span style="padding:4px 10px;background:#d1fae5;border-radius:999px;font-size:11px;">引用占位</span>
  </div>
</div>
"""

DEMO_VIZ_GO_KEGG = """
<table style="width:100%;border-collapse:collapse;font-family:system-ui;font-size:11px;">
  <thead><tr style="background:#f1f5f9;"><th style="padding:6px 8px;text-align:left;border:1px solid #e2e8f0;">Term</th>
  <th style="padding:6px 8px;text-align:left;border:1px solid #e2e8f0;">padj</th><th style="padding:6px 8px;text-align:left;border:1px solid #e2e8f0;">GeneRatio</th></tr></thead>
  <tbody>
  <tr><td style="padding:6px 8px;border:1px solid #e2e8f0;">GO:0006955 immune response</td><td style="padding:6px 8px;border:1px solid #e2e8f0;">2.1e-5</td><td style="padding:6px 8px;border:1px solid #e2e8f0;">12/180</td></tr>
  <tr><td style="padding:6px 8px;border:1px solid #e2e8f0;">KEGG:hsa04668 TNF signaling</td><td style="padding:6px 8px;border:1px solid #e2e8f0;">8.3e-4</td><td style="padding:6px 8px;border:1px solid #e2e8f0;">9/180</td></tr>
  </tbody>
</table>
"""

DEMO_VIZ_SC_CHECKLIST = """
<div style="font-family:system-ui;font-size:12px;padding:12px;background:#fff;border:1px solid #e2e8f0;border-radius:8px;">
  <label style="display:flex;gap:8px;margin-bottom:6px;"><input type="checkbox" checked disabled> 样本量 ≥ 3 生物学重复</label>
  <label style="display:flex;gap:8px;margin-bottom:6px;"><input type="checkbox" checked disabled> 批次效应评估计划</label>
  <label style="display:flex;gap:8px;margin-bottom:6px;"><input type="checkbox" disabled> 双胞检测工具选定</label>
  <label style="display:flex;gap:8px;"><input type="checkbox" disabled> 参考面板 / 质控阈值</label>
</div>
"""

DEMO_VIZ_ABSTRACT = """
<div style="font-family:Georgia,serif;font-size:12px;line-height:1.6;padding:14px;background:#fff;border:1px solid #e2e8f0;border-radius:8px;">
  <div style="font-size:11px;color:#64748b;margin-bottom:6px;font-family:system-ui;">精炼后 · 中英双语摘要</div>
  <p style="margin:0 0 8px;color:#1e293b;"><strong>EN:</strong> Single-cell RNA sequencing enables high-resolution profiling…</p>
  <p style="margin:0;color:#1e293b;"><strong>中:</strong> 单细胞 RNA 测序可实现高分辨率图谱构建…</p>
</div>
"""

DEMO_VIZ_CHEMBL = """
<table style="width:100%;border-collapse:collapse;font-family:system-ui;font-size:11px;">
  <thead><tr style="background:#1e40af;color:#fff;">
  <th style="padding:6px 8px;text-align:left;">CHEMBL ID</th><th style="padding:6px 8px;text-align:left;">名称</th><th style="padding:6px 8px;text-align:left;">类型</th></tr></thead>
  <tbody>
  <tr><td style="padding:6px 8px;border-bottom:1px solid #e2e8f0;">CHEMBL25</td><td style="padding:6px 8px;border-bottom:1px solid #e2e8f0;">Aspirin</td><td style="padding:6px 8px;border-bottom:1px solid #e2e8f0;">Small molecule</td></tr>
  <tr><td style="padding:6px 8px;border-bottom:1px solid #e2e8f0;">CHEMBL1200544</td><td style="padding:6px 8px;border-bottom:1px solid #e2e8f0;">Acetylsalicylic acid</td><td style="padding:6px 8px;border-bottom:1px solid #e2e8f0;">Drug</td></tr>
  </tbody>
</table>
"""

DEMO_VIZ_PUBMED = """
<table style="width:100%;border-collapse:collapse;font-family:system-ui;font-size:11px;">
  <thead><tr style="background:#f1f5f9;"><th style="padding:6px 8px;text-align:left;border:1px solid #e2e8f0;">PMID</th>
  <th style="padding:6px 8px;text-align:left;border:1px solid #e2e8f0;">Title</th><th style="padding:6px 8px;text-align:left;border:1px solid #e2e8f0;">Journal</th></tr></thead>
  <tbody>
  <tr><td style="padding:6px 8px;border-bottom:1px solid #e2e8f0;">38123456</td><td style="padding:6px 8px;border-bottom:1px solid #e2e8f0;">Single-cell atlas of T cell exhaustion…</td><td style="padding:6px 8px;border-bottom:1px solid #e2e8f0;">Nat Med</td></tr>
  <tr><td style="padding:6px 8px;border-bottom:1px solid #e2e8f0;">37890123</td><td style="padding:6px 8px;border-bottom:1px solid #e2e8f0;">scRNA-seq in immunotherapy response</td><td style="padding:6px 8px;border-bottom:1px solid #e2e8f0;">Cell</td></tr>
  </tbody>
</table>
"""

DEMO_VIZ_UNIPROT = """
<table style="width:100%;border-collapse:collapse;font-family:system-ui;font-size:11px;">
  <thead><tr style="background:#f1f5f9;"><th style="padding:6px 8px;text-align:left;border:1px solid #e2e8f0;">Entry</th>
  <th style="padding:6px 8px;text-align:left;border:1px solid #e2e8f0;">Protein</th><th style="padding:6px 8px;text-align:left;border:1px solid #e2e8f0;">Function</th></tr></thead>
  <tbody>
  <tr><td style="padding:6px 8px;border-bottom:1px solid #e2e8f0;">P04637</td><td style="padding:6px 8px;border-bottom:1px solid #e2e8f0;">Cellular tumor antigen p53</td><td style="padding:6px 8px;border-bottom:1px solid #e2e8f0;">Tumor suppressor</td></tr>
  </tbody>
</table>
"""

DEMO_VIZ_BEPIPRED = """
<div style="font-family:system-ui;font-size:11px;padding:8px;background:#fff;">
  <div style="font-weight:600;margin-bottom:8px;color:#334155;">BepiPred-3.0 Epitope Score（示意）</div>
  <div style="display:flex;align-items:flex-end;gap:2px;height:80px;border-bottom:1px solid #cbd5e1;padding-bottom:4px;">
    <div style="flex:1;background:linear-gradient(#6366f1,#a855f7);height:45%;border-radius:2px 2px 0 0;"></div>
    <div style="flex:1;background:linear-gradient(#6366f1,#eab308);height:72%;border-radius:2px 2px 0 0;"></div>
    <div style="flex:1;background:linear-gradient(#6366f1,#eab308);height:88%;border-radius:2px 2px 0 0;"></div>
    <div style="flex:1;background:linear-gradient(#6366f1,#a855f7);height:55%;border-radius:2px 2px 0 0;"></div>
    <div style="flex:1;background:linear-gradient(#6366f1,#6366f1);height:35%;border-radius:2px 2px 0 0;"></div>
  </div>
  <div style="margin-top:6px;color:#64748b;">— 阈值线（Top 20%）</div>
</div>
"""

DEMO_VIZ_SEQ_CONVERT = """
<div style="font-family:ui-monospace,monospace;font-size:10px;line-height:1.5;background:#f8fafc;padding:10px;border-radius:6px;">
  <div style="color:#64748b;margin-bottom:4px;">FASTA → GenBank</div>
  <div style="color:#1e40af;">LOCUS       demo_gene 42 bp ds-DNA linear UNK</div>
  <div>ORIGIN</div>
  <div>        1 atgcgtacgt tagctagcta gctagctagc tag</div>
  <div>//</div>
</div>
"""

DEMO_VIZ_MOL_CONVERT = """
<div style="font-family:ui-monospace,monospace;font-size:11px;padding:10px;background:#f8fafc;border-radius:6px;">
  <div style="color:#64748b;margin-bottom:4px;">SMILES → InChI</div>
  <div style="word-break:break-all;color:#1e293b;">InChI=1S/C9H8O4/c1-6(10)13-8-5-3-2-4-7(8)9(11)12/h2-5H,1H3,(H,11,12)</div>
</div>
"""

DEMO_VIZ_BLAST_HITS = """
<table style="width:100%;border-collapse:collapse;font-family:system-ui;font-size:11px;">
  <thead><tr style="background:#1e293b;color:#f8fafc;">
  <th style="padding:6px 8px;">Accession</th><th style="padding:6px 8px;">Description</th><th style="padding:6px 8px;">E-value</th></tr></thead>
  <tbody style="color:#334155;">
  <tr><td style="padding:5px 8px;border-bottom:1px solid #e2e8f0;">NM_001234</td><td style="padding:5px 8px;border-bottom:1px solid #e2e8f0;">Homo sapiens …</td><td style="padding:5px 8px;border-bottom:1px solid #e2e8f0;">2.3e-45</td></tr>
  <tr><td style="padding:5px 8px;border-bottom:1px solid #e2e8f0;">XM_567890</td><td style="padding:5px 8px;border-bottom:1px solid #e2e8f0;">Mus musculus …</td><td style="padding:5px 8px;border-bottom:1px solid #e2e8f0;">1.1e-12</td></tr>
  </tbody>
</table>
"""

DEMO_VIZ_RESUME = """
<div style="font-family:system-ui;font-size:12px;padding:12px;background:#fff;border:1px solid #e2e8f0;border-radius:8px;">
  <div style="font-weight:700;font-size:14px;">Jane Doe</div>
  <div style="color:#64748b;margin:4px 0 8px;">Senior Data Analyst · Healthcare</div>
  <div style="font-weight:600;margin-bottom:4px;">Summary</div>
  <div style="color:#475569;line-height:1.5;">5+ years SQL/Python · A/B testing · Tableau</div>
</div>
"""

DEMO_VIZ_COVER_LETTER = """
<div style="font-family:Georgia,serif;font-size:12px;line-height:1.6;padding:14px;background:#fff;border:1px solid #e2e8f0;border-radius:8px;color:#1e293b;">
  <p style="margin:0 0 8px;">Dear Editor,</p>
  <p style="margin:0 0 8px;">We submit our manuscript entitled <em>scRNA atlas under PD-1 blockade</em> for consideration as a Research Article in <strong>Nature Communications</strong>.</p>
  <p style="margin:0;color:#64748b;font-family:system-ui;font-size:11px;">（示意草稿 · 实际内容由技能生成）</p>
</div>
"""

DEMO_VIZ_SUBSTRUCTURE = """
<div style="font-family:system-ui;font-size:12px;padding:12px;background:#fff;border:1px solid #e2e8f0;border-radius:8px;">
  <div style="margin-bottom:8px;"><span style="padding:4px 8px;background:#dbeafe;border-radius:4px;">子结构: 苯环 c1ccccc1</span></div>
  <div><span style="padding:4px 8px;background:#d1fae5;border-radius:4px;">✓ c1ccccc1CCO 匹配</span></div>
</div>
"""

DEMO_VIZ_SMILES_CID = """
<table style="width:100%;border-collapse:collapse;font-family:system-ui;font-size:11px;">
  <thead><tr style="background:#f1f5f9;"><th style="padding:6px 8px;text-align:left;">CID</th><th style="padding:6px 8px;text-align:left;">Name</th></tr></thead>
  <tbody><tr><td style="padding:6px 8px;border-bottom:1px solid #e2e8f0;">2244</td><td style="padding:6px 8px;border-bottom:1px solid #e2e8f0;">Aspirin</td></tr></tbody>
</table>
"""
