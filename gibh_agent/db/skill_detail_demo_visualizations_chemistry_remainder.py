# -*- coding: utf-8 -*-
"""化学大类 · demo_visualization 批量具象化（Phase 1.3）。"""
from __future__ import annotations

import json
from typing import Dict

from gibh_agent.db._demo_mol_block_aspirin import DEMO_ASPIRIN_MOL_BLOCK


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


def _smiles_block(label: str, smiles: str, color: str = "#475569") -> str:
    return (
        f'<div style="margin-bottom:6px;"><span style="color:#64748b;">{label}</span> '
        f'<code style="font-size:10px;word-break:break-all;color:{color};">{smiles}</code></div>'
    )


def _mol3d_viewer_demo(mol_block: str, height: str = "260px") -> str:
    """技能广场抽屉 · 可交互 3D 预览（与运行时 omics-3d-viewer-mountpoint 同构）。"""
    payload = json.dumps(mol_block)
    return (
        '<div class="skill-vis-mol3d-block" data-skill-mol3d="1" data-mol-format="mol_block">'
        f'<div class="skill-vis-mol3d-payload" style="display:none" aria-hidden="true">{payload}</div>'
        '<p style="margin:0 0 8px;color:#64748b;font-size:11px;">'
        "鼠标拖拽旋转 · 滚轮缩放 · ETKDG + MMFF 3D 构象</p>"
        '<div class="omics-3d-viewer-mountpoint mol-viewer-container skill-vis-mol3d-viewer" '
        'data-molecule-format="mol_block" '
        f'style="width:100%;height:{height};min-height:{height};border:1px solid #e2e8f0;'
        'border-radius:8px;background:#f9fafb;position:relative;box-sizing:border-box;"></div>'
        "</div>"
    )


CHEMISTRY_VIZ: Dict[str, str] = {
    "sascore_analysis": _report(
        "SA Score · 合成可行性评估",
        _table(
            ["分子", "SMILES", "SA Score", "解读"],
            [
                ["阿司匹林", "CC(=O)Oc1ccccc1C(=O)O", "1.52", "较易合成"],
                ["咖啡因", "CN1C=NC2=C1C(=O)N(C(=O)N2C)C", "2.89", "中等复杂度"],
                ["紫杉醇 (示意)", "—", "8.7", "高度复杂 · 需多步合成"],
            ],
        )
        + '<p style="margin:8px 0 0;color:#334155;"><strong>说明：</strong> SA Score 越低通常表示合成路线越可行（1–10 量表，基于分子图复杂度启发式）。</p>',
        "* 示意数值；实际结果由 RDKit SA_Score 模块计算。",
    ),
    "pymol_analysis": _report(
        "PyMOL 残基互作分析 · 演示 α-螺旋片段",
        _table(
            ["互作类型", "残基对", "距离 (Å)", "链"],
            [
                ["氢键", "ALA 3 O — ALA 7 N", "2.9", "A—A"],
                ["疏水接触", "ALA 5 CB — ALA 9 CB", "4.1", "A—A"],
                ["盐桥", "—", "—", "无（演示肽无带电侧链）"],
            ],
        )
        + '<div style="margin-top:8px;padding:10px;background:#f8fafc;border-radius:6px;font-size:11px;color:#475569;">'
        "Cartoon 渲染 · 14×Ala α-螺旋 · 非共价接触阈值 4.0 Å</div>",
    ),
    "chem_openbabel": _report(
        "Open Babel 理化性质 · 阿司匹林",
        _smiles_block("输入 SMILES", "CC(=O)Oc1ccccc1C(=O)O")
        + _table(
            ["属性", "值", "单位"],
            [
                ["分子量 (MW)", "180.16", "Da"],
                ["LogP", "1.19", "—"],
                ["TPSA", "63.6", "Å²"],
                ["HBD / HBA", "1 / 4", "计数"],
                ["可旋转键", "3", "—"],
                ["形式电荷", "0", "—"],
            ],
        ),
    ),
    "chem_lipinski_five": _report(
        "Lipinski 五规则 · 分子分析工具",
        _smiles_block("查询", "CC(=O)Oc1ccccc1C(=O)O")
        + _table(
            ["规则", "实测", "阈值", "判定"],
            [
                ["分子量 MW", "180.16 Da", "≤ 500", "✓ 通过"],
                ["LogP", "1.19", "≤ 5", "✓ 通过"],
                ["氢键供体 HBD", "1", "≤ 5", "✓ 通过"],
                ["氢键受体 HBA", "4", "≤ 10", "✓ 通过"],
            ],
        )
        + '<p style="margin:8px 0 0;color:#059669;"><strong>违规数：0</strong> · 符合口服类药性经验门槛</p>',
    ),
    "lipinski_druglikeness": _report(
        "Lipinski 类药性快筛 · 阿司匹林",
        _smiles_block("SMILES", "CC(=O)Oc1ccccc1C(=O)O")
        + _table(
            ["描述符", "值", "Lipinski", "Veber"],
            [
                ["MW", "180.16", "✓", "—"],
                ["LogP", "1.19", "✓", "—"],
                ["HBD", "1", "✓", "✓ (≤10)"],
                ["HBA", "4", "✓", "—"],
                ["TPSA", "63.6 Å²", "—", "✓ (≤140)"],
                ["可旋转键", "3", "—", "✓ (≤10)"],
            ],
        )
        + '<p style="margin:8px 0 0;color:#334155;">模式：<strong>Lipinski 五规则快筛</strong>（非高通量库检索）</p>',
    ),
    "rdkit_mol_format_convert": _report(
        "分子格式转换 · SMILES → InChI / Mol",
        _smiles_block("输入 (SMILES)", "CC(=O)Oc1ccccc1C(=O)O", "#1e40af")
        + '<div style="font-family:ui-monospace,monospace;font-size:10px;background:#f8fafc;padding:10px;border-radius:6px;margin:8px 0;">'
        "<div style='color:#64748b;margin-bottom:4px;'>InChI</div>"
        "<div style='word-break:break-all;color:#1e293b;'>InChI=1S/C9H8O4/c1-6(10)13-8-5-3-2-4-7(8)9(11)12/h2-5H,1H3,(H,11,12)</div>"
        "<div style='color:#64748b;margin:8px 0 4px;'>Mol 块 (节选)</div>"
        "<div style='color:#475569;'>  6  6  0  0  0  0  0  0  0  0999 V2000<br/>"
        "    1.2124    0.8250    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0<br/>"
        "M  END</div></div>"
        + _table(
            ["校验", "状态"],
            [["Sanitize", "通过"], ["芳香性", "已感知"]],
        ),
    ),
    "rdkit_3d_mol_render": _report(
        "3D 分子结构渲染 · 阿司匹林",
        _smiles_block("输入 SMILES", "CC(=O)Oc1ccccc1C(=O)O")
        + _mol3d_viewer_demo(DEMO_ASPIRIN_MOL_BLOCK)
        + _table(
            ["指标", "值"],
            [
                ["原子数 (含 H)", "21"],
                ["重原子数", "13"],
                ["输出字段", "data.mol_block"],
                ["前端格式", "mol_block → 3Dmol sdf"],
            ],
        ),
        "* 预览 MOL 块与 launch_skill_demos 默认可运行示例同源；打开详情抽屉后由 3Dmol.js 实时渲染。",
    ),
    "chem_molecule_image": _report(
        "分子 2D 结构图 · PNG 输出",
        _smiles_block("SMILES", "CC(=O)Oc1ccccc1C(=O)O")
        + '<div style="height:100px;background:#fff;border:1px solid #e2e8f0;border-radius:8px;display:flex;align-items:center;justify-content:center;margin:8px 0;">'
        '<svg width="140" height="70" viewBox="0 0 140 70" xmlns="http://www.w3.org/2000/svg">'
        '<polygon points="70,10 95,25 95,55 70,70 45,55 45,25" fill="none" stroke="#2563eb" stroke-width="1.5"/>'
        '<circle cx="70" cy="40" r="8" fill="#dbeafe" stroke="#1e40af"/>'
        '<text x="70" y="44" text-anchor="middle" font-size="8" fill="#1e40af">2D</text></svg></div>'
        + _table(
            ["参数", "值"],
            [["分辨率", "400 × 400 px"], ["DPI", "300"], ["存储", "MinIO 对象 URL"]],
        ),
    ),
    "chem_molmass": _report(
        "化学式分子量 · Hill 记法",
        _table(
            ["化学式", "平均质量", "单同位素质量", "元素组成"],
            [
                ["C6H12O6", "180.156 g/mol", "180.0634 Da", "C 6 · H 12 · O 6"],
                ["H2SO4", "98.079 g/mol", "97.9673 Da", "H 2 · S 1 · O 4"],
                ["C8H10N4O2", "194.194 g/mol", "194.0804 Da", "C 8 · H 10 · N 4 · O 2"],
            ],
        ),
        "* 输入为化学式而非 SMILES；由结构求质量请用「分子量计算工具」。",
    ),
    "calc_molecular_weight": _report(
        "分子量计算 · RDKit Descriptors",
        _smiles_block("阿司匹林", "CC(=O)Oc1ccccc1C(=O)O")
        + _table(
            ["指标", "值"],
            [
                ["分子式", "C9H8O4"],
                ["平均分子量", "180.16 Da"],
                ["精确质量", "180.0423 Da"],
                ["重原子数", "13"],
            ],
        ),
    ),
    "chem_molecular_weight": _report(
        "批量分子量 · 多 SMILES 输入",
        _table(
            ["#", "SMILES", "MW (Da)", "分子式"],
            [
                ["1", "CCO", "46.07", "C2H6O"],
                ["2", "c1ccccc1", "78.11", "C6H6"],
                ["3", "CC(=O)Oc1ccccc1C(=O)O", "180.16", "C9H8O4"],
            ],
        ),
    ),
    "chem_element_query": _report(
        "元素周期表查询 · 金 (Au)",
        _table(
            ["属性", "值"],
            [
                ["符号 / 名称", "Au · Gold · 金"],
                ["原子序数", "79"],
                ["标准原子量", "196.966570 Da"],
                ["周期 / 族", "6 / 11 (IB)"],
                ["电子构型", "[Xe] 4f14 5d10 6s1"],
                ["类别", "过渡金属 · 贵金属"],
            ],
        ),
    ),
    "chem_tanimoto_similarity": _report(
        "Tanimoto 相似度 · Morgan 指纹 (2048-bit)",
        _table(
            ["分子", "SMILES"],
            [["A (阿司匹林)", "CC(=O)Oc1ccccc1C(=O)O"], ["B (水杨酸)", "OC(=O)c1ccccc1O"]],
        )
        + '<p style="margin:10px 0;font-size:22px;color:#1e40af;text-align:center;"><strong>0.68</strong></p>'
        + _table(
            ["参数", "值"],
            [["指纹", "Morgan radius=2"], ["位数", "2048"], ["解读", "中等结构相似 · 共享芳环+羧酸"]],
        ),
    ),
    "chem_tanimoto_matrix": _report(
        "Tanimoto 距离矩阵 · 3 分子队列",
        _table(
            ["", "乙醇", "苯", "乙酸"],
            [
                ["乙醇 CCO", "0.00", "0.82", "0.71"],
                ["苯 c1ccccc1", "0.82", "0.00", "0.79"],
                ["乙酸 CC(=O)O", "0.71", "0.79", "0.00"],
            ],
        )
        + '<p style="margin:8px 0 0;color:#64748b;font-size:11px;">距离 = 1 − Tanimoto 相似度 · 值越大结构差异越大</p>',
    ),
    "chem_morgan_fingerprint": _report(
        "Morgan 指纹摘要 · 乙醇",
        _smiles_block("SMILES", "CCO")
        + _table(
            ["字段", "值"],
            [
                ["指纹长度", "2048 bits"],
                ["置位 (on-bits)", "12"],
                ["半径", "2 (ECFP4 类)"],
                ["哈希预览", "0x1a3f… · 0x8c02… · 0x44b1…"],
            ],
        )
        + '<div style="margin-top:6px;font-family:ui-monospace;font-size:9px;color:#64748b;word-break:break-all;">'
        "011000010010… (2048-bit 向量节选)</div>",
    ),
    "chem_pattern_fingerprint": _report(
        "Pattern 指纹 · 咖啡因",
        _smiles_block("SMILES", "CN1C=NC2=C1C(=O)N(C(=O)N2C)C")
        + _table(
            ["指标", "值"],
            [
                ["向量长度", "2048"],
                ["置位数", "28"],
                ["命中模式 (节选)", "嘌呤双环 · 酰胺 C=O · N-甲基"],
            ],
        ),
    ),
    "chem_functional_groups": _report(
        "官能团识别 · 阿司匹林",
        _smiles_block("SMILES", "CC(=O)Oc1ccccc1C(=O)O")
        + _table(
            ["官能团", "计数", "位置提示"],
            [
                ["羧酸 (Carboxylic acid)", "1", "芳环对位"],
                ["酯 (Ester)", "1", "乙酰氧基"],
                ["芳香环 (Aromatic)", "1", "苯环"],
                ["酮/酰基 (Acetyl)", "1", "—"],
            ],
        ),
    ),
    "chem_kekulization": _report(
        "Kekulization · 芳香键展开",
        _table(
            ["阶段", "SMILES"],
            [
                ["输入 (小写芳香)", "c1ccccc1"],
                ["输出 (Kekulé)", "C1=CC=CC=C1"],
            ],
        )
        + '<p style="margin:8px 0 0;color:#334155;">显式单/双键交替形式，便于部分导出流程与键级可视化。</p>',
    ),
    "chem_aromaticity_perception": _report(
        "芳香性感知 · 凯库勒 → 芳香记法",
        _table(
            ["阶段", "SMILES"],
            [
                ["输入 (Kekulé)", "C1=CC=CC=C1"],
                ["输出 (芳香)", "c1ccccc1"],
            ],
        )
        + '<p style="margin:8px 0 0;color:#334155;">统一环内共轭体系的简洁记法，便于数据库检索与表示一致。</p>',
    ),
    "chem_pains_filter": _report(
        "PAINS 筛查 · 泛测定干扰片段",
        _smiles_block("待测", "CC(=O)Oc1ccccc1C(=O)O")
        + _table(
            ["过滤器", "命中", "片段"],
            [["PAINS_A", "否", "—"], ["PAINS_B", "否", "—"], ["PAINS_C", "否", "—"]],
        )
        + '<p style="margin:8px 0 0;color:#059669;"><strong>结论：</strong> 未检出已知 PAINS 子结构 · 苗头质量可接受</p>',
    ),
    "chem_brenk_filter": _report(
        "Brenk 不良片段筛查",
        _smiles_block("待测", "CC(=O)NS(=O)(=O)c1ccc(N)cc1")
        + _table(
            ["警示", "命中", "说明"],
            [
                ["磺胺类 (sulfonamide)", "是", "可能干扰某些生物测定"],
                ["芳香胺", "是", "代谢活化风险需关注"],
                ["Brenk 总计", "2 条", "建议进一步 SAR 优化"],
            ],
        ),
    ),
    "chem_bbb_assessment": _report(
        "血脑屏障透过 · 启发式评估",
        _smiles_block("咖啡因", "CN1C=NC2=C1C(=O)N(C(=O)N2C)C")
        + _table(
            ["描述符", "值", "BBB 经验"],
            [
                ["LogP", "−0.07", "亲水性偏高"],
                ["TPSA", "61.8 Å²", "< 90 有利"],
                ["MW", "194.19 Da", "< 450 有利"],
                ["HBD", "0", "有利"],
            ],
        )
        + '<p style="margin:8px 0 0;"><span style="padding:3px 8px;background:#fef3c7;border-radius:4px;color:#92400e;">'
        "中等 CNS 暴露风险 · 非 PBPK 结论</span></p>",
    ),
    "chem_gi_absorption": _report(
        "胃肠道吸收 · 经验规则评估",
        _table(
            ["分子", "SMILES", "TPSA", "LogP", "吸收倾向"],
            [
                ["异丁酸", "CC(C)C(=O)O", "37.3", "0.79", "高"],
                ["阿司匹林", "CC(=O)Oc1ccccc1C(=O)O", "63.6", "1.19", "高"],
                ["扑热息痛", "CC(=O)Nc1ccc(O)cc1", "49.3", "0.46", "高"],
            ],
        ),
        "* 基于 TPSA/LogP 等启发式，非 Caco-2 或体内实验数据。",
    ),
    "chembl_drug_search": _report(
        "ChEMBL 药物检索 · aspirin",
        _table(
            ["CHEMBL ID", "首选名称", "类型", "最大阶段"],
            [
                ["CHEMBL25", "ASPIRIN", "Small molecule", "Approved"],
                ["CHEMBL1200544", "Acetylsalicylic acid", "Small molecule", "Approved"],
                ["CHEMBL941", "Salicylic acid", "Small molecule", "Phase 4"],
            ],
        )
        + _smiles_block("代表结构", "CC(=O)Oc1ccccc1C(=O)O", "#1e40af"),
    ),
    "chembl_similar_molecules": _report(
        "ChEMBL 相似分子 · 查询阿司匹林",
        _smiles_block("查询 SMILES", "CC(=O)Oc1ccccc1C(=O)O")
        + _table(
            ["CHEMBL ID", "相似度", "名称"],
            [
                ["CHEMBL25", "1.00", "Aspirin"],
                ["CHEMBL941", "0.72", "Salicylic acid"],
                ["CHEMBL145", "0.58", "2-Hydroxybenzoic acid methyl ester"],
            ],
        ),
    ),
    "drug_similarity_search": _report(
        "高通量结构相似性检索 · PubChem/ChEMBL",
        _smiles_block("查询", "CC(=O)Oc1ccccc1C(=O)O")
        + _table(
            ["Rank", "CID / CHEMBL", "相似度", "名称"],
            [
                ["1", "2244", "1.00", "Aspirin"],
                ["2", "338", "0.85", "Salicylic acid"],
                ["3", "8369", "0.62", "Methyl salicylate"],
                ["4", "CHEMBL941", "0.72", "Salicylic acid (ChEMBL)"],
            ],
        )
        + '<p style="margin:8px 0 0;color:#64748b;">阈值 0.7 · Top 20 · 需联网 worker-pyskills</p>',
    ),
    "drug_id_crossref": _report(
        "药物标识符交叉检索 · imatinib",
        _table(
            ["数据库", "标识符", "名称"],
            [
                ["ChEMBL", "CHEMBL941", "Imatinib"],
                ["PubChem CID", "5291", "Imatinib"],
                ["DrugBank", "DB00619", "Imatinib"],
                ["UNII", "8A1O1M485", "Imatinib mesylate"],
            ],
        ),
    ),
    "refmet_metabolite_search": _report(
        "RefMet 代谢物检索 · glucose",
        _table(
            ["RefMet ID", "名称", "类", "Formula"],
            [
                ["RM0135891", "Glucose", "Monosaccharides", "C6H12O6"],
                ["RM0135892", "D-Glucose", "Monosaccharides", "C6H12O6"],
                ["RM0136012", "Glucose 6-phosphate", "Sugar phosphates", "C6H13O9P"],
            ],
        ),
    ),
    "opentargets_chembl_hierarchy": _report(
        "OpenTargets · ChEMBL 层级 · CHEMBL941",
        _table(
            ["关系", "CHEMBL ID", "名称", "SMILES (节选)"],
            [
                ["父分子", "CHEMBL1201583", "Imatinib (parent)", "Cc1ccc(NC(=O)c2ccc(C)c…"],
                ["目标", "CHEMBL941", "Imatinib", "Cc1ccc(NC(=O)c2ccc(C)c…"],
                ["盐型", "CHEMBL1201584", "Imatinib mesylate", "—"],
            ],
        ),
    ),
    "fda_drug_label_search": _report(
        "FDA 药品标签字段 · aspirin",
        _table(
            ["字段", "内容 (节选)"],
            [
                ["brand_name", "BAYER ASPIRIN"],
                ["generic_name", "ASPIRIN"],
                ["indications_and_usage", "Pain reliever / fever reducer…"],
                ["warnings", "Reye's syndrome · GI bleeding…"],
                ["dosage_and_administration", "325–650 mg q4–6h PRN…"],
            ],
        ),
    ),
    "eis_drt_analysis": _report(
        "EIS → DRT 弛豫时间分布 · 示意",
        _table(
            ["峰 #", "τ (s)", "强度", "过程归属 (示意)"],
            [
                ["1", "1.2×10⁻⁵", "0.42", "电荷转移 / 双电层"],
                ["2", "3.8×10⁻³", "0.28", "扩散 / 传质"],
                ["3", "0.15", "0.11", "界面弛豫"],
            ],
        )
        + '<div style="margin-top:8px;height:70px;background:linear-gradient(180deg,#eff6ff,#fff);border:1px solid #e2e8f0;border-radius:6px;display:flex;align-items:flex-end;padding:0 12px 8px;gap:4px;">'
        '<div style="width:18%;height:55%;background:#2563eb;border-radius:2px 2px 0 0;"></div>'
        '<div style="width:18%;height:38%;background:#3b82f6;border-radius:2px 2px 0 0;"></div>'
        '<div style="width:18%;height:22%;background:#93c5fd;border-radius:2px 2px 0 0;"></div></div>'
        + '<p style="margin:6px 0 0;color:#64748b;font-size:11px;">输入：frequency, Z_real, Z_imag CSV · 正则 λ=0.1</p>',
    ),
    "rdkit_substructure_search": _report(
        "子结构匹配 · 苯环筛选",
        _table(
            ["子结构", "目标 SMILES", "匹配"],
            [
                ["c1ccccc1 (苯环)", "c1ccccc1CCO", "✓ 命中"],
                ["c1ccccc1", "CCO", "✗ 未命中"],
                ["c1ccccc1", "CC(=O)Oc1ccccc1C(=O)O", "✓ 命中"],
            ],
        ),
    ),
    "smiles_to_cid": _report(
        "SMILES → PubChem CID",
        _smiles_block("输入", "CC(=O)Oc1ccccc1C(=O)O")
        + _table(
            ["CID", "IUPAC / 常用名", "分子式", "InChIKey (节选)"],
            [
                ["2244", "2-acetyloxybenzoic acid (Aspirin)", "C9H8O4", "BSYNRYMUTXBXSQ…"],
            ],
        )
        + '<p style="margin:8px 0 0;color:#64748b;">歧义 SMILES 可返回多条 CID 候选列表</p>',
    ),
}
