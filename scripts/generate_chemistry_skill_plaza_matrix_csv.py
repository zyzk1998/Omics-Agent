#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
从内置「技能广场 ↔ 磐石矩阵」行数据生成标准 CSV（UTF-8-SIG），Excel/WPS 可直接打开为表格。

用法:
  python scripts/generate_chemistry_skill_plaza_matrix_csv.py
  python scripts/generate_chemistry_skill_plaza_matrix_csv.py -o ~/Desktop/chemistry_matrix.csv

说明:
  - tool_chain_key 列会自动去掉前导「.」；若需与磐石 bundle 完全一致，请以 gibh_agent/data/scienceone_toolchain_tools.json 为准。
  - 第 5 行广场侧为「SMILES 分子量」；与第 19 行「化学式分子量」区分，第 19 行 key 使用 molWeight。
"""
from __future__ import annotations

import argparse
import csv
from pathlib import Path
from typing import Any

# 按你 WPS/截图中的 20 行整理（磐石矩阵_技能名称 / 磐石矩阵_技能简介 与广场一致时可与下表相同）
ROWS: list[dict[str, Any]] = [
    {
        "技能广场_技能名称": "分子假阳性片段检测工具",
        "技能广场_UI子分类标签": "药物发现",
        "技能广场_简介": "使用 PAINS filter 检查分子中是否含有可能导致假阳性的片段。",
        "磐石矩阵_tool_chain_key": "molecularFalsePositiveFragmentDetectionTool",
    },
    {
        "技能广场_技能名称": "DRTtools",
        "技能广场_UI子分类标签": "数据分析",
        "技能广场_简介": "用于处理基于物理时效分布方法的电化学阻抗谱分析工具。",
        "磐石矩阵_tool_chain_key": "DRTtools",
    },
    {
        "技能广场_技能名称": "CP2K",
        "技能广场_UI子分类标签": "计算",
        "技能广场_简介": "主要用于原子模拟的分子动力学软件包，基于密度泛函理论。",
        "磐石矩阵_tool_chain_key": "CP2K",
    },
    {
        "技能广场_技能名称": "Tanimoto 距离矩阵计算工具",
        "技能广场_UI子分类标签": "数据处理",
        "技能广场_简介": "基于分子指纹 (Morgan Fingerprint) 计算 Tanimoto 距离矩阵。",
        "磐石矩阵_tool_chain_key": "molecularTanimotoDistanceMatrixCalculationTool",
    },
    {
        "技能广场_技能名称": "分子量计算工具",
        "技能广场_UI子分类标签": "性质预测",
        "技能广场_简介": "根据 SMILES 表达式计算分子量。",
        "磐石矩阵_tool_chain_key": "molecularMolWeightCalculationTool",
    },
    {
        "技能广场_技能名称": "分子 BBB 评估工具",
        "技能广场_UI子分类标签": "性质预测",
        "技能广场_简介": "预测分子是否可通过血脑屏障 (BBB)。",
        "磐石矩阵_tool_chain_key": "molecularBBBAssessmentTool",
    },
    {
        "技能广场_技能名称": "Brenk filter 分子毒性检查工具",
        "技能广场_UI子分类标签": "药物发现",
        "技能广场_简介": "使用 Brenk filter 检查分子中是否含有毒性或不良片段。",
        "磐石矩阵_tool_chain_key": "molecularBrenkFilterToxicityCheckTool",
    },
    {
        "技能广场_技能名称": "分子 Morgan Fingerprint 生成工具",
        "技能广场_UI子分类标签": "特征提取",
        "技能广场_简介": "生成分子的 Morgan Fingerprint (也叫环状粘性结构)。",
        "磐石矩阵_tool_chain_key": "molecularMorganFingerprintGenerationTool",
    },
    {
        "技能广场_技能名称": "分子 Pattern Fingerprint 生成工具",
        "技能广场_UI子分类标签": "特征提取",
        "技能广场_简介": "生成分子的 Pattern Fingerprint。",
        "磐石矩阵_tool_chain_key": "molecularPatternFingerprintGenerationTool",
    },
    {
        "技能广场_技能名称": "分子芳香性矫正操作工具",
        "技能广场_UI子分类标签": "结构处理",
        "技能广场_简介": "对分子进行芳香性矫正操作。",
        "磐石矩阵_tool_chain_key": "molecularAromaticityPerceptionOperationTool",
    },
    {
        "技能广场_技能名称": "分子 Kekulization 转换工具",
        "技能广场_UI子分类标签": "结构处理",
        "技能广场_简介": "对分子进行 Kekulization 转换。",
        "磐石矩阵_tool_chain_key": "molecularKekulizationConversionTool",
    },
    {
        "技能广场_技能名称": "分子图像生成工具",
        "技能广场_UI子分类标签": "数据可视化",
        "技能广场_简介": "根据 SMILES 生成分子的图像 (PNG)。",
        "磐石矩阵_tool_chain_key": "molecularImageGenerationTool",
    },
    {
        "技能广场_技能名称": "Open Babel",
        "技能广场_UI子分类标签": "格式转换",
        "技能广场_简介": "支持分子文件的批量转换、分子读写。",
        "磐石矩阵_tool_chain_key": "openBabel",
    },
    {
        "技能广场_技能名称": "分子相似性评估工具",
        "技能广场_UI子分类标签": "特征提取",
        "技能广场_简介": "计算两个分子的 Tanimoto 相似度。",
        "磐石矩阵_tool_chain_key": "molecularSimilarityAssessmentTool",
    },
    {
        "技能广场_技能名称": "化学元素查询",
        "技能广场_UI子分类标签": "基础查询",
        "技能广场_简介": "查询特定化学元素的各种属性信息。",
        "磐石矩阵_tool_chain_key": "chemicalElementQueryTool",
    },
    {
        "技能广场_技能名称": "分子官能团识别工具",
        "技能广场_UI子分类标签": "结构识别",
        "技能广场_简介": "识别分子中的常见官能团 (functional groups)。",
        "磐石矩阵_tool_chain_key": "molecularFunctionalGroupIdentificationTool",
    },
    {
        "技能广场_技能名称": "分子胃肠道吸收能力评估工具",
        "技能广场_UI子分类标签": "性质预测",
        "技能广场_简介": "估算分子在胃肠道的吸收能力。",
        "磐石矩阵_tool_chain_key": "molecularGastrointestinalAbsorptionCapacityAssessmentTool",
    },
    {
        "技能广场_技能名称": "药物相似性评估工具",
        "技能广场_UI子分类标签": "性质预测",
        "技能广场_简介": "根据 Lipinski Rule of Five 检验分子是否具有药物潜力。",
        "磐石矩阵_tool_chain_key": "drugSimilarityAssessmentTool",
    },
    {
        "技能广场_技能名称": "分子量计算",
        "技能广场_UI子分类标签": "计算",
        "技能广场_简介": "根据化学式计算分子质量和元素组成。",
        "磐石矩阵_tool_chain_key": "molWeight",
    },
    {
        "技能广场_技能名称": "LAMMPS",
        "技能广场_UI子分类标签": "动力学模拟",
        "技能广场_简介": "开源分子动力学模拟软件，常用于模拟液体、固体或气态粒子集合。",
        "磐石矩阵_tool_chain_key": "lammps",
    },
]

FIELDNAMES = [
    "序号",
    "技能广场_主分类",
    "技能广场_技能名称",
    "技能广场_UI子分类标签",
    "技能广场_简介",
    "磐石矩阵_tool_chain_key",
    "磐石矩阵_技能名称",
    "磐石矩阵_子分类",
    "磐石矩阵_技能简介",
]

# 截图中统一写 chemistryTools；若需与 JSON 中 DRTtools/CP2K 的 toolNames 一致，可改为按 key 分支
DEFAULT_PANSHI_SUBCAT = "化学工具矩阵(chemistryTools)"


def build_output_rows(same_as_plaza_for_panshi: bool = True) -> list[dict[str, str]]:
    out: list[dict[str, str]] = []
    for i, r in enumerate(ROWS, start=1):
        pname = r["技能广场_技能名称"]
        pintro = r["技能广场_简介"]
        row = {
            "序号": str(i),
            "技能广场_主分类": "化学",
            "技能广场_技能名称": pname,
            "技能广场_UI子分类标签": r["技能广场_UI子分类标签"],
            "技能广场_简介": pintro,
            "磐石矩阵_tool_chain_key": r["磐石矩阵_tool_chain_key"],
            "磐石矩阵_技能名称": pname if same_as_plaza_for_panshi else pname,
            "磐石矩阵_子分类": DEFAULT_PANSHI_SUBCAT,
            "磐石矩阵_技能简介": pintro if same_as_plaza_for_panshi else pintro,
        }
        out.append(row)
    return out


def write_csv(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8-sig", newline="") as f:
        w = csv.DictWriter(f, fieldnames=FIELDNAMES, quoting=csv.QUOTE_MINIMAL)
        w.writeheader()
        w.writerows(rows)


def main() -> None:
    ap = argparse.ArgumentParser(description="生成技能广场-磐石矩阵对照 CSV")
    ap.add_argument(
        "-o",
        "--output",
        type=Path,
        default=Path(__file__).resolve().parents[1] / "docs" / "chemistry_skill_plaza_panshi_matrix.csv",
        help="输出 CSV 路径",
    )
    args = ap.parse_args()
    rows = build_output_rows()
    write_csv(args.output, rows)
    print(f"已写入 {args.output} ，共 {len(rows)} 行。")


if __name__ == "__main__":
    main()
