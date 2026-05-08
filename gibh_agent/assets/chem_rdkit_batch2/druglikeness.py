#!/usr/bin/env python3
"""
Druglikeness Checker - 药物相似性评估工具

基于 Lipinski's Rule of Five (五规则) 评估分子的药物相似性：
  - MW ≤ 500 Da
  - logP ≤ 5
  - HBD ≤ 5
  - HBA ≤ 10
违反 ≤ 1 条 → 具有药物相似性

Usage:
    python druglikeness.py --smiles "CCO"
    python druglikeness.py --smiles "CC(=O)OC1=CC=CC=C1C(=O)O" --format json
    python druglikeness.py --smiles "CC(C)CC1=CC=C(C=C1)C(C)C(=O)O" --format markdown
"""

import argparse
import json
import sys

try:
    from rdkit import Chem
    from rdkit.Chem import Descriptors
    HAS_RDKIT = True
except ImportError:
    HAS_RDKIT = False


# Lipinski Rule of Five 阈值
LIPINSKI_THRESHOLDS = {
    "MW": 500.0,
    "logP": 5.0,
    "HBD": 5,
    "HBA": 10
}


def check_druglikeness(smiles: str) -> dict:
    """
    基于 Lipinski 五规则检查分子的药物相似性。

    Args:
        smiles: SMILES 字符串

    Returns:
        dict: 包含各参数、违规数、药物相似性判断
    """
    result = {
        "smiles": smiles,
        "MW": None,
        "logP": None,
        "HBD": None,
        "HBA": None,
        "violations": [],
        "violation_count": None,
        "is_druglike": None,
        "Status": "success",
        "Message": "Druglikeness assessment completed"
    }

    if not HAS_RDKIT:
        result["Status"] = "error"
        result["Message"] = "RDKit not found. Please install: pip install rdkit"
        return result

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        result["Status"] = "error"
        result["Message"] = f"Invalid SMILES string: {smiles}"
        return result

    try:
        # 计算各项参数
        mw = Descriptors.MolWt(mol)
        logp = Descriptors.MolLogP(mol)
        hbd = Descriptors.NumHDonors(mol)
        hba = Descriptors.NumHAcceptors(mol)

        result["MW"] = round(mw, 2)
        result["logP"] = round(logp, 4)
        result["HBD"] = int(hbd)
        result["HBA"] = int(hba)

        # 检查违规
        violations = []
        if mw > 500:
            violations.append(f"MW ({round(mw, 2)}) > 500 Da")
        if logp > 5:
            violations.append(f"logP ({round(logp, 4)}) > 5")
        if hbd > 5:
            violations.append(f"HBD ({int(hbd)}) > 5")
        if hba > 10:
            violations.append(f"HBA ({int(hba)}) > 10")

        result["violations"] = violations
        result["violation_count"] = len(violations)
        # Lipinski规则：违反 ≤ 1 条视为具有药物相似性
        result["is_druglike"] = len(violations) <= 1

    except Exception as e:
        result["Status"] = "error"
        result["Message"] = f"Calculation error: {str(e)}"

    return result


def format_markdown(result: dict) -> str:
    """格式化为 Markdown"""
    if result["Status"] == "error":
        return f"❌ **Error:** {result['Message']}"

    mw = result["MW"]
    logp = result["logP"]
    hbd = result["HBD"]
    hba = result["HBA"]
    violations = result["violations"]
    is_druglike = result["is_druglike"]

    lines = [
        f"**Druglikeness Assessment (Lipinski's Rule of Five)**",
        "",
        f"* **SMILES:** `{result['smiles']}`",
        "",
        "| Parameter | Value | Threshold | Status |",
        "|-----------|-------|-----------|--------|",
    ]

    mw_ok = "✅" if mw <= 500 else "❌"
    logp_ok = "✅" if logp <= 5 else "❌"
    hbd_ok = "✅" if hbd <= 5 else "❌"
    hba_ok = "✅" if hba <= 10 else "❌"

    lines.append(f"| **MW** | {mw} Da | ≤ 500 | {mw_ok} |")
    lines.append(f"| **logP** | {logp} | ≤ 5 | {logp_ok} |")
    lines.append(f"| **HBD** | {hbd} | ≤ 5 | {hbd_ok} |")
    lines.append(f"| **HBA** | {hba} | ≤ 10 | {hba_ok} |")
    lines.append("")

    lines.append(f"**Violations:** {len(violations)}/4")
    if violations:
        lines.append("")
        for v in violations:
            lines.append(f"- ❌ {v}")
    else:
        lines.append("✅ 无违规")

    lines.append("")
    if is_druglike:
        lines.append(f"**判定:** ✅ **具有药物相似性**（违反规则 ≤ 1 条）")
        if not violations:
            lines.append("该分子完全符合 Lipinski 五规则，是典型的类药分子。")
    else:
        lines.append(f"**判定:** ❌ **不符合药物相似性**（违反规则 > 1 条）")
        lines.append("该分子可能面临口服吸收差、溶解度低或药代动力学问题。")

    lines.append("")
    return "\n".join(lines)


def main():
    parser = argparse.ArgumentParser(
        description='Druglikeness Checker - Lipinski Rule of Five assessment'
    )
    parser.add_argument('--smiles', '-s', required=True,
                        help='Input SMILES string')
    parser.add_argument('--format', '-f', choices=['json', 'markdown'],
                        default='json', help='Output format')

    args = parser.parse_args()

    result = check_druglikeness(args.smiles)

    if args.format == 'json':
        output = json.dumps(result, indent=2, ensure_ascii=False)
    else:
        output = format_markdown(result)

    print(output)
    return 0 if result["Status"] == "success" else 1


if __name__ == '__main__':
    sys.exit(main())
