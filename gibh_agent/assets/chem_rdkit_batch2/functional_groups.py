#!/usr/bin/env python3
"""
Functional Groups Identifier - 分子官能团识别工具

识别给定SMILES分子中的常见官能团。

Usage:
    python functional_groups.py --smiles "CC(=O)O"
    python functional_groups.py --smiles "CC(=O)OC1=CC=CC=C1C(=O)O" --format json
    python functional_groups.py --smiles "CCN(CC)CC" --format markdown
"""

import argparse
import json
import sys

try:
    from rdkit import Chem
    from rdkit.Chem import AllChem
    HAS_RDKIT = True
except ImportError:
    HAS_RDKIT = False


# 官能团 SMARTS 定义
# 按优先级排序（更具体的模式在前，避免重复检测）
FUNCTIONAL_GROUPS = [
    # 羧酸及其衍生物（最优先）
    ("Carboxylic acid", "[CX3](=O)[OX2H1]"),
    ("Ester", "[CX3](=O)[OX2H0]"),
    ("Amide", "[CX3](=O)[NX3]"),
    ("Acid chloride", "[CX3](=O)[Cl]"),
    ("Anhydride", "[CX3](=O)O[CX3](=O)"),
    
    # 羰基化合物（排除已识别的羧酸衍生物）
    ("Aldehyde", "[CX3H1](=O)[#6]"),
    # 酮：羰基碳连接两个碳，且不连接O/N/Cl/S（排除羧酸衍生物）
    ("Ketone", "[$([#6X3](=O)([#6])[#6]);!$([#6X3](=O)[O,N,Cl,S])]"),

    # 醇和酚（排除羧酸中的OH）
    ("Phenol", "[OX2H1]c1ccccc1"),
    ("Alcohol", "[$([OX2H1]);!$([OX2H1]C(=O));!$([OX2H1]c1ccccc1)]"),

    # 胺类（排除酰胺中的N）
    ("Primary amine", "[$([NX3H2]);!$([NX3H2]C=O)]"),
    ("Secondary amine", "[$([NX3H1]);!$([NX3H1]C=O)]"),
    ("Tertiary amine", "[$([NX3H0]);!$([NX3H0]C=O)]"),
    ("Quaternary ammonium", "[NX4+]"),
    ("Imine", "[NX2]=[CX3]"),
    ("Nitrile", "[CX1]#[NX1]"),
    ("Hydrazine", "[NX3][NX3]"),
    ("Azide", "[NX1]~[NX1]~[NX1]"),
    ("Nitro", "[$([NX3](=O)=O),$([NX3+](=O)[O-])]"),

    # 醚和缩醛
    ("Epoxide", "C1OC1"),
    ("Ether", "[$([OX2]([#6])[#6]);!$([OX2]C(=O))]"),
    ("Acetal", "[CX4]([OX2])([OX2])"),

    # 含硫官能团
    ("Thiol", "[SX2H1]"),
    ("Thioether", "[SX2]([#6])[#6]"),
    ("Disulfide", "[SX2][SX2]"),
    ("Sulfonic acid", "[SX4](=O)(=O)[OX2H1]"),
    ("Sulfone", "[SX4](=O)(=O)([#6])[#6]"),
    ("Sulfoxide", "[SX3](=O)([#6])[#6]"),

    # 含磷官能团
    ("Phosphate", "[PX4](=O)([OX2])([OX2])([OX2])"),
    ("Phosphonate", "[PX4](=O)([OX2])([OX2])[#6]"),

    # 卤代烃
    ("Fluoride", "[F][#6]"),
    ("Chloride", "[Cl][#6]"),
    ("Bromide", "[Br][#6]"),
    ("Iodide", "[I][#6]"),
    ("Perfluoroalkyl", "[F][CX4]([F])([F])"),

    # 不饱和结构
    ("Alkene", "[CX3]=[CX3]"),
    ("Alkyne", "[CX2]#[CX2]"),
    ("Allene", "[CX3]=C=[CX3]"),

    # 芳香结构
    ("Aromatic ring", "c1ccccc1"),
    ("Heteroaromatic", "[a;R]1[n,o,s][a;R][a;R][a;R][a;R]1"),

    # 其他
    ("Carbonate", "[OX2][CX3](=O)[OX2]"),
    ("Carbamate", "[NX3][CX3](=O)[OX2]"),
    ("Urea", "[NX3][CX3](=O)[NX3]"),
    ("Isocyanate", "[NX2]=[CX2]=[OX1]"),
    ("Isothiocyanate", "[NX2]=[CX2]=[SX1]"),
    ("Carbodiimide", "[NX2]=[CX2]=[NX2]"),
    ("Oxime", "[CX3]=[NX2][OX2H1]"),
    ("Hydrazone", "[CX3]=[NX2][NX3]"),
    ("Peroxide", "[OX2][OX2]"),
    ("Hydroperoxide", "[OX2H1][OX2]"),
    ("Boronic acid", "[BX3]([OX2H1])[OX2H1]"),
    ("Silanol", "[SiX4]([OX2H1])"),
    ("Siloxane", "[SiX4]([OX2][Si])"),
]


def identify_functional_groups(smiles: str) -> dict:
    """
    识别分子中的官能团。

    Args:
        smiles: SMILES字符串

    Returns:
        dict: 包含官能团列表、SMILES、统计信息等
    """
    result = {
        "smiles": smiles,
        "functional_groups": [],
        "count": 0,
        "Status": "success",
        "Message": "Functional groups identified successfully"
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

    # 对每个官能团进行 SMARTS 匹配
    found_groups = set()
    for group_name, smarts in FUNCTIONAL_GROUPS:
        try:
            pattern = Chem.MolFromSmarts(smarts)
            if pattern is None:
                continue
            matches = mol.GetSubstructMatches(pattern)
            if len(matches) > 0:
                found_groups.add(group_name)
        except Exception:
            continue

    # 排序并去重
    result["functional_groups"] = sorted(list(found_groups))
    result["count"] = len(found_groups)

    return result


def format_markdown(result: dict) -> str:
    """格式化为 Markdown"""
    if result["Status"] == "error":
        return f"❌ **Error:** {result['Message']}"

    lines = [
        f"**Functional Group Analysis**",
        "",
        f"* **SMILES:** `{result['smiles']}`",
        f"* **Total groups found:** {result['count']}",
        "",
        "**Identified Functional Groups:**",
        ""
    ]

    if result["functional_groups"]:
        for i, group in enumerate(result["functional_groups"], 1):
            lines.append(f"{i}. {group}")
    else:
        lines.append("_No common functional groups identified._")

    lines.append("")
    return "\n".join(lines)


def main():
    parser = argparse.ArgumentParser(
        description='Functional Groups Identifier - Identify functional groups in a molecule'
    )
    parser.add_argument('--smiles', '-s', required=True,
                        help='Input SMILES string')
    parser.add_argument('--format', '-f', choices=['json', 'markdown', 'compact'],
                        default='json', help='Output format')

    args = parser.parse_args()

    result = identify_functional_groups(args.smiles)

    if args.format == 'json':
        output = json.dumps(result, indent=2, ensure_ascii=False)
    elif args.format == 'markdown':
        output = format_markdown(result)
    else:  # compact
        if result["Status"] == "error":
            output = f"Error: {result['Message']}"
        else:
            groups = ", ".join(result["functional_groups"]) if result["functional_groups"] else "None"
            output = f"{result['smiles']} -> {groups}"

    print(output)
    return 0 if result["Status"] == "success" else 1


if __name__ == '__main__':
    sys.exit(main())
