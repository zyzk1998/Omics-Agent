#!/usr/bin/env python3
"""
Aromaticity-Perception - 分子芳香性感知工具

对化学分子的SMILES表达式执行芳香性感知操作，将凯库勒式（Kekulé form）
转换为芳香式（Aromatic form），识别分子结构中的芳香环体系。

Usage:
    python aromaticity_perception.py --smiles "C1=CC=CC=C1"
    python aromaticity_perception.py --smiles "C1=CC=C(C=C1)O" --output result.md
"""

import argparse
import sys


def aromaticity_perception(smiles: str) -> dict:
    """
    对输入的SMILES进行芳香性感知操作。

    Args:
        smiles: 输入的SMILES字符串（可为Kekulé式或混合形式）

    Returns:
        dict: 包含输入SMILES、芳香SMILES、状态等信息的字典
    """
    try:
        from rdkit import Chem
    except ImportError:
        return {
            "InputSMILES": smiles,
            "AromaticSMILES": "",
            "Status": "error",
            "Message": "RDKit not found. Please install it: pip install rdkit-pypi"
        }

    # 解析SMILES（RDKit会自动进行芳香性感知）
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return {
            "InputSMILES": smiles,
            "AromaticSMILES": "",
            "Status": "error",
            "Message": f"Invalid SMILES string: {smiles}"
        }

    try:
        # 输出芳香式SMILES
        aromatic_smiles = Chem.MolToSmiles(mol)

        return {
            "InputSMILES": smiles,
            "AromaticSMILES": aromatic_smiles,
            "Status": "success",
            "Message": "Aromaticity perception completed successfully"
        }

    except Exception as e:
        return {
            "InputSMILES": smiles,
            "AromaticSMILES": "",
            "Status": "error",
            "Message": f"Error during aromaticity perception: {str(e)}"
        }


def format_markdown(result: dict) -> str:
    """将结果格式化为Markdown字符串。"""
    if result["Status"] == "error":
        return f"""**Aromaticity Perception Result:**
* **Input SMILES:** `{result['InputSMILES']}`
* **Status:** ❌ Error
* **Message:** {result['Message']}
"""

    return f"""**Aromaticity Perception Result:**
* **Input SMILES:** `{result['InputSMILES']}`
* **Aromatic SMILES:** `{result['AromaticSMILES']}`
* **Status:** ✅ {result['Message']}
"""


def main():
    parser = argparse.ArgumentParser(
        description='Aromaticity-Perception: Convert Kekulé SMILES to Aromatic form'
    )
    parser.add_argument('--smiles', '-s', required=True,
                        help='Input SMILES string (e.g., "C1=CC=CC=C1")')
    parser.add_argument('--output', '-o', default=None,
                        help='Output file path (optional)')
    parser.add_argument('--format', '-f', default='markdown',
                        choices=['markdown', 'json'],
                        help='Output format (default: markdown)')

    args = parser.parse_args()

    # 执行芳香性感知
    result = aromaticity_perception(args.smiles)

    # 格式化输出
    if args.format == 'json':
        import json
        output = json.dumps(result, indent=2, ensure_ascii=False)
    else:
        output = format_markdown(result)

    # 输出到文件或控制台
    if args.output:
        with open(args.output, 'w', encoding='utf-8') as f:
            f.write(output)
        print(f"Result saved to: {args.output}")
    else:
        print(output)

    # 如果出错，返回非零退出码
    if result['Status'] == 'error':
        sys.exit(1)

    return 0


if __name__ == '__main__':
    sys.exit(main())
