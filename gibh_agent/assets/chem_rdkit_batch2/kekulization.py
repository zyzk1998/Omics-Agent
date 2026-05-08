#!/usr/bin/env python3
"""
Kekulization - 分子凯库勒化转换工具

对化学分子的SMILES表达式执行Kekulization转换，
将芳香式（小写字母表示）转换为凯库勒式（显式单双键交替）。

Usage:
    python kekulization.py --smiles "c1ccccc1"
    python kekulization.py --smiles "c1ccc(cc1)O" --output result.md
"""

import argparse
import sys


def kekulization(smiles: str) -> dict:
    """
    对输入的SMILES进行Kekulization转换。

    Args:
        smiles: 输入的SMILES字符串（可为芳香式或凯库勒式）

    Returns:
        dict: 包含输入SMILES、凯库勒SMILES、状态等信息的字典
    """
    try:
        from rdkit import Chem
    except ImportError:
        return {
            "InputSMILES": smiles,
            "KekulizedSMILES": "",
            "Status": "error",
            "Message": "RDKit not found. Please install it: pip install rdkit-pypi"
        }

    # 解析SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return {
            "InputSMILES": smiles,
            "KekulizedSMILES": "",
            "Status": "error",
            "Message": f"Invalid SMILES string: {smiles}"
        }

    try:
        # Kekulize分子（将芳香键转为单双键交替）
        Chem.Kekulize(mol, clearAromaticFlags=True)
        
        # 输出凯库勒式SMILES
        kekulized_smiles = Chem.MolToSmiles(mol, kekuleSmiles=True)

        return {
            "InputSMILES": smiles,
            "KekulizedSMILES": kekulized_smiles,
            "Status": "success",
            "Message": "Kekulization completed successfully"
        }

    except Exception as e:
        # 如果Kekulization失败（可能已经是凯库勒式），尝试直接输出
        try:
            kekulized_smiles = Chem.MolToSmiles(mol, kekuleSmiles=True)
            return {
                "InputSMILES": smiles,
                "KekulizedSMILES": kekulized_smiles,
                "Status": "success",
                "Message": "Kekulization completed (molecule may already be in Kekulé form)"
            }
        except:
            return {
                "InputSMILES": smiles,
                "KekulizedSMILES": "",
                "Status": "error",
                "Message": f"Error during Kekulization: {str(e)}"
            }


def format_markdown(result: dict) -> str:
    """将结果格式化为Markdown字符串。"""
    if result["Status"] == "error":
        return f"""**Kekulization Result:**
* **Input SMILES:** `{result['InputSMILES']}`
* **Status:** ❌ Error
* **Message:** {result['Message']}
"""

    return f"""**Kekulization Result:**
* **Input SMILES:** `{result['InputSMILES']}`
* **Kekulized SMILES:** `{result['KekulizedSMILES']}`
* **Status:** ✅ {result['Message']}
"""


def main():
    parser = argparse.ArgumentParser(
        description='Kekulization: Convert Aromatic SMILES to Kekulé form'
    )
    parser.add_argument('--smiles', '-s', required=True,
                        help='Input SMILES string (e.g., "c1ccccc1")')
    parser.add_argument('--output', '-o', default=None,
                        help='Output file path (optional)')
    parser.add_argument('--format', '-f', default='markdown',
                        choices=['markdown', 'json'],
                        help='Output format (default: markdown)')

    args = parser.parse_args()

    # 执行Kekulization
    result = kekulization(args.smiles)

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
