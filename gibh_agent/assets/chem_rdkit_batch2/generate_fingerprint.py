#!/usr/bin/env python3
"""
Pattern-Fingerprint Generator

从分子的SMILES表达式生成模式指纹（Pattern Fingerprint），
并提供详细的比特向量信息。

Usage:
    python generate_fingerprint.py --smiles "CC(=O)O"
    python generate_fingerprint.py --smiles "CC(=O)O" --output result.json
    python generate_fingerprint.py --smiles "CC(=O)O" --format binary
"""

import argparse
import json
import sys


def generate_pattern_fingerprint(smiles: str, n_bits: int = 2048) -> dict:
    """
    从SMILES字符串生成分子的Pattern Fingerprint。

    Args:
        smiles: 分子的SMILES表达式
        n_bits: 指纹长度（默认2048）

    Returns:
        dict: 包含NumBits, NumOnBits, BitVector, Binary等信息的字典
    """
    try:
        from rdkit import Chem
        from rdkit.Chem import PatternFingerprint
    except ImportError:
        return {
            "NumBits": 0,
            "NumOnBits": 0,
            "BitVector": "",
            "Binary": "",
            "SMILES": smiles,
            "Status": "error",
            "Message": "RDKit not found. Please install it: pip install rdkit-pypi"
        }

    # 解析SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return {
            "NumBits": 0,
            "NumOnBits": 0,
            "BitVector": "",
            "Binary": "",
            "SMILES": smiles,
            "Status": "error",
            "Message": f"Invalid SMILES string: {smiles}"
        }

    try:
        # 生成Pattern Fingerprint
        fp = PatternFingerprint(mol, fpSize=n_bits)

        # 获取指纹信息
        num_bits = fp.GetNumBits()
        num_on_bits = fp.GetNumOnBits()

        # 生成BitVector（列表字符串形式）
        bit_vector = [1 if fp.GetBit(i) else 0 for i in range(num_bits)]
        bit_vector_str = str(bit_vector)

        # 生成Binary（二进制字符串形式）
        binary_str = ''.join(str(b) for b in bit_vector)

        return {
            "NumBits": num_bits,
            "NumOnBits": num_on_bits,
            "BitVector": bit_vector_str,
            "Binary": binary_str,
            "SMILES": smiles,
            "Status": "success",
            "Message": "Pattern fingerprint generated successfully"
        }

    except Exception as e:
        return {
            "NumBits": 0,
            "NumOnBits": 0,
            "BitVector": "",
            "Binary": "",
            "SMILES": smiles,
            "Status": "error",
            "Message": f"Error generating fingerprint: {str(e)}"
        }


def main():
    parser = argparse.ArgumentParser(
        description='Generate Pattern Fingerprint from SMILES string'
    )
    parser.add_argument('--smiles', '-s', required=True,
                        help='Input SMILES string (e.g., "CC(=O)O")')
    parser.add_argument('--output', '-o', default=None,
                        help='Output JSON file path (optional)')
    parser.add_argument('--format', '-f', default='json',
                        choices=['json', 'binary', 'vector'],
                        help='Output format (default: json)')
    parser.add_argument('--bits', '-b', type=int, default=2048,
                        help='Fingerprint size in bits (default: 2048)')

    args = parser.parse_args()

    # 生成指纹
    result = generate_pattern_fingerprint(args.smiles, n_bits=args.bits)

    # 根据格式输出
    if args.format == 'binary':
        output = result.get('Binary', '')
    elif args.format == 'vector':
        output = result.get('BitVector', '')
    else:
        output = json.dumps(result, indent=2, ensure_ascii=False)

    # 输出到文件或控制台
    if args.output:
        with open(args.output, 'w', encoding='utf-8') as f:
            f.write(output)
        print(f"Result saved to: {args.output}")
    else:
        print(output)

    # 如果出错，返回非零退出码
    if result.get('Status') == 'error':
        sys.exit(1)

    return 0


if __name__ == '__main__':
    sys.exit(main())
