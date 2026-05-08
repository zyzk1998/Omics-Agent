#!/usr/bin/env python3
"""
Tanimoto Distance Matrix Calculator

计算一组分子之间的两两Tanimoto距离矩阵，基于Morgan指纹。

Usage:
    python tanimoto_matrix.py --smiles "CCO" "c1ccccc1" "CC(=O)O"
    python tanimoto_matrix.py --file molecules.txt
    python tanimoto_matrix.py --smiles "CCO" "CCN" "CCCN" --output matrix.json
"""

import argparse
import json
import sys
from typing import List, Tuple

try:
    from rdkit import Chem
    from rdkit.Chem import AllChem, DataStructs
    HAS_RDKIT = True
except ImportError:
    HAS_RDKIT = False
    print("Error: RDKit is required. Install with: pip install rdkit", file=sys.stderr)
    sys.exit(1)


def smiles_to_morgan_fingerprint(smiles: str, radius: int = 2, n_bits: int = 2048) -> DataStructs.ExplicitBitVect:
    """
    将SMILES字符串转换为Morgan指纹。
    
    Args:
        smiles: SMILES字符串
        radius: Morgan指纹半径（默认2）
        n_bits: 指纹位数（默认2048）
    
    Returns:
        Morgan指纹位向量
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Invalid SMILES string: {smiles}")
    
    fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius=radius, nBits=n_bits)
    return fp


def calculate_tanimoto_similarity(fp1: DataStructs.ExplicitBitVect, 
                                   fp2: DataStructs.ExplicitBitVect) -> float:
    """
    计算两个指纹之间的Tanimoto相似度。
    
    Tanimoto Similarity = (A ∩ B) / (A ∪ B)
    范围：0（完全不同）到1（完全相同）
    
    Args:
        fp1: 第一个指纹
        fp2: 第二个指纹
    
    Returns:
        Tanimoto相似度 (0-1)
    """
    return DataStructs.TanimotoSimilarity(fp1, fp2)


def calculate_tanimoto_distance_matrix(smiles_list: List[str], 
                                        radius: int = 2, 
                                        n_bits: int = 2048) -> Tuple[List[List[float]], List[str]]:
    """
    计算SMILES列表的Tanimoto距离矩阵。
    
    距离 = 1 - 相似度
    距离范围：0（结构相同）到1（结构完全不同）
    
    Args:
        smiles_list: SMILES字符串列表
        radius: Morgan指纹半径
        n_bits: 指纹位数
    
    Returns:
        (距离矩阵, 有效的SMILES列表)
    """
    # 生成指纹
    fingerprints = []
    valid_smiles = []
    
    for smiles in smiles_list:
        try:
            fp = smiles_to_morgan_fingerprint(smiles, radius, n_bits)
            fingerprints.append(fp)
            valid_smiles.append(smiles)
        except ValueError as e:
            print(f"Warning: {e}", file=sys.stderr)
    
    if len(fingerprints) == 0:
        raise ValueError("No valid SMILES strings provided")
    
    n = len(fingerprints)
    
    # 计算距离矩阵
    distance_matrix = [[0.0 for _ in range(n)] for _ in range(n)]
    
    for i in range(n):
        for j in range(i, n):
            if i == j:
                distance_matrix[i][j] = 0.0
            else:
                similarity = calculate_tanimoto_similarity(fingerprints[i], fingerprints[j])
                distance = 1.0 - similarity
                distance_matrix[i][j] = distance
                distance_matrix[j][i] = distance  # 对称矩阵
    
    return distance_matrix, valid_smiles


def format_matrix_markdown(smiles_list: List[str], 
                           distance_matrix: List[List[float]]) -> str:
    """
    将距离矩阵格式化为Markdown表格。
    
    Args:
        smiles_list: SMILES字符串列表
        distance_matrix: 距离矩阵
    
    Returns:
        Markdown格式字符串
    """
    n = len(smiles_list)
    
    # 构建分子索引表
    lines = []
    lines.append("**Molecule Index:**")
    lines.append("")
    lines.append("| Index | SMILES |")
    lines.append("|-------|--------|")
    for i, smiles in enumerate(smiles_list):
        lines.append(f"| {i} | `{smiles}` |")
    lines.append("")
    
    # 构建距离矩阵表
    lines.append("**Tanimoto Distance Matrix:**")
    lines.append("")
    
    # 表头
    header = "| | " + " | ".join(str(i) for i in range(n)) + " |"
    lines.append(header)
    lines.append("|" + "|".join(["-------"] * (n + 1)) + "|")
    
    # 数据行
    for i in range(n):
        row_values = [f"{distance_matrix[i][j]:.3f}" for j in range(n)]
        row = f"| **{i}** | " + " | ".join(row_values) + " |"
        lines.append(row)
    
    lines.append("")
    
    return "\n".join(lines)


def read_smiles_from_file(filepath: str) -> List[str]:
    """
    从文件读取SMILES列表。
    
    每行一个SMILES字符串，忽略空行和以#开头的注释行。
    
    Args:
        filepath: 文件路径
    
    Returns:
        SMILES字符串列表
    """
    smiles_list = []
    with open(filepath, 'r') as f:
        for line in f:
            line = line.strip()
            if line and not line.startswith('#'):
                smiles_list.append(line)
    return smiles_list


def main():
    parser = argparse.ArgumentParser(
        description='Calculate Tanimoto Distance Matrix from SMILES using Morgan Fingerprints'
    )
    parser.add_argument('--smiles', '-s', nargs='+', 
                        help='SMILES strings to analyze')
    parser.add_argument('--file', '-f', 
                        help='File containing SMILES strings (one per line)')
    parser.add_argument('--radius', '-r', type=int, default=2,
                        help='Morgan fingerprint radius (default: 2)')
    parser.add_argument('--n-bits', '-b', type=int, default=2048,
                        help='Number of bits in fingerprint (default: 2048)')
    parser.add_argument('--output', '-o', 
                        help='Output JSON file (optional)')
    parser.add_argument('--format', choices=['markdown', 'json', 'csv'], default='markdown',
                        help='Output format (default: markdown)')
    
    args = parser.parse_args()
    
    # 获取SMILES列表
    if args.smiles:
        smiles_list = args.smiles
    elif args.file:
        try:
            smiles_list = read_smiles_from_file(args.file)
        except FileNotFoundError:
            print(f"Error: File not found: {args.file}", file=sys.stderr)
            sys.exit(1)
    else:
        print("Error: Either --smiles or --file must be provided", file=sys.stderr)
        sys.exit(1)
    
    if len(smiles_list) < 2:
        print("Error: At least 2 SMILES strings are required", file=sys.stderr)
        sys.exit(1)
    
    if args.format != "json":
        print(f"Calculating Tanimoto distance matrix for {len(smiles_list)} molecules...")
        print(f"Morgan fingerprint: radius={args.radius}, n_bits={args.n_bits}")
        print()
    
    # 计算距离矩阵
    try:
        distance_matrix, valid_smiles = calculate_tanimoto_distance_matrix(
            smiles_list, 
            radius=args.radius, 
            n_bits=args.n_bits
        )
    except ValueError as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)
    
    # 输出结果
    if args.format == 'markdown':
        output = format_matrix_markdown(valid_smiles, distance_matrix)
        print(output)
    elif args.format == 'json':
        result = {
            "smiles_list": valid_smiles,
            "distance_matrix": distance_matrix,
            "parameters": {
                "radius": args.radius,
                "n_bits": args.n_bits
            }
        }
        output = json.dumps(result, indent=2)
        print(output)
    elif args.format == 'csv':
        # CSV格式
        n = len(valid_smiles)
        lines = []
        lines.append("," + ",".join(str(i) for i in range(n)))
        for i in range(n):
            row = f"{i}," + ",".join(f"{distance_matrix[i][j]:.6f}" for j in range(n))
            lines.append(row)
        output = "\n".join(lines)
        print(output)
    
    # 保存到文件（进度信息勿写入 stdout，以免破坏 JSON 管道）
    if args.output:
        with open(args.output, 'w', encoding='utf-8') as f:
            f.write(output)
        print(f"Output saved to: {args.output}", file=sys.stderr)
    
    return 0


if __name__ == '__main__':
    sys.exit(main())
