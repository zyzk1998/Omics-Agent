#!/usr/bin/env python3
"""
Molecular-Image-Generation - 分子结构图像生成工具

根据SMILES表达式生成分子的二维结构图像（PNG格式），
并返回图像文件路径。

Usage:
    python generate_molecule_image.py --smiles "CC(=O)O"
    python generate_molecule_image.py --smiles "CC(=O)O" --output /tmp/molecule.png
    python generate_molecule_image.py --smiles "CC(=O)O" --size 800 --format json
"""

import argparse
import json
import os
import sys
import hashlib


def generate_molecule_image(smiles: str, output_path: str = None, 
                            size: int = 400, dpi: int = 150) -> dict:
    """
    根据SMILES字符串生成分子结构图像。

    Args:
        smiles: 输入的SMILES字符串
        output_path: 输出文件路径（可选，默认自动生成）
        size: 图像尺寸（像素，默认400x400）
        dpi: 图像DPI（默认150）

    Returns:
        dict: 包含smiles、image_path、status等信息的字典
    """
    try:
        from rdkit import Chem
        from rdkit.Chem import Draw
    except ImportError:
        return {
            "smiles": smiles,
            "molecule_image_url": "",
            "Status": "error",
            "Message": "RDKit not found. Please install it: pip install rdkit-pypi"
        }

    # 解析SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return {
            "smiles": smiles,
            "molecule_image_url": "",
            "Status": "error",
            "Message": f"Invalid SMILES string: {smiles}"
        }

    try:
        # 生成图像
        img = Draw.MolToImage(mol, size=(size, size), kekulize=True, 
                             fitImage=True, imageType='PNG')
        
        # 计算文件名（基于SMILES的哈希）
        smiles_hash = hashlib.md5(smiles.encode()).hexdigest()[:8]
        
        # 确定输出路径
        if output_path is None:
            # 默认输出到 /tmp/molecule_images/
            output_dir = "/tmp/molecule_images"
            os.makedirs(output_dir, exist_ok=True)
            output_path = os.path.join(output_dir, f"mol_{smiles_hash}.png")
        else:
            # 确保目录存在
            output_dir = os.path.dirname(output_path)
            if output_dir:
                os.makedirs(output_dir, exist_ok=True)
        
        # 保存图像
        img.save(output_path, "PNG", dpi=(dpi, dpi))
        
        # 获取绝对路径
        abs_path = os.path.abspath(output_path)
        
        return {
            "smiles": smiles,
            "molecule_image_url": f"file://{abs_path}",
            "image_path": abs_path,
            "image_size": f"{size}x{size}",
            "Status": "success",
            "Message": "Molecule image generated successfully"
        }

    except Exception as e:
        return {
            "smiles": smiles,
            "molecule_image_url": "",
            "Status": "error",
            "Message": f"Error generating image: {str(e)}"
        }


def main():
    parser = argparse.ArgumentParser(
        description='Molecular-Image-Generation: Generate 2D molecule structure image from SMILES'
    )
    parser.add_argument('--smiles', '-s', required=True,
                        help='Input SMILES string (e.g., "CC(=O)O")')
    parser.add_argument('--output', '-o', default=None,
                        help='Output image file path (optional)')
    parser.add_argument('--size', type=int, default=400,
                        help='Image size in pixels (default: 400)')
    parser.add_argument('--dpi', type=int, default=150,
                        help='Image DPI (default: 150)')
    parser.add_argument('--format', '-f', default='json',
                        choices=['json', 'markdown'],
                        help='Output format (default: json)')

    args = parser.parse_args()

    # 执行图像生成
    result = generate_molecule_image(
        smiles=args.smiles,
        output_path=args.output,
        size=args.size,
        dpi=args.dpi
    )

    # 格式化输出
    if args.format == 'json':
        output = json.dumps(result, indent=2, ensure_ascii=False)
    else:
        if result['Status'] == 'error':
            output = f"""**Molecular Image Generation Result:**
* **Input SMILES:** `{result['smiles']}`
* **Status:** ❌ Error
* **Message:** {result['Message']}
"""
        else:
            output = f"""**Molecular Image Generation Result:**
* **Input SMILES:** `{result['smiles']}`
* **Image Path:** `{result['image_path']}`
* **Image Size:** {result['image_size']}
* **Status:** ✅ {result['Message']}
"""

    print(output)

    # 如果出错，返回非零退出码
    if result['Status'] == 'error':
        sys.exit(1)

    return 0


if __name__ == '__main__':
    sys.exit(main())
