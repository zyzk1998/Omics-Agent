#!/usr/bin/env python3
"""
Open Babel Skill - 分子文件格式转换与结构处理工具

支持超过110种化学文件格式的互转、结构预处理、3D构象生成及基础性质计算。

Usage:
    python openbabel_converter.py --input "CCO" --in-format smi --out-format mol --output ethanol.mol
    python openbabel_converter.py --input aspirin.sdf --out-format pdb --add-hydrogens --gen3d
    python openbabel_converter.py --input "CC(=O)OC1=CC=CC=C1C(=O)O" --in-format smi --properties
"""

import argparse
import json
import os
import subprocess
import sys
import tempfile
from pathlib import Path


def run_obabel(input_data: str, input_format: str = None, output_format: str = None,
               output_path: str = None, add_hydrogens: bool = False, 
               remove_hydrogens: bool = False, gen3d: bool = False,
               optimize: bool = False, properties: bool = False) -> dict:
    """
    调用 Open Babel 进行分子格式转换和处理。
    
    Args:
        input_data: 输入分子数据（文件路径或SMILES/InChI字符串）
        input_format: 输入格式（如 smi, mol, sdf, pdb, xyz, inchi）
        output_format: 输出格式（如 mol, pdb, xyz, sdf, smiles, inchi, svg, png）
        output_path: 输出文件路径（可选）
        add_hydrogens: 是否添加氢原子
        remove_hydrogens: 是否移除氢原子
        gen3d: 是否生成3D构象
        optimize: 是否进行几何优化
        properties: 是否计算分子性质
    
    Returns:
        dict: 包含转换结果、文件路径、性质信息等的字典
    """
    
    result = {
        "input": input_data,
        "input_format": input_format,
        "output_format": output_format,
        "Status": "success",
        "Message": "Conversion successful",
        "output_path": None,
        "output_content": None,
        "properties": {}
    }
    
    # 判断 input_data 是文件路径还是字符串
    is_file = os.path.isfile(input_data)
    
    try:
        # 构建 obabel 命令
        cmd = ["obabel"]
        
        if is_file:
            cmd.extend(["-i", input_format if input_format else guess_format(input_data)])
            cmd.append(input_data)
            
            # 判断是否为二进制输出
            is_binary_output = output_format in ["png", "jpg", "jpeg"]
            
            if output_path and is_binary_output:
                # 二进制格式直接输出到文件，不需要捕获stdout
                cmd.extend(["-O", output_path])
                process = subprocess.run(
                    cmd,
                    capture_output=True,
                    text=True,
                    timeout=60
                )
                if process.returncode != 0:
                    result["Status"] = "error"
                    result["Message"] = f"Open Babel error: {process.stderr.strip()}"
                    return result
                output_content = ""
            else:
                # 文本格式，捕获stdout
                process = subprocess.run(
                    cmd,
                    capture_output=True,
                    text=True,
                    timeout=60
                )
                if process.returncode != 0:
                    result["Status"] = "error"
                    result["Message"] = f"Open Babel error: {process.stderr.strip()}"
                    return result
                output_content = process.stdout.strip()
        else:
            # 从字符串通过stdin读取
            input_fmt = input_format if input_format else "smi"
            
            # 判断是否为二进制输出
            is_binary_output = output_format in ["png", "jpg", "jpeg"]
            
            if output_path and is_binary_output:
                # 二进制格式直接输出到文件
                cmd = ["obabel", "-i", input_fmt, "-O", output_path]
            else:
                cmd = ["obabel", "-i", input_fmt, "-o", output_format if output_format else "smi"]
            
            if add_hydrogens:
                cmd.append("-h")
            if remove_hydrogens:
                cmd.append("-d")
            if gen3d:
                cmd.append("--gen3d")
            if optimize:
                cmd.append("--minimize")
            
            if output_path and is_binary_output:
                # 二进制格式不需要text=True（stdout只有消息）
                process = subprocess.run(
                    cmd,
                    input=input_data + "\n",
                    capture_output=True,
                    text=True,
                    timeout=60
                )
                if process.returncode != 0:
                    result["Status"] = "error"
                    result["Message"] = f"Open Babel error: {process.stderr.strip()}"
                    return result
                output_content = ""
            else:
                process = subprocess.run(
                    cmd,
                    input=input_data + "\n",
                    capture_output=True,
                    text=True,
                    timeout=60
                )
                if process.returncode != 0:
                    result["Status"] = "error"
                    result["Message"] = f"Open Babel error: {process.stderr.strip()}"
                    return result
                output_content = process.stdout.strip()
        
        # 如果请求了性质计算
        if properties:
            result["properties"] = calculate_properties(input_data, input_format)
        
        # 如果指定了输出文件路径
        if output_path:
            if output_format in ["png", "jpg", "jpeg"]:
                # 二进制格式已通过 -O 直接输出到文件
                if os.path.exists(output_path):
                    result["output_path"] = os.path.abspath(output_path)
                    result["output_content"] = None
                else:
                    result["Status"] = "error"
                    result["Message"] = f"Failed to generate binary output file"
            else:
                with open(output_path, 'w', encoding='utf-8') as f:
                    f.write(output_content)
                result["output_path"] = os.path.abspath(output_path)
                result["output_content"] = None
        else:
            result["output_content"] = output_content
        
        return result
        
    except subprocess.TimeoutExpired:
        result["Status"] = "error"
        result["Message"] = "Open Babel execution timed out (60s)"
        return result
    except FileNotFoundError:
        result["Status"] = "error"
        result["Message"] = "Open Babel (obabel) not found. Please install: apt-get install openbabel"
        return result
    except Exception as e:
        result["Status"] = "error"
        result["Message"] = f"Error: {str(e)}"
        return result


def guess_format(filepath: str) -> str:
    """根据文件扩展名猜测格式"""
    ext = Path(filepath).suffix.lower()
    format_map = {
        '.mol': 'mol',
        '.sdf': 'sdf',
        '.pdb': 'pdb',
        '.xyz': 'xyz',
        '.cml': 'cml',
        '.smi': 'smi',
        '.smiles': 'smi',
        '.inchi': 'inchi',
        '.cdx': 'cdx',
        '.mol2': 'mol2',
        '.fasta': 'fasta',
        '.pqr': 'pqr',
        '.mmcif': 'mmcif',
        '.cif': 'cif',
    }
    return format_map.get(ext, 'mol')


def _parse_obabel_append_line(line: str, n_props: int) -> list:
    """
    解析 ``obabel -o smi --append ...`` 的首行输出。

    常见两种形态：
    1) 多列制表：SMILES 与每个性质各占一列（须按列对齐）；旧实现误用 ``parts[-1]`` 会把最后一列当成全部性质。
    2) 单列续行：仅一个制表符，第二段内以空格分隔多个性质。
    """
    line = (line or "").strip()
    if not line or n_props <= 0:
        return []
    if "\t" not in line:
        toks = line.split()
        return toks[1 : 1 + n_props] if len(toks) > 1 else []

    parts = line.split("\t")
    if len(parts) >= n_props + 1:
        return [parts[i + 1].strip() for i in range(n_props)]
    if len(parts) == 2:
        return parts[1].split()
    toks = line.split()
    return toks[1 : 1 + n_props] if len(toks) > 1 else []


def calculate_properties(input_data: str, input_format: str = None) -> dict:
    """计算分子基础性质（依赖系统 ``obabel`` 的 ``--append`` 描述符列）。"""
    props = {}
    
    try:
        is_file = os.path.isfile(input_data)
        fmt = input_format if input_format else (guess_format(input_data) if is_file else "smi")
        
        # 使用 obabel --append 计算性质
        # 支持的性质: molwt, logP, TPSA, HBA1, HBD, nrot 等
        properties_to_calc = ["molwt", "logP", "TPSA", "HBA1", "HBD", "nrot"]
        append_str = " ".join(properties_to_calc)
        
        if is_file:
            cmd = ["obabel", "-i", fmt, input_data, "-o", "smi", "--append", append_str]
        else:
            cmd = ["obabel", "-i", fmt, "-o", "smi", "--append", append_str]
        
        process = subprocess.run(
            cmd,
            input=(input_data + "\n") if not is_file else None,
            capture_output=True,
            text=True,
            timeout=30
        )
        
        if process.returncode == 0:
            output = process.stdout.strip()
            if output:
                first_line = output.splitlines()[0]
                vals = _parse_obabel_append_line(first_line, len(properties_to_calc))
                for i, prop_name in enumerate(properties_to_calc):
                    if i < len(vals):
                        props[prop_name] = vals[i]
        
    except FileNotFoundError:
        props["error"] = "obabel not found"
    except Exception as e:
        props["error"] = str(e)
    
    return props


def list_supported_formats() -> list:
    """列出 Open Babel 支持的格式"""
    try:
        process = subprocess.run(["obabel", "-L", "formats"], 
                                capture_output=True, text=True, timeout=10)
        if process.returncode == 0:
            formats = [line.strip() for line in process.stdout.split('\n') if line.strip()]
            return formats
    except:
        pass
    return []


def batch_convert(input_file: str, output_dir: str, output_format: str, 
                   options: dict = None) -> dict:
    """批量转换文件中的多个分子"""
    result = {
        "input_file": input_file,
        "output_dir": output_dir,
        "output_format": output_format,
        "Status": "success",
        "Message": "Batch conversion completed",
        "output_files": []
    }
    
    try:
        os.makedirs(output_dir, exist_ok=True)
        input_format = guess_format(input_file)
        
        # 使用 -m 选项进行多分子拆分
        cmd = [
            "obabel",
            "-i", input_format,
            input_file,
            "-o", output_format,
            "-m"  # 多个输出文件
        ]
        
        if options:
            if options.get('add_hydrogens'):
                cmd.append("-h")
            if options.get('gen3d'):
                cmd.append("--gen3d")
        
        # 设置输出目录前缀
        base_name = Path(input_file).stem
        output_pattern = os.path.join(output_dir, f"{base_name}_#." + 
                                     (output_format if len(output_format) <= 3 else output_format[:3]))
        cmd.extend(["-O", output_pattern])
        
        process = subprocess.run(cmd, capture_output=True, text=True, timeout=120)
        
        # 查找生成的文件
        for f in os.listdir(output_dir):
            if f.startswith(f"{base_name}_"):
                result["output_files"].append(os.path.join(output_dir, f))
        
        if not result["output_files"]:
            result["Message"] = "No output files generated"
            result["Status"] = "warning"
            
    except Exception as e:
        result["Status"] = "error"
        result["Message"] = str(e)
    
    return result


def main():
    parser = argparse.ArgumentParser(
        description='Open Babel Skill: Molecular format conversion and processing'
    )
    
    # 输入
    parser.add_argument('--input', '-i', required=True,
                        help='Input molecule (file path or SMILES/InChI string)')
    parser.add_argument('--in-format', '-if',
                        help='Input format (e.g., smi, mol, sdf, pdb, xyz, inchi)')
    
    # 输出
    parser.add_argument('--out-format', '-of',
                        help='Output format (e.g., mol, pdb, xyz, sdf, smiles, inchi, svg, png)')
    parser.add_argument('--output', '-o',
                        help='Output file path (optional)')
    
    # 处理选项
    parser.add_argument('--add-hydrogens', '-h_plus', action='store_true',
                        help='Add hydrogen atoms')
    parser.add_argument('--remove-hydrogens', action='store_true',
                        help='Remove hydrogen atoms')
    parser.add_argument('--gen3d', action='store_true',
                        help='Generate 3D coordinates')
    parser.add_argument('--optimize', action='store_true',
                        help='Optimize geometry')
    parser.add_argument('--properties', '-p', action='store_true',
                        help='Calculate molecular properties')
    
    # 批量转换
    parser.add_argument('--batch', action='store_true',
                        help='Batch mode for multi-molecule files')
    parser.add_argument('--output-dir',
                        help='Output directory for batch conversion')
    
    # 其他
    parser.add_argument('--list-formats', action='store_true',
                        help='List supported formats')
    parser.add_argument('--format', '-f', default='json',
                        choices=['json', 'markdown', 'raw'],
                        help='Output format (default: json)')
    
    args = parser.parse_args()
    
    # 列出支持的格式
    if args.list_formats:
        formats = list_supported_formats()
        if args.format == 'json':
            print(json.dumps({"formats": formats}, indent=2, ensure_ascii=False))
        else:
            print("Supported formats:")
            for f in formats:
                print(f"  - {f}")
        return 0
    
    # 批量转换
    if args.batch:
        if not args.output_dir:
            print("Error: --output-dir required for batch mode")
            return 1
        options = {
            'add_hydrogens': args.add_hydrogens,
            'gen3d': args.gen3d
        }
        result = batch_convert(args.input, args.output_dir, args.out_format, options)
    else:
        # 单分子转换
        result = run_obabel(
            input_data=args.input,
            input_format=args.in_format,
            output_format=args.out_format,
            output_path=args.output,
            add_hydrogens=args.add_hydrogens,
            remove_hydrogens=args.remove_hydrogens,
            gen3d=args.gen3d,
            optimize=args.optimize,
            properties=args.properties
        )
    
    # 格式化输出
    if args.format == 'json':
        output = json.dumps(result, indent=2, ensure_ascii=False)
    elif args.format == 'markdown':
        lines = ["**Open Babel Conversion Result:**"]
        lines.append(f"* **Input:** `{result['input']}`")
        lines.append(f"* **Status:** {'✅' if result['Status'] == 'success' else '❌'} {result['Message']}")
        if result.get('output_path'):
            lines.append(f"* **Output File:** `{result['output_path']}`")
        if result.get('output_content') and len(result['output_content']) < 1000:
            lines.append(f"* **Output:**")
            lines.append(f"```")
            lines.append(result['output_content'])
            lines.append(f"```")
        if result.get('properties'):
            lines.append(f"* **Properties:**")
            for k, v in result['properties'].items():
                lines.append(f"  - {k}: {v}")
        if result.get('output_files'):
            lines.append(f"* **Output Files:** {len(result['output_files'])} generated")
            for f in result['output_files'][:5]:
                lines.append(f"  - `{f}`")
        output = "\n".join(lines)
    else:
        # raw: 只输出内容或错误信息
        if result['Status'] == 'success' and result.get('output_content'):
            output = result['output_content']
        else:
            output = result['Message']
    
    print(output)
    
    return 0 if result['Status'] == 'success' else 1


if __name__ == '__main__':
    sys.exit(main())
