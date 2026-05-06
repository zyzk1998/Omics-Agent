#!/usr/bin/env python3
"""
MolMass Calculator - 分子质量计算工具

基于 molmass 库计算分子质量、元素组成和同位素分布谱。

Usage:
    python molmass_cli.py --formula C6H12O6
    python molmass_cli.py --formula C8H10N4O2 --format markdown
    python molmass_cli.py --formula H2SO4 --spectrum
"""

import argparse
import json
import sys
from pathlib import Path

# 将 molmass_lib 加入路径
sys.path.insert(0, str(Path(__file__).parent))
from molmass_lib.molmass import Formula


def calculate_molmass(formula_str: str, include_spectrum: bool = False) -> dict:
    """
    计算分子质量相关信息。

    Args:
        formula_str: 化学式字符串，如 "C6H12O6"
        include_spectrum: 是否包含同位素分布谱

    Returns:
        dict: 包含分子量、元素组成等信息
    """
    result = {
        "formula": formula_str,
        "Status": "success",
        "Message": "Calculation completed"
    }

    try:
        f = Formula(formula_str)

        # 基本信息
        result["hill_notation"] = f.formula
        result["empirical_formula"] = f.empirical
        result["atom_count"] = f.atoms
        result["charge"] = f.charge

        # 质量
        result["average_mass"] = round(f.mass, 6)
        result["nominal_mass"] = int(f.nominal_mass)
        result["monoisotopic_mass"] = round(f.monoisotopic_mass, 6)

        # 元素组成
        comp = f.composition()
        composition_list = []
        for symbol, item in comp.items():
            composition_list.append({
                "element": symbol,
                "count": int(item.count),
                "relative_mass": round(float(item.mass), 6),
                "fraction": round(float(item.fraction), 6)
            })
        result["elemental_composition"] = composition_list

        # 同位素分布谱
        if include_spectrum:
            spec = f.spectrum(min_intensity=0.001)
            spectrum_list = []
            for massnumber, entry in spec.items():
                spectrum_list.append({
                    "mass_number": int(massnumber),
                    "relative_mass": round(float(entry.mass), 6),
                    "fraction": round(float(entry.fraction), 6),
                    "intensity_percent": round(float(entry.intensity), 4)
                })
            result["isotope_spectrum"] = spectrum_list

    except Exception as e:
        result["Status"] = "error"
        result["Message"] = str(e)

    return result


def format_markdown(result: dict) -> str:
    """格式化为 Markdown"""
    if result["Status"] == "error":
        return f"❌ **Error:** {result['Message']}"

    lines = [
        f"**Molecular Mass Analysis**",
        "",
        f"* **Formula:** `{result['formula']}`",
        f"* **Hill Notation:** `{result['hill_notation']}`",
        f"* **Empirical Formula:** `{result['empirical_formula']}`",
        f"* **Atom Count:** {result['atom_count']}",
        "",
        "**Mass Data:**",
        "",
        "| Property | Value | Unit |",
        "|----------|-------|------|",
        f"| **Average Mass** | {result['average_mass']} | g/mol |",
        f"| **Nominal Mass** | {result['nominal_mass']} | Da |",
        f"| **Monoisotopic Mass** | {result['monoisotopic_mass']} | Da |",
        "",
        "**Elemental Composition:**",
        "",
        "| Element | Count | Relative Mass | Fraction | Mass % |",
        "|---------|-------|---------------|----------|--------|",
    ]

    for comp in result["elemental_composition"]:
        mass_percent = round(comp["fraction"] * 100, 2)
        lines.append(
            f"| {comp['element']} | {comp['count']} | {comp['relative_mass']} | "
            f"{comp['fraction']} | {mass_percent}% |"
        )

    # 同位素谱
    if "isotope_spectrum" in result and result["isotope_spectrum"]:
        lines.append("")
        lines.append("**Isotope Distribution Spectrum:**")
        lines.append("")
        lines.append("| Mass Number | Relative Mass | Fraction | Intensity % |")
        lines.append("|-------------|---------------|----------|-------------|")
        for peak in result["isotope_spectrum"]:
            lines.append(
                f"| {peak['mass_number']} | {peak['relative_mass']} | "
                f"{peak['fraction']} | {peak['intensity_percent']}% |"
            )

    lines.append("")
    return "\n".join(lines)


def main():
    parser = argparse.ArgumentParser(
        description='MolMass Calculator - Molecular mass and composition analysis'
    )
    parser.add_argument('--formula', '-f', required=True,
                        help='Chemical formula (e.g., C6H12O6, C8H10N4O2)')
    parser.add_argument('--format', '-fmt', choices=['json', 'markdown'],
                        default='json', help='Output format')
    parser.add_argument('--spectrum', '-s', action='store_true',
                        help='Include isotope distribution spectrum')

    args = parser.parse_args()

    result = calculate_molmass(args.formula, include_spectrum=args.spectrum)

    if args.format == 'json':
        output = json.dumps(result, indent=2, ensure_ascii=False)
    else:
        output = format_markdown(result)

    print(output)
    return 0 if result["Status"] == "success" else 1


if __name__ == '__main__':
    sys.exit(main())
