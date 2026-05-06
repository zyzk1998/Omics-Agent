#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
GI Absorption Estimator - 胃肠道吸收能力估算（RDKit）

由技能包 gi-absorption.zip 同步；编排层经子进程调用以隔离 RDKit。
支持 --output-dir 时写出 chem_summary.json + structure.png（与 chem_rdkit 右栏一致）。

Usage:
    python gi_absorption.py --smiles "CC(C)C(=O)O" [--output-dir DIR]
    python gi_absorption.py --file /path/to/first_line_smiles.txt --output-dir DIR
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Any, Dict, Optional

try:
    from rdkit import Chem
    from rdkit.Chem import Descriptors, Draw, rdMolDescriptors

    HAS_RDKIT = True
except ImportError:
    HAS_RDKIT = False


def _read_first_smiles(path: Path) -> str:
    raw = path.read_text(encoding="utf-8", errors="replace")
    for line in raw.splitlines():
        s = line.strip()
        if not s or s.startswith("#"):
            continue
        return s.split()[0].strip()
    raise ValueError(f"文件中未找到 SMILES: {path}")


def _mol_png(mol: Chem.Mol, out_path: Path, size: tuple[int, int] = (420, 420)) -> None:
    img = Draw.MolToImage(mol, size=size)
    img.save(str(out_path))


def estimate_gi_absorption(smiles: str) -> dict:
    """
    估算分子的胃肠道吸收能力（TPSA + logP 经验阈值）。

    Returns:
        dict: logP、TPSA、is_high_gi_absorption、Status、Message
    """
    result: Dict[str, Any] = {
        "smiles": smiles,
        "logP": None,
        "TPSA": None,
        "is_high_gi_absorption": None,
        "Status": "success",
        "Message": "Estimation completed",
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
        tpsa = float(rdMolDescriptors.CalcTPSA(mol))
        logp = float(Descriptors.MolLogP(mol))

        result["TPSA"] = round(tpsa, 4)
        result["logP"] = round(logp, 4)
        result["is_high_gi_absorption"] = (tpsa <= 140) and (logp <= 5)

    except Exception as e:  # noqa: BLE001
        result["Status"] = "error"
        result["Message"] = f"Calculation error: {str(e)}"

    return result


def format_markdown(result: dict) -> str:
    """格式化为 Markdown（与上游 SKILL.md 示例一致）。"""
    if result["Status"] == "error":
        return f"❌ **Error:** {result['Message']}"

    tpsa = result["TPSA"]
    logp = result["logP"]
    is_high = result["is_high_gi_absorption"]

    if tpsa <= 60:
        tpsa_interpretation = "极低，渗透性优秀"
    elif tpsa <= 90:
        tpsa_interpretation = "较低，渗透性良好"
    elif tpsa <= 140:
        tpsa_interpretation = "中等，被动扩散尚可"
    else:
        tpsa_interpretation = "过高，渗透性差（可能依赖主动转运）"

    if logp <= 0:
        logp_interpretation = "极低，亲水性过强"
    elif logp <= 1:
        logp_interpretation = "较低，偏亲水"
    elif logp <= 3:
        logp_interpretation = "适中，平衡良好"
    elif logp <= 5:
        logp_interpretation = "中等偏高，仍可被动扩散"
    else:
        logp_interpretation = "过高，亲脂性过强（可能蓄积/结合蛋白）"

    lines = [
        "**GI Absorption Estimation**",
        "",
        f"* **SMILES:** `{result['smiles']}`",
        "",
        "| Property | Value | Threshold | Status |",
        "|----------|-------|-----------|--------|",
    ]

    tpsa_status = "✅" if tpsa <= 140 else "❌"
    logp_status = "✅" if logp <= 5 else "❌"
    overall = "✅ High absorption probability" if is_high else "⚠️ Low absorption probability"

    lines.append(
        f"| **TPSA** | {tpsa} Å² | ≤ 140 | {tpsa_status} {tpsa_interpretation} |"
    )
    lines.append(f"| **logP** | {logp} | ≤ 5 | {logp_status} {logp_interpretation} |")
    lines.append("")
    lines.append(f"**Overall:** {overall}")
    lines.append("")

    if is_high:
        lines.append(
            "该分子同时满足 TPSA ≤ 140 且 logP ≤ 5，推测具有**高胃肠道被动吸收可能性**。"
        )
    else:
        if tpsa > 140:
            lines.append("TPSA 过高，分子极性太强，可能难以通过被动扩散穿过细胞膜。")
        if logp > 5:
            lines.append(
                "logP 过高，分子亲脂性太强，可能在脂肪组织中过度蓄积或与血浆蛋白过度结合。"
            )
        lines.append("该分子**不满足**高吸收的经验规则，口服生物利用度可能较低。")

    lines.append("")
    return "\n".join(lines)


def main() -> int:
    parser = argparse.ArgumentParser(
        description="GI Absorption Estimator - Predict gastrointestinal absorption from SMILES"
    )
    parser.add_argument("--smiles", "-s", default="", help="Input SMILES string")
    parser.add_argument(
        "--file",
        "-F",
        default="",
        help="文本文件首行 SMILES（与 --smiles 二选一）",
    )
    parser.add_argument(
        "--format",
        "-f",
        choices=["json", "markdown"],
        default="json",
        help="Output format（stdout）",
    )
    parser.add_argument(
        "--output-dir",
        default="",
        help="若指定则写入 chem_summary.json、structure.png（成功时）",
    )

    args = parser.parse_args()

    smi = (args.smiles or "").strip()
    if (args.file or "").strip():
        fpth = Path(args.file.strip()).expanduser()
        if not fpth.is_file():
            print(json.dumps({"Status": "error", "Message": f"文件不存在: {fpth}"}, ensure_ascii=False), file=sys.stderr)
            return 1
        smi = _read_first_smiles(fpth)
    if not smi:
        print(
            json.dumps(
                {"Status": "error", "Message": "请提供 --smiles 或 --file"},
                ensure_ascii=False,
            ),
            file=sys.stderr,
        )
        return 1

    result = estimate_gi_absorption(smi)

    out_dir: Optional[Path] = None
    if (args.output_dir or "").strip():
        out_dir = Path(args.output_dir.strip()).resolve()
        out_dir.mkdir(parents=True, exist_ok=True)
        mol: Any = None
        if HAS_RDKIT and result.get("Status") == "success":
            mol = Chem.MolFromSmiles(smi)
        if mol is not None:
            can_smi = Chem.MolToSmiles(mol)
        else:
            can_smi = smi
        inner = {k: v for k, v in result.items() if k not in ("Status", "Message")}
        summary = {
            "tool": "gi_absorption",
            "canonical_smiles": can_smi,
            "result": inner,
            "status_line": result.get("Status"),
            "message_line": result.get("Message"),
        }
        (out_dir / "chem_summary.json").write_text(
            json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8"
        )
        if mol is not None:
            try:
                _mol_png(mol, out_dir / "structure.png")
            except Exception as exc:  # noqa: BLE001
                print(f"WARN: 2D 结构图写入失败: {exc}", file=sys.stderr)

    if args.format == "json":
        output = json.dumps(result, indent=2, ensure_ascii=False)
    else:
        output = format_markdown(result)

    print(output)
    return 0 if result["Status"] == "success" else 1


if __name__ == "__main__":
    sys.exit(main())
