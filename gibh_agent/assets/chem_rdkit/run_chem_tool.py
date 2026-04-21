#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
RDKit 化学原子工具统一子进程入口（BBB / PAINS / Brenk / Morgan / MolWt / Tanimoto）。
随仓库分发；由 gibh_agent.tools.chem_rdkit_tools 以 subprocess 调用，隔离 RDKit 崩溃。

用法:
  python run_chem_tool.py --tool <bbb|pains|brenk|morgan_fp|mol_weight|tanimoto> \\
      [--smiles S] [--file PATH] [--smiles2 S] --output-dir DIR

环境:
  CHEM_STAGING_SCRIPT 若指向替代脚本则不应来到此处（由上层选用）。
"""
from __future__ import annotations

import argparse
import json
import sys
import traceback
from pathlib import Path
import hashlib
from typing import Any, Dict, List, Tuple

# ---------------------------------------------------------------------------
# RDKit（须在解释器环境中可用；API 镜像通过 rdkit-pypi 安装）
# ---------------------------------------------------------------------------
from rdkit import Chem
from rdkit.Chem import AllChem, DataStructs, Descriptors, Draw
from rdkit.Chem.FilterCatalog import FilterCatalog, FilterCatalogParams


def _stderr(msg: str) -> None:
    print(msg, file=sys.stderr)


def _filter_catalog_hit_description(h: Any) -> str:
    """兼容不同 RDKit 版本：GetMatches 可能返回含 filterMatch 的对象，或直接为 FilterCatalogEntry。"""
    fm = getattr(h, "filterMatch", None)
    if fm is not None and hasattr(fm, "GetDescription"):
        return str(fm.GetDescription())
    if hasattr(h, "GetDescription"):
        return str(h.GetDescription())
    return str(h)


def _read_smiles_from_file(path: Path) -> str:
    raw = path.read_text(encoding="utf-8", errors="replace")
    for line in raw.splitlines():
        s = line.strip()
        if not s or s.startswith("#"):
            continue
        return s.split()[0].strip()
    raise ValueError(f"文件中未找到 SMILES: {path}")


def _read_two_smiles(path: Path) -> Tuple[str, str]:
    raw = path.read_text(encoding="utf-8", errors="replace")
    lines = [ln.strip() for ln in raw.splitlines() if ln.strip() and not ln.strip().startswith("#")]
    if len(lines) >= 2:
        a = lines[0].split()[0].strip()
        b = lines[1].split()[0].strip()
        return a, b
    raise ValueError("需要至少两行 SMILES（或单列 CSV 前两行）")


def _mol_png(mol: Chem.Mol, out_path: Path, size: Tuple[int, int] = (420, 420)) -> None:
    img = Draw.MolToImage(mol, size=size)
    img.save(str(out_path))


def _bbb_assessment(mol: Chem.Mol) -> Dict[str, Any]:
    """启发式 BBB 渗透倾向（Descriptors；非临床 PBPK 模型，仅供筛选参考）。"""
    mw = float(Descriptors.MolWt(mol))
    logp = float(Descriptors.MolLogP(mol))
    tpsa = float(Descriptors.TPSA(mol))
    hbd = int(Descriptors.NumHDonors(mol))
    # 文献常见 CNS 渗透有利区（多规则综合的简化打分）
    score = 0.0
    reasons: List[str] = []
    if mw <= 500:
        score += 0.25
    else:
        reasons.append(f"MW 偏大 ({mw:.1f} > 500)")
    if logp <= 5 and logp >= -0.4:
        score += 0.25
    else:
        reasons.append(f"logP 异常 ({logp:.2f})")
    if tpsa <= 90:
        score += 0.35
    else:
        reasons.append(f"TPSA 偏高 ({tpsa:.1f} > 90 Å²)")
    if hbd <= 3:
        score += 0.15
    else:
        reasons.append(f"HBD 偏多 ({hbd})")
    # 可旋转键略惩罚
    nrb = int(Descriptors.NumRotatableBonds(mol))
    if nrb <= 10:
        score += max(0.0, 0.1 - 0.01 * max(0, nrb - 7))
    else:
        reasons.append(f"可旋转键较多 ({nrb})")
    label = "可能利于 CNS 暴露（启发式）" if score >= 0.65 else "BBB 透过风险不确定或偏低（启发式）"
    return {
        "MW": mw,
        "MolLogP": logp,
        "TPSA": tpsa,
        "HBD": hbd,
        "NumRotatableBonds": nrb,
        "bbb_heuristic_score": round(min(1.0, score), 4),
        "summary_label": label,
        "notes": reasons or ["未触发明显红旗（启发式）"],
        "disclaimer": "此为基于理化描述符的筛选启发式，非实验或 PBPK 预测。",
    }


def _pains(mol: Chem.Mol) -> Dict[str, Any]:
    params = FilterCatalogParams()
    params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS_A)
    params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS_B)
    params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS_C)
    fc = FilterCatalog(params)
    hits = fc.GetMatches(mol)
    detail = [{"name": _filter_catalog_hit_description(h)} for h in hits]
    return {
        "has_pains_hit": len(hits) > 0,
        "hit_count": len(hits),
        "hits": detail[:50],
    }


def _brenk(mol: Chem.Mol) -> Dict[str, Any]:
    hits: List[Dict[str, str]] = []
    # 新版 RDKit 含 BRENK 目录时优先使用
    brenk_cat = getattr(FilterCatalogParams.FilterCatalogs, "BRENK", None)
    if brenk_cat is not None:
        params = FilterCatalogParams()
        params.AddCatalog(brenk_cat)
        fc = FilterCatalog(params)
        for h in fc.GetMatches(mol):
            hits.append({"description": _filter_catalog_hit_description(h)})
    else:
        # 极简 Brenk-like  SMARTS 子集（兜底）
        fallback = [
            ("迈克尔受体", "[CX3:1]=[CX3]-[CX4](-[OX2])-[OX2]"),
            ("氮芥类", "[NX3;H2,H1][CX4]"),
            ("硫芥类", "[SX2][CX4][Cl]"),
            ("醛基反应头", "[CX3H1](=O)"),
            ("叠氮化物", "[N;X2-]=[N+]=[N-]"),
            ("环氧化物", "[OX2r3]"),
        ]
        for label, sma in fallback:
            q = Chem.MolFromSmarts(sma)
            if q and mol.HasSubstructMatch(q):
                hits.append({"description": f"{label} ({sma})"})
    return {
        "has_brenk_alert": len(hits) > 0,
        "hit_count": len(hits),
        "hits": hits[:50],
        "catalog": "FilterCatalog.BRENK" if brenk_cat is not None else "fallback_smarts",
    }


def _morgan_fp(mol: Chem.Mol, radius: int = 2, n_bits: int = 2048) -> Dict[str, Any]:
    fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius, nBits=n_bits)
    fp_bin = fp.ToBinary()
    fp_digest = hashlib.sha256(fp_bin).hexdigest()[:48]
    return {
        "radius": radius,
        "n_bits": n_bits,
        "fp_sha256_preview": fp_digest,
        "density": float(fp.GetNumOnBits() / max(1, n_bits)),
        "num_on_bits": int(fp.GetNumOnBits()),
    }


def _tanimoto(smiles_a: str, smiles_b: str) -> Dict[str, Any]:
    m1 = Chem.MolFromSmiles(smiles_a)
    m2 = Chem.MolFromSmiles(smiles_b)
    if not m1 or not m2:
        raise ValueError("无法解析 SMILES 对")
    fp1 = AllChem.GetMorganFingerprintAsBitVect(m1, 2, nBits=2048)
    fp2 = AllChem.GetMorganFingerprintAsBitVect(m2, 2, nBits=2048)
    sim = DataStructs.TanimotoSimilarity(fp1, fp2)
    dist = 1.0 - float(sim)
    return {
        "smiles_a": Chem.MolToSmiles(m1),
        "smiles_b": Chem.MolToSmiles(m2),
        "tanimoto_similarity": round(float(sim), 6),
        "tanimoto_distance": round(float(dist), 6),
        "fingerprint": "Morgan r=2, 2048-bit",
    }


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--tool",
        required=True,
        choices=("bbb", "pains", "brenk", "morgan_fp", "mol_weight", "tanimoto"),
    )
    ap.add_argument("--smiles", default="")
    ap.add_argument("--smiles2", default="")
    ap.add_argument("--file", default="")
    ap.add_argument("--output-dir", required=True)
    args = ap.parse_args()
    out = Path(args.output_dir).resolve()
    out.mkdir(parents=True, exist_ok=True)

    try:
        smi = (args.smiles or "").strip()
        result: Dict[str, Any]

        if args.tool == "tanimoto":
            if (args.file or "").strip():
                fpth = Path(args.file.strip()).expanduser()
                if not fpth.is_file():
                    raise ValueError(f"文件不存在: {fpth}")
                sa, sb = _read_two_smiles(fpth)
            else:
                sa = (args.smiles or "").strip()
                sb = (args.smiles2 or "").strip()
                if "|" in sa and not sb:
                    parts = sa.split("|", 1)
                    sa, sb = parts[0].strip(), parts[1].strip()
            if not sa or not sb:
                raise ValueError(
                    "Tanimoto 需要两条 SMILES：--smiles 与 --smiles2，或 --file 两行，或 --smiles 内用 | 分隔"
                )
            result = _tanimoto(sa, sb)
            m1 = Chem.MolFromSmiles(result["smiles_a"])
            m2 = Chem.MolFromSmiles(result["smiles_b"])
            if m1:
                _mol_png(m1, out / "structure_a.png")
            if m2:
                _mol_png(m2, out / "structure_b.png")
            summary: Dict[str, Any] = {"tool": "tanimoto", "result": result}
        else:
            if (args.file or "").strip():
                fpth = Path(args.file.strip()).expanduser()
                if not fpth.is_file():
                    raise ValueError(f"文件不存在: {fpth}")
                smi = _read_smiles_from_file(fpth)
            if not smi:
                raise ValueError("请提供 --smiles 或 --file（含首行 SMILES）")
            mol = Chem.MolFromSmiles(smi)
            if mol is None:
                raise ValueError(f"无法解析 SMILES: {smi!r}")
            if args.tool == "bbb":
                result = _bbb_assessment(mol)
            elif args.tool == "pains":
                result = _pains(mol)
            elif args.tool == "brenk":
                result = _brenk(mol)
            elif args.tool == "morgan_fp":
                result = _morgan_fp(mol)
            elif args.tool == "mol_weight":
                result = {
                    "exact_molecular_weight": float(round(Descriptors.ExactMolWt(mol), 6)),
                    "average_molecular_weight": float(round(Descriptors.MolWt(mol), 6)),
                }
            else:
                raise ValueError(f"未知工具: {args.tool}")
            _mol_png(mol, out / "structure.png")
            summary = {
                "tool": args.tool,
                "canonical_smiles": Chem.MolToSmiles(mol),
                "result": result,
            }

        js = out / "chem_summary.json"
        js.write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")
        _stderr(f"OK wrote {js}")
        return 0
    except Exception:
        traceback.print_exc(file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
