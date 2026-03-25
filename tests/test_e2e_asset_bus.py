#!/usr/bin/env python3
"""
实弹 E2E：资产总线 → 体检路径 → WorkflowExecutor 参数注入
运行: python tests/test_e2e_asset_bus.py
"""
from __future__ import annotations

import logging
import os
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

# 仓库根目录
ROOT = Path(__file__).resolve().parent.parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

logging.basicConfig(
    level=logging.INFO,
    format="%(levelname)s %(name)s: %(message)s",
)
log = logging.getLogger("test_e2e_asset_bus")


def _find_csv_under(base: Path, limit: int = 50) -> Optional[Path]:
    if not base.is_dir():
        return None
    count = 0
    for p in base.rglob("*.csv"):
        if p.is_file() and "human_cachexia" in p.name.lower():
            return p
        count += 1
        if count > limit:
            break
    for p in base.rglob("*.csv"):
        if p.is_file():
            return p
    return None


def _find_nii_pair_under(base: Path) -> Optional[Tuple[Path, Path]]:
    if not base.is_dir():
        return None
    masks = [p for p in base.rglob("*.nii.gz") if p.is_file() and "mask" in p.name.lower()]
    for m in masks:
        parent = m.parent
        stem = m.name.lower().replace("_mask.nii.gz", "").replace("mask.nii.gz", "").strip("_")
        for p in parent.glob("*.nii.gz"):
            if p == m:
                continue
            if "mask" not in p.name.lower():
                return (p.resolve(), m.resolve())
    # 任意两个 nii.gz，一名含 mask
    all_nii = [p for p in base.rglob("*.nii.gz") if p.is_file()]
    masks = [p for p in all_nii if "mask" in p.name.lower() or "seg" in p.name.lower()]
    imgs = [p for p in all_nii if p not in masks]
    if masks and imgs:
        return (imgs[0].resolve(), masks[0].resolve())
    return None


def _find_10x_dir_under(base: Path) -> Optional[Path]:
    if not base.is_dir():
        return None
    for mtx in base.rglob("matrix.mtx"):
        if mtx.is_file():
            parent = mtx.parent
            names = {c.name.lower() for c in parent.iterdir() if c.is_file()}
            if any("barcodes" in n for n in names) and any(
                "features" in n or "genes" in n for n in names
            ):
                return parent.resolve()
    return None


def ensure_test_artifacts() -> Tuple[Path, Dict[str, Any]]:
    """
    寻找真实文件；找不到则在 /tmp/mock_omics 生成可解析的最小结构。
    返回 (upload_root, meta)，meta 含绝对路径键。
    """
    candidates = [
        Path("/app/uploads"),
        ROOT / "tests" / "data",
        ROOT / "uploads",
    ]
    csv_p: Optional[Path] = None
    pair: Optional[Tuple[Path, Path]] = None
    tenx: Optional[Path] = None

    for c in candidates:
        if csv_p is None:
            csv_p = _find_csv_under(c)
        if pair is None:
            pair = _find_nii_pair_under(c)
        if tenx is None:
            tenx = _find_10x_dir_under(c)

    mock_root = Path(os.environ.get("MOCK_OMICS_DIR", "/tmp/mock_omics"))
    need_mock = csv_p is None or pair is None or tenx is None

    if need_mock:
        mock_root.mkdir(parents=True, exist_ok=True)
        log.info("部分样本缺失，写入 Mock 目录: %s", mock_root)

    if csv_p is None:
        csv_p = mock_root / "synthetic_metabo.csv"
        csv_p.write_text("sample_id,group,metab_1,metab_2\nS1,Control,1.0,2.0\nS2,Treat,1.1,2.1\n", encoding="utf-8")
        log.info("已生成代谢组 CSV Mock: %s", csv_p)

    if pair is None:
        rad = mock_root / "radiomics_case"
        rad.mkdir(parents=True, exist_ok=True)
        img = rad / "case_t1w.nii.gz"
        msk = rad / "case_t1w_mask.nii.gz"
        img.write_bytes(b"")
        msk.write_bytes(b"")
        pair = (img.resolve(), msk.resolve())
        log.info("已生成影像 Mock: %s + %s", pair[0], pair[1])

    if tenx is None:
        td = mock_root / "tenx_sample" / "filtered_gene_bc_matrices" / "mm10"
        td.mkdir(parents=True, exist_ok=True)
        (td / "matrix.mtx").write_text("", encoding="utf-8")
        (td / "barcodes.tsv").write_text("AAAC\n", encoding="utf-8")
        (td / "features.tsv").write_text("gene\tname\nG1\tg1\n", encoding="utf-8")
        tenx = td.resolve()
        log.info("已生成 10x Mock 目录: %s", tenx)

    if need_mock:
        upload_root = mock_root.resolve()
    else:
        paths = [csv_p.resolve(), pair[0], pair[1], tenx]
        try:
            upload_root = Path(os.path.commonpath([str(p) for p in paths]))
        except ValueError:
            upload_root = csv_p.parent.resolve()
        if not upload_root.is_dir():
            upload_root = csv_p.parent.resolve()

    return upload_root, {
        "csv": csv_p.resolve(),
        "image_nii": pair[0],
        "mask_nii": pair[1],
        "tenx_dir": tenx,
    }


def orchestrator_style_path_for_inspect(
    upload_dir: Path,
    files: List[Dict[str, str]],
) -> str:
    """与 orchestrator 多文件分支一致的资产优先级（简化版，仅本测试使用）。"""
    from gibh_agent.core.file_inspector import OmicsAssetManager

    resolved_paths: List[Path] = []
    for f in files:
        p = f.get("path") or f.get("file_path") or ""
        if not p:
            continue
        po = Path(p)
        if not po.is_absolute():
            po = (upload_dir / po).resolve()
        else:
            po = po.resolve()
        if po.exists():
            resolved_paths.append(po)

    if len(resolved_paths) <= 1:
        p0 = resolved_paths[0] if resolved_paths else Path(files[0]["path"])
        if not p0.is_absolute():
            p0 = (upload_dir / p0).resolve()
        return str(p0.resolve())

    unique: List[str] = []
    seen = set()
    for rp in resolved_paths:
        s = str(rp)
        if s not in seen:
            seen.add(s)
            unique.append(s)

    assets = OmicsAssetManager().classify_assets(unique)
    chosen: Optional[str] = None
    for a in assets:
        if a.asset_type == "10x_bundle":
            chosen = a.primary_path
            break
    if not chosen:
        for a in assets:
            if a.asset_type == "radiomics_pair":
                chosen = a.primary_path
                break
    if not chosen:
        for a in assets:
            if a.asset_type in ("metabolomics_csv", "h5ad_single"):
                chosen = a.primary_path
                break
    if not chosen:
        for a in assets:
            if a.asset_type == "radiomics_image":
                chosen = a.primary_path
                break
    if not chosen and unique:
        chosen = unique[0]

    p_ins = Path(chosen or unique[0])
    if not p_ins.is_absolute():
        p_ins = (upload_dir / p_ins).resolve()
    else:
        p_ins = p_ins.resolve()
    return str(p_ins)


def rel_under(root: Path, path: Path) -> str:
    try:
        return str(path.resolve().relative_to(root.resolve()))
    except ValueError:
        return str(path.resolve())


def main() -> None:
    # 先刷 stdout，避免 logging 走 stderr 时视觉上跑到横幅前面
    print("=" * 72, flush=True)
    print(" E2E 资产总线实弹演习", flush=True)
    print("=" * 72, flush=True)

    upload_root, meta = ensure_test_artifacts()
    print("\n[弹药库] upload_root =", upload_root, flush=True)
    print("  CSV        :", meta["csv"], flush=True)
    print("  影像       :", meta["image_nii"], flush=True)
    print("  Mask       :", meta["mask_nii"], flush=True)
    print("  10x 目录   :", meta["tenx_dir"], flush=True)

    from gibh_agent.core.asset_manager import OmicsAssetManager
    from gibh_agent.core.file_inspector import FileInspector
    from gibh_agent.core.executor import WorkflowExecutor

    # ------------------------------------------------------------------ 场景 A
    print("\n" + "=" * 72, flush=True)
    print(" 场景 A：单文件代谢组", flush=True)
    print("=" * 72, flush=True)
    csv_rel = rel_under(upload_root, meta["csv"])
    files_a = [{"path": csv_rel}]
    abs_csv = str(meta["csv"])
    assets_a = OmicsAssetManager().classify_assets([abs_csv])
    print("\n[步骤1] classify_assets([csv]) ->", [a.asset_type for a in assets_a])
    assert any(a.asset_type == "metabolomics_csv" for a in assets_a), "应产生 metabolomics_csv 资产"
    metab_asset = next(a for a in assets_a if a.asset_type == "metabolomics_csv")

    path_insp_a = orchestrator_style_path_for_inspect(upload_root, files_a)
    print("[步骤2] path_for_inspect =", path_insp_a)
    fi = FileInspector(str(upload_root))
    meta_ins_a = fi.inspect_file(path_insp_a)
    print("[步骤2] FileInspector.inspect_file -> status =", meta_ins_a.get("status"))
    print("        keys:", list(meta_ins_a.keys())[:12], "...")

    ex = WorkflowExecutor(
        output_dir=str(upload_root / "e2e_results"),
        upload_dir=str(upload_root),
    )
    params_a: Dict[str, Any] = {}
    ok = ex._inject_paths_from_asset(
        params_a, "Metabolomics", metab_asset, "file_path"
    )
    print("[步骤3] _inject_paths_from_asset(Metabolomics) done =", ok)
    print("        params =", params_a)
    assert "file_path" in params_a and params_a["file_path"], "应注入 file_path"
    print("✅ 场景 A 通过")

    # ------------------------------------------------------------------ 场景 B
    print("\n" + "=" * 72, flush=True)
    print(" 场景 B：多文件影像组（图 + mask）", flush=True)
    print("=" * 72, flush=True)
    img_rel = rel_under(upload_root, meta["image_nii"])
    msk_rel = rel_under(upload_root, meta["mask_nii"])
    files_b = [{"path": img_rel}, {"path": msk_rel}]
    abs_pair = [str(meta["image_nii"]), str(meta["mask_nii"])]
    assets_b = OmicsAssetManager().classify_assets(abs_pair)
    print("\n[步骤1] classify_assets([img, mask]) ->", [(a.asset_type, a.primary_path) for a in assets_b])
    assert any(a.asset_type == "radiomics_pair" for a in assets_b), "应产生 radiomics_pair 资产"
    rad_asset = next(a for a in assets_b if a.asset_type == "radiomics_pair")

    path_insp_b = orchestrator_style_path_for_inspect(upload_root, files_b)
    print("[步骤2] path_for_inspect =", path_insp_b)
    assert Path(path_insp_b).resolve() == meta["image_nii"].resolve(), "体检路径应为主图"
    meta_ins_b = fi.inspect_file(path_insp_b)
    print("[步骤2] FileInspector.inspect_file -> status =", meta_ins_b.get("status"))
    print("        file_type:", meta_ins_b.get("file_type"), "preview snippet:", str(meta_ins_b.get("preview"))[:120])

    params_b: Dict[str, Any] = {}
    ok_b = ex._inject_paths_from_asset(
        params_b, "Radiomics", rad_asset, "image_path"
    )
    print("[步骤3] _inject_paths_from_asset(Radiomics, 对应工具 extract_radiomics_features) done =", ok_b)
    print("        params =", params_b)
    assert params_b.get("image_path"), "应注入 image_path"
    assert params_b.get("mask_path"), "应注入 mask_path"
    print("✅ 场景 B 通过")

    # ------------------------------------------------------------------ 场景 C（可选）
    print("\n" + "=" * 72, flush=True)
    print(" 场景 C：10x 多文件 → 目录级资产", flush=True)
    print("=" * 72, flush=True)
    tenx = meta["tenx_dir"]
    mtx = tenx / "matrix.mtx"
    bc = tenx / "barcodes.tsv"
    ft = tenx / "features.tsv"
    rels = [rel_under(upload_root, p) for p in (mtx, bc, ft)]
    files_c = [{"path": r} for r in rels]
    abs_three = [str(mtx.resolve()), str(bc.resolve()), str(ft.resolve())]
    assets_c = OmicsAssetManager().classify_assets(abs_three)
    print("classify_assets(10x x3) ->", [(a.asset_type, a.primary_path) for a in assets_c])
    assert any(a.asset_type == "10x_bundle" for a in assets_c), "应产生 10x_bundle"
    path_insp_c = orchestrator_style_path_for_inspect(upload_root, files_c)
    print("path_for_inspect =", path_insp_c)
    assert Path(path_insp_c).is_dir(), "10x 体检应为目录"
    meta_ins_c = fi.inspect_file(path_insp_c)
    print("FileInspector.inspect_file(10x dir) -> status =", meta_ins_c.get("status"))
    tx_asset = next(a for a in assets_c if a.asset_type == "10x_bundle")
    params_c: Dict[str, Any] = {}
    ex._inject_paths_from_asset(params_c, "scRNA-seq", tx_asset, "adata_path")
    print("_inject_paths_from_asset(scRNA-seq) params =", params_c)
    assert params_c.get("adata_path") and Path(params_c["adata_path"]).is_dir()
    print("✅ 场景 C 通过")

    print("\n" + "=" * 72, flush=True)
    print(" 战报：全链路畅通，实弹演习结束", flush=True)
    print("=" * 72, flush=True)


if __name__ == "__main__":
    main()
