#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
将「单细胞时空动力学」沉浸式 Demo 持久化资产导出为本地可按类型浏览的目录，并可选打 ZIP。

数据来源：
1. 扫描 demo_execution_stream.txt（及可选 demo_planning_stream.txt）中出现的
   /assets/demo/spatiotemporal_tutorial/<文件名> 与容器路径（/app/results/... 等，需 --results-root）。
2. 从 SSE 的 data: JSON 中提取最终分析报告（report_data.report）为 Markdown。

用法（**必须在仓库根目录执行**，结果只写入 ``test_data/``）::

  python3 scripts/export_spatiotemporal_demo_bundle.py
  python3 scripts/export_spatiotemporal_demo_bundle.py --zip

目录结构（均在 ``<仓库根>/test_data/`` 下）::

  test_data/files/images/   — .png .jpg .jpeg .gif .webp .svg
  test_data/files/tabular/  — .csv
  test_data/files/archives/ — .zip
  test_data/files/h5ad/     — .h5ad
  test_data/text/           — 提取的报告 .md、剧本 .txt、manifest.json

使用 ``--zip`` 时，压缩包路径为 ``test_data/spatiotemporal_demo_bundle.zip``（与上列目录同级，不另建其它输出目录）。

若部署在服务器上，也可在本机用 curl 拉取（需替换主机与端口）::

  curl -fLo umap_time.png https://<host>/assets/demo/spatiotemporal_tutorial/umap_time.png
"""
from __future__ import annotations

import argparse
import json
import re
import shutil
import sys
import zipfile
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Set, Tuple

# 与 scripts/build_demo_assets.py 对齐，并扩展常见后缀
_EXT = r"\.(?:png|jpg|jpeg|gif|webp|csv|json|h5ad|md|txt|zip|pdf|svg)"
DOCKER_PATH_RES = [
    re.compile(rf"/app/uploads/results/run_\d{{8}}_\d{{6}}/[\w./-]+{_EXT}", re.I),
    re.compile(rf"/app/results/run_\d{{8}}_\d{{6}}/[\w./-]+{_EXT}", re.I),
    re.compile(rf"/results/run_\d{{8}}_\d{{6}}/[\w./-]+{_EXT}", re.I),
    re.compile(rf"/app/test_data/[\w./-]+{_EXT}", re.I),
]
STATIC_ASSET_RE = re.compile(
    rf"/assets/demo/spatiotemporal_tutorial/([\w./-]+{_EXT})", re.I
)


def collect_paths_from_text(text: str) -> Set[str]:
    found: Set[str] = set()
    for rx in DOCKER_PATH_RES:
        found.update(rx.findall(text))
    for m in STATIC_ASSET_RE.finditer(text):
        found.add(m.group(0))
    return found


def resolve_host_path(
    docker_path: str,
    uploads_root: Optional[Path],
    results_root: Optional[Path],
    test_data_root: Optional[Path],
) -> Optional[Path]:
    if docker_path.startswith("/app/uploads/"):
        if not uploads_root:
            return None
        rel = docker_path[len("/app/uploads/") :].lstrip("/")
        return uploads_root / rel
    if docker_path.startswith("/app/results/"):
        if not results_root:
            return None
        rel = docker_path[len("/app/results/") :].lstrip("/")
        return results_root / rel
    if docker_path.startswith("/results/"):
        if not results_root:
            return None
        rel = docker_path[len("/results/") :].lstrip("/")
        return results_root / rel
    if docker_path.startswith("/app/test_data/"):
        if not test_data_root:
            return None
        rel = docker_path[len("/app/test_data/") :].lstrip("/")
        return test_data_root / rel
    return None


def static_url_to_repo_path(url_path: str, demo_dir: Path) -> Optional[Path]:
    prefix = "/assets/demo/spatiotemporal_tutorial/"
    if not url_path.startswith(prefix):
        return None
    name = url_path[len(prefix) :].lstrip("/")
    if ".." in name or name.startswith("/"):
        return None
    p = demo_dir / name
    return p if p.is_file() else None


def type_subdir(ext: str) -> str:
    e = ext.lower().lstrip(".")
    if e in ("png", "jpg", "jpeg", "gif", "webp", "svg"):
        return "images"
    if e == "csv":
        return "tabular"
    if e == "zip":
        return "archives"
    if e == "h5ad":
        return "h5ad"
    if e in ("md", "txt"):
        return "text"
    if e == "json":
        return "json"
    if e == "pdf":
        return "pdf"
    return "other"


def iter_sse_json_payloads(sse_text: str) -> Iterable[Dict[str, Any]]:
    for line in sse_text.splitlines():
        if not line.startswith("data:"):
            continue
        raw = line[5:].strip()
        if not raw or raw == "[DONE]":
            continue
        try:
            yield json.loads(raw)
        except json.JSONDecodeError:
            continue


def extract_final_report_md(sse_text: str) -> Optional[str]:
    """取最后一次出现的 report_data.report（通常为最长正文）。"""
    best: Optional[str] = None
    best_len = 0
    for payload in iter_sse_json_payloads(sse_text):
        rd = payload.get("report_data")
        if not isinstance(rd, dict):
            continue
        r = rd.get("report")
        if isinstance(r, str) and len(r) > best_len:
            best = r
            best_len = len(r)
    return best


def copy_unique(
    src: Path,
    dest_dir: Path,
    used_names: Set[str],
) -> Tuple[Path, str]:
    name = src.name
    if name not in used_names:
        used_names.add(name)
        dest_name = name
    else:
        stem, suf = src.stem, src.suffix
        i = 1
        while True:
            cand = f"{stem}__dup{i}{suf}"
            if cand not in used_names:
                used_names.add(cand)
                dest_name = cand
                break
            i += 1
    sub = type_subdir(src.suffix)
    out_sub = dest_dir / sub
    out_sub.mkdir(parents=True, exist_ok=True)
    dest = out_sub / dest_name
    shutil.copy2(src, dest)
    return dest, sub


def main() -> int:
    ap = argparse.ArgumentParser(description="Export spatiotemporal demo assets to typed folders + zip")
    ap.add_argument(
        "--demo-dir",
        type=Path,
        default=None,
        help="教程资源目录，默认 <repo>/services/nginx/html/assets/demo/spatiotemporal_tutorial",
    )
    ap.add_argument(
        "--streams",
        type=str,
        default="demo_execution_stream.txt,demo_planning_stream.txt",
        help="逗号分隔，相对于 demo-dir 的 SSE 剧本文件名",
    )
    ap.add_argument(
        "--out",
        type=Path,
        default=None,
        help=argparse.SUPPRESS,
    )
    ap.add_argument("--uploads-root", type=Path, default=None, help="宿主机路径，对应容器 /app/uploads")
    ap.add_argument(
        "--results-root",
        type=Path,
        default=None,
        help="宿主机路径，对应容器 /app/results（默认 <repo>/data/results）",
    )
    ap.add_argument("--test-data-root", type=Path, default=None, help="宿主机路径，对应 /app/test_data")
    ap.add_argument(
        "--zip",
        action="store_true",
        help="在输出目录内生成 spatiotemporal_demo_bundle.zip（与 files/、text/ 同属一个结果文件夹）",
    )
    ap.add_argument(
        "--no-report-md",
        action="store_true",
        help="不写入从 SSE 解析的 final_analysis_report.md",
    )
    args = ap.parse_args()
    repo = Path(__file__).resolve().parent.parent
    demo_dir = args.demo_dir or (repo / "services/nginx/html/assets/demo/spatiotemporal_tutorial")
    demo_dir = demo_dir.resolve()
    if not demo_dir.is_dir():
        print(f"ERROR: demo-dir 不存在: {demo_dir}", file=sys.stderr)
        return 1

    results_root = args.results_root or (repo / "data/results")
    out_root = (args.out or (repo / "test_data")).expanduser().resolve()
    files_root = out_root / "files"
    text_root = out_root / "text"
    files_root.mkdir(parents=True, exist_ok=True)
    text_root.mkdir(parents=True, exist_ok=True)

    stream_names = [s.strip() for s in args.streams.split(",") if s.strip()]
    combined = ""
    for sn in stream_names:
        p = demo_dir / sn
        if not p.is_file():
            print(f"WARN: 跳过不存在的剧本: {p}", file=sys.stderr)
            continue
        combined += p.read_text(encoding="utf-8", errors="replace") + "\n"

    if not combined.strip():
        print("ERROR: 未读取到任何 SSE 剧本内容", file=sys.stderr)
        return 1

    paths = collect_paths_from_text(combined)
    used: Set[str] = set()
    manifest: List[Dict[str, Any]] = []
    missing: List[str] = []

    for url_path in sorted(paths):
        src: Optional[Path] = None
        if url_path.startswith("/assets/demo/spatiotemporal_tutorial/"):
            src = static_url_to_repo_path(url_path, demo_dir)
        else:
            src = resolve_host_path(
                url_path,
                args.uploads_root,
                results_root,
                args.test_data_root,
            )
        if src is None or not src.is_file():
            missing.append(url_path)
            continue
        dest, sub = copy_unique(src, files_root, used)
        manifest.append(
            {
                "source_reference": url_path,
                "host_path": str(src),
                "exported_to": str(dest.relative_to(out_root)),
                "type_folder": sub,
            }
        )

    if not args.no_report_md:
        report = extract_final_report_md(combined)
        if report:
            rp = text_root / "final_analysis_report.md"
            rp.write_text(report, encoding="utf-8")
            manifest.append(
                {
                    "source_reference": "sse:report_data.report",
                    "exported_to": str(rp.relative_to(out_root)),
                    "type_folder": "text",
                }
            )

    for sn in stream_names:
        sp = demo_dir / sn
        if sp.is_file():
            shutil.copy2(sp, text_root / sn)
    blueprint = demo_dir / "demo_blueprint.md"
    if blueprint.is_file():
        shutil.copy2(blueprint, text_root / "demo_blueprint.md")

    man_path = text_root / "manifest.json"
    man_path.write_text(
        json.dumps(
            {
                "demo_dir": str(demo_dir),
                "exported_files": len([m for m in manifest if "host_path" in m]),
                "missing_references": missing,
                "entries": manifest,
            },
            ensure_ascii=False,
            indent=2,
        ),
        encoding="utf-8",
    )

    readme = text_root / "README_EXPORT.txt"
    readme.write_text(
        f"本目录由 scripts/export_spatiotemporal_demo_bundle.py 生成。\n"
        f"导出根目录（本说明所在目录的上一级）: {out_root}\n"
        "files/ 下按类型分子文件夹：images、tabular、archives、h5ad、text、json、pdf、other。\n"
        "text/final_analysis_report.md 为从 demo 执行流式剧本中解析的最终解读报告（若存在）。\n"
        "manifest.json 列出每个引用路径与导出相对路径；missing_references 为剧本中出现但未在本地解析到的路径。\n"
        "若使用 --zip，压缩包在同根目录下 spatiotemporal_demo_bundle.zip。\n",
        encoding="utf-8",
    )

    print(f"OK 导出目录: {out_root}")
    print(f"  已复制数据文件: {len([m for m in manifest if 'host_path' in m])} 个")
    if missing:
        print(f"  缺失引用（需调整 --results-root 或补齐 demo 目录）: {len(missing)} 个")
        for u in missing[:15]:
            print(f"    - {u}")
        if len(missing) > 15:
            print(f"    ... 另有 {len(missing) - 15} 条，见 manifest.json")

    if args.zip:
        zip_path = out_root / "spatiotemporal_demo_bundle.zip"
        if zip_path.is_file():
            zip_path.unlink()
        to_zip = [
            f
            for f in out_root.rglob("*")
            if f.is_file() and f.resolve() != zip_path.resolve()
        ]
        with zipfile.ZipFile(zip_path, "w", zipfile.ZIP_DEFLATED) as zf:
            for f in to_zip:
                zf.write(f, f.relative_to(out_root))
        print(f"OK ZIP（已置于输出目录内）: {zip_path}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
