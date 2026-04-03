#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
从真实 SSE 捕获文本生成沉浸式 Demo 静态资产：
- 扫描 raw_exec.txt 中的结果路径（支持 /app/results/、/results/、/app/uploads/results/、/app/test_data/）
- 复制文件到前端目录，并将原文路径替换为 /assets/demo/spatiotemporal_tutorial/<文件名>
- 输出 demo_execution_stream.txt / demo_planning_stream.txt

用法（仓库根目录）：
  python3 scripts/build_demo_assets.py \\
    --raw-exec services/nginx/html/assets/demo/spatiotemporal_tutorial/raw_exec.txt \\
    --raw-plan services/nginx/html/assets/demo/spatiotemporal_tutorial/raw_plan.txt \\
    --results-root data/results \\
    --test-data-root test_data

与 docker-compose 一致：容器内 /app/results → 宿主机 ./data/results；/app/test_data → ./test_data
"""
from __future__ import annotations

import argparse
import re
import shutil
import sys
from pathlib import Path

STATIC_PREFIX = "/assets/demo/spatiotemporal_tutorial"

# 严格匹配带后缀的真实文件路径，避免 JSON/Markdown 粘连（如 .h5ad。中文）导致误截断
_EXT = r"\.(?:png|jpg|jpeg|gif|webp|csv|json|h5ad|md|txt|zip)"
PATH_RES = [
    re.compile(rf"/app/uploads/results/run_\d{{8}}_\d{{6}}/[\w./-]+{_EXT}", re.I),
    re.compile(rf"/app/results/run_\d{{8}}_\d{{6}}/[\w./-]+{_EXT}", re.I),
    re.compile(rf"/results/run_\d{{8}}_\d{{6}}/[\w./-]+{_EXT}", re.I),
    re.compile(rf"/app/test_data/[\w./-]+{_EXT}", re.I),
]


def _unique_dest_name(src: Path, used: set[str]) -> str:
    name = src.name
    if name not in used:
        used.add(name)
        return name
    stem = src.stem
    suf = src.suffix
    i = 1
    while True:
        cand = f"{stem}__{i}{suf}"
        if cand not in used:
            used.add(cand)
            return cand
        i += 1


def collect_paths(text: str) -> list[str]:
    found: list[str] = []
    seen: set[str] = set()
    for rx in PATH_RES:
        for m in rx.findall(text):
            if m not in seen:
                seen.add(m)
                found.append(m)
    return found


def resolve_host_path(docker_path: str, uploads_root: Path | None, results_root: Path | None, test_data_root: Path | None) -> Path | None:
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


def main() -> int:
    ap = argparse.ArgumentParser(description="Build immersive demo assets under services/nginx/html/assets/demo/")
    ap.add_argument("--raw-exec", type=Path, required=True, help="原始执行阶段 SSE 文本")
    ap.add_argument("--raw-plan", type=Path, required=True, help="原始规划阶段 SSE 文本")
    ap.add_argument(
        "--uploads-root",
        type=Path,
        default=None,
        help="宿主机目录：对应容器内 /app/uploads",
    )
    ap.add_argument(
        "--results-root",
        type=Path,
        default=None,
        help="宿主机目录：对应容器内 /app/results（默认 <repo>/data/results）",
    )
    ap.add_argument(
        "--test-data-root",
        type=Path,
        default=None,
        help="宿主机目录：对应容器内 /app/test_data（默认 <repo>/test_data）",
    )
    ap.add_argument(
        "--out-dir",
        type=Path,
        default=None,
        help="输出目录，默认 services/nginx/html/assets/demo/spatiotemporal_tutorial",
    )
    args = ap.parse_args()
    repo = Path(__file__).resolve().parent.parent
    out_dir = args.out_dir or (repo / "services/nginx/html/assets/demo/spatiotemporal_tutorial")
    out_dir.mkdir(parents=True, exist_ok=True)

    results_root = args.results_root or (repo / "data/results")
    test_data_root = args.test_data_root or (repo / "test_data")
    uploads_root = args.uploads_root

    if not args.raw_exec.is_file():
        print(f"ERROR: --raw-exec 不是文件: {args.raw_exec}", file=sys.stderr)
        return 1
    if not args.raw_plan.is_file():
        print(f"ERROR: --raw-plan 不是文件: {args.raw_plan}", file=sys.stderr)
        return 1

    raw_exec = args.raw_exec.read_text(encoding="utf-8", errors="replace")
    raw_plan = args.raw_plan.read_text(encoding="utf-8", errors="replace")

    paths = collect_paths(raw_exec)
    mapping: dict[str, str] = {}
    used_names: set[str] = set()

    for docker_path in paths:
        src = resolve_host_path(docker_path, uploads_root, results_root, test_data_root)
        if src is None:
            print(f"WARN: 无法解析宿主机路径（缺少对应 --*-root）: {docker_path}", file=sys.stderr)
            continue
        if not src.is_file():
            print(f"WARN: 源文件不存在，跳过: {src} (来自 {docker_path})", file=sys.stderr)
            continue
        dest_name = _unique_dest_name(src, used_names)
        dest = out_dir / dest_name
        shutil.copy2(src, dest)
        mapping[docker_path] = f"{STATIC_PREFIX}/{dest_name}"
        print(f"COPY {src} -> {dest}")

    def rewrite(text: str) -> str:
        out = text
        for old, new in sorted(mapping.items(), key=lambda x: -len(x[0])):
            out = out.replace(old, new)
        return out

    plan_out = rewrite(raw_plan)
    exec_out = rewrite(raw_exec)

    (out_dir / "demo_planning_stream.txt").write_text(plan_out, encoding="utf-8")
    (out_dir / "demo_execution_stream.txt").write_text(exec_out, encoding="utf-8")
    print(f"WROTE {out_dir / 'demo_planning_stream.txt'}")
    print(f"WROTE {out_dir / 'demo_execution_stream.txt'}")
    print(f"映射条目数: {len(mapping)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
