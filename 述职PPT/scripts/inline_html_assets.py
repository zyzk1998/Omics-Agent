#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
将 index.html 中所有「本地」<img src="..."> 转为 data URI，便于单文件分发。
不处理 <script src="https:...">。

默认读取 ../index.html，写入 ../index.self-contained.html

用法:
  cd 述职PPT && python3 scripts/inline_html_assets.py
  cd 述职PPT && python3 scripts/inline_html_assets.py --inplace   # 覆盖 index.html（慎用）
"""
from __future__ import annotations

import argparse
import base64
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
REPO_ROOT = ROOT.parent

IMG_SRC = re.compile(r'(<img\b[^>]*\bsrc=")([^"]+)(")', re.IGNORECASE | re.DOTALL)


def resolve_src(raw: str) -> Path | None:
    s = raw.strip()
    if s.startswith(("http://", "https://", "data:")):
        return None
    if s.startswith("../"):
        p = (REPO_ROOT / s[3:]).resolve()
    else:
        p = (ROOT / s).resolve()
    try:
        p.relative_to(REPO_ROOT)
    except ValueError:
        print("skip escape:", p, file=sys.stderr)
        return None
    return p if p.is_file() else None


def mime_for(path: Path) -> str:
    ext = path.suffix.lower()
    if ext == ".svg":
        return "image/svg+xml"
    if ext == ".png":
        return "image/png"
    if ext in (".jpg", ".jpeg"):
        return "image/jpeg"
    if ext == ".webp":
        return "image/webp"
    return "application/octet-stream"


def inline_html(html: str) -> str:
    def repl(m: re.Match[str]) -> str:
        pre, src, post = m.group(1), m.group(2), m.group(3)
        path = resolve_src(src)
        if path is None:
            return m.group(0)
        b64 = base64.b64encode(path.read_bytes()).decode("ascii")
        mt = mime_for(path)
        return f'{pre}data:{mt};base64,{b64}{post}'

    return IMG_SRC.sub(repl, html)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--inplace", action="store_true", help="覆盖 index.html")
    args = ap.parse_args()
    src_path = ROOT / "index.html"
    html = src_path.read_text(encoding="utf-8")
    out = inline_html(html)
    if args.inplace:
        out_path = src_path
    else:
        out_path = ROOT / "index.self-contained.html"
    out_path.write_text(out, encoding="utf-8")
    print("written:", out_path, "size_kb=", round(out_path.stat().st_size / 1024, 1))


if __name__ == "__main__":
    main()
