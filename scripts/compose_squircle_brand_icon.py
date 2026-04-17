#!/usr/bin/env python3
"""
从「浅灰方底 + 无穷符号」类 PNG 生成品牌 app-icon（1024²）：
  - 边缘 flood-fill 识别旧底 → 透明；
  - 与边缘不连通的同色孔洞 → 透明（∞ 形内旧灰底）；
  - 透明画布上绘制纯白圆角矩形「卡片」；
  - 将抠出的主体居中叠放。

与早期「仅无穷符号 + 透明底」一致之处：画布四角保持透明，不把整张图压成不透明白方块。

用法:
  python3 scripts/compose_squircle_brand_icon.py /path/to/flat-source.png
可选:
  --out PATH   默认写入 gibh-desktop-app/app-icon.png
  --radius R   圆角占边长比例，默认 0.21
"""
from __future__ import annotations

import argparse
import os
import sys
from collections import deque

import numpy as np
from PIL import Image, ImageDraw

_REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
_DEFAULT_OUT = os.path.join(_REPO_ROOT, "gibh-desktop-app", "app-icon.png")


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("source", help="带浅灰方底的源 PNG（建议正方形）")
    p.add_argument(
        "--out",
        default=_DEFAULT_OUT,
        help="输出路径（默认仓库内 gibh-desktop-app/app-icon.png）",
    )
    p.add_argument("--radius", type=float, default=0.21, help="圆角半径 / 边长")
    p.add_argument("--size", type=int, default=1024, help="输出边长")
    return p.parse_args()


def is_bg(rgb, ref, tol_l1: int = 45) -> bool:
    return int(abs(rgb[0] - ref[0]) + abs(rgb[1] - ref[1]) + abs(rgb[2] - ref[2])) <= tol_l1


def flood_edge_mask(rgba, ref):
    h, w = rgba.shape[0], rgba.shape[1]
    rgb = rgba[:, :, :3].astype(int)
    visited = np.zeros((h, w), dtype=bool)
    q: deque[tuple[int, int]] = deque()
    for x in range(w):
        for y in (0, h - 1):
            if is_bg(tuple(rgb[y, x]), ref):
                visited[y, x] = True
                q.append((x, y))
    for y in range(h):
        for x in (0, w - 1):
            if not visited[y, x] and is_bg(tuple(rgb[y, x]), ref):
                visited[y, x] = True
                q.append((x, y))
    while q:
        x, y = q.popleft()
        for dx, dy in ((1, 0), (-1, 0), (0, 1), (0, -1)):
            nx, ny = x + dx, y + dy
            if 0 <= nx < w and 0 <= ny < h and not visited[ny, nx] and is_bg(tuple(rgb[ny, nx]), ref):
                visited[ny, nx] = True
                q.append((nx, ny))
    return visited


def main() -> int:
    args = parse_args()
    try:
        src_im = Image.open(args.source).convert("RGBA")
    except OSError as e:
        print("无法读取源图:", e, file=sys.stderr)
        return 1

    src = np.array(src_im)
    h0, w0 = src.shape[:2]
    if w0 != h0:
        print("警告: 源图非正方形，将先居中裁成正方形", file=sys.stderr)
        s = min(w0, h0)
        x0, y0 = (w0 - s) // 2, (h0 - s) // 2
        src = src[y0 : y0 + s, x0 : x0 + s]
        h0 = w0 = s

    ref = np.mean(
        [src[0, 0, :3], src[w0 - 1, 0, :3], src[0, h0 - 1, :3], src[w0 - 1, h0 - 1, :3]],
        axis=0,
    )
    visited = flood_edge_mask(src, ref.astype(int))
    rgb = src[:, :, :3].astype(int)
    inner_gray = (~visited) & (
        np.abs(rgb - ref.astype(int)).sum(axis=2) <= 45
    )

    layer = src.copy()
    layer[:, :, 3] = 255
    layer[visited, 3] = 0
    layer[inner_gray, 3] = 0

    alpha = layer[:, :, 3]
    ys, xs = np.where(alpha > 16)
    if ys.size == 0:
        print("未识别到前景主体，请检查源图背景色是否与四角一致。", file=sys.stderr)
        return 1
    y0, y1 = ys.min(), ys.max()
    x0, x1 = xs.min(), xs.max()
    fg = Image.fromarray(layer[y0 : y1 + 1, x0 : x1 + 1], mode="RGBA")

    W = H = args.size
    corner_r = int(round(min(W, H) * args.radius))
    plate = Image.new("RGBA", (W, H), (0, 0, 0, 0))
    dr = ImageDraw.Draw(plate, "RGBA")
    dr.rounded_rectangle([0, 0, W - 1, H - 1], radius=corner_r, fill=(255, 255, 255, 255))

    pad = 0.10
    inner = int(min(W, H) * (1 - 2 * pad))
    fw, fh = fg.size
    scale = inner / max(fw, fh)
    nw, nh = int(round(fw * scale)), int(round(fh * scale))
    fg_r = fg.resize((nw, nh), Image.Resampling.LANCZOS)

    canvas = Image.new("RGBA", (W, H), (0, 0, 0, 0))
    canvas.alpha_composite(plate)
    ox = (W - nw) // 2
    oy = (H - nh) // 2
    canvas.alpha_composite(fg_r, (ox, oy))

    out_path = args.out
    canvas.save(out_path, format="PNG", optimize=True)
    print("已写入", out_path)
    return 0


if __name__ == "__main__":
    sys.exit(main())
