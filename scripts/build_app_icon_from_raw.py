#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""从带灰底的客户端原图生成透明底、白圆角底板、主体约 80% 占位的 app 图标。"""
from __future__ import annotations

import argparse
import math
import os
import sys
from collections import deque

import numpy as np
from PIL import Image, ImageDraw

def _edge_flood_mask(rgb: np.ndarray, ref: np.ndarray, tol_l1: int = 45) -> np.ndarray:
    h, w = rgb.shape[0], rgb.shape[1]

    def is_bg(rgbv: np.ndarray) -> bool:
        return int(np.abs(rgbv - ref).sum()) <= tol_l1

    visited = np.zeros((h, w), dtype=bool)
    q: deque[tuple[int, int]] = deque()
    ref = ref.astype(np.int64)
    for x in range(w):
        for y in (0, h - 1):
            if is_bg(rgb[y, x]):
                visited[y, x] = True
                q.append((x, y))
    for y in range(h):
        for x in (0, w - 1):
            if not visited[y, x] and is_bg(rgb[y, x]):
                visited[y, x] = True
                q.append((x, y))
    while q:
        x, y = q.popleft()
        for dx, dy in ((1, 0), (-1, 0), (0, 1), (0, -1)):
            nx, ny = x + dx, y + dy
            if 0 <= nx < w and 0 <= ny < h and not visited[ny, nx] and is_bg(rgb[ny, nx]):
                visited[ny, nx] = True
                q.append((nx, ny))
    return visited


def _inner_gray_holes_mask(
    rgb: np.ndarray, ref: np.ndarray, visited_edge: np.ndarray, tol_l1: int = 45
) -> np.ndarray:
    """与边缘不连通的、与参考色相近的灰底孔洞（∞ 内部等）→ 透明。"""
    refi = ref.astype(np.int64)
    d = np.abs(rgb.astype(np.int64) - refi).sum(axis=2)
    inner = (~visited_edge) & (d <= tol_l1)
    return inner


def _crop_content_rgba(im: Image.Image) -> tuple[Image.Image, tuple[int, int, int, int]]:
    """去边连灰底与内孔灰底，裁剪到白卡+符号外接。"""
    rgba = np.array(im.convert("RGBA"))
    h, w = rgba.shape[:2]
    rgb = rgba[:, :, :3]
    # 角点：array 为 (行 y, 列 x) = (高, 宽)
    corners = np.array(
        [rgb[0, 0], rgb[0, w - 1], rgb[h - 1, 0], rgb[h - 1, w - 1]], dtype=np.float64
    )
    ref = corners.mean(axis=0)
    visited = _edge_flood_mask(rgb, ref)
    inner_gray = _inner_gray_holes_mask(rgb, ref, visited)
    keep = (~visited) & (~inner_gray)
    ys, xs = np.where(keep)
    if ys.size == 0:
        raise RuntimeError("未检测到主体：请检查原图是否与四角背景色一致。")
    y0, y1 = int(ys.min()), int(ys.max())
    x0, x1 = int(xs.min()), int(xs.max())
    crop = rgba[y0 : y1 + 1, x0 : x1 + 1].copy()
    crop[visited[y0 : y1 + 1, x0 : x1 + 1], 3] = 0
    crop[inner_gray[y0 : y1 + 1, x0 : x1 + 1], 3] = 0
    patch = Image.fromarray(crop, mode="RGBA")
    return patch, (x0, y0, x1, y1)


def _extract_symbol_rgba(cropped: Image.Image, white_thr: int = 248) -> Image.Image:
    """从白底卡片区分离莫比乌斯彩色主体为独立 RGBA。"""
    a = np.array(cropped.convert("RGBA"))
    rgb = a[:, :, :3].astype(np.int64)
    # 近白纸张 → 透明
    paper = (rgb[:, :, 0] >= white_thr) & (rgb[:, :, 1] >= white_thr) & (rgb[:, :, 2] >= white_thr)
    sym = a.copy()
    sym[paper, 3] = 0
    return Image.fromarray(sym, mode="RGBA")


def _rounded_rectangle_mask(size: int, radius: float) -> Image.Image:
    """抗锯齿圆角矩形 Alpha 蒙版。"""
    r = Image.new("L", (size, size), 0)
    dr = ImageDraw.Draw(r)
    rr = int(round(radius))
    dr.rounded_rectangle([0, 0, size - 1, size - 1], radius=rr, fill=255)
    return r


def _compose_final(
    symbol_rgba: Image.Image,
    canvas_size: int,
    radius_ratio: float,
    symbol_area_ratio: float,
) -> Image.Image:
    """透明底 + 白圆角底板 + 居中缩放符号（面积约 symbol_area_ratio）。"""
    sym = symbol_rgba.copy()
    bbox = sym.getbbox()
    if bbox is None:
        raise RuntimeError("符号层全透明。")
    sym = sym.crop(bbox)
    sw, sh = sym.size
    sym_area = float(sw * sh)
    plate_area = float(canvas_size * canvas_size)
    # 目标：符号像素包围盒面积 ≈ symbol_area_ratio * 白板面积（白板近似满幅圆角方）
    k_area = math.sqrt(symbol_area_ratio * plate_area / max(sym_area, 1.0))
    # 保险：不超过画布可放下的最大均匀缩放
    k_fit = min(canvas_size / sw, canvas_size / sh)
    scale = min(k_area, k_fit * 0.98)
    nw = max(1, int(round(sw * scale)))
    nh = max(1, int(round(sh * scale)))
    sym_r = sym.resize((nw, nh), Image.Resampling.LANCZOS)

    canvas = Image.new("RGBA", (canvas_size, canvas_size), (0, 0, 0, 0))
    radius = canvas_size * radius_ratio
    mask = _rounded_rectangle_mask(canvas_size, radius)
    plate = Image.new("RGBA", (canvas_size, canvas_size), (0, 0, 0, 0))
    white = Image.new("RGBA", (canvas_size, canvas_size), (255, 255, 255, 255))
    plate.paste(white, mask=mask)

    ox = (canvas_size - nw) // 2
    oy = (canvas_size - nh) // 2
    plate_rgba = Image.new("RGBA", (canvas_size, canvas_size), (0, 0, 0, 0))
    plate_rgba.paste(plate, (0, 0))
    plate_rgba.paste(sym_r, (ox, oy), sym_r)
    return plate_rgba


def _save_ico(png1024: Image.Image, ico_path: str, sizes: tuple[int, ...]) -> None:
    rgba = png1024.convert("RGBA")
    imgs = [rgba.resize((s, s), Image.Resampling.LANCZOS) for s in sizes]
    imgs[0].save(
        ico_path,
        format="ICO",
        append_images=imgs[1:] if len(imgs) > 1 else [],
    )


def main() -> int:
    ap = argparse.ArgumentParser(description="从原图生成 app_icon.png / app_icon.ico")
    ap.add_argument("source", help="原图 PNG 路径（带灰底+白区+图案）")
    ap.add_argument(
        "--out-dir",
        default="",
        help="输出目录（默认：仓库内 gibh-desktop-app/，即本脚本上级的该子目录）",
    )
    ap.add_argument("--size", type=int, default=1024, help="输出正方形边长，默认 1024")
    ap.add_argument("--radius-ratio", type=float, default=0.225, help="圆角半径 / 边长")
    ap.add_argument(
        "--symbol-area-ratio",
        type=float,
        default=0.80,
        help="符号包围盒相对整板目标面积占比（约 0.80）",
    )
    ap.add_argument(
        "--also-app-icon",
        action="store_true",
        help="同时写入 app-icon.png / build/icon.ico（electron-builder 用名）",
    )
    args = ap.parse_args()

    src = os.path.abspath(args.source)
    if not os.path.isfile(src):
        print("找不到原图:", src, file=sys.stderr)
        return 1

    root = os.path.dirname(os.path.abspath(__file__))
    repo = os.path.dirname(root)
    out_dir = args.out_dir or os.path.join(repo, "gibh-desktop-app")
    os.makedirs(out_dir, exist_ok=True)
    build_dir = os.path.join(out_dir, "build")
    os.makedirs(build_dir, exist_ok=True)

    im0 = Image.open(src).convert("RGBA")
    cropped, _bb = _crop_content_rgba(im0)
    symbol = _extract_symbol_rgba(cropped)
    final = _compose_final(symbol, args.size, args.radius_ratio, args.symbol_area_ratio)

    png_out = os.path.join(out_dir, "app_icon.png")
    final.save(png_out, format="PNG", optimize=True)

    ico_sizes = (256, 48, 32, 16)
    ico_out = os.path.join(out_dir, "app_icon.ico")
    _save_ico(final, ico_out, ico_sizes)

    print("已写入:", png_out)
    print("已写入:", ico_out)

    if args.also_app_icon:
        app_png = os.path.join(out_dir, "app-icon.png")
        final.save(app_png, format="PNG", optimize=True)
        _save_ico(final, os.path.join(build_dir, "icon.ico"), ico_sizes)
        print("已写入:", app_png)
        print("已写入:", os.path.join(build_dir, "icon.ico"))

    return 0


if __name__ == "__main__":
    sys.exit(main())
