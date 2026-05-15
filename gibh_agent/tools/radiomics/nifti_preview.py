"""
NIfTI / 医学影像：中心切片 PNG → Base64 data URI，供步骤 summary/markdown 内嵌预览。

优先 nibabel 读入体数据；失败时回退 SimpleITK（与 radiomics/io 一致）。
主进程内仅生成小图（默认单张轴位），符合 TaaS 外轻量可视化约定。
"""
from __future__ import annotations

import base64
import io
import logging
from pathlib import Path
from typing import Optional

import numpy as np

logger = logging.getLogger(__name__)


def _volume_array_from_nifti(path: str) -> Optional[np.ndarray]:
    """返回 3D float32 数组 (z, y, x) 或 (x, y, z) 统一为可切片形式。"""
    p = Path(path)
    if not p.is_file():
        return None
    try:
        import nibabel as nib  # type: ignore

        img = nib.load(str(p))
        data = np.asanyarray(img.dataobj, dtype=np.float32)
        if data.ndim > 3:
            data = data[..., 0].astype(np.float32, copy=False)
        if data.ndim != 3:
            return None
        # nib 常为 (x, y, z)，转为 (z, y, x) 以便轴位取中间层 = 沿 z
        return np.transpose(data, (2, 1, 0))
    except Exception as e:
        logger.debug("nibabel load failed, try SimpleITK: %s", e)
    try:
        import SimpleITK as sitk  # type: ignore

        img = sitk.ReadImage(str(p))
        return sitk.GetArrayFromImage(img).astype(np.float32, copy=False)
    except Exception as e2:
        logger.warning("SimpleITK slice preview failed: %s", e2)
        return None


def build_radiomics_nifti_preview_markdown(
    image_path: str,
    *,
    title: str = "中心轴位预览（Axial）",
    figsize: tuple[float, float] = (4.6, 4.2),
    dpi: int = 96,
) -> str:
    """
    截取体数据中间层，matplotlib(Agg) → PNG → data URI，嵌入 Markdown。

    失败时返回空字符串（调用方仍可展示纯文本摘要）。
    """
    vol = _volume_array_from_nifti(image_path)
    if vol is None or vol.size == 0:
        return ""
    zc = int(vol.shape[0] // 2)
    sl = np.asarray(vol[zc, :, :], dtype=np.float32)
    # 稳健归一化到 0–1 便于显示
    lo, hi = float(np.nanmin(sl)), float(np.nanmax(sl))
    if not np.isfinite(lo) or not np.isfinite(hi) or hi <= lo:
        sl_disp = np.zeros_like(sl)
    else:
        sl_disp = np.clip((sl - lo) / (hi - lo), 0.0, 1.0)
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError as exc:
        logger.debug("matplotlib unavailable for NIfTI preview: %s", exc)
        return ""

    try:
        fig, ax = plt.subplots(figsize=figsize, dpi=dpi)
        ax.imshow(np.rot90(sl_disp, k=1), cmap="gray", aspect="auto")
        ax.set_title(title, fontsize=10)
        ax.axis("off")
        buf = io.BytesIO()
        fig.savefig(buf, format="png", bbox_inches="tight", pad_inches=0.08)
        plt.close(fig)
        b64 = base64.standard_b64encode(buf.getvalue()).decode("ascii")
        return f"#### {title}\n\n![](data:image/png;base64,{b64})\n\n"
    except Exception as e:
        logger.warning("matplotlib NIfTI preview render failed: %s", e)
        return ""
