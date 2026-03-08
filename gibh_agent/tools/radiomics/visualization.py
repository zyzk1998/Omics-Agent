"""
Radiomics visualization: export mid-slice PNG for UI (browsers cannot show .nii.gz directly).
When output_path is not provided, saves under RESULTS_DIR so Nginx can serve at /results/.
"""
import logging
import os
from pathlib import Path
from typing import Dict, Any, Optional

from ...core.tool_registry import registry

logger = logging.getLogger(__name__)


@registry.register(
    name="radiomics_plot_mid_slice",
    description="Extract the middle slice of a 3D volume, overlay mask if provided, and save as PNG for display in the UI.",
    category="Radiomics",
    output_type="json",
)
def plot_mid_slice(
    image_path: str,
    output_path: Optional[str] = None,
    mask_path: Optional[str] = None,
) -> Dict[str, Any]:
    """
    Export middle slice of 3D NIfTI (and optional mask overlay) as PNG for UI.

    Args:
        image_path: Path to .nii or .nii.gz image.
        output_path: Optional path for output PNG; if None, saves next to image as radiomics_preview.png.
        mask_path: Optional mask; will be overlaid in red.

    Returns:
        Dict with status, output_path, and optional error.
    """
    try:
        import SimpleITK as sitk
        import numpy as np
    except ImportError as e:
        logger.warning("SimpleITK/numpy not installed: %s", e)
        return {"status": "error", "error": "SimpleITK and numpy are required."}

    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError as e:
        logger.warning("matplotlib not installed: %s", e)
        return {"status": "error", "error": "matplotlib is required for PNG export."}

    p = Path(image_path)
    if not p.exists():
        return {"status": "error", "error": f"Image not found: {image_path}"}

    try:
        img = sitk.ReadImage(str(p))
        arr = sitk.GetArrayFromImage(img)
        if arr.ndim != 3:
            return {"status": "error", "error": f"Expected 3D volume, got ndim={arr.ndim}"}

        mid_idx = arr.shape[0] // 2
        slice_2d = arr[mid_idx, :, :]

        if output_path and str(output_path).strip():
            out_p = Path(output_path)
        else:
            results_dir = os.environ.get("RESULTS_DIR", "").strip()
            if results_dir and Path(results_dir).exists():
                out_p = Path(results_dir) / "radiomics_preview.png"
            else:
                out_p = p.parent / "radiomics_preview.png"
        if out_p.suffix.lower() not in (".png", ".jpg", ".jpeg"):
            out_p = out_p.parent / (out_p.name + ".png")
        out_p.parent.mkdir(parents=True, exist_ok=True)

        fig, ax = plt.subplots(1, 1, figsize=(6, 6))
        ax.imshow(slice_2d, cmap="gray", aspect="equal")
        if mask_path and Path(mask_path).exists():
            mask = sitk.ReadImage(mask_path)
            mask_arr = sitk.GetArrayFromImage(mask)
            if mask_arr.shape == arr.shape:
                mask_slice = mask_arr[mid_idx, :, :]
                ax.contour(mask_slice, levels=[0.5], colors=["red"], linewidths=1)
            else:
                logger.warning("Mask shape %s != image shape %s", mask_arr.shape, arr.shape)
        ax.set_title("Mid-slice preview")
        ax.axis("off")
        plt.tight_layout()
        plt.savefig(str(out_p), dpi=100, bbox_inches="tight")
        plt.close()

        return {
            "status": "success",
            "output_path": str(out_p.resolve()),
            "image_path": image_path,
            "slice_index": mid_idx,
        }
    except Exception as e:
        logger.exception("plot_mid_slice failed: %s", e)
        return {"status": "error", "error": str(e), "image_path": image_path}
