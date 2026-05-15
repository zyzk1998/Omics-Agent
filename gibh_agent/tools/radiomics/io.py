"""
Radiomics I/O: load medical images (NIfTI, DICOM) via SimpleITK.
"""
import logging
from pathlib import Path
from typing import Dict, Any, Optional

from ...core.tool_registry import registry
from .nifti_preview import build_radiomics_nifti_preview_markdown

logger = logging.getLogger(__name__)


@registry.register(
    name="radiomics_load_medical_image",
    description="Load a medical image (NIfTI .nii/.nii.gz or DICOM) and optionally a mask. Returns metadata: size, spacing, origin.",
    category="Radiomics",
    output_type="json",
)
def load_medical_image(
    image_path: str,
    mask_path: Optional[str] = None,
) -> Dict[str, Any]:
    """
    Load medical image using SimpleITK. Return metadata for downstream tools.

    Args:
        image_path: Path to image (.nii, .nii.gz, or DICOM).
        mask_path: Optional path to segmentation mask (same format).

    Returns:
        Dict with status, size, spacing, origin, and optional error.
    """
    try:
        import SimpleITK as sitk
    except ImportError as e:
        logger.warning("SimpleITK not installed: %s", e)
        return {
            "status": "error",
            "error": "Critical Dependency Missing: SimpleITK is required for Radiomics. Please rebuild Docker or run: pip install SimpleITK",
        }

    p = Path(image_path)
    if not p.exists():
        return {"status": "error", "error": f"Image not found: {image_path}"}

    try:
        img = sitk.ReadImage(str(p))
        size = list(img.GetSize())
        spacing = list(img.GetSpacing())
        origin = list(img.GetOrigin())
        direction = list(img.GetDirection())

        out = {
            "status": "success",
            "image_path": str(p.resolve()),
            "size": size,
            "spacing": spacing,
            "origin": origin,
            "direction": direction,
        }

        if mask_path:
            mp = Path(mask_path)
            if mp.exists():
                mask = sitk.ReadImage(str(mp))
                out["mask_size"] = list(mask.GetSize())
                out["mask_path"] = str(mp.resolve())
            else:
                out["mask_path"] = None
                out["mask_warning"] = f"Mask not found: {mask_path}"
        suf = p.suffix.lower()
        _is_vol = suf in (".nii", ".gz") or ".nii" in p.name.lower() or suf == ".dcm"
        _preview_md = ""
        if _is_vol:
            _preview_md = build_radiomics_nifti_preview_markdown(str(p.resolve()), title="中心轴位预览（Axial）")
        _meta_lines = (
            f"- **文件**：`{p.name}`\n"
            f"- **体素尺寸**：{size}\n"
            f"- **间距 (mm)**：{spacing}\n"
        )
        out["markdown"] = (_preview_md + "### 影像元数据\n" + _meta_lines) if _preview_md else ("### 影像元数据\n" + _meta_lines)
        out["summary"] = out["markdown"]
        out["message"] = f"已从 `{p.name}` 加载医学影像（维度 {size}）。"
        return out
    except Exception as e:
        logger.exception("load_medical_image failed: %s", e)
        return {"status": "error", "error": str(e), "image_path": image_path}
