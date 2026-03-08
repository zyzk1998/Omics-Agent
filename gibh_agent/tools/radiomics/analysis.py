"""
Radiomics analysis: extract texture and shape features via PyRadiomics.
"""
import logging
from pathlib import Path
from typing import Dict, Any, Optional

from ...core.tool_registry import registry

logger = logging.getLogger(__name__)


def _create_full_image_mask(image_path: str, out_mask_path: Optional[Path] = None) -> Optional[Path]:
    """
    Create a full-image binary mask (all voxels = 1) for PyRadiomics when no ROI is provided.
    Returns path to the written mask file, or None on failure.
    """
    try:
        import SimpleITK as sitk
        import tempfile
        img = sitk.ReadImage(image_path)
        arr = sitk.GetArrayFromImage(img)
        mask_arr = (arr != 0).astype("uint8")  # Non-zero as ROI; or use ones_like for full volume
        if mask_arr.sum() == 0:
            mask_arr = (arr * 0 + 1).astype("uint8")  # full volume
        mask_img = sitk.GetImageFromArray(mask_arr)
        mask_img.CopyInformation(img)
        if out_mask_path is None:
            fd, p = tempfile.mkstemp(suffix=".nii.gz", prefix="radiomics_fullmask_")
            import os
            os.close(fd)
            out_mask_path = Path(p)
        sitk.WriteImage(mask_img, str(out_mask_path))
        return out_mask_path
    except Exception as e:
        logger.warning("Create full-image mask failed: %s", e)
        return None


@registry.register(
    name="radiomics_extract_features",
    description="Extract radiomics features (shape, first-order, texture) from a medical image and optional mask using PyRadiomics. If no mask provided, uses a full-image mask. Returns path to CSV of features.",
    category="Radiomics",
    output_type="json",
)
def extract_radiomics_features(
    image_path: str,
    mask_path: Optional[str] = None,
    output_path: Optional[str] = None,
    config_path: Optional[str] = None,
) -> Dict[str, Any]:
    """
    Extract radiomics features with PyRadiomics. Uses default settings if no YAML config.
    If mask_path is None or missing, generates a full-image mask so the tool completes gracefully.

    Args:
        image_path: Path to image (.nii, .nii.gz).
        mask_path: Optional path to segmentation mask; if None or file missing, use full-image mask.
        output_path: Optional path for CSV output; default next to image.
        config_path: Optional PyRadiomics YAML config; if not provided, use defaults.

    Returns:
        Dict with status, csv_path, n_features, and optional error.
    """
    try:
        import SimpleITK as sitk  # noqa: F401
    except ImportError as e:
        logger.warning("SimpleITK not installed: %s", e)
        return {
            "status": "error",
            "error": "Critical Dependency Missing: SimpleITK is required for Radiomics. Please rebuild Docker or run: pip install SimpleITK",
        }
    try:
        from radiomics import featureextractor
    except ImportError as e:
        logger.warning("pyradiomics import failed: %s", e)
        err_detail = str(e).strip() or "unknown"
        return {
            "status": "error",
            "error": f"Critical Dependency Missing: pyradiomics is required for Radiomics. Import error: {err_detail}. Run in api-server container: pip show pyradiomics. If missing there, restart api-server after no-cache rebuild.",
        }

    img_p = Path(image_path)
    if not img_p.exists():
        return {"status": "error", "error": f"Image not found: {image_path}"}

    mask_to_use: Optional[Path] = None
    temp_mask: Optional[Path] = None
    if mask_path and Path(mask_path).exists():
        mask_to_use = Path(mask_path)
    else:
        if mask_path:
            logger.info("Mask not found or empty; generating full-image mask for %s", image_path)
        temp_mask = _create_full_image_mask(image_path)
        if temp_mask is None:
            return {"status": "error", "error": "No mask provided and could not generate full-image mask. Please provide a valid mask_path.", "image_path": image_path}
        mask_to_use = temp_mask

    try:
        if config_path and Path(config_path).exists():
            extractor = featureextractor.RadiomicsFeatureExtractor(config_path)
        else:
            extractor = featureextractor.RadiomicsFeatureExtractor()
            extractor.disableAllFeatures()
            extractor.enableFeatureClassByName("shape", True)
            extractor.enableFeatureClassByName("firstorder", True)
            extractor.enableFeatureClassByName("glcm", True)
            extractor.enableFeatureClassByName("glrlm", True)
            extractor.enableFeatureClassByName("gldm", True)
            extractor.enableFeatureClassByName("ngtdm", True)
            extractor.enableFeatureClassByName("glszm", True)

        result = extractor.execute(str(img_p), str(mask_to_use))
        if not result:
            return {"status": "error", "error": "PyRadiomics returned no features.", "image_path": image_path, "mask_path": str(mask_to_use)}

        out_dir = img_p.parent
        out_name = img_p.stem.replace(".nii", "") + "_radiomics_features.csv"
        csv_path = Path(output_path) if output_path else (out_dir / out_name)
        csv_path = csv_path.resolve()
        if csv_path.suffix.lower() != ".csv":
            csv_path = csv_path.parent / (csv_path.name + ".csv")
        csv_path.parent.mkdir(parents=True, exist_ok=True)

        import csv
        with open(csv_path, "w", newline="", encoding="utf-8") as f:
            writer = csv.writer(f)
            writer.writerow(["feature", "value"])
            for k, v in sorted(result.items()):
                if not k.startswith("diagnostics"):
                    writer.writerow([k, v])

        n_features = sum(1 for k in result if not k.startswith("diagnostics"))
        out_result = {
            "status": "success",
            "csv_path": str(csv_path),
            "n_features": n_features,
            "image_path": image_path,
            "mask_path": str(mask_to_use),
        }
        if temp_mask and temp_mask.exists():
            try:
                temp_mask.unlink()
            except Exception:
                pass
        return out_result
    except Exception as e:
        logger.exception("extract_radiomics_features failed: %s", e)
        if temp_mask and temp_mask.exists():
            try:
                temp_mask.unlink()
            except Exception:
                pass
        return {"status": "error", "error": str(e), "image_path": image_path, "mask_path": mask_path or ""}
