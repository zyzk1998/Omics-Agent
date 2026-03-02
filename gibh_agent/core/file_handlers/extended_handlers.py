"""
Extended file handlers: Spatial (Visium), Radiomics, Genomics.
Priority-based so Visium is detected before generic 10x.
SpatialVisiumHandler: can_handle only when spatial/ + valid matrix (.h5 or .mtx).
"""
import logging
import os
from pathlib import Path
from typing import Dict, Any, Optional, Tuple

from ..file_inspector import BaseFileHandler

logger = logging.getLogger(__name__)

try:
    from .structure_normalizer import pair_radiomics_files
except ImportError:
    pair_radiomics_files = None  # type: ignore

# Standard Visium matrix names (10x industry standard)
PREFERRED_MATRIX_H5 = "filtered_feature_bc_matrix.h5"
MATRIX_H5_SUFFIX = "_feature_bc_matrix.h5"
SPATIAL_IMAGE_NAMES = ("tissue_hires_image.png", "tissue_lowres_image.png", "tissue_positions_image.png")


def _is_visium_matrix_h5(name: str) -> bool:
    """True if filename is a 10x-style matrix .h5."""
    n = name.lower()
    return n == PREFERRED_MATRIX_H5 or n == "raw_feature_bc_matrix.h5" or n.endswith(MATRIX_H5_SUFFIX)


def _has_spatial_image(spatial_dir: Path) -> bool:
    """True if spatial/ contains a known tissue image."""
    if not spatial_dir.is_dir():
        return False
    return any((spatial_dir / img).is_file() for img in SPATIAL_IMAGE_NAMES)


# ---------------------------------------------------------------------------
# Spatial: 10x Visium (The Sensor)
# - Priority=20 (registered in file_inspector.py; higher than generic 10x handler).
# - can_handle: True IF AND ONLY IF directory contains spatial/ AND valid matrix (.h5 OR .mtx).
# - Metadata: file_type="visium", matrix_format="h5"|"mtx".
# ---------------------------------------------------------------------------

class SpatialVisiumHandler(BaseFileHandler):
    """
    10x Visium 目录检查器。
    条件：目录下存在 spatial/ 子目录，且存在有效矩阵（.h5 或 .mtx）。
    与 scRNA-seq 的 10x 区分：Visium 一定有 spatial/。
    元数据：file_type=visium, matrix_format=h5|mtx, has_image=True|False。
    """

    def _find_visium_root_and_matrix(
        self, root_path: Path
    ) -> Optional[Tuple[Path, str, Optional[Path]]]:
        """
        查找同时包含 spatial/ 与有效矩阵的目录。
        返回 (visium_root, matrix_format, matrix_path) 或 None。
        matrix_format: "h5" | "mtx"
        matrix_path: Path to the matrix file (for .h5, may be filtered or raw).
        """
        if not root_path.is_dir():
            return None
        try:
            for root, dirs, files in os.walk(root_path):
                p = Path(root)
                has_spatial = "spatial" in dirs and (p / "spatial").is_dir()
                if not has_spatial:
                    continue
                # Prefer .h5 (industry standard); accept any *_feature_bc_matrix.h5 or matrix.mtx
                matrix_h5 = None
                for f in files:
                    if _is_visium_matrix_h5(f):
                        matrix_h5 = p / f
                        break
                if matrix_h5 is not None:
                    return (p, "h5", matrix_h5)
                if "matrix.mtx" in files:
                    return (p, "mtx", p / "matrix.mtx")
                if "matrix.mtx.gz" in files:
                    return (p, "mtx", p / "matrix.mtx.gz")
        except Exception as e:
            logger.debug("SpatialVisiumHandler walk failed: %s", e)
        return None

    def can_handle(self, path: Path) -> bool:
        """Return True IF AND ONLY IF directory contains spatial/ AND a valid matrix (.h5 OR .mtx)."""
        if not path.is_dir():
            return False
        return self._find_visium_root_and_matrix(path) is not None

    def inspect(self, path: Path) -> Dict[str, Any]:
        result = self._find_visium_root_and_matrix(path)
        if not result:
            return {
                "status": "error",
                "error": "Not a valid Visium directory (missing spatial/ or valid matrix .h5/.mtx).",
                "file_type": "directory",
            }
        visium_root, matrix_format, matrix_path = result
        spatial_dir = visium_root / "spatial"
        has_image = _has_spatial_image(spatial_dir)

        n_obs, n_vars = 0, 0
        try:
            if matrix_format == "h5" and matrix_path and matrix_path.is_file():
                import h5py
                with h5py.File(matrix_path, "r") as f:
                    if "matrix" in f and "shape" in f["matrix"]:
                        sh = f["matrix"]["shape"][:]
                        if len(sh) >= 2:
                            n_vars, n_obs = int(sh[0]), int(sh[1])
            if (n_obs, n_vars) == (0, 0) and matrix_format == "mtx" and matrix_path and matrix_path.is_file():
                import gzip
                open_fn = gzip.open if str(matrix_path).endswith(".gz") else open
                with open_fn(matrix_path, "rt", encoding="utf-8") as f:
                    for line in f:
                        line = line.strip()
                        if line.startswith("%") or line.startswith("#"):
                            continue
                        parts = line.split()
                        if len(parts) >= 2:
                            try:
                                n_obs, n_vars = int(parts[0]), int(parts[1])
                                break
                            except ValueError:
                                pass
        except Exception as e:
            logger.debug("SpatialVisiumHandler shape detection failed: %s", e)
        if n_obs == 0 and n_vars == 0:
            n_obs = n_vars = -1

        # Prefer passing standard name for workflow; loader uses filtered_feature_bc_matrix.h5 by default
        matrix_file_name = matrix_path.name if matrix_path else (
            "filtered_feature_bc_matrix.h5" if matrix_format == "h5" else "matrix.mtx"
        )

        return {
            "status": "success",
            "file_path": str(path.resolve()),
            "real_data_path": str(visium_root.resolve()),
            "file_type": "visium",
            "domain": "Spatial",
            "matrix_format": matrix_format,
            "matrix_file": matrix_file_name,
            "has_image": has_image,
            "shape": {"rows": n_obs, "cols": n_vars},
            "n_obs": n_obs,
            "n_vars": n_vars,
            "n_samples": n_obs,
            "n_features": n_vars,
            "columns": None,
            "head": {
                "markdown": (
                    f"10x Visium 空间转录组\n- 数据根: {visium_root}\n"
                    f"- 矩阵格式: {matrix_format}\n- 含 spatial/ 与表达矩阵"
                    + ("\n- 含组织图像" if has_image else "")
                ),
                "json": {
                    "domain": "Spatial",
                    "file_type": "visium",
                    "matrix_format": matrix_format,
                    "has_image": has_image,
                },
            },
        }


# ---------------------------------------------------------------------------
# Radiomics: medical imaging (.nii, .nii.gz, .dcm) — Priority=10 (set in file_inspector)
# ---------------------------------------------------------------------------

def _read_radiomics_header(image_path: Path) -> tuple:
    """Read image header (size, spacing, origin, direction). Returns (size, spacing, origin, direction)."""
    size, spacing, origin, direction = [], [], [], []
    try:
        import SimpleITK as sitk
        img = sitk.ReadImage(str(image_path))
        size = list(img.GetSize())
        spacing = list(img.GetSpacing())
        origin = list(img.GetOrigin())
        direction = list(img.GetDirection())
    except Exception as e:
        logger.debug("RadiomicsHandler header read failed: %s", e)
    return (size, spacing, origin, direction)


class RadiomicsHandler(BaseFileHandler):
    """
    Radiomics / 影像组学：.nii, .nii.gz, .dcm。
    file_type="medical_image", domain="Radiomics".
    When given a file, scans its directory for a paired mask; when given a directory,
    finds image + mask and returns file_path (image) and mask_path (or None).
    """

    def can_handle(self, path: Path) -> bool:
        if path.is_file():
            name = path.name.lower()
            return name.endswith(".nii.gz") or name.endswith(".nii") or name.endswith(".dcm")
        if path.is_dir() and pair_radiomics_files:
            image_path, _ = pair_radiomics_files(path)
            return image_path is not None
        return False

    def inspect(self, path: Path) -> Dict[str, Any]:
        # Resolve image + mask: single file -> pair from same dir; directory -> pair from dir
        if path.is_file():
            search_dir = path.parent
            image_path = path
            mask_path = None
            if pair_radiomics_files:
                img_cand, mask_cand = pair_radiomics_files(search_dir)
                if img_cand is not None:
                    image_path = img_cand
                if mask_cand is not None:
                    mask_path = mask_cand
        else:
            image_path = mask_path = None
            if pair_radiomics_files:
                image_path, mask_path = pair_radiomics_files(path)
            if not image_path:
                return {
                    "status": "error",
                    "error": "未在目录中发现影像文件（.nii/.nii.gz/.dcm）",
                    "file_type": "unknown",
                    "domain": "Radiomics",
                }
        size, spacing, origin, direction = _read_radiomics_header(image_path)
        abs_image = str(image_path.resolve())
        abs_mask = str(mask_path.resolve()) if mask_path else None
        head_md = (
            f"医学影像（Radiomics）\n- 影像: {image_path.name}\n- 尺寸: {size or 'N/A'}\n- 间距: {spacing or 'N/A'}"
            + (f"\n- 掩膜: {mask_path.name}" if mask_path else "\n- 掩膜: 无")
        )
        return {
            "status": "success",
            "file_path": abs_image,
            "mask_path": abs_mask,
            "file_type": "medical_image",
            "domain": "Radiomics",
            "modality": "Radiomics",
            "shape": {"size": size, "spacing": spacing, "origin": origin, "direction": direction},
            "head": {
                "markdown": head_md,
                "json": {"domain": "Radiomics", "file_type": "medical_image", "mask_path": abs_mask},
            },
        }


class GenomicsHandler(BaseFileHandler):
    """
    基因组学：.vcf。
    占位实现，后续扩展。
    """

    def can_handle(self, path: Path) -> bool:
        if path.is_file():
            return path.name.lower().endswith(".vcf") or path.name.lower().endswith(".vcf.gz")
        return False

    def inspect(self, path: Path) -> Dict[str, Any]:
        return {
            "status": "success",
            "file_path": str(path.resolve()),
            "file_type": "vcf",
            "domain": "Genomics",
            "shape": {},
            "head": {"markdown": "VCF 变异文件（Genomics 占位）", "json": {}},
        }
