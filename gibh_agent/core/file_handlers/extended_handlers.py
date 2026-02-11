"""
Extended file handlers: Spatial (Visium), Radiomics, Genomics.
Priority-based so Visium is detected before generic 10x.
"""
import logging
import os
from pathlib import Path
from typing import Dict, Any, Optional

from ..file_inspector import BaseFileHandler

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Spatial: 10x Visium (must run BEFORE generic TenXDirectoryHandler)
# ---------------------------------------------------------------------------

class SpatialVisiumHandler(BaseFileHandler):
    """
    10x Visium 目录检查器。
    条件：目录下存在 spatial/ 子目录，且存在 filtered_feature_bc_matrix.h5 或 matrix.mtx。
    与 scRNA-seq 的 10x 区分：Visium 一定有 spatial/。
    """

    def _find_visium_root(self, root_path: Path) -> Optional[Path]:
        """
        查找同时包含 spatial/ 与 (filtered_feature_bc_matrix.h5 或 matrix.mtx) 的目录。
        支持 Space Ranger 的 outs/ 或顶层结构。
        """
        if not root_path.is_dir():
            return None
        try:
            for root, dirs, files in os.walk(root_path):
                p = Path(root)
                has_spatial = "spatial" in dirs and (p / "spatial").is_dir()
                has_matrix_h5 = "filtered_feature_bc_matrix.h5" in files
                has_mtx = "matrix.mtx" in files or "matrix.mtx.gz" in files
                if has_spatial and (has_matrix_h5 or has_mtx):
                    return p
        except Exception as e:
            logger.debug("SpatialVisiumHandler walk failed: %s", e)
        return None

    def can_handle(self, path: Path) -> bool:
        if not path.is_dir():
            return False
        return self._find_visium_root(path) is not None

    def inspect(self, path: Path) -> Dict[str, Any]:
        visium_root = self._find_visium_root(path)
        if not visium_root:
            return {
                "status": "error",
                "error": "Not a valid Visium directory (missing spatial/ or matrix).",
                "file_type": "directory",
            }
        n_obs, n_vars = 0, 0
        try:
            h5 = visium_root / "filtered_feature_bc_matrix.h5"
            if h5.is_file():
                import h5py
                with h5py.File(h5, "r") as f:
                    if "matrix" in f and "shape" in f["matrix"]:
                        sh = f["matrix"]["shape"][:]
                        if len(sh) >= 2:
                            n_vars, n_obs = int(sh[0]), int(sh[1])  # shape is (genes, barcodes)
            if (n_obs, n_vars) == (0, 0):
                mtx = visium_root / "matrix.mtx"
                if not mtx.is_file():
                    mtx = visium_root / "matrix.mtx.gz"
                if mtx.is_file():
                    import gzip
                    open_fn = gzip.open if str(mtx).endswith(".gz") else open
                    with open_fn(mtx, "rt", encoding="utf-8") as f:
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
        return {
            "status": "success",
            "file_path": str(path.resolve()),
            "real_data_path": str(visium_root.resolve()),
            "file_type": "visium",
            "domain": "Spatial",
            "shape": {"rows": n_obs, "cols": n_vars},
            "n_obs": n_obs,
            "n_vars": n_vars,
            "n_samples": n_obs,
            "n_features": n_vars,
            "columns": None,
            "head": {
                "markdown": f"10x Visium 空间转录组\n- 数据根: {visium_root}\n- 含 spatial/ 与表达矩阵",
                "json": {"domain": "Spatial", "file_type": "visium"},
            },
        }


# ---------------------------------------------------------------------------
# Stubs for future modalities (do not register yet)
# ---------------------------------------------------------------------------

class RadiomicsHandler(BaseFileHandler):
    """
    Radiomics / 影像组学：.nii, .nii.gz。
    占位实现，后续扩展。
    """

    def can_handle(self, path: Path) -> bool:
        if path.is_file():
            name = path.name.lower()
            return name.endswith(".nii.gz") or name.endswith(".nii")
        return False

    def inspect(self, path: Path) -> Dict[str, Any]:
        return {
            "status": "success",
            "file_path": str(path.resolve()),
            "file_type": "nifti",
            "domain": "Radiomics",
            "shape": {},
            "head": {"markdown": "NIfTI 影像（Radiomics 占位）", "json": {}},
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
