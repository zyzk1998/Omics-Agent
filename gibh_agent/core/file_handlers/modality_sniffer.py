"""
Universal Ingestion Pipeline — Phase 2: Modality Sniffer

After unpacking, determines "what is this?" using heuristic rules (priority order).
Returns dominant modality and the path(s) to use for inspection/execution.
"""
import logging
from pathlib import Path
from typing import Any, List, Optional, Tuple

logger = logging.getLogger(__name__)


def _has_spatial_dir_and_matrix(root: Path) -> bool:
    """True if root contains spatial/ and a valid matrix (.h5 or matrix.mtx)."""
    spatial = root / "spatial"
    if not spatial.is_dir():
        return False
    for p in root.iterdir():
        if p.is_file():
            name = p.name.lower()
            if name.endswith(".h5") and "feature" in name and "matrix" in name:
                return True
            if name in ("matrix.mtx", "matrix.mtx.gz"):
                return True
    return False


def _find_radiomics(root: Path) -> Tuple[Optional[Path], Optional[Path]]:
    """Find image and mask in root (and one level of subdirs). Returns (image_path, mask_path)."""
    try:
        from .structure_normalizer import pair_radiomics_files
    except ImportError:
        pair_radiomics_files = None
    if not pair_radiomics_files:
        # Fallback: any .nii.gz/.nii/.dcm
        exts = (".nii.gz", ".nii", ".dcm")
        candidates = [p for p in root.iterdir() if p.is_file() and p.name.lower().endswith(exts)]
        if not candidates:
            for d in root.iterdir():
                if d.is_dir():
                    candidates = [p for p in d.iterdir() if p.is_file() and p.name.lower().endswith(exts)]
                    if candidates:
                        break
        if not candidates:
            return (None, None)
        return (candidates[0], None)
    image_path, mask_path = pair_radiomics_files(root)
    if image_path is not None:
        return (image_path, mask_path)
    for d in root.iterdir():
        if d.is_dir():
            image_path, mask_path = pair_radiomics_files(d)
            if image_path is not None:
                return (image_path, mask_path)
    return (None, None)


def _has_scrna_10x(root: Path) -> bool:
    """True if root contains matrix.mtx and barcodes.tsv (and optionally features/genes.tsv)."""
    files = {p.name for p in root.iterdir() if p.is_file()}
    if "matrix.mtx" in files and "barcodes.tsv" in files:
        return True
    for d in root.iterdir():
        if d.is_dir():
            sub = {p.name for p in d.iterdir() if p.is_file()}
            if "matrix.mtx" in sub and "barcodes.tsv" in sub:
                return True
    return False


def _find_vcf(root: Path) -> Optional[Path]:
    """Return first .vcf or .vcf.gz in root or one level down."""
    for p in root.iterdir():
        if p.is_file() and p.name.lower().endswith((".vcf", ".vcf.gz")):
            return p
    for d in root.iterdir():
        if d.is_dir():
            for p in d.iterdir():
                if p.is_file() and p.name.lower().endswith((".vcf", ".vcf.gz")):
                    return p
    return None


def _find_tabular(root: Path) -> Optional[Path]:
    """Return first .csv or .xlsx in root or one level down."""
    for p in root.iterdir():
        if p.is_file() and p.name.lower().endswith((".csv", ".xlsx")):
            return p
    for d in root.iterdir():
        if d.is_dir():
            for p in d.iterdir():
                if p.is_file() and p.name.lower().endswith((".csv", ".xlsx")):
                    return p
    return None


def detect_dominant_modality(session_dir: Path) -> Tuple[str, Any]:
    """
    Heuristic: decide dominant modality and return path(s) for inspection/execution.

    Priority:
    1. Spatial: spatial/ + .h5 or .mtx -> ("Spatial", dir_path)
    2. Radiomics: .nii/.nii.gz/.dcm -> ("Radiomics", image_path, mask_path|None)
    3. scRNA-seq: matrix.mtx + barcodes.tsv -> ("scRNA-seq", dir_path)
    4. Genomics: .vcf -> ("Genomics", vcf_path)
    5. Tabular: .csv/.xlsx -> ("Tabular", csv_path)
    6. Else -> ("unknown", None)

    Returns:
        (modality, payload) where payload is:
        - Spatial: Path (directory)
        - Radiomics: (image_path: Path, mask_path: Optional[Path])
        - scRNA-seq: Path (directory)
        - Genomics: Path (file)
        - Tabular: Path (file)
        - unknown: None
    """
    if not session_dir.is_dir():
        return ("unknown", None)
    root = session_dir.resolve()

    # 1. Spatial (root or one level of subdirs, so Visium inside extracted_content/SampleName works)
    if _has_spatial_dir_and_matrix(root):
        logger.info("ModalitySniffer: detected Spatial (spatial/ + matrix) at root")
        return ("Spatial", root)
    for sub in root.iterdir():
        if sub.is_dir() and _has_spatial_dir_and_matrix(sub):
            logger.info("ModalitySniffer: detected Spatial in subdir %s", sub.name)
            return ("Spatial", sub)

    # 2. Radiomics
    image_path, mask_path = _find_radiomics(root)
    if image_path is not None:
        logger.info("ModalitySniffer: detected Radiomics image=%s mask=%s", image_path.name, mask_path.name if mask_path else None)
        return ("Radiomics", (image_path, mask_path))

    # 3. scRNA-seq (10x)
    if _has_scrna_10x(root):
        logger.info("ModalitySniffer: detected scRNA-seq (matrix.mtx + barcodes)")
        return ("scRNA-seq", root)

    # 4. Genomics
    vcf_path = _find_vcf(root)
    if vcf_path is not None:
        logger.info("ModalitySniffer: detected Genomics vcf=%s", vcf_path.name)
        return ("Genomics", vcf_path)

    # 5. Tabular
    tab_path = _find_tabular(root)
    if tab_path is not None:
        logger.info("ModalitySniffer: detected Tabular file=%s", tab_path.name)
        return ("Tabular", tab_path)

    return ("unknown", None)


def paths_for_response(
    modality: str,
    payload: Any,
    upload_dir: Path,
) -> List[Path]:
    """
    Convert (modality, payload) into a list of absolute paths to return to the API.
    Used for path rewriting after unpack + sniffer.
    """
    if payload is None:
        return []
    upload_dir = upload_dir.resolve()
    out: List[Path] = []
    if modality == "Spatial":
        out.append(payload)
    elif modality == "Radiomics":
        image_path, mask_path = payload
        out.append(image_path)
        if mask_path:
            out.append(mask_path)
    elif modality in ("scRNA-seq", "Genomics", "Tabular"):
        out.append(payload)
    return [p.resolve() for p in out if p is not None and (p.is_file() or p.is_dir())]
