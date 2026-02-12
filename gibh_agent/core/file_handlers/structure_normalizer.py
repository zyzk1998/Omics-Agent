"""
SpatialStructureNormalizer (The Assembler): reorganize loose uploads into 10x Visium directory structure.

Target state (enforced before analysis):
  session_dir/
    ├── filtered_feature_bc_matrix.h5   (expression data; sibling of spatial/)
    └── spatial/
         ├── tissue_hires_image.png
         └── scalefactors_json.json

Logic (Phase 1 - Complete Body):
  1. Scan for spatial.tar.gz / spatial.zip; extract to a folder named spatial/.
  2. Scan for *.h5 (e.g. filtered_feature_bc_matrix.h5 or raw_...h5).
  3. Action: move any .h5 from subdirs to session_dir root so it is a *sibling* of spatial/.
  4. Ensure filtered_feature_bc_matrix.h5 exists (copy from raw/other if needed).
  Result: root/spatial/ and root/matrix.h5 (industry standard).
"""
import logging
import shutil
import tarfile
import zipfile
from pathlib import Path
from typing import List, Optional, Tuple

logger = logging.getLogger(__name__)

# Radiomics: extensions and naming for image vs mask
RADIOMICS_EXTENSIONS = (".nii.gz", ".nii", ".dcm")
RADIOMICS_MASK_KEYWORDS = ("mask", "label", "seg", "segmentation", "roi")

# Archive extensions and naming patterns for "spatial" pack
SPATIAL_ARCHIVE_SUFFIXES = (".tar.gz", ".tgz", ".tar", ".zip")
SPATIAL_NAME_PATTERNS = ("spatial", "spatial_images", "spatial_metadata")

# Standard 10x Visium matrix filename (industry standard)
PREFERRED_MATRIX_H5 = "filtered_feature_bc_matrix.h5"
RAW_MATRIX_H5 = "raw_feature_bc_matrix.h5"
MATRIX_H5_SUFFIX = "_feature_bc_matrix.h5"


def _is_spatial_archive(path: Path) -> bool:
    """True if path looks like a spatial archive (name + extension)."""
    if not path.is_file():
        return False
    name_lower = path.name.lower()
    if name_lower.endswith(".spatial.tar.gz") or name_lower == "spatial.tar.gz":
        return True
    if name_lower.endswith(".spatial.zip") or name_lower == "spatial.zip":
        return True
    if any(name_lower.endswith(s) for s in SPATIAL_ARCHIVE_SUFFIXES) and any(
        p in name_lower for p in SPATIAL_NAME_PATTERNS
    ):
        return True
    return False


def _is_visium_matrix_h5(path: Path) -> bool:
    """True if path is a 10x-style matrix .h5 (filtered, raw, or any *feature_bc_matrix.h5)."""
    if not path.is_file() or path.suffix.lower() != ".h5":
        return False
    name_lower = path.name.lower()
    return (
        name_lower == PREFERRED_MATRIX_H5
        or name_lower == RAW_MATRIX_H5
        or name_lower.endswith(MATRIX_H5_SUFFIX)
    )


def _copy_h5_from_parent(session_dir: Path) -> bool:
    """
    If session_dir has spatial/ or spatial-like content but no matrix .h5,
    copy any *_feature_bc_matrix.h5 from session_dir.parent into session_dir.
    (Typical case: user uploads spatial.tar.gz + .h5; unpack puts tar in extracted_content/X, .h5 stays in parent.)
    Returns True if a copy was made.
    """
    root_h5 = [f for f in session_dir.iterdir() if f.is_file() and _is_visium_matrix_h5(f)]
    if root_h5:
        return False
    has_spatial = (session_dir / "spatial").is_dir()
    if not has_spatial:
        for f in session_dir.iterdir():
            if f.is_file() and ("tissue_positions" in f.name.lower() or "scalefactors" in f.name.lower()):
                has_spatial = True
                break
    if not has_spatial:
        return False
    # 先查 parent，再查 parent.parent（分体上传时 .h5 在 user_dir，解压目录在 extracted_content/X）
    for parent in (session_dir.parent, session_dir.parent.parent):
        if not parent or not parent.is_dir():
            continue
        for f in parent.iterdir():
            if f.is_file() and _is_visium_matrix_h5(f):
                dest = session_dir / f.name
                if dest.exists():
                    return False
                try:
                    shutil.copy2(str(f), str(dest))
                    logger.info(
                        "SpatialStructureNormalizer: copied %s from %s into %s",
                        f.name,
                        parent.name,
                        session_dir.name,
                    )
                    return True
                except Exception as e:
                    logger.warning("SpatialStructureNormalizer: copy .h5 from parent failed: %s", e)
                    return False
    return False


def _move_h5_to_root(session_dir: Path) -> bool:
    """
    Move any Visium matrix .h5 from subdirs to session_dir root so it is a *sibling* of spatial/.
    Only moves if the .h5 is not already in session_dir root. Returns True if any move was done.
    """
    root_h5 = [f for f in session_dir.iterdir() if f.is_file() and _is_visium_matrix_h5(f)]
    if root_h5:
        return False  # already at least one .h5 at root
    for child in session_dir.iterdir():
        if not child.is_dir():
            continue
        for f in child.iterdir():
            if f.is_file() and _is_visium_matrix_h5(f):
                dest = session_dir / f.name
                if dest == f:
                    continue
                if dest.exists():
                    dest = session_dir / f"{f.stem}_moved{f.suffix}"
                try:
                    shutil.move(str(f), str(dest))
                    logger.info(
                        "SpatialStructureNormalizer: moved %s to root (sibling of spatial/)",
                        f.name,
                    )
                    return True
                except Exception as e:
                    logger.warning("SpatialStructureNormalizer: move .h5 to root failed: %s", e)
    return False


def _ensure_preferred_matrix_h5(session_dir: Path) -> bool:
    """
    Ensure session_dir contains filtered_feature_bc_matrix.h5.
    If only raw_feature_bc_matrix.h5 or another *_feature_bc_matrix.h5 exists, copy it to
    filtered_feature_bc_matrix.h5 so downstream always sees the standard name.
    Returns True if a change was made.
    """
    preferred = session_dir / PREFERRED_MATRIX_H5
    if preferred.is_file():
        return False

    candidates: List[Path] = []
    for f in session_dir.iterdir():
        if f.is_file() and _is_visium_matrix_h5(f):
            if f.name == RAW_MATRIX_H5:
                candidates.insert(0, f)  # prefer raw first if no filtered
            elif PREFERRED_MATRIX_H5 not in [p.name for p in candidates]:
                candidates.append(f)
    if not candidates:
        return False

    source = candidates[0]
    try:
        shutil.copy2(str(source), str(preferred))
        logger.info(
            "SpatialStructureNormalizer: ensured %s from %s",
            PREFERRED_MATRIX_H5,
            source.name,
        )
        return True
    except Exception as e:
        logger.warning("SpatialStructureNormalizer: copy to %s failed: %s", preferred, e)
        return False


def _list_archive_top_level(archive_path: Path) -> Optional[List[str]]:
    """List top-level names in archive (single folder vs loose files)."""
    name_lower = archive_path.name.lower()
    try:
        if name_lower.endswith(".zip"):
            with zipfile.ZipFile(archive_path, "r") as z:
                names = set()
                for i in z.namelist():
                    top = i.split("/")[0].split("\\")[0]
                    if top and top != ".":
                        names.add(top)
                return list(names)
        if name_lower.endswith((".tar.gz", ".tgz", ".tar")):
            with tarfile.open(archive_path, "r:*") as t:
                names = [m.name.split("/")[0].split("\\")[0] for m in t.getmembers() if m.name]
                return list(dict.fromkeys(names))
    except Exception as e:
        logger.warning("Could not list archive %s: %s", archive_path, e)
    return None


def _extract_zip(zip_path: Path, dest_dir: Path) -> bool:
    try:
        with zipfile.ZipFile(zip_path, "r") as z:
            z.extractall(dest_dir)
        return True
    except Exception as e:
        logger.warning("Zip extract failed %s: %s", zip_path, e)
        return False


def _extract_tar(tar_path: Path, dest_dir: Path) -> bool:
    try:
        with tarfile.open(tar_path, "r:*") as t:
            t.extractall(dest_dir)
        return True
    except Exception as e:
        logger.warning("Tar extract failed %s: %s", tar_path, e)
        return False


def _is_radiomics_file(path: Path) -> bool:
    """True if path is a radiomics-capable file (.nii.gz, .nii, .dcm)."""
    if not path.is_file():
        return False
    return path.name.lower().endswith(RADIOMICS_EXTENSIONS)


def _is_radiomics_mask_name(name: str) -> bool:
    """True if filename suggests a segmentation mask (mask/label/seg)."""
    n = name.lower()
    return any(kw in n for kw in RADIOMICS_MASK_KEYWORDS)


def pair_radiomics_files(session_dir: Path) -> Tuple[Optional[Path], Optional[Path]]:
    """
    Identify Image and Mask files in a directory (Radiomics pairing).
    - Image: usually larger or lacks "mask"/"label"/"seg" in name.
    - Mask: contains "mask", "label", "seg", "segmentation", "roi" in name.

    Returns:
        (image_path, mask_path). mask_path may be None if no mask found.
    """
    if not session_dir.is_dir():
        return (None, None)
    candidates = [p for p in session_dir.iterdir() if _is_radiomics_file(p)]
    if not candidates:
        return (None, None)
    images = [p for p in candidates if not _is_radiomics_mask_name(p.name)]
    masks = [p for p in candidates if _is_radiomics_mask_name(p.name)]
    # Prefer image without mask/label in name; fallback to larger file as image
    if images:
        image_path = images[0]
        if len(images) > 1:
            try:
                image_path = max(images, key=lambda p: p.stat().st_size)
            except OSError:
                image_path = images[0]
    else:
        image_path = masks[0] if masks else None  # single file: treat as image
        masks = masks[1:] if len(masks) > 1 else []
    mask_path = masks[0] if masks else None
    return (image_path, mask_path)


def normalize_session_directory(session_dir: Path) -> bool:
    """
    Reorganize session_dir to the standard Visium layout:
      session_dir/
        ├── filtered_feature_bc_matrix.h5
        └── spatial/
             ├── tissue_hires_image.png
             └── scalefactors_json.json

    - Detects spatial.tar.gz, spatial.zip, or archives with "spatial" in the name; extracts into spatial/.
    - Detects loose .h5 matrix (filtered/raw/*_feature_bc_matrix.h5); ensures filtered_feature_bc_matrix.h5 exists (copy from raw/other if needed).
    - .h5 remains adjacent to spatial/ (both in session_dir).

    Args:
        session_dir: Directory that may contain loose files (e.g. matrix .h5 + spatial.tar.gz).

    Returns:
        True if any change was made, False otherwise.
    """
    if not session_dir.is_dir():
        return False

    session_dir = session_dir.resolve()
    changed = False

    # 0) If we have spatial/ (or spatial-like files) but no .h5, copy from parent (split upload: tar + .h5)
    if _copy_h5_from_parent(session_dir):
        changed = True

    spatial_dir = session_dir / "spatial"
    # 0.5) 若根目录有 tissue_positions_list.csv 等但无 spatial/ 子目录，创建 spatial/ 并移入（部分 tar 解压结构）
    if not spatial_dir.is_dir():
        spatial_like = [
            f for f in session_dir.iterdir()
            if f.is_file()
            and (
                "tissue_positions" in f.name.lower()
                or "scalefactors" in f.name.lower()
                or f.suffix.lower() in (".json", ".csv")
                and "spatial" in f.name.lower()
            )
        ]
        if spatial_like:
            spatial_dir.mkdir(parents=True, exist_ok=True)
            for f in spatial_like:
                try:
                    shutil.move(str(f), str(spatial_dir / f.name))
                    changed = True
                except Exception as e:
                    logger.warning("SpatialStructureNormalizer: move %s into spatial/ failed: %s", f.name, e)
            if changed:
                logger.info("SpatialStructureNormalizer: 已从根目录创建 spatial/ 并移入 %s 个文件", len(spatial_like))

    # 1) Extract spatial archive if present
    if not (spatial_dir.is_dir() and any(spatial_dir.iterdir())):
        candidates: List[Path] = [
            f for f in session_dir.iterdir()
            if f.is_file() and _is_spatial_archive(f)
        ]
        if candidates:
            archive_path = candidates[0]
            if len(candidates) > 1:
                logger.info(
                    "SpatialStructureNormalizer: multiple spatial archives, using %s",
                    archive_path.name,
                )
            top_level = _list_archive_top_level(archive_path)
            if top_level:
                extract_to = session_dir / f"_spatial_extract_{archive_path.stem}"
                extract_to.mkdir(parents=True, exist_ok=True)
                try:
                    if archive_path.suffix.lower() == ".zip":
                        ok = _extract_zip(archive_path, extract_to)
                    else:
                        ok = _extract_tar(archive_path, extract_to)
                    if ok:
                        subdirs = [p for p in extract_to.iterdir() if p.is_dir()]
                        files = [p for p in extract_to.iterdir() if p.is_file()]
                        if len(subdirs) == 1 and subdirs[0].name.lower() == "spatial":
                            src_spatial = subdirs[0]
                            if spatial_dir.exists():
                                shutil.rmtree(spatial_dir, ignore_errors=True)
                            shutil.move(str(src_spatial), str(spatial_dir))
                            changed = True
                            logger.info(
                                "SpatialStructureNormalizer: moved extracted spatial/ into %s",
                                session_dir,
                            )
                        elif subdirs or files:
                            if spatial_dir.exists():
                                shutil.rmtree(spatial_dir, ignore_errors=True)
                            spatial_dir.mkdir(parents=True, exist_ok=True)
                            for p in subdirs + files:
                                try:
                                    dest = spatial_dir / p.name
                                    if dest.exists():
                                        dest = spatial_dir / f"{p.stem}_1{p.suffix}"
                                    shutil.move(str(p), str(dest))
                                    changed = True
                                except Exception as e:
                                    logger.warning(
                                        "SpatialStructureNormalizer: move %s into spatial/ failed: %s",
                                        p,
                                        e,
                                    )
                            if changed:
                                logger.info(
                                    "SpatialStructureNormalizer: created spatial/ and moved loose content in %s",
                                    session_dir,
                                )
                    if extract_to.exists():
                        shutil.rmtree(extract_to, ignore_errors=True)
                except Exception as e:
                    logger.warning("SpatialStructureNormalizer: extract failed: %s", e)
                    if extract_to.exists():
                        shutil.rmtree(extract_to, ignore_errors=True)

    # 2) Move any .h5 from subdirs to root so it is a *sibling* of spatial/
    if _move_h5_to_root(session_dir):
        changed = True

    # 3) Ensure standard matrix name: filtered_feature_bc_matrix.h5 (copy from raw/other .h5 if needed)
    if _ensure_preferred_matrix_h5(session_dir):
        changed = True

    return changed
