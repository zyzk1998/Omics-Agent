"""
SpatialStructureNormalizer: reorganize loose uploads into 10x Visium directory structure.

Runs after file upload so that spatial.tar.gz / spatial.zip and .h5 matrix
uploaded together become a valid Visium directory (spatial/ + matrix) for FileInspector.
"""
import logging
import shutil
import tarfile
import zipfile
from pathlib import Path
from typing import List, Optional

logger = logging.getLogger(__name__)

# Archive extensions and naming patterns for "spatial" pack
SPATIAL_ARCHIVE_SUFFIXES = (".tar.gz", ".tgz", ".tar", ".zip")
SPATIAL_NAME_PATTERNS = ("spatial", "spatial_images", "spatial_metadata")


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


def normalize_session_directory(session_dir: Path) -> bool:
    """
    Reorganize session_dir so it has a spatial/ folder and (optionally) .h5 matrix side-by-side.
    - Detects spatial.tar.gz, spatial.zip, or archives with "spatial" in the name.
    - Extracts them; if result is a single folder named "spatial", keeps it; if loose files, creates spatial/ and moves them in.
    - Does nothing if structure is already correct or no spatial archive is found.

    Args:
        session_dir: Directory that may contain loose files (e.g. matrix.h5 + spatial.tar.gz).

    Returns:
        True if any change was made, False otherwise.
    """
    if not session_dir.is_dir():
        return False

    session_dir = session_dir.resolve()
    changed = False

    # 1) Already correct: spatial/ exists and has content
    spatial_dir = session_dir / "spatial"
    if spatial_dir.is_dir() and any(spatial_dir.iterdir()):
        logger.debug("SpatialStructureNormalizer: spatial/ already present in %s", session_dir)
        return False

    # 2) Find spatial archive(s); process first match
    candidates: List[Path] = []
    for f in session_dir.iterdir():
        if f.is_file() and _is_spatial_archive(f):
            candidates.append(f)

    if not candidates:
        logger.debug("SpatialStructureNormalizer: no spatial archive in %s", session_dir)
        return False

    archive_path = candidates[0]
    if len(candidates) > 1:
        logger.info("SpatialStructureNormalizer: multiple spatial archives, using %s", archive_path.name)

    top_level = _list_archive_top_level(archive_path)
    if not top_level:
        return False

    # 3) Extract to a temp dir under session_dir to avoid collisions
    extract_to = session_dir / f"_spatial_extract_{archive_path.stem}"
    extract_to.mkdir(parents=True, exist_ok=True)
    try:
        if archive_path.suffix.lower() == ".zip":
            ok = _extract_zip(archive_path, extract_to)
        else:
            ok = _extract_tar(archive_path, extract_to)
        if not ok:
            return False
    except Exception as e:
        logger.warning("SpatialStructureNormalizer: extract failed: %s", e)
        if extract_to.exists():
            shutil.rmtree(extract_to, ignore_errors=True)
        return False

    # 4) Smart handling: one folder named "spatial" vs loose files
    subdirs = [p for p in extract_to.iterdir() if p.is_dir()]
    files = [p for p in extract_to.iterdir() if p.is_file()]

    if len(subdirs) == 1 and subdirs[0].name.lower() == "spatial":
        # Single top-level "spatial" folder: move it to session_dir/spatial
        src_spatial = subdirs[0]
        if spatial_dir.exists():
            shutil.rmtree(spatial_dir, ignore_errors=True)
        try:
            shutil.move(str(src_spatial), str(spatial_dir))
            changed = True
        except Exception as e:
            logger.warning("SpatialStructureNormalizer: move spatial folder failed: %s", e)
        logger.info("SpatialStructureNormalizer: moved extracted spatial/ into %s", session_dir)
    elif subdirs or files:
        # Loose folders or files: create spatial/ and move all into it
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
                logger.warning("SpatialStructureNormalizer: move %s into spatial/ failed: %s", p, e)
        if changed:
            logger.info("SpatialStructureNormalizer: created spatial/ and moved loose content in %s", session_dir)
    else:
        logger.debug("SpatialStructureNormalizer: extract empty or no top-level items")

    # Cleanup temp extract dir
    if extract_to.exists():
        shutil.rmtree(extract_to, ignore_errors=True)

    return changed
