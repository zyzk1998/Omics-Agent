"""
Universal Ingestion Pipeline — Phase 1: Universal Unpacker

Handles .zip, .tar, .tar.gz, .tgz (and optionally .7z) for ANY modality.
Extracts archives, flattens single-folder extractions, and optionally recurses.
Does not interpret modality; that is left to the Modality Sniffer.
"""
import logging
import shutil
import tarfile
import zipfile
from pathlib import Path
from typing import List, Optional, Tuple

logger = logging.getLogger(__name__)

EXTRACTED_DIR_NAME = "extracted_content"
ARCHIVE_EXTENSIONS = (".zip", ".tar", ".tgz", ".tar.gz")
# .gz alone is ambiguous (could be gzipped nifti); only .tar.gz / .tgz are treated as archives


def _is_archive(path: Path) -> bool:
    """True if path is a known archive type."""
    if not path.is_file():
        return False
    name = path.name.lower()
    return any(name.endswith(ext) for ext in ARCHIVE_EXTENSIONS)


def _list_archive_top_level(archive_path: Path) -> Optional[List[str]]:
    """List top-level entry names in archive (files and dirs)."""
    name_lower = archive_path.name.lower()
    try:
        if name_lower.endswith(".zip"):
            with zipfile.ZipFile(archive_path, "r") as z:
                names = set()
                for i in z.namelist():
                    top = i.split("/")[0].split("\\")[0].rstrip("/")
                    if top and top != ".":
                        names.add(top)
                return list(names)
        if name_lower.endswith((".tar.gz", ".tgz", ".tar")):
            with tarfile.open(archive_path, "r:*") as t:
                names = [m.name.split("/")[0].split("\\")[0].rstrip("/") for m in t.getmembers() if m.name]
                return list(dict.fromkeys(names))
    except Exception as e:
        logger.warning("UniversalNormalizer: list archive failed %s: %s", archive_path.name, e)
    return None


def _extract_zip(archive_path: Path, dest_dir: Path) -> bool:
    try:
        with zipfile.ZipFile(archive_path, "r") as z:
            z.extractall(dest_dir)
        return True
    except Exception as e:
        logger.warning("UniversalNormalizer: zip extract failed %s: %s", archive_path.name, e)
        return False


def _extract_tar(archive_path: Path, dest_dir: Path) -> bool:
    try:
        with tarfile.open(archive_path, "r:*") as t:
            t.extractall(dest_dir)
        return True
    except Exception as e:
        logger.warning("UniversalNormalizer: tar extract failed %s: %s", archive_path.name, e)
        return False


def _extract_archive(archive_path: Path, dest_dir: Path) -> bool:
    """Extract archive to dest_dir. Returns True on success."""
    dest_dir.mkdir(parents=True, exist_ok=True)
    name_lower = archive_path.name.lower()
    if name_lower.endswith(".zip"):
        return _extract_zip(archive_path, dest_dir)
    if name_lower.endswith((".tar.gz", ".tgz", ".tar")):
        return _extract_tar(archive_path, dest_dir)
    return False


def _flatten_single_folder(target: Path) -> bool:
    """
    If target contains exactly one directory and no other files, move that directory's
    contents up to target and remove the wrapper. Returns True if flattening was done.
    """
    if not target.is_dir():
        return False
    items = [p for p in target.iterdir()]
    dirs = [p for p in items if p.is_dir()]
    files = [p for p in items if p.is_file()]
    if len(dirs) != 1 or len(files) != 0:
        return False
    single = dirs[0]
    moved = False
    for p in single.iterdir():
        dest = target / p.name
        if dest.exists():
            dest = target / f"{p.stem}_moved{p.suffix}" if p.is_file() else target / f"{p.name}_moved"
        try:
            shutil.move(str(p), str(dest))
            moved = True
        except Exception as e:
            logger.warning("UniversalNormalizer: flatten move failed %s: %s", p.name, e)
    try:
        if single.exists():
            shutil.rmtree(single, ignore_errors=True)
    except Exception as e:
        logger.debug("UniversalNormalizer: rmtree wrapper failed: %s", e)
    return moved


def _unpack_one(archive_path: Path, session_dir: Path, remove_archive: bool) -> Optional[Path]:
    """
    Extract one archive into session_dir/EXTRACTED_DIR_NAME/<stem>/.
    Apply flattening. Optionally remove archive.
    Returns the effective root path (where data lives after extract), or None.
    """
    stem = archive_path.stem
    if archive_path.name.lower().endswith(".tar.gz"):
        stem = archive_path.name[:-7]  # avoid .tar
    extract_to = session_dir / EXTRACTED_DIR_NAME / stem
    if extract_to.exists():
        shutil.rmtree(extract_to, ignore_errors=True)
    if not _extract_archive(archive_path, extract_to):
        return None
    # Flatten: if extract_to has only one subdir, move contents up
    _flatten_single_folder(extract_to)
    if remove_archive:
        try:
            archive_path.unlink()
        except OSError as e:
            logger.warning("UniversalNormalizer: remove archive failed: %s", e)
    return extract_to


def normalize_session_directory(
    session_dir: Path,
    remove_archive: bool = False,
    max_nested: int = 1,
) -> Tuple[bool, Optional[Path]]:
    """
    Universal unpack: find archives in session_dir, extract, flatten, and optionally recurse.

    Logic:
    1. Iterate session_dir. If an archive is found, extract to extracted_content/<stem>/.
    2. If the extracted content is a single folder, move its contents up (flatten).
    3. If max_nested > 0, scan extracted dirs for nested archives and unpack one level.
    4. Optionally remove the original archive.

    Returns:
        (changed, effective_root): changed=True if any extraction was done;
        effective_root = path to the directory that contains the actual data (session_dir
        or session_dir/extracted_content/<stem>), or None if nothing was unpacked.
    """
    if not session_dir.is_dir():
        return (False, None)
    session_dir = session_dir.resolve()
    archives = [p for p in session_dir.iterdir() if p.is_file() and _is_archive(p)]
    if not archives:
        return (False, None)

    # Use first archive (e.g. single zip upload)
    archive_path = archives[0]
    effective_root = _unpack_one(archive_path, session_dir, remove_archive)
    if not effective_root:
        return (False, None)

    # Optional: handle nested archives inside effective_root (one level)
    if max_nested > 0:
        nested = [p for p in effective_root.iterdir() if p.is_file() and _is_archive(p)]
        if nested:
            sub_root = _unpack_one(nested[0], effective_root, remove_archive=False)
            if sub_root:
                effective_root = sub_root

    logger.info("UniversalNormalizer: unpacked %s -> %s", archive_path.name, effective_root)
    return (True, effective_root)
