"""
Atomic Epigenomics tools (pybedtools).
Peak intersection and merging. Requires pybedtools>=0.9.0 and bedtools installed (apt-get install bedtools).
"""
import logging
from pathlib import Path
from typing import Dict, Any, Optional

from ..core.tool_registry import registry

logger = logging.getLogger(__name__)


@registry.register(
    name="intersect_peaks",
    description="Find overlapping regions between two BED files. Returns a BED file path of intersections.",
    category="Epigenomics",
    output_type="file_path",
)
def intersect_peaks(
    bed_file_a: str,
    bed_file_b: str,
    output_path: Optional[str] = None,
) -> Dict[str, Any]:
    """
    Return overlapping regions between two BED files.

    Args:
        bed_file_a: Path to first BED file.
        bed_file_b: Path to second BED file.
        output_path: Optional path to write result BED; if not set, a temp file is used.

    Returns:
        Dict with status, output_path (BED of intersections), n_overlaps, and optional error.
    """
    try:
        import pybedtools
    except ImportError as e:
        logger.warning("pybedtools not installed: %s", e)
        return {"status": "error", "error": "pybedtools is required. Install with: pip install pybedtools>=0.9.0. Also install bedtools: apt-get install bedtools"}
    a, b = Path(bed_file_a), Path(bed_file_b)
    if not a.is_file():
        return {"status": "error", "error": f"BED file not found: {bed_file_a}"}
    if not b.is_file():
        return {"status": "error", "error": f"BED file not found: {bed_file_b}"}
    try:
        A = pybedtools.BedTool(str(a))
        B = pybedtools.BedTool(str(b))
        result = A.intersect(B)
        if output_path:
            out = Path(output_path)
            out.parent.mkdir(parents=True, exist_ok=True)
            result.saveas(str(out))
            out_str = str(out.resolve())
        else:
            out_str = result.fn
        n = sum(1 for _ in open(out_str)) if Path(out_str).exists() else 0
        return {
            "status": "success",
            "output_path": out_str,
            "n_overlaps": n,
            "bed_file_a": str(a.resolve()),
            "bed_file_b": str(b.resolve()),
        }
    except Exception as e:
        logger.exception("intersect_peaks failed: %s", e)
        return {"status": "error", "error": str(e)}


@registry.register(
    name="merge_peaks",
    description="Merge nearby or overlapping peaks in a BED file. Returns a BED file path of merged regions.",
    category="Epigenomics",
    output_type="file_path",
)
def merge_peaks(
    bed_file: str,
    output_path: Optional[str] = None,
    distance: int = 0,
) -> Dict[str, Any]:
    """
    Merge overlapping or nearby peaks in a BED file.

    Args:
        bed_file: Path to BED file.
        output_path: Optional path to write merged BED; if not set, a temp file is used.
        distance: Merge peaks within this distance (default 0 = only overlapping).

    Returns:
        Dict with status, output_path (merged BED), n_merged, and optional error.
    """
    try:
        import pybedtools
    except ImportError as e:
        logger.warning("pybedtools not installed: %s", e)
        return {"status": "error", "error": "pybedtools is required. Install with: pip install pybedtools>=0.9.0. Also install bedtools: apt-get install bedtools"}
    p = Path(bed_file)
    if not p.is_file():
        return {"status": "error", "error": f"BED file not found: {bed_file}"}
    try:
        bt = pybedtools.BedTool(str(p))
        merged = bt.merge(d=distance)
        if output_path:
            out = Path(output_path)
            out.parent.mkdir(parents=True, exist_ok=True)
            merged.saveas(str(out))
            out_str = str(out.resolve())
        else:
            out_str = merged.fn
        n = sum(1 for _ in open(out_str)) if Path(out_str).exists() else 0
        return {
            "status": "success",
            "output_path": out_str,
            "n_merged": n,
            "bed_file": str(p.resolve()),
        }
    except Exception as e:
        logger.exception("merge_peaks failed: %s", e)
        return {"status": "error", "error": str(e)}
