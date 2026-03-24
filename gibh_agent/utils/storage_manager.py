# -*- coding: utf-8 -*-
"""
存储管理：扫描大文件、按策略自动清理重型产物，保留轻量报告（图表/表格/Markdown/JSON）。
"""
from __future__ import annotations

import logging
import os
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Literal, Optional, Sequence, Tuple

logger = logging.getLogger(__name__)

SortKey = Literal["mtime_desc", "size_desc", "mtime_asc"]

# 永不删除：历史快照、报告、图表、配置等
_LIGHT_EXTENSIONS = frozenset(
    {
        ".csv",
        ".tsv",
        ".md",
        ".png",
        ".json",
        ".svg",
        ".jpg",
        ".jpeg",
        ".gif",
        ".webp",
        ".html",
        ".htm",
        ".txt",
        ".pdf",
        ".yaml",
        ".yml",
        ".log",
    }
)

# 明确视为「重型」的数据扩展（按最长后缀优先匹配）
_HEAVY_ENDSWITH: Tuple[str, ...] = (
    ".fastq.gz",
    ".fq.gz",
    ".tar.gz",
    ".nii.gz",
    ".mtx.gz",
    ".h5ad",
    ".loom",
    ".bam",
    ".bai",
    ".cram",
    ".mtx",
    ".h5",
    ".fastq",
    ".fq",
    ".zip",
    ".tar",
    ".zarr",
    ".nii",
    ".dcm",
    ".rds",
    ".parquet",
    ".feather",
    ".hdf5",
)

# 未列入「轻量」且体积极大时，视为可清理冷数据（避免未知扩展撑盘）
_LARGE_ORPHAN_BYTES = 200 * 1024 * 1024


def _is_protected_path(path: Path) -> bool:
    """轻量报告与文本：一律保留，不参与删除。"""
    name = path.name.lower()
    if name.endswith(".nii.gz"):
        return False
    suf = path.suffix.lower()
    return suf in _LIGHT_EXTENSIONS


def _is_heavy_file(path: Path, stat_size: int) -> bool:
    """
    可清理的重型候选：匹配重型后缀，或「非保护 + 超大 orphan」。
    核心原则：protected（.csv/.md/.png/.json 等）永远不删。
    """
    if _is_protected_path(path):
        return False
    nl = path.name.lower()
    for suf in _HEAVY_ENDSWITH:
        if nl.endswith(suf):
            return True
    if stat_size >= _LARGE_ORPHAN_BYTES:
        return True
    return False


@dataclass
class CleanupResult:
    deleted_count: int = 0
    freed_bytes: int = 0
    errors: List[str] = field(default_factory=list)
    deleted_paths_sample: List[str] = field(default_factory=list)

    def to_dict(self, sample_limit: int = 50) -> Dict[str, Any]:
        freed_mb = self.freed_bytes / (1024 * 1024)
        freed_gb = self.freed_bytes / (1024 * 1024 * 1024)
        return {
            "deleted_count": self.deleted_count,
            "freed_bytes": self.freed_bytes,
            "freed_mb": round(freed_mb, 4),
            "freed_gb": round(freed_gb, 6),
            "errors": self.errors,
            "deleted_paths_sample": self.deleted_paths_sample[:sample_limit],
        }


class StorageManager:
    """扫描 RESULTS / UPLOAD，并按策略清理过期重型文件。"""

    def __init__(
        self,
        results_dir: Optional[str] = None,
        upload_dir: Optional[str] = None,
    ) -> None:
        self.results_dir = Path(results_dir or os.getenv("RESULTS_DIR", "/app/results")).resolve()
        self.upload_dir = Path(upload_dir or os.getenv("UPLOAD_DIR", "/app/uploads")).resolve()

    def _roots(self) -> List[Path]:
        roots = []
        for p in (self.results_dir, self.upload_dir):
            if p.is_dir():
                roots.append(p)
        return roots

    def scan_large_files(
        self,
        sort_by: SortKey = "mtime_desc",
        max_entries: Optional[int] = None,
    ) -> List[Dict[str, Any]]:
        """
        扫描 RESULTS_DIR 与 UPLOAD_DIR 下所有普通文件，返回大小与 mtime。
        """
        rows: List[Dict[str, Any]] = []
        for root in self._roots():
            for dirpath, _, filenames in os.walk(root, followlinks=False):
                base = Path(dirpath)
                for fn in filenames:
                    fp = base / fn
                    try:
                        if not fp.is_file():
                            continue
                        st = fp.stat()
                    except OSError:
                        continue
                    rows.append(
                        {
                            "path": str(fp),
                            "size_bytes": st.st_size,
                            "mtime": st.st_mtime,
                            "protected": _is_protected_path(fp),
                            "heavy": _is_heavy_file(fp, st.st_size),
                        }
                    )

        if sort_by == "size_desc":
            rows.sort(key=lambda r: (-r["size_bytes"], -r["mtime"]))
        elif sort_by == "mtime_asc":
            rows.sort(key=lambda r: (r["mtime"], -r["size_bytes"]))
        else:
            rows.sort(key=lambda r: (-r["mtime"], -r["size_bytes"]))

        if max_entries is not None:
            rows = rows[: max(0, max_entries)]
        return rows

    def _total_bytes(self) -> int:
        total = 0
        for root in self._roots():
            for dirpath, _, filenames in os.walk(root, followlinks=False):
                base = Path(dirpath)
                for fn in filenames:
                    fp = base / fn
                    try:
                        if fp.is_file():
                            total += fp.stat().st_size
                    except OSError:
                        pass
        return total

    def _iter_heavy_files(self) -> List[Tuple[Path, int, float]]:
        """所有「重型且非保护」文件 (path, size, mtime)。"""
        found: List[Tuple[Path, int, float]] = []
        for root in self._roots():
            for dirpath, _, filenames in os.walk(root, followlinks=False):
                base = Path(dirpath)
                for fn in filenames:
                    fp = base / fn
                    try:
                        if not fp.is_file():
                            continue
                        st = fp.stat()
                    except OSError:
                        continue
                    if _is_protected_path(fp):
                        continue
                    if not _is_heavy_file(fp, st.st_size):
                        continue
                    found.append((fp, st.st_size, st.st_mtime))
        return found

    def _delete_one(
        self,
        path: Path,
        size: int,
        result: CleanupResult,
        dry_run: bool,
    ) -> bool:
        try:
            if dry_run:
                result.deleted_count += 1
                result.freed_bytes += size
                if len(result.deleted_paths_sample) < 50:
                    result.deleted_paths_sample.append(str(path))
                return True
            path.unlink()
            result.deleted_count += 1
            result.freed_bytes += size
            if len(result.deleted_paths_sample) < 50:
                result.deleted_paths_sample.append(str(path))
            return True
        except OSError as e:
            result.errors.append(f"{path}: {e}")
            return False

    def auto_cleanup(
        self,
        max_days: int = 7,
        max_size_gb: float = 100.0,
        *,
        deep: bool = False,
        dry_run: bool = False,
    ) -> Dict[str, Any]:
        """
        自动清理：
        - 仅删除「重型」且「修改时间早于 max_days」的文件；轻量报告扩展名永不删除。
        - 若两目录合计仍超过 max_size_gb，按「最旧重型优先」继续删除直至低于阈值（仍不删保护类型）。

        deep=True 时有效保留期缩短为 min(max_days, 3)。
        """
        import time

        eff_days = max(1, min(max_days, 3) if deep else max(1, max_days))
        cutoff_ts = time.time() - eff_days * 86400
        limit_bytes = max(1.0, float(max_size_gb)) * (1024**3)

        result = CleanupResult()
        heavy = self._iter_heavy_files()
        # 阶段 1：过期重型文件（mtime < cutoff）
        age_candidates = [(p, s, m) for p, s, m in heavy if m < cutoff_ts]
        age_candidates.sort(key=lambda x: (x[2], -x[1]))
        for p, s, _m in age_candidates:
            self._delete_one(p, s, result, dry_run)

        # 阶段 2：总占用仍超上限 → 删最旧重型（即使未过 max_days）
        if dry_run:
            total = self._total_bytes() - result.freed_bytes
        else:
            total = self._total_bytes()
        safety = 0
        max_quota_passes = 500000
        while total > limit_bytes and safety < max_quota_passes:
            safety += 1
            heavy2 = self._iter_heavy_files()
            if not heavy2:
                break
            heavy2.sort(key=lambda x: (x[2], -x[1]))
            p, s, _m = heavy2[0]
            if not p.exists():
                continue
            self._delete_one(p, s, result, dry_run)
            total -= s

        logger.info(
            "[StorageManager] cleanup dry_run=%s deep=%s deleted=%s freed_mb=%.2f errors=%s",
            dry_run,
            deep,
            result.deleted_count,
            result.freed_bytes / (1024 * 1024),
            len(result.errors),
        )
        out = result.to_dict()
        out["max_days_effective"] = eff_days
        out["dry_run"] = dry_run
        out["deep"] = deep
        out["projected_total_bytes_after"] = int(max(0, total))
        return out
