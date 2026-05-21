"""
文件资产来源规范化：Composer 草稿 / 上传响应 / 工作流复用统一 source 词汇表。
"""
from __future__ import annotations

from typing import Any, Dict, Optional

SOURCE_HPC = "hpc_asset"
SOURCE_DATA = "data_asset"


def normalize_file_source(raw: Any, *, path: str = "") -> str:
    """
    将前端/侧栏遗留的 hpc|local|files 与显式 hpc_asset|data_asset 归一为法定值。
    规则：超算内部读取 → hpc_asset；其余上传/数据资产选取 → data_asset。
    """
    if raw is not None:
        s = str(raw).strip().lower()
        if s in (SOURCE_HPC, "hpc", "hpc-open", "compute", "supercompute"):
            return SOURCE_HPC
        if s in (SOURCE_DATA, "data_asset", "data", "local", "files", "upload", "asset"):
            return SOURCE_DATA
    p = (path or "").strip().lower()
    if "/hpc/" in p or p.startswith("hpc://") or "hpc-open" in p:
        return SOURCE_HPC
    return SOURCE_DATA


def enrich_file_record(
    *,
    path: str,
    name: str = "",
    source: Any = None,
    size: Optional[int] = None,
    word_count: Optional[int] = None,
) -> Dict[str, Any]:
    p = (path or "").strip()
    n = (name or "").strip() or (p.split("/")[-1] if p else "")
    src = normalize_file_source(source, path=p)
    row: Dict[str, Any] = {
        "name": n,
        "file_name": n,
        "path": p,
        "file_path": p,
        "source": src,
        "type": _infer_type_from_name(n or p),
    }
    if size is not None:
        try:
            row["size"] = int(size)
        except (TypeError, ValueError):
            pass
    if word_count is not None:
        try:
            row["word_count"] = int(word_count)
        except (TypeError, ValueError):
            pass
    return row


def _infer_type_from_name(name: str) -> str:
    n = (name or "").lower()
    if n.endswith((".fastq", ".fq")) or n.endswith((".fastq.gz", ".fq.gz")):
        return "fastq"
    if ".vcf" in n:
        return "vcf"
    if n.endswith(".bam"):
        return "bam"
    if n.endswith(".csv"):
        return "csv"
    return "file"


def normalize_attachment_list(attachments: Any) -> list:
    if not isinstance(attachments, list):
        return []
    out = []
    for item in attachments:
        if not isinstance(item, dict):
            continue
        p = str(item.get("path") or item.get("file_path") or "").strip()
        n = str(item.get("name") or item.get("file_name") or "").strip()
        if not p and not n:
            continue
        out.append(
            enrich_file_record(
                path=p,
                name=n,
                source=item.get("source"),
                size=item.get("size"),
                word_count=item.get("word_count"),
            )
        )
    return out
