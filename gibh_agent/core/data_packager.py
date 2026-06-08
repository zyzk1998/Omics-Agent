# -*- coding: utf-8 -*-
"""
三轨合一入库前 · 最终版数据归档打包（Artifacts Archive）。

将专家报告、HITL 标注导出、工作流步骤摘要、精选输出文件打成标准 tar.gz，
避免把零散中间 log 全量上传。
"""
from __future__ import annotations

import json
import logging
import os
import tarfile
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, Optional

logger = logging.getLogger(__name__)

_ARCHIVE_ROOT = Path(
    os.getenv("GIBH_ARTIFACTS_ARCHIVE_DIR", "data/artifacts_archive")
).expanduser()

# 归档时允许纳入的扩展名（轻量报告 + 结构化标注 + 关键图表）
_ARCHIVE_EXTENSIONS = frozenset(
    {
        ".md", ".json", ".csv", ".tsv", ".png", ".jpg", ".jpeg", ".svg", ".pdf",
        ".html", ".yaml", ".yml", ".txt",
    }
)
_MAX_FILE_BYTES = int(os.getenv("GIBH_ARCHIVE_MAX_FILE_BYTES", str(50 * 1024 * 1024)))


def _safe_read_text(path: Path, limit: int = 512_000) -> str:
    try:
        raw = path.read_text(encoding="utf-8", errors="replace")
        return raw[:limit]
    except OSError:
        return ""


def _collect_output_files(output_dir: Optional[str], max_files: int = 40) -> List[Dict[str, str]]:
    if not output_dir:
        return []
    root = Path(output_dir).expanduser()
    if not root.is_dir():
        return []
    picked: List[Dict[str, str]] = []
    for p in sorted(root.rglob("*")):
        if not p.is_file():
            continue
        if p.suffix.lower() not in _ARCHIVE_EXTENSIONS and not p.name.endswith(".nii.gz"):
            continue
        try:
            if p.stat().st_size > _MAX_FILE_BYTES:
                continue
        except OSError:
            continue
        picked.append({"relative_path": str(p.relative_to(root)), "absolute_path": str(p.resolve())})
        if len(picked) >= max_files:
            break
    return picked


def build_artifacts_archive(
    *,
    session_id: str,
    owner_id: str,
    expert_report_markdown: str = "",
    steps_details: Optional[List[Dict[str, Any]]] = None,
    hitl_annotations: Any = None,
    hitl_meta: Optional[Dict[str, Any]] = None,
    output_dir: Optional[str] = None,
    skip_hitl: bool = False,
) -> Dict[str, Any]:
    """
    构建标准归档包并返回 manifest + archive_path。

    目录结构（tar 内）:
        manifest.json
        expert_report_final.md
        hitl/annotations.json   (若有)
        workflow/steps_summary.json
        outputs/...              (精选文件)
    """
    sid = str(session_id or "unknown").strip()
    oid = str(owner_id or "anonymous").strip()
    ts = datetime.utcnow().strftime("%Y%m%d_%H%M%S")
    work = _ARCHIVE_ROOT / oid / sid / f"bundle_{ts}"
    work.mkdir(parents=True, exist_ok=True)

    manifest: Dict[str, Any] = {
        "session_id": sid,
        "owner_id": oid,
        "created_at": datetime.utcnow().isoformat() + "Z",
        "skip_hitl": bool(skip_hitl),
        "hitl_meta": hitl_meta or {},
        "output_dir": output_dir,
        "files": [],
    }

    report_md = (expert_report_markdown or "").strip()
    if not report_md:
        report_md = "# 专家分析报告\n\n（归档时尚未生成最终版报告正文）\n"
    report_path = work / "expert_report_final.md"
    report_path.write_text(report_md, encoding="utf-8")
    manifest["files"].append("expert_report_final.md")

    steps_path = work / "workflow" / "steps_summary.json"
    steps_path.parent.mkdir(parents=True, exist_ok=True)
    steps_summary = steps_details or []
    steps_path.write_text(
        json.dumps(steps_summary, ensure_ascii=False, indent=2, default=str),
        encoding="utf-8",
    )
    manifest["files"].append("workflow/steps_summary.json")

    if hitl_annotations is not None and not skip_hitl:
        hitl_dir = work / "hitl"
        hitl_dir.mkdir(parents=True, exist_ok=True)
        ann_path = hitl_dir / "annotations.json"
        ann_path.write_text(
            json.dumps(hitl_annotations, ensure_ascii=False, indent=2, default=str),
            encoding="utf-8",
        )
        manifest["files"].append("hitl/annotations.json")

    out_pick = _collect_output_files(output_dir)
    for item in out_pick:
        rel = item["relative_path"]
        src = Path(item["absolute_path"])
        dest = work / "outputs" / rel
        dest.parent.mkdir(parents=True, exist_ok=True)
        try:
            dest.write_bytes(src.read_bytes())
            manifest["files"].append(f"outputs/{rel}")
        except OSError as exc:
            logger.warning("[DataPackager] skip file %s: %s", src, exc)

    manifest_path = work / "manifest.json"
    manifest_path.write_text(json.dumps(manifest, ensure_ascii=False, indent=2), encoding="utf-8")

    archive_path = work.parent / f"{work.name}.tar.gz"
    with tarfile.open(archive_path, "w:gz") as tar:
        tar.add(work, arcname=work.name)

    logger.info(
        "[DataPackager] archive ready session=%s path=%s files=%s",
        sid,
        archive_path,
        len(manifest.get("files") or []),
    )
    return {
        "status": "success",
        "archive_path": str(archive_path.resolve()),
        "manifest": manifest,
        "bundle_dir": str(work.resolve()),
    }
