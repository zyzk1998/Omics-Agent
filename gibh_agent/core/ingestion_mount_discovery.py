# -*- coding: utf-8 -*-
"""Docker 容器内可写挂载目录自动探测（一键入库 Smart Mount Discovery）。"""
from __future__ import annotations

import os
from pathlib import Path
from typing import Any, Dict, List, Optional


def _candidate_mount_paths() -> List[str]:
    """按优先级返回候选容器内路径（去重保序）。"""
    ordered: List[str] = []
    env_default = str(os.getenv("DEFAULT_INGESTION_DIR", "") or "").strip()
    if env_default:
        ordered.append(env_default)
    ordered.extend(
        [
            "/app/data/omics_output",
            "/app/data",
            "/app/workspace",
            "/data/omics_output",
            str(Path(os.getenv("RESULTS_DIR", "/app/results")).expanduser()),
        ]
    )
    seen: set[str] = set()
    out: List[str] = []
    for raw in ordered:
        p = str(Path(raw).expanduser())
        if not p or p in seen:
            continue
        seen.add(p)
        out.append(p)
    return out


def probe_mount_path(path: str) -> Optional[Dict[str, Any]]:
    """探测单一路径是否可作为入库落盘目录。"""
    raw = str(path or "").strip()
    if not raw:
        return None
    p = Path(raw).expanduser()
    try:
        resolved = p.resolve()
    except OSError:
        return None
    if not resolved.is_dir():
        return None
    if not os.access(resolved, os.W_OK):
        return None
    return {
        "path": str(resolved),
        "writable": True,
        "source": "env" if raw == os.getenv("DEFAULT_INGESTION_DIR", "").strip() else "probe",
    }


def discover_ingestion_mount_paths() -> List[Dict[str, Any]]:
    """返回所有可用挂载点（存在且可写）。"""
    found: List[Dict[str, Any]] = []
    for candidate in _candidate_mount_paths():
        info = probe_mount_path(candidate)
        if info:
            found.append(info)
    return found


def get_default_ingestion_mount_path() -> Optional[str]:
    """返回系统推荐默认入库路径（第一个可用挂载点）。"""
    paths = discover_ingestion_mount_paths()
    if not paths:
        return None
    return str(paths[0]["path"])


def validate_ingestion_mount_path(path: str) -> Optional[str]:
    """
    校验前端/localStorage 传入的路径是否为容器内可写目录。
    返回规范化绝对路径；无效则 None。
    """
    info = probe_mount_path(path)
    if not info:
        return None
    return str(info["path"])
