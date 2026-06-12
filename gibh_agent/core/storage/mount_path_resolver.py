# -*- coding: utf-8 -*-
"""挂载路径解析：区分容器可写卷、临时缓存与宿主机/Sidecar 路径。"""
from __future__ import annotations

import os
from pathlib import Path
from typing import Any, Dict, Optional

from gibh_agent.core.ingestion_mount_discovery import (
    get_default_ingestion_mount_path,
    probe_mount_path,
    validate_ingestion_mount_path,
)
from gibh_agent.core.workspace_context import get_workspace_context
from gibh_agent.utils.path_resolver import is_windows_abs_path

_EPHEMERAL_PREFIXES = (
    "/app/results",
    "/app/uploads",
    "/results",
    "/uploads",
)


def is_ephemeral_container_path(path: str) -> bool:
    """临时运行时缓存，非用户持久挂载卷。"""
    norm = str(path or "").replace("\\", "/").rstrip("/").lower()
    if not norm:
        return False
    for prefix in _EPHEMERAL_PREFIXES:
        p = prefix.lower()
        if norm == p or norm.startswith(p + "/"):
            return True
    return False


def is_host_side_mount_path(path: str) -> bool:
    """无法在 API 容器内直接 open/write 的路径（Windows 盘符等）。"""
    raw = str(path or "").strip()
    if not raw:
        return False
    if is_windows_abs_path(raw):
        return True
    info = probe_mount_path(raw)
    return info is None


def normalize_workspace_mount_root(path: str) -> str:
    """
    将工作区路径规范为挂载根目录。
    若传入 ``.../rawdata/result``（Sidecar 标准子目录），回退到 ``.../rawdata``。
    """
    s = str(path or "").strip()
    if not s:
        return s
    use_backslash = "\\" in s and not s.replace("\\", "/").startswith("/")
    norm = s.replace("\\", "/").rstrip("/")
    lower = norm.lower()
    if lower.endswith("/result"):
        parent = norm[: -len("/result")]
        if parent:
            if use_backslash and len(parent) >= 3 and parent[1] == ":":
                return parent[0] + ":" + parent[2:].replace("/", "\\")
            return parent
    return s


def resolve_container_writable_mount(
    raw_path: Optional[str],
    *,
    workspace_context: Optional[Dict[str, Any]] = None,
    allow_default: bool = True,
) -> Optional[str]:
    """
    解析可在 API 容器内直接读写的挂载根路径。
    拒绝 /app/results 等临时目录；Windows 盘符路径返回 None（须走 Sidecar）。
    当 raw_path 显式为宿主机路径时，禁止回退到容器 default（避免货不对板）。
    """
    ctx = workspace_context if workspace_context is not None else get_workspace_context()
    explicit = str(raw_path or "").strip()

    if explicit and is_windows_abs_path(explicit):
        return None

    candidates: list[str] = []
    for item in (explicit or None, ctx.get("container_workspace_path"), ctx.get("workspace_path")):
        s = str(item or "").strip()
        if s and s not in candidates:
            candidates.append(s)

    for raw in candidates:
        if is_windows_abs_path(raw):
            continue
        if is_host_side_mount_path(raw) and not is_windows_abs_path(raw):
            continue
        info = probe_mount_path(raw)
        if not info:
            continue
        resolved = str(info["path"])
        if is_ephemeral_container_path(resolved):
            continue
        return resolved

    if allow_default and not explicit:
        default = get_default_ingestion_mount_path()
        if default and not is_ephemeral_container_path(default):
            return default
    return None


def resolve_host_workspace_path(
    raw_path: Optional[str] = None,
    *,
    workspace_context: Optional[Dict[str, Any]] = None,
) -> Optional[str]:
    """返回宿主机/Sidecar 可访问的工作区根路径（用于 Sidecar 落盘）。"""
    ctx = workspace_context if workspace_context is not None else get_workspace_context()
    for item in (
        raw_path,
        ctx.get("host_workspace_path"),
        ctx.get("workspace_path"),
    ):
        s = str(item or "").strip()
        if not s:
            continue
        if is_windows_abs_path(s):
            return normalize_workspace_mount_root(s)
        if is_host_side_mount_path(s) and probe_mount_path(s) is None:
            return normalize_workspace_mount_root(s)
    host = str(ctx.get("host_workspace_path") or "").strip()
    if host:
        return normalize_workspace_mount_root(host)
    return None


def resolve_ingestion_mount_target(
    mount_path_override: Optional[str],
    *,
    workspace_context: Optional[Dict[str, Any]] = None,
    local_workspace_mounted: bool = False,
) -> tuple[Optional[str], Optional[str], Optional[str]]:
    """
    解析入库/落盘目标。

    Returns:
        (mount_type, path, error_message)
        mount_type: ``local_sidecar`` | ``local_volume`` | None
    """
    raw = str(mount_path_override or "").strip()
    ctx = workspace_context if workspace_context is not None else get_workspace_context()
    prefer_sidecar = local_workspace_mounted or bool(ctx.get("host_workspace_path"))

    if raw:
        if is_windows_abs_path(raw) or (
            is_host_side_mount_path(raw) and not validate_ingestion_mount_path(raw)
        ):
            host = normalize_workspace_mount_root(raw)
            if host:
                return "local_sidecar", host, None
            return None, None, f"无法解析本地挂载路径: {raw}"

        validated = validate_ingestion_mount_path(raw)
        if not validated:
            validated = resolve_container_writable_mount(raw, workspace_context=ctx, allow_default=False)
        if validated:
            return "local_volume", validated, None

        if prefer_sidecar:
            host = resolve_host_workspace_path(raw, workspace_context=ctx)
            if host:
                return "local_sidecar", host, None

        return None, None, (
            f"挂载路径不可用（须为容器内可写目录，或已启动 Local Sidecar 的宿主机路径）: {raw}"
        )

    if prefer_sidecar:
        host = resolve_host_workspace_path(workspace_context=ctx)
        if host:
            return "local_sidecar", host, None
        return None, None, "已声明本地项目挂载，但未获取到宿主机工作区路径；请重新挂载本地目录。"

    container = resolve_container_writable_mount(None, workspace_context=ctx, allow_default=True)
    if container:
        return "local_volume", container, None

    host = resolve_host_workspace_path(workspace_context=ctx)
    if host:
        return "local_sidecar", host, None

    return None, None, "尚未配置业务数据库挂载，且未能自动探测到可用路径。"
