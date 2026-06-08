# -*- coding: utf-8 -*-
"""Label Studio 浏览器可达基址解析（同源反代 /label-studio/ 优先）。"""
from __future__ import annotations

import os


def is_same_origin_proxy_enabled() -> bool:
    raw = os.getenv("LABEL_STUDIO_SAME_ORIGIN_PROXY", "1").strip().lower()
    return raw not in ("0", "false", "no", "off")


def _is_docker_internal_web_url(url: str) -> bool:
    u = url.lower()
    return (
        "host.docker.internal" in u
        or "api-server" in u
        or u.startswith("http://nginx")
        or u.startswith("https://nginx")
    )


def resolve_ls_public_base_url() -> str:
    """
    供浏览器 iframe 打开的 LS 基址。
    默认同源反代：相对路径 /label-studio（与 Omics Agent 同 host:port）。
    显式 LABEL_STUDIO_PUBLIC_URL 非空时优先。
    """
    explicit = (os.getenv("LABEL_STUDIO_PUBLIC_URL") or "").strip().rstrip("/")
    if explicit:
        return explicit
    if is_same_origin_proxy_enabled():
        return "/label-studio"
    return (os.getenv("LABEL_STUDIO_URL") or "http://127.0.0.1:8082").rstrip("/")


def resolve_ls_browser_project_url(project_id: int) -> str:
    """供浏览器 iframe / 新窗口打开的 LS 项目 URL。"""
    base = resolve_ls_public_base_url()
    suffix = f"/projects/{int(project_id)}/data"
    if base.startswith("/"):
        return f"{base.rstrip('/')}{suffix}"
    return f"{base}{suffix}"
