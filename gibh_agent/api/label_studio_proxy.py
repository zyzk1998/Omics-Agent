# -*- coding: utf-8 -*-
"""Label Studio 同源反代：/label-studio/* → LABEL_STUDIO_URL（容器内 label-studio:8080）。"""
from __future__ import annotations

import os
from typing import Iterable
from urllib.parse import parse_qs, urlencode, urlsplit, urlunsplit

import httpx
from fastapi import APIRouter, Request, Response

_PROXY_PUBLIC_PREFIX = "/label-studio"
_LS_REDIRECT_PREFIXES = (
    "/user/",
    "/projects/",
    "/static/",
    "/react-app/",
    "/label-studio-frontend/",
    "/data/",
    "/tasks/",
    "/api/",
)

router = APIRouter(include_in_schema=False)

_HOP_BY_HOP = frozenset(
    {
        "connection",
        "keep-alive",
        "proxy-authenticate",
        "proxy-authorization",
        "te",
        "trailers",
        "transfer-encoding",
        "upgrade",
        "host",
        "content-length",
    }
)

_UPSTREAM_BASE = (os.getenv("LABEL_STUDIO_URL") or "http://127.0.0.1:8082").rstrip("/")
_PROXY_TIMEOUT = float(os.getenv("LABEL_STUDIO_PROXY_TIMEOUT_SEC", "600"))
# 上游 LS 监听根路径 /；浏览器侧 /label-studio/ 由反代剥前缀。设非空则保留子路径（特殊部署）。
_UPSTREAM_SUBPATH = (os.getenv("LABEL_STUDIO_UPSTREAM_SUBPATH") or "").strip().rstrip("/")


def _upstream_url(path: str, query: str) -> str:
    sub = _UPSTREAM_SUBPATH
    if sub:
        target = f"{_UPSTREAM_BASE}{sub}/{path}" if path else f"{_UPSTREAM_BASE}{sub}/"
    else:
        target = f"{_UPSTREAM_BASE}/{path}" if path else f"{_UPSTREAM_BASE}/"
    if query:
        target = f"{target}?{query}"
    return target


def _forwarding_headers(request: Request) -> dict[str, str]:
    out: dict[str, str] = {}
    for key, value in request.headers.items():
        if key.lower() in _HOP_BY_HOP:
            continue
        out[key] = value
    # 让 LS Django 识别外部访问 host / 端口（CSRF_TRUSTED_ORIGINS 对齐）
    out["X-Forwarded-Host"] = request.headers.get("x-forwarded-host") or request.headers.get("host", "")
    out["X-Forwarded-Proto"] = request.headers.get("x-forwarded-proto") or request.url.scheme
    out["X-Forwarded-Port"] = str(request.url.port or (443 if request.url.scheme == "https" else 80))
    return out


def _prefix_ls_public_path(path: str) -> str:
    raw = str(path or "").strip()
    if not raw.startswith("/"):
        return raw
    if raw.startswith(_PROXY_PUBLIC_PREFIX + "/") or raw == _PROXY_PUBLIC_PREFIX:
        return raw
    if any(raw.startswith(p) for p in _LS_REDIRECT_PREFIXES):
        return f"{_PROXY_PUBLIC_PREFIX}{raw}"
    return raw


def _rewrite_next_query(query: str) -> str:
    if not query:
        return query
    qs = parse_qs(query, keep_blank_values=True)
    if "next" not in qs:
        return query
    rewritten: list[str] = []
    for val in qs["next"]:
        if isinstance(val, str) and val.startswith("/") and not val.startswith(_PROXY_PUBLIC_PREFIX):
            rewritten.append(_prefix_ls_public_path(val))
        else:
            rewritten.append(val)
    qs["next"] = rewritten
    return urlencode(qs, doseq=True)


def rewrite_ls_proxy_location(location: str, request: Request) -> str:
    """
    将 LS 上游相对重定向（如 /user/login/?next=/projects/1/data/）
    改写为同源 /label-studio/...，避免 iframe 落到 FastAPI 404。
    """
    loc = str(location or "").strip()
    if not loc:
        return loc

    parts = urlsplit(loc)
    if parts.scheme and parts.netloc:
        req_host = (request.headers.get("host") or "").split(":")[0]
        loc_host = parts.hostname or ""
        if loc_host and req_host and loc_host != req_host:
            return loc
        path = _prefix_ls_public_path(parts.path or "/")
        query = _rewrite_next_query(parts.query)
        return urlunsplit((parts.scheme, parts.netloc, path, query, parts.fragment))

    if loc.startswith("/"):
        path = _prefix_ls_public_path(urlsplit(loc).path or loc)
        query = _rewrite_next_query(urlsplit(loc).query)
        return urlunsplit(("", "", path, query, ""))

    return loc


def _response_headers(
    items: Iterable[tuple[str, str]],
    *,
    request: Request | None = None,
) -> dict[str, str]:
    out: dict[str, str] = {}
    for key, value in items:
        lk = key.lower()
        if lk in _HOP_BY_HOP:
            continue
        if lk == "location" and request is not None:
            out[key] = rewrite_ls_proxy_location(value, request)
            continue
        out[key] = value
    return out


@router.api_route("/label-studio", methods=["GET", "POST", "PUT", "PATCH", "DELETE", "OPTIONS", "HEAD"])
@router.api_route("/label-studio/{full_path:path}", methods=["GET", "POST", "PUT", "PATCH", "DELETE", "OPTIONS", "HEAD"])
async def proxy_label_studio(request: Request, full_path: str = "") -> Response:
    path = full_path.lstrip("/")
    target = _upstream_url(path, request.url.query)

    body = await request.body()
    headers = _forwarding_headers(request)

    async with httpx.AsyncClient(follow_redirects=False, timeout=_PROXY_TIMEOUT) as client:
        upstream = await client.request(
            request.method,
            target,
            headers=headers,
            content=body if body else None,
        )

    return Response(
        content=upstream.content,
        status_code=upstream.status_code,
        headers=_response_headers(upstream.headers.items(), request=request),
        media_type=upstream.headers.get("content-type"),
    )
