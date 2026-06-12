# -*- coding: utf-8 -*-
"""API 容器调用宿主机 Local Sidecar 执行落盘/检测。"""
from __future__ import annotations

import base64
import logging
import os
from pathlib import Path
from typing import Any, Dict, List, Optional

import httpx

logger = logging.getLogger(__name__)


def sidecar_base_url() -> str:
    return os.getenv("OMICS_LOCAL_SIDECAR_URL", "http://host.docker.internal:8019").rstrip("/")


def _post_json(path: str, payload: Dict[str, Any], *, timeout: float = 120.0) -> Dict[str, Any]:
    url = f"{sidecar_base_url()}{path}"
    try:
        with httpx.Client(timeout=timeout) as client:
            resp = client.post(url, json=payload)
            body = resp.json() if resp.content else {}
            if resp.status_code >= 400:
                detail = body.get("detail") or body.get("message") or resp.text
                raise RuntimeError(str(detail))
            return body if isinstance(body, dict) else {"status": "success", "data": body}
    except httpx.RequestError as exc:
        logger.warning("[SidecarClient] request failed %s: %s", url, exc)
        raise RuntimeError(f"Sidecar 不可达 ({sidecar_base_url()})") from exc


def sidecar_silent_deploy(
    *,
    session_id: str,
    session_title: Optional[str],
    host_mount_path: str,
    file_specs: List[Dict[str, str]],
    session_folder: Optional[str] = None,
    folder_timestamp: Optional[str] = None,
) -> Dict[str, Any]:
    """
    将容器内已读出的文件内容经 Sidecar 写入宿主机挂载目录。
    file_specs: [{ "dest_name": "sft_corpus.json", "content_b64": "..." }, ...]
    """
    payload = {
        "session_id": session_id,
        "session_title": session_title or session_id,
        "host_mount_path": host_mount_path,
        "files": file_specs,
    }
    if session_folder:
        payload["session_folder"] = session_folder
    if folder_timestamp:
        payload["folder_timestamp"] = folder_timestamp
    return _post_json("/api/workspace/silent_deploy", payload)


def sidecar_check_local_result(
    *,
    host_mount_path: str,
    session_id: str,
    session_title: Optional[str] = None,
    session_folder: Optional[str] = None,
) -> Dict[str, Any]:
    payload = {
        "host_mount_path": host_mount_path,
        "session_id": session_id,
        "session_title": session_title or session_id,
    }
    if session_folder:
        payload["session_folder"] = session_folder
    try:
        return _post_json("/api/workspace/local_result_check", payload, timeout=30.0)
    except RuntimeError as exc:
        return {
            "has_local_result": False,
            "mount_path": host_mount_path,
            "message": str(exc),
        }


def read_file_as_b64(path: Path) -> str:
    return base64.b64encode(path.read_bytes()).decode("ascii")
