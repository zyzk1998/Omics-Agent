# -*- coding: utf-8 -*-
"""
用户数据库挂载配置持久化（Phase 1 · 内存 + JSON 文件，按 owner_id 隔离）。

后续 Phase 可迁移至 ORM 表；当前不改动核心工作流生命周期。
"""
from __future__ import annotations

import json
import logging
import os
import socket
from pathlib import Path
from typing import Any, Dict, Literal, Optional
from urllib.parse import urlparse

import httpx

logger = logging.getLogger(__name__)

MountType = Literal["local_volume", "hpc_slurm", "api_url"]

_SETTINGS_DIR = Path(
    os.getenv("GIBH_USER_DATABASE_SETTINGS_DIR", "data/user_database_settings")
).expanduser()
_CACHE: Dict[str, Dict[str, Any]] = {}


def _settings_path(owner_id: str) -> Path:
    safe = "".join(c if c.isalnum() or c in "-_" else "_" for c in str(owner_id))
    return _SETTINGS_DIR / f"{safe}.json"


def get_database_mount_config(owner_id: str) -> Dict[str, Any]:
    """读取用户数据库挂载配置；不存在则返回默认空配置。"""
    oid = str(owner_id or "").strip()
    if not oid:
        return _default_config()

    if oid in _CACHE:
        return dict(_CACHE[oid])

    path = _settings_path(oid)
    if path.is_file():
        try:
            data = json.loads(path.read_text(encoding="utf-8"))
            if isinstance(data, dict):
                merged = _default_config()
                merged.update(data)
                _CACHE[oid] = merged
                return dict(merged)
        except (OSError, json.JSONDecodeError) as exc:
            logger.warning("读取用户数据库配置失败 owner=%s: %s", oid, exc)

    return _default_config()


def save_database_mount_config(owner_id: str, config: Dict[str, Any]) -> Dict[str, Any]:
    """校验并保存用户数据库挂载配置。"""
    oid = str(owner_id or "").strip()
    if not oid:
        raise ValueError("owner_id 不能为空")

    normalized = _normalize_config(config)
    _SETTINGS_DIR.mkdir(parents=True, exist_ok=True)
    path = _settings_path(oid)
    path.write_text(json.dumps(normalized, ensure_ascii=False, indent=2), encoding="utf-8")
    _CACHE[oid] = dict(normalized)
    logger.info("用户数据库挂载配置已保存: owner=%s mount_type=%s", oid, normalized.get("mount_type"))
    return dict(normalized)


def probe_database_mount_connection(config: Dict[str, Any]) -> Dict[str, Any]:
    """
    模拟/基础连通性测试，供前端「关联测试」按钮调用。

    - local_volume：检查 mount_path 是否存在且为目录
    - hpc_slurm：TCP 探测 host:port（默认 22）
    - api_url：HTTP HEAD/GET 探测 endpoint（可选 Bearer Token）
    """
    normalized = _normalize_config(config)
    mount_type = normalized["mount_type"]

    if mount_type == "local_volume":
        return _test_local_volume(normalized.get("local_volume") or {})
    if mount_type == "hpc_slurm":
        return _test_hpc_slurm(normalized.get("hpc_slurm") or {})
    if mount_type == "api_url":
        return _test_api_url(normalized.get("api_url") or {})

    return {"status": "error", "message": f"未知 mount_type: {mount_type}"}


def _default_config() -> Dict[str, Any]:
    return {
        "mount_type": "local_volume",
        "is_auto_ingestion_enabled": False,
        "local_volume": {"mount_path": ""},
        "hpc_slurm": {"host": "127.0.0.1", "port": 22, "username": "", "base_path": ""},
        "api_url": {"endpoint": "", "token": ""},
    }


def _normalize_config(raw: Dict[str, Any]) -> Dict[str, Any]:
    base = _default_config()
    if not isinstance(raw, dict):
        return base

    mount_type = str(raw.get("mount_type") or "local_volume").strip()
    if mount_type not in ("local_volume", "hpc_slurm", "api_url"):
        raise ValueError("mount_type 必须是 local_volume | hpc_slurm | api_url")

    base["mount_type"] = mount_type
    base["is_auto_ingestion_enabled"] = bool(raw.get("is_auto_ingestion_enabled", False))

    lv = raw.get("local_volume") if isinstance(raw.get("local_volume"), dict) else {}
    base["local_volume"] = {
        "mount_path": str((lv or {}).get("mount_path") or "").strip(),
    }

    hpc = raw.get("hpc_slurm") if isinstance(raw.get("hpc_slurm"), dict) else {}
    port_raw = (hpc or {}).get("port", 22)
    try:
        port = int(port_raw)
    except (TypeError, ValueError):
        port = 22
    base["hpc_slurm"] = {
        "host": str((hpc or {}).get("host") or "127.0.0.1").strip(),
        "port": port,
        "username": str((hpc or {}).get("username") or "").strip(),
        "base_path": str((hpc or {}).get("base_path") or "").strip(),
    }

    api = raw.get("api_url") if isinstance(raw.get("api_url"), dict) else {}
    base["api_url"] = {
        "endpoint": str((api or {}).get("endpoint") or "").strip(),
        "token": str((api or {}).get("token") or "").strip(),
    }
    return base


def _test_local_volume(cfg: Dict[str, Any]) -> Dict[str, Any]:
    mount_path = str(cfg.get("mount_path") or "").strip()
    if not mount_path:
        return {"status": "error", "message": "local_volume.mount_path 不能为空"}
    path = Path(mount_path).expanduser()
    if not path.exists():
        return {"status": "error", "message": f"路径不存在: {path}"}
    if not path.is_dir():
        return {"status": "error", "message": f"路径不是目录: {path}"}
    return {
        "status": "success",
        "message": f"本地目录可访问: {path}",
        "resolved_path": str(path.resolve()),
    }


def _test_hpc_slurm(cfg: Dict[str, Any]) -> Dict[str, Any]:
    host = str(cfg.get("host") or "").strip()
    if not host:
        return {"status": "error", "message": "hpc_slurm.host 不能为空"}
    port = int(cfg.get("port") or 22)
    try:
        with socket.create_connection((host, port), timeout=5.0):
            pass
    except OSError as exc:
        return {
            "status": "error",
            "message": f"无法连接 HPC/Slurm 主机 {host}:{port} — {exc}",
        }
    return {
        "status": "success",
        "message": f"HPC/Slurm 主机 TCP 可达: {host}:{port}（Slurm 作业提交尚未在本阶段验证）",
        "host": host,
        "port": port,
    }


def _test_api_url(cfg: Dict[str, Any]) -> Dict[str, Any]:
    endpoint = str(cfg.get("endpoint") or "").strip()
    if not endpoint:
        return {"status": "error", "message": "api_url.endpoint 不能为空"}
    parsed = urlparse(endpoint)
    if parsed.scheme not in ("http", "https"):
        return {"status": "error", "message": "api_url.endpoint 须为 http(s) URL"}

    headers: Dict[str, str] = {}
    token = str(cfg.get("token") or "").strip()
    if token:
        headers["Authorization"] = f"Bearer {token}"

    try:
        with httpx.Client(timeout=10.0, follow_redirects=True) as client:
            resp = client.get(endpoint, headers=headers)
    except httpx.RequestError as exc:
        return {"status": "error", "message": f"API 端点不可达: {exc}"}

    if resp.status_code >= 400:
        return {
            "status": "error",
            "message": f"API 返回 HTTP {resp.status_code}",
            "http_status": resp.status_code,
        }
    return {
        "status": "success",
        "message": f"API 端点可达 (HTTP {resp.status_code})",
        "http_status": resp.status_code,
    }
