# -*- coding: utf-8 -*-
"""
用户数据库挂载配置 API（Phase 1 · 三轨合一入库前置配置中心）。

路由前缀由 server.py 挂载为 /api/settings。
"""
from __future__ import annotations

import logging
from typing import Any, Dict, Literal, Optional

from fastapi import APIRouter, Depends
from pydantic import BaseModel, Field

from gibh_agent.core.deps import get_current_owner_id
from gibh_agent.core.user_database_settings_store import (
    get_database_mount_config,
    probe_database_mount_connection,
    save_database_mount_config,
)

logger = logging.getLogger(__name__)

router = APIRouter(prefix="/api/settings", tags=["UserSettings"])


class LocalVolumeConfig(BaseModel):
    mount_path: str = Field(default="", description="本地挂载目录绝对路径")


class HpcSlurmConfig(BaseModel):
    host: str = Field(default="127.0.0.1", description="HPC/Slurm 登录节点 IP 或主机名")
    port: int = Field(default=22, ge=1, le=65535, description="SSH 端口")
    username: str = Field(default="", description="SSH 用户名（暂存）")
    base_path: str = Field(default="", description="远端数据根路径（暂存）")


class ApiUrlConfig(BaseModel):
    endpoint: str = Field(default="", description="REST API 基址或入库 Endpoint")
    token: str = Field(default="", description="Bearer Token（暂存，响应时可脱敏）")


class DatabaseMountSettingsBody(BaseModel):
    """
    PUT /api/settings/database 请求体。

    mount_type 决定生效的配置块；其余块可一并提交以便前端表单保留草稿。
    """

    mount_type: Literal["local_volume", "hpc_slurm", "api_url"] = Field(
        default="local_volume",
        description="挂载类型：本地卷 | HPC/Slurm | API/URL",
    )
    is_auto_ingestion_enabled: bool = Field(
        default=False,
        description="是否在 pipeline 完成后自动触发三轨合一入库",
    )
    local_volume: LocalVolumeConfig = Field(default_factory=LocalVolumeConfig)
    hpc_slurm: HpcSlurmConfig = Field(default_factory=HpcSlurmConfig)
    api_url: ApiUrlConfig = Field(default_factory=ApiUrlConfig)


class DatabaseMountTestBody(DatabaseMountSettingsBody):
    """POST /api/settings/database/test-connection 与 PUT 结构相同，仅做连通性探测。"""


def _mask_token(config: Dict[str, Any]) -> Dict[str, Any]:
    """GET 响应中脱敏 token。"""
    out = dict(config)
    api = dict(out.get("api_url") or {})
    token = str(api.get("token") or "")
    if token:
        api["token"] = token[:4] + "***" if len(token) > 4 else "***"
    out["api_url"] = api
    return out


@router.get("/database")
def get_database_settings(
    owner_id: str = Depends(get_current_owner_id),
) -> Dict[str, Any]:
    """读取当前用户的数据库挂载配置。"""
    cfg = get_database_mount_config(owner_id)
    return {"status": "success", "config": _mask_token(cfg)}


@router.put("/database")
def put_database_settings(
    body: DatabaseMountSettingsBody,
    owner_id: str = Depends(get_current_owner_id),
) -> Dict[str, Any]:
    """保存当前用户的数据库挂载配置。"""
    try:
        saved = save_database_mount_config(owner_id, body.model_dump())
    except ValueError as exc:
        return {"status": "error", "message": str(exc)}
    return {
        "status": "success",
        "message": "数据库挂载配置已保存",
        "config": _mask_token(saved),
    }


@router.post("/database/test-connection")
def post_database_test_connection(
    body: DatabaseMountTestBody,
    owner_id: str = Depends(get_current_owner_id),
) -> Dict[str, Any]:
    """关联测试：按 mount_type 执行基础连通性/路径校验。"""
    result = probe_database_mount_connection(body.model_dump())
    return {
        "status": result.get("status", "error"),
        "message": result.get("message", ""),
        "detail": {k: v for k, v in result.items() if k not in ("status", "message")},
    }
