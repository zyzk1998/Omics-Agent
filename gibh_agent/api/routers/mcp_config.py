# -*- coding: utf-8 -*-
"""HPC MCP 配置 API：重连与状态查询。"""
from __future__ import annotations

import logging

from fastapi import APIRouter
from pydantic import BaseModel, Field

from gibh_agent.mcp_client import hpc_mcp_manager

logger = logging.getLogger(__name__)

router = APIRouter(tags=["config"])


class MCPConfigBody(BaseModel):
    url: str = Field(default="http://127.0.0.1:8001/mcp", description="MCP Server URL")
    transport: str = Field(default="streamableHttp", description="streamableHttp | sse")


@router.get("/mcp/status")
async def get_mcp_status():
    return await hpc_mcp_manager.status_dict_async()


@router.post("/mcp")
async def post_mcp_config(body: MCPConfigBody):
    try:
        await hpc_mcp_manager.connect(body.url.strip(), body.transport.strip())
        return {
            "status": "success",
            "message": "已调度后台断开并重连（不阻塞主线程）",
            "config": {"url": body.url.strip().rstrip("/"), "transport": body.transport},
        }
    except Exception as e:  # noqa: BLE001
        logger.warning("POST /api/config/mcp 失败（仍返回 200 侧信息）: %s", e, exc_info=True)
        return {
            "status": "error",
            "message": str(e),
        }
