# -*- coding: utf-8 -*-
"""IM Gateway：接收钉钉 / 企业微信 / 飞书 / Slack Webhook，异步转发至主 API。"""
from __future__ import annotations

import logging
from typing import Any, Dict

from fastapi import FastAPI, Request

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

app = FastAPI(title="GIBH IM Gateway", version="0.1.0")


@app.get("/ping")
async def ping() -> Dict[str, str]:
    """健康检查（容器 healthcheck / 运维探活）。"""
    return {"status": "ok"}


@app.get("/health")
async def health() -> Dict[str, str]:
    return {"status": "ok"}


@app.post("/api/webhook/dingtalk")
async def webhook_dingtalk(request: Request) -> Dict[str, str]:
    """
    钉钉事件回调占位。
    正式实现：校验 Signature → 入队 Celery / 后台任务 → 调用 AgentOrchestrator。
    """
    body: Any = {}
    try:
        body = await request.json()
    except Exception:  # noqa: BLE001
        body = {}
    logger.info("[im-gateway] dingtalk webhook placeholder, keys=%s", list(body.keys()) if isinstance(body, dict) else type(body).__name__)
    return {"status": "ok", "message": "placeholder"}


@app.post("/api/webhook/wecom")
async def webhook_wecom() -> Dict[str, str]:
    """企业微信回调占位。"""
    return {"status": "ok", "message": "placeholder"}


@app.post("/api/webhook/feishu")
async def webhook_feishu() -> Dict[str, str]:
    """飞书事件回调占位。"""
    return {"status": "ok", "message": "placeholder"}


@app.post("/api/webhook/slack")
async def webhook_slack() -> Dict[str, str]:
    """Slack Events API 占位。"""
    return {"status": "ok", "message": "placeholder"}
