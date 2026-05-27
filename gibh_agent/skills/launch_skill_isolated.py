# -*- coding: utf-8 -*-
"""首发技能隔离执行：api-server 经 HTTP 委托 launch-skills 微服务（TaaS 轻量栈）。"""
from __future__ import annotations

import logging
import os
from typing import Any, Dict, FrozenSet

import requests

from gibh_agent.skills.launch_skill_demos import LAUNCH_SKILL_DEMO_ARGS

logger = logging.getLogger(__name__)

LAUNCH_CORE_TOOL_IDS: FrozenSet[str] = frozenset(LAUNCH_SKILL_DEMO_ARGS.keys())


def launch_skills_in_worker() -> bool:
    return os.getenv("LAUNCH_SKILLS_IN_WORKER", "").strip() in ("1", "true", "yes")


def launch_skills_base_url() -> str:
    return (os.getenv("LAUNCH_SKILLS_BASE_URL") or "").strip().rstrip("/")


def should_delegate_to_launch_worker(skill_id: str) -> bool:
    if launch_skills_in_worker():
        return False
    base = launch_skills_base_url()
    if not base:
        return False
    return (skill_id or "").strip() in LAUNCH_CORE_TOOL_IDS


def delegate_launch_skill_execute(skill_id: str, kwargs: Dict[str, Any]) -> Dict[str, Any]:
    base = launch_skills_base_url()
    url = f"{base}/api/launch/run"
    timeout_s = float(os.getenv("LAUNCH_SKILLS_HTTP_TIMEOUT", "600"))
    payload = {"tool_id": skill_id, "kwargs": kwargs or {}}
    try:
        resp = requests.post(url, json=payload, timeout=timeout_s)
        resp.raise_for_status()
        data = resp.json()
        if isinstance(data, dict):
            return data
        return {"status": "error", "message": f"launch-skills 返回非对象: {type(data).__name__}"}
    except requests.exceptions.ConnectionError:
        return {
            "status": "error",
            "message": (
                f"无法连接 launch-skills（{base}）。请执行: "
                "docker compose up -d launch-skills"
            ),
        }
    except requests.exceptions.Timeout:
        return {"status": "error", "message": f"launch-skills 请求超时（>{timeout_s}s）"}
    except requests.exceptions.HTTPError as exc:
        body = (getattr(exc.response, "text", None) or "")[:800]
        return {"status": "error", "message": f"launch-skills HTTP 错误: {exc}; {body}"}
    except Exception as exc:
        logger.exception("delegate_launch_skill_execute failed")
        return {"status": "error", "message": str(exc)}
