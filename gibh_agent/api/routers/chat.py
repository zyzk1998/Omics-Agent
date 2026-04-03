# -*- coding: utf-8 -*-
"""
聊天 API 分流：compute_scheduler 开启时由 server 主动引流至 hpc_orchestrator，不实例化 AgentOrchestrator。

本模块不重复注册 /api/chat（入口仍在 server.py），仅提供判定与常量，满足「路由在 api/routers/chat 中集中声明」的架构约束。
"""
from __future__ import annotations

import json
from typing import Any, List

COMPUTE_SCHEDULER_MCP_KEY = "compute_scheduler"


def normalize_enabled_mcps(raw: Any) -> List[str]:
    if raw is None:
        return []
    if isinstance(raw, str):
        try:
            parsed = json.loads(raw)
            if isinstance(parsed, list):
                raw = parsed
            else:
                s = raw.strip()
                return [s] if s else []
        except Exception:
            s = raw.strip()
            return [s] if s else []
    if not isinstance(raw, list):
        return []
    return [str(x).strip() for x in raw if x is not None and str(x).strip()]


def should_use_hpc_isolated_route(enabled_mcps: Any) -> bool:
    """为 True 时 chat SSE 应走 gibh_agent.core.hpc_orchestrator.stream_hpc_direct_chat。"""
    return COMPUTE_SCHEDULER_MCP_KEY in set(normalize_enabled_mcps(enabled_mcps))
