# -*- coding: utf-8 -*-
"""
聊天 API 路由常量与 enabled_mcps 归一化。

HPC/Open（compute_scheduler）请求与普聊一致：一律经 server.py → AgentOrchestrator；启用超算时由 server 侧强制
thinking_mode=deep（DeepReActRunner + 动态 MCP schema），不再使用 hpc_orchestrator 直连短路。

本模块不重复注册 /api/chat（入口仍在 server.py）。
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
    """
    已废弃：曾用于将 chat SSE 引流至 hpc_orchestrator（绕过 LLM）。
    现恒为 False；保留函数以免外部脚本 import 报错。
    """
    return False
