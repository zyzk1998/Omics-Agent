# -*- coding: utf-8 -*-
"""HITL 会话 ↔ Label Studio 项目映射（供 Webhook 反查 session_id）。"""
from __future__ import annotations

import json
import logging
import os
from pathlib import Path
from typing import Any, Dict, Optional

logger = logging.getLogger(__name__)

_REGISTRY_DIR = Path(
    os.getenv("GIBH_HITL_REGISTRY_DIR", "data/hitl_session_registry")
).expanduser()


def _path(session_id: str) -> Path:
    safe = "".join(c if c.isalnum() or c in "-_" else "_" for c in str(session_id))
    return _REGISTRY_DIR / f"{safe}.json"


def register_hitl_session(session_id: str, payload: Dict[str, Any]) -> None:
    sid = str(session_id or "").strip()
    if not sid:
        return
    _REGISTRY_DIR.mkdir(parents=True, exist_ok=True)
    data = dict(payload or {})
    data["session_id"] = sid
    _path(sid).write_text(json.dumps(data, ensure_ascii=False, indent=2), encoding="utf-8")
    pid = data.get("project_id")
    if pid is not None:
        (_REGISTRY_DIR / f"project_{int(pid)}.json").write_text(
            json.dumps({"session_id": sid, "project_id": int(pid)}, ensure_ascii=False),
            encoding="utf-8",
        )
    logger.info("[HITL Registry] registered session=%s project_id=%s", sid, pid)


def get_hitl_session(session_id: str) -> Optional[Dict[str, Any]]:
    sid = str(session_id or "").strip()
    if not sid:
        return None
    p = _path(sid)
    if not p.is_file():
        return None
    try:
        return json.loads(p.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return None


def resolve_session_by_project(project_id: int) -> Optional[str]:
    p = _REGISTRY_DIR / f"project_{int(project_id)}.json"
    if not p.is_file():
        return None
    try:
        data = json.loads(p.read_text(encoding="utf-8"))
        return str(data.get("session_id") or "").strip() or None
    except (OSError, json.JSONDecodeError):
        return None
