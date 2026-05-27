# -*- coding: utf-8 -*-
"""动态技能包 — 由 skill_registry 在启动时扫描本目录下的 ``skill_*.py``。"""
from __future__ import annotations

import logging

logger = logging.getLogger(__name__)


def bootstrap_skills(*, force_reload: bool = False) -> dict:
    from gibh_agent.core.skill_registry import discover_and_register_skills

    return discover_and_register_skills(force_reload=force_reload)
