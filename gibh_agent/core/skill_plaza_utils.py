# -*- coding: utf-8 -*-
"""
技能广场：不落库的「已实现」推断与排序（仅基于 prompt_template 快车道暗号）。
"""
from __future__ import annotations

from datetime import datetime
from typing import Any, Dict, List, Optional

__all__ = [
    "infer_skill_implemented_from_prompt",
    "is_prompt_placeholder_only",
    "is_multimodal_skill_hidden_from_plaza",
    "skill_payload_created_ts",
    "sort_skill_payloads_inplace",
]

# 曾用于隐藏未接链路的子类；三大组学快车道已接 Mock 全 DAG，橱窗与 GET /api/skills 均展示
MULTIMODAL_PLAZA_HIDDEN_SUBCATEGORIES = frozenset()


def infer_skill_implemented_from_prompt(prompt_template: Optional[str]) -> bool:
    """含快车道路由暗号则视为已对接可执行链路（用于置顶与标签）。"""
    p = prompt_template or ""
    return "[Skill_Route:" in p or "[Omics_Route:" in p


def is_multimodal_skill_hidden_from_plaza(
    main_category: Optional[str],
    sub_category: Optional[str],
) -> bool:
    """
    「多模态组学」下若 sub_category 命中 MULTIMODAL_PLAZA_HIDDEN_SUBCATEGORIES 则 API 侧不展示。
    集合默认可为空；需要灰度隐藏某子类时在此维护子类名，并同步前端 `filterDeferredMultimodalSkillsForPlaza`。
    """
    if (main_category or "").strip() != "多模态组学":
        return False
    return (sub_category or "").strip() in MULTIMODAL_PLAZA_HIDDEN_SUBCATEGORIES


def is_prompt_placeholder_only(prompt_template: Optional[str]) -> bool:
    """
    与 seed_skills.PLACEHOLDER_PROMPT 全文一致时视为「仅占位、未落地」。
    用于橱窗/API 列表过滤；**不修改数据库**，仅控制展示。
    """
    try:
        from gibh_agent.db.seed_skills import PLACEHOLDER_PROMPT
    except Exception:
        return False
    a = " ".join((prompt_template or "").split())
    b = " ".join(PLACEHOLDER_PROMPT.split())
    return len(a) > 0 and a == b


def skill_payload_created_ts(iso_created: Optional[str]) -> float:
    if not iso_created:
        return 0.0
    try:
        return datetime.fromisoformat(iso_created.replace("Z", "+00:00")).timestamp()
    except Exception:
        return 0.0


def sort_skill_payloads_inplace(payloads: List[Dict[str, Any]]) -> None:
    """已实现（快车道）在前；其次 created_at 降序。"""
    payloads.sort(
        key=lambda x: (
            not x.get("is_implemented"),
            -skill_payload_created_ts(x.get("created_at")),
        )
    )
