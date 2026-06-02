# -*- coding: utf-8 -*-
"""技能广场 detailed_spec 附加逻辑（API 与 server 共用）。"""
from __future__ import annotations

from typing import Any, Dict, Optional

from gibh_agent.db.skill_detailed_specs import enrich_skill_plaza_item


def build_skill_plaza_payload(
    *,
    skill_id: Any,
    name: str,
    description: Optional[str],
    main_category: Optional[str],
    sub_category: Optional[str],
    prompt_template: Optional[str],
    is_implemented: bool,
    author_id: Optional[str],
    created_at: Optional[str],
    saved: bool,
    db_detailed_spec: Any = None,
    extra: Optional[Dict[str, Any]] = None,
) -> Dict[str, Any]:
    """构造 GET /api/skills 单条响应并附加 detailed_spec。"""
    item: Dict[str, Any] = {
        "id": skill_id,
        "name": name,
        "description": description,
        "main_category": main_category,
        "sub_category": sub_category,
        "prompt_template": prompt_template,
        "is_implemented": is_implemented,
        "author_id": author_id,
        "created_at": created_at,
        "saved": saved,
    }
    if extra:
        item.update(extra)
    return enrich_skill_plaza_item(item, db_detailed_spec=db_detailed_spec)
