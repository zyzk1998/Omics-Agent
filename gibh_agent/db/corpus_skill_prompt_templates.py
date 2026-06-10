# -*- coding: utf-8 -*-
"""
科学语料数据加工技能 prompt_template（与 seed_skills / patch_corpus_skill 同源）。
"""
from __future__ import annotations

from gibh_agent.db.corpus_skill_user_guide import (
    CORPUS_SKILL_DESCRIPTION_SHORT,
    CORPUS_SKILL_PROMPT_TEMPLATE,
    CORPUS_SKILL_USER_GUIDE,
)

CORPUS_SKILL_NAME = "科学语料数据加工"

CORPUS_SKILL_SEED: dict[str, str] = {
    "name": CORPUS_SKILL_NAME,
    "main_category": "特色科研流程",
    "sub_category": "语料加工",
    "description": CORPUS_SKILL_DESCRIPTION_SHORT,
}
