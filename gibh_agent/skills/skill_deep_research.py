# -*- coding: utf-8 -*-
"""深度调研 — Prompt-Engineered 软技能。"""
from __future__ import annotations

from typing import Any, Dict

from gibh_agent.skills._prompt_skill_terminal import execute_terminal_spec_from_kwargs
from gibh_agent.skills.base_skill import BaseSkill
from gibh_agent.skills.launch_skill_demos import apply_launch_demo_defaults

_SPEC_ID = "deep_research"


class DeepResearchSkill(BaseSkill):
    __abstractskill__ = False

    skill_id = "deep_research"
    display_name = "深度调研"
    description = (
        "多视角深度调研：执行摘要、关键发现、分主题分析、共识/争议与编号引用来源，"
        "输出结构化 Markdown 报告。"
    )
    category = "其他技能"
    sub_category = "信息检索"
    aliases = ["deep research", "调研", "文献综述", "深度研究"]
    __dependencies__: list[str] = []

    def execute(
        self,
        user_request: str = "",
        context: str = "",
        **kwargs: Any,
    ) -> Dict[str, Any]:
        filled = apply_launch_demo_defaults(
            self.skill_id,
            {
                "user_request": (user_request or kwargs.get("topic") or "").strip(),
                "context": (context or kwargs.get("sources") or "").strip(),
            },
        )
        return execute_terminal_spec_from_kwargs(_SPEC_ID, filled, kwargs, max_tokens=8192)
