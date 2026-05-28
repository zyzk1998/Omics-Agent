# -*- coding: utf-8 -*-
"""学术摘要精炼 — Prompt-Engineered 软技能。"""
from __future__ import annotations

from typing import Any, Dict

from gibh_agent.skills._prompt_skill_terminal import execute_terminal_spec_from_kwargs
from gibh_agent.skills.base_skill import BaseSkill
from gibh_agent.skills.launch_skill_demos import apply_launch_demo_defaults

_SPEC_ID = "academic_abstract_refiner"


class AcademicAbstractRefinerSkill(BaseSkill):
    __abstractskill__ = False

    skill_id = "academic_abstract_refiner"
    display_name = "学术摘要精炼"
    description = (
        "将长文草稿精炼为 SCI 风格无小标题单段摘要，输出中英双语 Summary_Report Markdown。"
    )
    category = "其他技能"
    sub_category = "文本处理"
    aliases = ["摘要", "abstract", "SCI摘要", "双语摘要"]
    __dependencies__: list[str] = []

    def execute(
        self,
        user_request: str = "",
        context: str = "",
        file_path: str = "",
        **kwargs: Any,
    ) -> Dict[str, Any]:
        filled = apply_launch_demo_defaults(
            self.skill_id,
            {
                "user_request": (user_request or kwargs.get("source_text") or "").strip(),
                "context": (context or "").strip(),
                "file_path": (file_path or "").strip(),
            },
        )
        return execute_terminal_spec_from_kwargs(_SPEC_ID, filled, kwargs, max_tokens=6144)
