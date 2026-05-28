# -*- coding: utf-8 -*-
"""学术会议海报 — Prompt-Engineered 软技能。"""
from __future__ import annotations

from typing import Any, Dict

from gibh_agent.skills._prompt_skill_terminal import execute_terminal_spec_from_kwargs
from gibh_agent.skills.base_skill import BaseSkill
from gibh_agent.skills.launch_skill_demos import apply_launch_demo_defaults

_SPEC_ID = "academic_poster_generator"


class AcademicPosterGeneratorSkill(BaseSkill):
    __abstractskill__ = False

    skill_id = "academic_poster_generator"
    display_name = "学术会议海报生成"
    description = (
        "从论文 PDF 或摘要生成会议海报 Markdown 故事板：分节要点、配图方案与质量检查清单。"
    )
    category = "其他技能"
    sub_category = "文本处理"
    aliases = ["学术海报", "poster", "会议海报", "beamerposter"]
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
                "user_request": (user_request or "").strip(),
                "context": (context or "").strip(),
                "file_path": (file_path or kwargs.get("pdf_path") or "").strip(),
            },
        )
        return execute_terminal_spec_from_kwargs(_SPEC_ID, filled, kwargs, max_tokens=8192)
