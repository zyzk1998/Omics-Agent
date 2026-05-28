# -*- coding: utf-8 -*-
"""定制化简历 — 终结态 JSON + 多轮记忆闭环。"""
from __future__ import annotations

from typing import Any, Dict

from gibh_agent.skills._prompt_skill_terminal import execute_terminal_spec_from_kwargs
from gibh_agent.skills.base_skill import BaseSkill
from gibh_agent.skills.launch_skill_demos import apply_launch_demo_defaults

_SPEC_ID = "tailored_resume"


class TailoredResumeSkill(BaseSkill):
    __abstractskill__ = False

    skill_id = "tailored_resume"
    display_name = "定制化简历生成"
    description = "根据职位描述（JD）与候选人背景生成 ATS 友好、关键词优化的 Markdown 简历及求职建议。"
    category = "其他技能"
    sub_category = "文本处理"
    aliases = ["简历", "resume", "CV", "求职简历"]
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
                "user_request": (user_request or kwargs.get("job_description") or "").strip(),
                "context": (context or kwargs.get("background") or "").strip(),
                "file_path": (file_path or "").strip(),
            },
        )
        return execute_terminal_spec_from_kwargs(_SPEC_ID, filled, kwargs, max_tokens=8192)
