# -*- coding: utf-8 -*-
"""邮件管理助手 — Prompt-Engineered 软技能（仅生成文稿，不访问邮箱）。"""
from __future__ import annotations

from typing import Any, Dict

from gibh_agent.skills._prompt_skill_terminal import execute_terminal_spec_from_kwargs
from gibh_agent.skills.base_skill import BaseSkill
from gibh_agent.skills.launch_skill_demos import apply_launch_demo_defaults

_SPEC_ID = "email_manager"


class EmailManagerSkill(BaseSkill):
    __abstractskill__ = False

    skill_id = "email_manager"
    display_name = "邮件管理助手"
    description = (
        "撰写/润色商务邮件、跟进信、冷启动 outreach、主题行与多语气回复选项；"
        "不连接邮箱，仅生成可复制文稿。"
    )
    category = "其他技能"
    sub_category = "文本处理"
    aliases = ["邮件", "email", "写邮件", "follow up email"]
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
                "user_request": (user_request or kwargs.get("email_context") or "").strip(),
                "context": (context or "").strip(),
            },
        )
        return execute_terminal_spec_from_kwargs(_SPEC_ID, filled, kwargs, max_tokens=4096)
