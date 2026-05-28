# -*- coding: utf-8 -*-
"""周报撰写助手 — 终结态 JSON + 多轮记忆闭环。"""
from __future__ import annotations

from typing import Any, Dict

from gibh_agent.skills._prompt_skill_terminal import execute_terminal_spec_from_kwargs
from gibh_agent.skills.base_skill import BaseSkill
from gibh_agent.skills.launch_skill_demos import apply_launch_demo_defaults

_SPEC_ID = "weekly_report_writer"


class WeeklyReportWriterSkill(BaseSkill):
    __abstractskill__ = False

    skill_id = "weekly_report_writer"
    display_name = "周报撰写助手"
    description = (
        "按日期范围与笔记库日志合成周报/双周报：个人状态同步 + 团队对外同步，"
        "支持待办继承与 [[Wiki Link]] 溯源。"
    )
    category = "其他技能"
    sub_category = "文本处理"
    aliases = ["周报", "双周报", "weekly report", "工作汇报"]
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
                "user_request": (user_request or kwargs.get("analysis_request") or "").strip(),
                "context": (context or kwargs.get("context") or "").strip(),
                "file_path": (file_path or kwargs.get("vault_path") or "").strip(),
            },
        )
        return execute_terminal_spec_from_kwargs(_SPEC_ID, filled, kwargs, max_tokens=8192)
