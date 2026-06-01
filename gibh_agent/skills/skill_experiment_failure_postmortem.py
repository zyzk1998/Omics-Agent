# -*- coding: utf-8 -*-
"""实验失败复盘报告 — 纯 Prompt → Markdown。"""
from __future__ import annotations

from typing import Any, Dict

from gibh_agent.skills._array1_batch2_system import EXPERIMENT_FAILURE_POSTMORTEM
from gibh_agent.skills._prompt_skill_markdown import execute_markdown_prompt_skill
from gibh_agent.skills.base_skill import BaseSkill


class ExperimentFailurePostmortemSkill(BaseSkill):
    __abstractskill__ = False

    skill_id = "experiment_failure_postmortem"
    display_name = "实验失败复盘报告"
    description = "将失败现象整理为组会可用复盘稿：原因树、排查与重复实验建议。"
    category = "生物医药"
    sub_category = "文本处理"
    aliases = ['实验失败', '复盘', '排障']
    required_parameters = ['实验目标、预期结果、实际现象与已尝试排查']
    output_type = "markdown"
    __dependencies__: list[str] = []

    def execute(
        self,
        user_request: str = "",
        context: str = "",
        **kwargs: Any,
    ) -> Dict[str, Any]:
        return execute_markdown_prompt_skill(
            skill_id=self.skill_id,
            system_prompt=EXPERIMENT_FAILURE_POSTMORTEM,
            user_request=(user_request or "").strip(),
            context=(context or "").strip(),
            deliver_message="实验失败复盘报告已生成，请在右侧工作台查看",
        )
