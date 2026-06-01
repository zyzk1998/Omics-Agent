# -*- coding: utf-8 -*-
"""统计方法选择建议书 — 纯 Prompt → Markdown。"""
from __future__ import annotations

from typing import Any, Dict

from gibh_agent.skills._array1_batch2_system import STATS_METHOD_ADVISOR
from gibh_agent.skills._prompt_skill_markdown import execute_markdown_prompt_skill
from gibh_agent.skills.base_skill import BaseSkill


class StatsMethodAdvisorSkill(BaseSkill):
    __abstractskill__ = False

    skill_id = "stats_method_advisor"
    display_name = "统计方法选择建议书"
    description = "按设计类型推荐检验方法、前提假设与报告规范提醒。"
    category = "生物医药"
    sub_category = "文本处理"
    aliases = ['统计方法', '假设检验', '样本量']
    required_parameters = ['研究设计、分组、结局变量类型与重复测量结构']
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
            system_prompt=STATS_METHOD_ADVISOR,
            user_request=(user_request or "").strip(),
            context=(context or "").strip(),
            deliver_message="统计方法建议书已生成，请在右侧工作台查看",
        )
