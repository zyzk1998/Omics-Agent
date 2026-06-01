# -*- coding: utf-8 -*-
"""蛋白质功能假说推演 — 纯 Prompt → Markdown。"""
from __future__ import annotations

from typing import Any, Dict

from gibh_agent.skills._array1_batch2_system import PROTEIN_FUNCTION_HYPOTHESIS
from gibh_agent.skills._prompt_skill_markdown import execute_markdown_prompt_skill
from gibh_agent.skills.base_skill import BaseSkill


class ProteinFunctionHypothesisSkill(BaseSkill):
    __abstractskill__ = False

    skill_id = "protein_function_hypothesis"
    display_name = "蛋白质功能假说推演"
    description = "从用户给定结构域/修饰/定位注释推演功能假说与验证路线（含 Mermaid）。"
    category = "生物医药"
    sub_category = "文本处理"
    aliases = ['蛋白功能', '结构域', '功能注释']
    required_parameters = ['蛋白 ID、结构域、修饰、亚细胞定位等注释（粘贴至 context）']
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
            system_prompt=PROTEIN_FUNCTION_HYPOTHESIS,
            user_request=(user_request or "").strip(),
            context=(context or "").strip(),
            deliver_message="蛋白功能假说报告已生成，请在右侧工作台查看",
        )
