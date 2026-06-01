# -*- coding: utf-8 -*-
"""变异临床意义解读草案（ACMG 风格） — 纯 Prompt → Markdown。"""
from __future__ import annotations

from typing import Any, Dict

from gibh_agent.skills._array1_batch2_system import ACMG_VARIANT_INTERPRETATION
from gibh_agent.skills._prompt_skill_markdown import execute_markdown_prompt_skill
from gibh_agent.skills.base_skill import BaseSkill


class AcmgVariantInterpretationSkill(BaseSkill):
    __abstractskill__ = False

    skill_id = "acmg_variant_interpretation"
    display_name = "变异临床意义解读草案（ACMG 风格）"
    description = "根据变异位点、基因与人群频率文字生成 ACMG 风格证据条目草稿（须标注需数据库核实）。"
    category = "生物医药"
    sub_category = "文本处理"
    aliases = ['ACMG', '变异解读', '临床意义', '胚系变异', 'VUS']
    required_parameters = ['变异位点、基因名、转录本、人群频率与功能预测（粘贴至 context）']
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
            system_prompt=ACMG_VARIANT_INTERPRETATION,
            user_request=(user_request or "").strip(),
            context=(context or "").strip(),
            deliver_message="变异解读草稿已生成，请在右侧工作台查看",
        )
