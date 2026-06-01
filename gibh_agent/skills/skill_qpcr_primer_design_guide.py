# -*- coding: utf-8 -*-
"""引物/qPCR 实验设计说明 — 纯 Prompt → Markdown。"""
from __future__ import annotations

from typing import Any, Dict

from gibh_agent.skills._array1_batch2_system import QPCR_PRIMER_DESIGN_GUIDE
from gibh_agent.skills._prompt_skill_markdown import execute_markdown_prompt_skill
from gibh_agent.skills.base_skill import BaseSkill


class QpcrPrimerDesignGuideSkill(BaseSkill):
    __abstractskill__ = False

    skill_id = "qpcr_primer_design_guide"
    display_name = "引物/qPCR 实验设计说明"
    description = "生成引物设计原则、对照设置与扩增/qPCR 参数建议（序列由用户提供）。"
    category = "生物医药"
    sub_category = "文本处理"
    aliases = ['qPCR', '引物设计', 'RT-qPCR', '内参基因']
    required_parameters = ['扩增目标、样本类型与已有引物序列（如有）']
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
            system_prompt=QPCR_PRIMER_DESIGN_GUIDE,
            user_request=(user_request or "").strip(),
            context=(context or "").strip(),
            deliver_message="qPCR 实验设计说明已生成，请在右侧工作台查看",
        )
