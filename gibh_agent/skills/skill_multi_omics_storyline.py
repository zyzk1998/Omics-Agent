# -*- coding: utf-8 -*-
"""多组学整合分析故事线 — 纯 Prompt → Markdown。"""
from __future__ import annotations

from typing import Any, Dict

from gibh_agent.skills._array1_batch2_system import MULTI_OMICS_STORYLINE
from gibh_agent.skills._prompt_skill_markdown import execute_markdown_prompt_skill
from gibh_agent.skills.base_skill import BaseSkill


class MultiOmicsStorylineSkill(BaseSkill):
    __abstractskill__ = False

    skill_id = "multi_omics_storyline"
    display_name = "多组学整合分析故事线"
    description = "撰写基因组→转录组→蛋白组逻辑链条与图表编排建议（含 Mermaid，无真实运算）。"
    category = "生物医药"
    sub_category = "文本处理"
    aliases = ['多组学', '故事线', 'Figure编排']
    required_parameters = ['各组学层级的关键发现文字摘要']
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
            system_prompt=MULTI_OMICS_STORYLINE,
            user_request=(user_request or "").strip(),
            context=(context or "").strip(),
            deliver_message="多组学故事线已生成，请在右侧工作台查看",
        )
