# -*- coding: utf-8 -*-
"""文献矩阵精读笔记 — 纯 Prompt → Markdown。"""
from __future__ import annotations

from typing import Any, Dict

from gibh_agent.skills._array1_batch2_system import LITERATURE_MATRIX_NOTES
from gibh_agent.skills._prompt_skill_markdown import execute_markdown_prompt_skill
from gibh_agent.skills.base_skill import BaseSkill


class LiteratureMatrixNotesSkill(BaseSkill):
    __abstractskill__ = False

    skill_id = "literature_matrix_notes"
    display_name = "文献矩阵精读笔记"
    description = "3–10 篇摘要→对比维度表（设计、样本、结论、局限）。"
    category = "生物医药"
    sub_category = "文本处理"
    aliases = ['文献矩阵', '精读笔记', '综述']
    required_parameters = ['多篇文献摘要或笔记（粘贴至 context）']
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
            system_prompt=LITERATURE_MATRIX_NOTES,
            user_request=(user_request or "").strip(),
            context=(context or "").strip(),
            deliver_message="文献矩阵笔记已生成，请在右侧工作台查看",
        )
