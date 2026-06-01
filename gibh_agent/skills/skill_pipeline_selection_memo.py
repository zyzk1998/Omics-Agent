# -*- coding: utf-8 -*-
"""生物信息 Pipeline 选型备忘录 — 纯 Prompt → Markdown。"""
from __future__ import annotations

from typing import Any, Dict

from gibh_agent.skills._array1_batch2_system import PIPELINE_SELECTION_MEMO
from gibh_agent.skills._prompt_skill_markdown import execute_markdown_prompt_skill
from gibh_agent.skills.base_skill import BaseSkill


class PipelineSelectionMemoSkill(BaseSkill):
    __abstractskill__ = False

    skill_id = "pipeline_selection_memo"
    display_name = "生物信息 Pipeline 选型备忘录"
    description = "针对 RNA-seq/WGS/蛋白组等目标对比 2–3 套工具链优缺点（不执行计算）。"
    category = "生物医药"
    sub_category = "文本处理"
    aliases = ['流程选型', 'pipeline', 'RNA-seq流程', '分析方案']
    required_parameters = ['分析目标、数据类型与约束（样本量、预算、时间）']
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
            system_prompt=PIPELINE_SELECTION_MEMO,
            user_request=(user_request or "").strip(),
            context=(context or "").strip(),
            deliver_message="Pipeline 选型备忘录已生成，请在右侧工作台查看",
        )
