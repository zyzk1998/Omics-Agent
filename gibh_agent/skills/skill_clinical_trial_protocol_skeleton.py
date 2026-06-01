# -*- coding: utf-8 -*-
"""临床注册方案骨架生成（CONSORT 导向） — 纯 Prompt → Markdown。"""
from __future__ import annotations

from typing import Any, Dict

from gibh_agent.skills._array1_batch2_system import CLINICAL_TRIAL_PROTOCOL_SKELETON
from gibh_agent.skills._prompt_skill_markdown import execute_markdown_prompt_skill
from gibh_agent.skills.base_skill import BaseSkill


class ClinicalTrialProtocolSkeletonSkill(BaseSkill):
    __abstractskill__ = False

    skill_id = "clinical_trial_protocol_skeleton"
    display_name = "临床注册方案骨架生成（CONSORT 导向）"
    description = "输出试验目的、入排标准、主要/次要终点与统计假设框架骨架。"
    category = "生物医药"
    sub_category = "文本处理"
    aliases = ['临床试验', 'CONSORT', '方案骨架', '注册研究']
    required_parameters = ['干预类型、目标人群与主要终点描述']
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
            system_prompt=CLINICAL_TRIAL_PROTOCOL_SKELETON,
            user_request=(user_request or "").strip(),
            context=(context or "").strip(),
            deliver_message="临床方案骨架已生成，请在右侧工作台查看",
        )
