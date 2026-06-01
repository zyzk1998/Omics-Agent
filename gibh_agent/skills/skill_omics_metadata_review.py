# -*- coding: utf-8 -*-
"""组学样本 Metadata 规范审查 — 纯 Prompt → Markdown。"""
from __future__ import annotations

from typing import Any, Dict

from gibh_agent.skills._array1_batch2_system import OMICS_METADATA_REVIEW
from gibh_agent.skills._prompt_skill_markdown import execute_markdown_prompt_skill
from gibh_agent.skills.base_skill import BaseSkill


class OmicsMetadataReviewSkill(BaseSkill):
    __abstractskill__ = False

    skill_id = "omics_metadata_review"
    display_name = "组学样本 Metadata 规范审查"
    description = "审查样本表字段完整性、分组一致性与批次混杂风险。"
    category = "生物医药"
    sub_category = "文本处理"
    aliases = ['metadata', '样本表', '临床信息', '批次效应']
    required_parameters = ['样本 metadata 表或字段列表（粘贴至 context）']
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
            system_prompt=OMICS_METADATA_REVIEW,
            user_request=(user_request or "").strip(),
            context=(context or "").strip(),
            deliver_message="Metadata 审查报告已生成，请在右侧工作台查看",
        )
