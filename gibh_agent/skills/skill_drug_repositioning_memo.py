# -*- coding: utf-8 -*-
"""药物重定位假说备忘录 — 纯 Prompt → Markdown。"""
from __future__ import annotations

from typing import Any, Dict

from gibh_agent.skills._array1_batch2_system import DRUG_REPOSITIONING_MEMO
from gibh_agent.skills._prompt_skill_markdown import execute_markdown_prompt_skill
from gibh_agent.skills.base_skill import BaseSkill


class DrugRepositioningMemoSkill(BaseSkill):
    __abstractskill__ = False

    skill_id = "drug_repositioning_memo"
    display_name = "药物重定位假说备忘录"
    description = "疾病—靶点—药物文字→可验证假说与文献检索关键词（不调用 OpenTargets）。"
    category = "生物医药"
    sub_category = "文本处理"
    aliases = ['药物重定位', '老药新用', '靶点假说']
    required_parameters = ['疾病、靶点、候选药物与已知证据（文字）']
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
            system_prompt=DRUG_REPOSITIONING_MEMO,
            user_request=(user_request or "").strip(),
            context=(context or "").strip(),
            deliver_message="药物重定位备忘录已生成，请在右侧工作台查看",
        )
