# -*- coding: utf-8 -*-
"""伦理与知情同意要点清单 — 纯 Prompt → Markdown。"""
from __future__ import annotations

from typing import Any, Dict

from gibh_agent.skills._array1_batch2_system import ETHICS_CONSENT_CHECKLIST
from gibh_agent.skills._prompt_skill_markdown import execute_markdown_prompt_skill
from gibh_agent.skills.base_skill import BaseSkill


class EthicsConsentChecklistSkill(BaseSkill):
    __abstractskill__ = False

    skill_id = "ethics_consent_checklist"
    display_name = "伦理与知情同意要点清单"
    description = "按干预类型输出伦理审查 Q&A 与知情同意书章节 Checklist。"
    category = "生物医药"
    sub_category = "文本处理"
    aliases = ['伦理审查', '知情同意', 'IRB']
    required_parameters = ['研究类型、干预风险与人群特征']
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
            system_prompt=ETHICS_CONSENT_CHECKLIST,
            user_request=(user_request or "").strip(),
            context=(context or "").strip(),
            deliver_message="伦理与知情同意清单已生成，请在右侧工作台查看",
        )
