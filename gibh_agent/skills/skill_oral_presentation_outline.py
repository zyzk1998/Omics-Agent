# -*- coding: utf-8 -*-
"""会议口头报告讲稿大纲 — Pydantic 结构化 JSON（复用 PPT 大纲组件）。"""
from __future__ import annotations

import logging
from typing import Any, Dict

from gibh_agent.skills._array1_batch2_system import ORAL_PRESENTATION_OUTLINE
from gibh_agent.skills._prompt_skill_completeness import PPT_AND_VISUAL_COMPLETENESS_RULES
from gibh_agent.skills._prompt_skill_llm import build_user_prompt, llm_chat_structured
from gibh_agent.skills._prompt_skill_schemas import PPTOutline
from gibh_agent.skills.base_skill import BaseSkill
from gibh_agent.skills.launch_skill_demos import apply_launch_demo_defaults

logger = logging.getLogger(__name__)

_SYSTEM = ORAL_PRESENTATION_OUTLINE + "\n\n" + PPT_AND_VISUAL_COMPLETENESS_RULES


class OralPresentationOutlineSkill(BaseSkill):
    __abstractskill__ = False

    skill_id = "oral_presentation_outline"
    display_name = "会议口头报告讲稿大纲"
    description = (
        "将论文/课题结构转为 10–15 分钟口述提纲（分页 JSON：标题、要点、建议时长），"
        "适合组会口头报告与答辩演练。"
    )
    category = "生物医药"
    sub_category = "文本处理"
    aliases = ["口头报告", "讲稿大纲", "答辩提纲", "组会口述", "oral presentation"]
    required_parameters = [
        "报告主题、听众类型与时长（写入 user_request）；论文结构或摘要可放 context",
    ]
    output_type = "json"
    __dependencies__: list[str] = []

    def execute(
        self,
        user_request: str = "",
        context: str = "",
        **kwargs: Any,
    ) -> Dict[str, Any]:
        filled = apply_launch_demo_defaults(
            self.skill_id,
            {
                "user_request": (user_request or "").strip(),
                "context": (context or "").strip(),
            },
        )
        try:
            user = build_user_prompt(
                str(filled.get("user_request") or ""),
                str(filled.get("context") or "") or None,
                allow_demo=True,
            )
        except ValueError as exc:
            return {"status": "error", "message": str(exc)}

        try:
            outline_model = llm_chat_structured(
                _SYSTEM,
                user,
                PPTOutline,
                max_tokens=6144,
            )
            outline = outline_model.to_frontend_dict()
        except Exception as exc:
            logger.exception("oral_presentation_outline LLM 失败")
            return {"status": "error", "message": f"生成讲稿大纲失败：{exc}"}

        return {
            "status": "success",
            "message": f"已生成 {outline['total_pages']} 页口头报告讲稿大纲",
            "data": {"schema": "oral_presentation_outline.v1", "ppt_outline": outline},
            "ppt_outline": outline,
            "total_pages": outline["total_pages"],
            "slides": outline["slides"],
            "title": outline.get("title") or "",
        }
