# -*- coding: utf-8 -*-
"""科研汇报 PPT 大纲 — Pydantic 强制结构化 JSON。"""
from __future__ import annotations

import logging
from typing import Any, Dict

from gibh_agent.skills._prompt_skill_completeness import PPT_AND_VISUAL_COMPLETENESS_RULES
from gibh_agent.skills._prompt_skill_llm import build_user_prompt, llm_chat_structured
from gibh_agent.skills._prompt_skill_schemas import PPTOutline
from gibh_agent.skills.base_skill import BaseSkill
from gibh_agent.skills.launch_skill_demos import apply_launch_demo_defaults

logger = logging.getLogger(__name__)

_PPT_SYSTEM = """你是顶级生物医学科研咨询顾问，擅长将复杂研究整理为组会/答辩级幻灯片大纲。

""" + PPT_AND_VISUAL_COMPLETENESS_RULES + """

输出 JSON 须符合 Schema：
- theme：汇报总标题
- slides：数组，每项含 title、content_bullets（3–6 条要点字符串）
- total_pages：等于 slides 长度

内容须覆盖：封面/研究背景/科学问题或假设/方法要点（若有）/核心发现（可分多页）/结论与展望。
缺失细节在 content_bullets 中用 `[待补充：…]` 占位，保持页结构完整。
禁止 Markdown 围栏与任何非 JSON 文字。"""


class PptOutlineSkill(BaseSkill):
    __abstractskill__ = False

    skill_id = "ppt_outline"
    display_name = "科研汇报 PPT 大纲生成"
    description = (
        "调用大模型将研究主题与结论整理为分页 PPT 大纲（JSON：页数、标题、要点），"
        "适合组会汇报与幻灯片制作。"
    )
    category = "其他技能"
    sub_category = "文本处理"
    aliases = ["PPT大纲", "汇报大纲", "幻灯片大纲", "presentation outline"]
    required_parameters = []
    tool_chain_key = ""
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
                "user_request": (user_request or kwargs.get("analysis_request") or "").strip(),
                "context": (context or kwargs.get("context") or "").strip(),
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
                _PPT_SYSTEM,
                user,
                PPTOutline,
                max_tokens=6144,
            )
            outline = outline_model.to_frontend_dict()
        except Exception as exc:
            logger.exception("ppt_outline LLM 失败")
            return {
                "status": "error",
                "message": f"生成 PPT 大纲失败：{exc}",
            }

        data_payload = {
            "schema": "ppt_outline.v1",
            "ppt_outline": outline,
        }
        return {
            "status": "success",
            "message": f"已生成 {outline['total_pages']} 页科研汇报大纲",
            "data": data_payload,
            "ppt_outline": outline,
            "total_pages": outline["total_pages"],
            "slides": outline["slides"],
            "title": outline.get("title") or "",
        }
