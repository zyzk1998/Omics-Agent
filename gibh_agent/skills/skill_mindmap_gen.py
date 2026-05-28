# -*- coding: utf-8 -*-
"""分析逻辑思维导图 — Prompt-Engineered 软技能（纯 LLM → Mermaid）。"""
from __future__ import annotations

import logging
from typing import Any, Dict

from gibh_agent.skills._prompt_skill_completeness import (
    MERMAID_QUALITY_RULES,
    PPT_AND_VISUAL_COMPLETENESS_RULES,
)
from gibh_agent.skills._prompt_skill_llm import (
    build_user_prompt,
    llm_chat_text,
    normalize_mermaid_code,
)
from gibh_agent.skills.base_skill import BaseSkill
from gibh_agent.skills.launch_skill_demos import apply_launch_demo_defaults

logger = logging.getLogger(__name__)

_MINDMAP_SYSTEM = """你是生物医学科研逻辑梳理专家。请从用户材料中提取核心逻辑树，输出用于 Mermaid.js 渲染的思维导图。

""" + PPT_AND_VISUAL_COMPLETENESS_RULES + "\n" + MERMAID_QUALITY_RULES + """

输出要求（严格遵守）：
1. 仅输出一段 Mermaid 源码正文：不要使用 Markdown 代码围栏（不要 ```），不要任何解释文字。
2. 优先使用 flowchart TD 或 graph TD；节点 ID 用英文字母或字母数字（如 A、B1），节点标签用中文放在方括号内，例如 A[研究背景] --> B[核心发现]。
3. 体现因果、并列与层级（根因→机制→表型→验证等），节点数建议 8–24 个，避免过深嵌套。
4. 示例形态（勿照抄内容，仅参考语法）：flowchart TD; R[主题]-->M1[分支一]; R-->M2[分支二]; M1-->L1[细节]"""


class MindmapGenSkill(BaseSkill):
    __abstractskill__ = False

    """
    将用户输入整理为 Mermaid 思维导图语法，供前端 Mermaid.js 渲染。

    参数:
        user_request: 希望导图涵盖的主题、逻辑链或关键词（必填）。
        context: 可选补充信息（如分析结论、基因列表）。
    """

    skill_id = "mindmap_gen"
    display_name = "分析逻辑思维导图生成"
    description = (
        "调用大模型提取分析逻辑树，输出 Mermaid.js 语法，"
        "在工作台渲染为可交互思维导图。"
    )
    category = "其他技能"
    sub_category = "文本处理"
    aliases = ["思维导图", "逻辑导图", "mindmap", "Mermaid导图"]
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
            raw = llm_chat_text(_MINDMAP_SYSTEM, user, max_tokens=4096, temperature=0.35)
            mermaid_code = normalize_mermaid_code(raw)
        except Exception as exc:
            logger.exception("mindmap_gen LLM 失败")
            return {
                "status": "error",
                "message": f"生成思维导图失败：{exc}",
            }

        return {
            "status": "success",
            "message": "思维导图已生成，见下方可视化",
            "mermaid_code": mermaid_code,
        }
