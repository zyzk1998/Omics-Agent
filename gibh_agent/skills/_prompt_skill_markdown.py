# -*- coding: utf-8 -*-
"""Prompt 软技能 · 直接 Markdown 交付（无 JSON 围栏，工作台 safeMarkedParse 渲染）。"""
from __future__ import annotations

import logging
import re
from typing import Any, Dict, Optional

from gibh_agent.skills._prompt_skill_completeness import COMPLETENESS_AND_PLACEHOLDER_RULES
from gibh_agent.skills._prompt_skill_llm import build_user_prompt, llm_chat_text
from gibh_agent.skills.launch_skill_demos import apply_launch_demo_defaults

logger = logging.getLogger(__name__)

_DIRECT_MARKDOWN_OUTPUT_RULES = """
【输出铁律 — 最高优先级，覆盖一切礼貌用语习惯】
1. **禁止**任何前言、寒暄、元话语（如「好的，我将为您…」「以下是…」「希望对您有帮助」）。
2. **禁止**用 Markdown 代码围栏（```）包裹全文；第一行必须是正文标题（`#` 或 `##`）。
3. 用户未提供实质输入时：文首用 `> **【Demo 演示】**` 引用块标注，并输出高质量生物医学范例；表内数字须视为示意，勿冒充用户真实结果。
4. 用户已粘贴表格/列表时：**仅解读给定数值**；禁止编造未出现的基因、p 值、FC、样本量。
5. 次要缺失信息用 `[待补充：…]` 占位，**禁止**「请补充」「待确认项」章节。
""" + COMPLETENESS_AND_PLACEHOLDER_RULES

_CHATTER_PREFIX_RE = re.compile(
    r"^(?:\s*(?:好的|当然|没问题|Sure|Certainly|Here is|以下是)[^\n#]*\n+)+",
    re.IGNORECASE | re.MULTILINE,
)


def strip_markdown_preamble(text: str) -> str:
    """去掉模型偶发的寒暄前缀，避免污染 marked 渲染。"""
    s = (text or "").strip()
    if not s:
        return s
    for _ in range(3):
        m = _CHATTER_PREFIX_RE.match(s)
        if not m:
            break
        s = s[m.end() :].lstrip()
    if s.startswith("```"):
        from gibh_agent.skills._prompt_skill_llm import strip_code_fence

        s = strip_code_fence(s)
    return s.strip()


def execute_markdown_prompt_skill(
    *,
    skill_id: str,
    system_prompt: str,
    user_request: str = "",
    context: str = "",
    deliver_message: str = "内容已生成，请在右侧工作台查看",
    max_tokens: int = 6144,
    temperature: float = 0.32,
) -> Dict[str, Any]:
    """
    纯 Prompt 软技能统一出口：返回 phase=deliver + markdown 供 SkillAgent 右栏渲染。
    """
    filled = apply_launch_demo_defaults(
        skill_id,
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

    full_system = (system_prompt.strip() + "\n\n" + _DIRECT_MARKDOWN_OUTPUT_RULES).strip()
    try:
        raw_md = llm_chat_text(
            full_system,
            user,
            max_tokens=max_tokens,
            temperature=temperature,
        )
        markdown = strip_markdown_preamble(raw_md)
        if not markdown:
            raise RuntimeError("LLM 返回空 Markdown")
    except Exception as exc:
        logger.exception("%s LLM 失败", skill_id)
        return {"status": "error", "message": f"生成失败：{exc}"}

    payload = {
        "schema": f"{skill_id}.markdown.v1",
        "markdown": markdown,
    }
    return {
        "status": "success",
        "message": deliver_message,
        "phase": "deliver",
        "markdown": markdown,
        "data": payload,
    }
