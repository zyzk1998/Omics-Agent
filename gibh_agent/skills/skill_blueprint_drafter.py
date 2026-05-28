# -*- coding: utf-8 -*-
"""工程蓝图 HTML 制图 — Prompt-Engineered 软技能（扁平工程示意图）。"""
from __future__ import annotations

import logging
import re
from typing import Any, Dict

from gibh_agent.skills._prompt_skill_llm import build_user_prompt, llm_chat_text, strip_code_fence
from gibh_agent.skills._prompt_skill_spec import load_prompt_spec, run_prompt_spec_skill
from gibh_agent.skills.base_skill import BaseSkill
from gibh_agent.skills.launch_skill_demos import apply_launch_demo_defaults

logger = logging.getLogger(__name__)

_SPEC_ID = "blueprint_drafter"
_HTML_SUFFIX = (
    "\n【输出附加约束】仅返回完整 HTML 文档正文，以 <!DOCTYPE html> 开头，"
    "禁止 Markdown 代码围栏与解释性前后文。"
)


def _normalize_blueprint_html(raw: str) -> str:
    text = strip_code_fence(raw or "").strip()
    if text.lower().startswith("html"):
        text = text.split("\n", 1)[-1].strip()
    match = re.search(r"(<!DOCTYPE[\s\S]*</html>)", text, re.IGNORECASE)
    if match:
        return match.group(1).strip()
    if "<html" in text.lower():
        if not text.lower().startswith("<!doctype"):
            text = "<!DOCTYPE html>\n" + text
        return text
    raise ValueError("LLM 未返回完整 HTML 文档（需含 <!DOCTYPE html> 或 <html>）")


class BlueprintDrafterSkill(BaseSkill):
    __abstractskill__ = False

    skill_id = "blueprint_drafter"
    display_name = "工程蓝图制图"
    description = (
        "生成扁平化、高数据墨水比的工程/架构示意图 HTML（系统字体、无 CDN），"
        "适合技术规格说明与流程蓝图。"
    )
    category = "其他技能"
    sub_category = "数据可视化"
    aliases = ["drafter", "蓝图", "架构图", "工程示意图", "blueprint"]
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
                "user_request": (user_request or kwargs.get("diagram_topic") or "").strip(),
                "context": (context or "").strip(),
            },
        )
        spec = load_prompt_spec(_SPEC_ID)
        if not spec:
            return run_prompt_spec_skill(_SPEC_ID, user_request=str(filled.get("user_request") or ""))

        try:
            user = build_user_prompt(
                str(filled.get("user_request") or ""),
                str(filled.get("context") or "") or None,
            )
        except ValueError as exc:
            return {"status": "error", "message": str(exc)}

        system = spec[:100_000] + _HTML_SUFFIX
        try:
            raw = llm_chat_text(system, user, max_tokens=8192, temperature=0.25)
            html = _normalize_blueprint_html(raw)
        except Exception as exc:
            logger.exception("blueprint_drafter LLM 失败")
            return {"status": "error", "message": f"生成工程蓝图失败：{exc}"}

        title_m = re.search(r"<title>([^<]*)</title>", html, re.IGNORECASE)
        title = (title_m.group(1).strip() if title_m else "") or "工程蓝图"

        return {
            "status": "success",
            "message": f"已生成工程蓝图 HTML：{title}",
            "html_content": html,
            "markdown": f"### 工程蓝图已生成\n\n**{title}** — 见下方 HTML 预览。",
        }
