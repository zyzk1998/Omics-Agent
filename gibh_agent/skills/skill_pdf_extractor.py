# -*- coding: utf-8 -*-
"""PDF 内容提取 — Prompt-Engineered 软技能（可选 pypdf 读 file_path）。"""
from __future__ import annotations

from typing import Any, Dict

from gibh_agent.skills._prompt_skill_terminal import execute_terminal_spec_from_kwargs
from gibh_agent.skills.base_skill import BaseSkill
from gibh_agent.skills.launch_skill_demos import apply_launch_demo_defaults

_SPEC_ID = "pdf_extractor"


class PdfExtractorSkill(BaseSkill):
    __abstractskill__ = False

    skill_id = "pdf_extractor"
    display_name = "PDF 内容提取助手"
    description = (
        "从 PDF 提取文本与结构摘要，输出 Markdown 大纲与表格说明；"
        "提供 file_path 时尝试本地解析（需 pypdf）。"
    )
    category = "其他技能"
    sub_category = "文本处理"
    aliases = ["PDF", "pdf提取", "pdfplumber", "文档提取"]
    __dependencies__: list[str] = ["pip:pypdf"]
    required_parameters = []

    def execute(
        self,
        user_request: str = "",
        context: str = "",
        file_path: str = "",
        **kwargs: Any,
    ) -> Dict[str, Any]:
        filled = apply_launch_demo_defaults(
            self.skill_id,
            {
                "user_request": (user_request or "").strip(),
                "context": (context or "").strip(),
                "file_path": (file_path or kwargs.get("pdf_path") or "").strip(),
            },
        )
        extra = ""
        pages = (kwargs.get("pages") or "").strip()
        if pages:
            extra = f"仅关注 PDF 页码范围：{pages}"
        return execute_terminal_spec_from_kwargs(
            _SPEC_ID,
            filled,
            {**kwargs, "extra_system": extra} if extra else kwargs,
            max_tokens=8192,
        )
