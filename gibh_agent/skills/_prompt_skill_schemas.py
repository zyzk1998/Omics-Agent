# -*- coding: utf-8 -*-
"""Prompt 软技能 · Pydantic 结构化输出模型（与前端 renderSkillResultVisually 对齐）。"""
from __future__ import annotations

from typing import List, Literal, Optional

from pydantic import BaseModel, Field, field_validator


class Slide(BaseModel):
    """单页幻灯片大纲。"""

    title: str = ""
    content_bullets: List[str] = Field(default_factory=list, alias="bullets")
    page: Optional[int] = None

    model_config = {"populate_by_name": True}

    @field_validator("title", mode="before")
    @classmethod
    def _coerce_title(cls, v: object) -> str:
        s = str(v or "").strip()
        return s if s else "未命名页面"

    @field_validator("content_bullets", mode="before")
    @classmethod
    def _coerce_bullets(cls, v: object) -> List[str]:
        if v is None:
            return []
        if isinstance(v, str):
            return [ln.strip() for ln in v.split("\n") if ln.strip()]
        if isinstance(v, list):
            return [str(x).strip() for x in v if str(x).strip()]
        return [str(v).strip()] if str(v).strip() else []


class PPTOutline(BaseModel):
    """科研汇报 PPT 大纲（序列化后供前端卡片渲染）。"""

    theme: str = Field(default="", alias="title")
    slides: List[Slide] = Field(min_length=1)
    total_pages: Optional[int] = None

    model_config = {"populate_by_name": True}

    def to_frontend_dict(self) -> dict:
        slides_out: List[dict] = []
        for idx, slide in enumerate(self.slides, start=1):
            bullets = list(slide.content_bullets or [])
            page_no = slide.page if slide.page is not None else idx
            slides_out.append(
                {
                    "page": int(page_no),
                    "title": (slide.title or f"第 {idx} 页").strip(),
                    "bullets": bullets[:12],
                }
            )
        total = self.total_pages if self.total_pages is not None else len(slides_out)
        return {
            "title": (self.theme or "").strip(),
            "total_pages": int(total),
            "slides": slides_out,
        }


class PromptSkillTerminalResponse(BaseModel):
    """周报/简历等软技能 · 终结态 JSON（禁止无限追问）。"""

    action: Literal["deliver", "missing_params"]
    markdown: str = ""
    missing_params: List[str] = Field(default_factory=list)
    summary: str = ""

    @field_validator("missing_params", mode="before")
    @classmethod
    def _coerce_missing(cls, v: object) -> List[str]:
        if v is None:
            return []
        if isinstance(v, str):
            return [v.strip()] if v.strip() else []
        if isinstance(v, list):
            return [str(x).strip() for x in v if str(x).strip()][:5]
        return []
