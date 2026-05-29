# -*- coding: utf-8 -*-
"""Prompt-Engineered 软技能共用 LLM 调用（非 BaseSkill，不注册 ToolRegistry）。"""
from __future__ import annotations

import json
import logging
import re
from typing import Any, Dict, Optional, Type, TypeVar

from pydantic import BaseModel, ValidationError

from gibh_agent.core.llm_client import LLMClientFactory
from gibh_agent.core.llm_cloud_providers import get_default_chat_model

logger = logging.getLogger(__name__)

TModel = TypeVar("TModel", bound=BaseModel)


def strip_code_fence(content: str) -> str:
    s = (content or "").strip()
    if not s.startswith("```"):
        return s
    lines = s.split("\n")
    if len(lines) >= 2 and lines[0].strip().startswith("```"):
        lines = lines[1:]
    if lines and lines[-1].strip() == "```":
        lines = lines[:-1]
    return "\n".join(lines).strip()


def llm_chat_text(system: str, user: str, *, max_tokens: int = 4096, temperature: float = 0.3) -> str:
    client = LLMClientFactory.create_for_model(get_default_chat_model())
    completion = client.chat(
        messages=[
            {"role": "system", "content": system},
            {"role": "user", "content": user},
        ],
        temperature=temperature,
        max_tokens=max_tokens,
    )
    choices = getattr(completion, "choices", None) or []
    if not choices:
        raise RuntimeError("LLM 返回空 choices")
    msg = choices[0].message
    text = (getattr(msg, "content", None) or "").strip()
    if not text:
        raise RuntimeError("LLM 正文为空")
    return strip_code_fence(text)


def llm_chat_structured(
    system: str,
    user: str,
    model_cls: Type[TModel],
    *,
    max_tokens: int = 4096,
    temperature: float = 0.2,
) -> TModel:
    """调用 LLM 并按 Pydantic 模型校验（Prompt 内嵌 JSON Schema）。"""
    schema = model_cls.model_json_schema()
    schema_hint = (
        "\n\n【JSON Schema 强制约束】\n"
        "仅输出一个 JSON 对象，须严格符合以下 Schema（字段名不可自造）：\n"
        f"{json.dumps(schema, ensure_ascii=False)}\n"
    )
    raw = llm_chat_json(
        system + schema_hint,
        user,
        max_tokens=max_tokens,
        temperature=temperature,
    )
    if "data" in raw and isinstance(raw["data"], dict):
        merged = {**raw, **raw["data"]}
        raw = merged
    try:
        return model_cls.model_validate(raw)
    except ValidationError as exc:
        raise ValueError(f"LLM 输出未通过 Schema 校验: {exc}") from exc


def llm_chat_json(
    system: str,
    user: str,
    *,
    max_tokens: int = 4096,
    temperature: float = 0.25,
) -> Dict[str, Any]:
    raw = llm_chat_text(system, user, max_tokens=max_tokens, temperature=temperature)
    try:
        parsed = json.loads(raw)
    except json.JSONDecodeError as exc:
        match = re.search(r"\{[\s\S]*\}", raw)
        if not match:
            raise ValueError(f"LLM 未返回合法 JSON: {exc}") from exc
        parsed = json.loads(match.group(0))
    if not isinstance(parsed, dict):
        raise ValueError("LLM JSON 根节点须为对象")
    return parsed


_MARKDOWN_FIELD_ALIASES = (
    "markdown",
    "markmarkdown",
    "mark_down",
    "markDown",
    "md",
    "content",
    "body",
    "report_markdown",
    "report",
    "text",
    "正文",
)


def _extract_markdown_from_mapping(raw: Dict[str, Any]) -> str:
    for key in _MARKDOWN_FIELD_ALIASES:
        val = raw.get(key)
        if isinstance(val, str) and val.strip():
            return val.strip()
    return ""


def normalize_terminal_response_payload(parsed: Dict[str, Any]) -> Dict[str, Any]:
    """将 LLM 终结态 JSON 规范为固定 markdown 字段（兼容错别字键名）。"""
    if not isinstance(parsed, dict):
        raise ValueError("终结态 JSON 根节点须为对象")
    root = dict(parsed)
    merged: Dict[str, Any] = dict(root)
    if isinstance(root.get("data"), dict):
        merged.update(root["data"])

    md = _extract_markdown_from_mapping(merged)
    if not md:
        for val in merged.values():
            if not isinstance(val, str):
                continue
            snippet = val.strip()
            if not snippet.startswith("{"):
                continue
            try:
                inner = json.loads(snippet)
            except json.JSONDecodeError:
                continue
            if isinstance(inner, dict):
                md = _extract_markdown_from_mapping(inner)
                if md:
                    break

    if md:
        merged["markdown"] = md

    action = str(merged.get("action") or merged.get("phase") or "").strip().lower()
    if action not in ("deliver", "missing_params"):
        if merged.get("markdown") and not merged.get("missing_params"):
            action = "deliver"
        else:
            action = "missing_params"
    merged["action"] = action
    return merged


def normalize_ppt_outline_payload(parsed: Dict[str, Any]) -> Dict[str, Any]:
    """将 LLM 输出规范为前端可渲染的 ppt_outline 结构。"""
    root = parsed
    if "ppt_outline" in root and isinstance(root["ppt_outline"], dict):
        root = root["ppt_outline"]
    slides_in = root.get("slides")
    if not isinstance(slides_in, list) or not slides_in:
        raise ValueError("缺少 slides 数组")

    slides: list[Dict[str, Any]] = []
    for idx, item in enumerate(slides_in, start=1):
        if not isinstance(item, dict):
            continue
        title = (
            item.get("title")
            or item.get("slide_title")
            or item.get("heading")
            or f"第 {idx} 页"
        )
        bullets = (
            item.get("bullets")
            or item.get("key_points")
            or item.get("points")
            or item.get("要点")
            or []
        )
        if isinstance(bullets, str):
            bullets = [b.strip() for b in bullets.split("\n") if b.strip()]
        if not isinstance(bullets, list):
            bullets = [str(bullets)]
        bullets = [str(b).strip() for b in bullets if str(b).strip()]
        page_no = item.get("page") or item.get("page_number") or idx
        slides.append(
            {
                "page": int(page_no) if str(page_no).isdigit() else idx,
                "title": str(title).strip(),
                "bullets": bullets[:12],
            }
        )
    if not slides:
        raise ValueError("slides 无有效页")

    total = root.get("total_pages") or root.get("page_count") or len(slides)
    deck_title = (root.get("title") or root.get("deck_title") or root.get("主题") or "").strip()
    return {
        "title": deck_title,
        "total_pages": int(total) if str(total).isdigit() else len(slides),
        "slides": slides,
    }


def normalize_mermaid_code(raw: str) -> str:
    code = strip_code_fence(raw or "")
    if code.lower().startswith("mermaid"):
        code = code.split("\n", 1)[-1].strip()
    code = code.strip()
    if not code:
        raise ValueError("Mermaid 代码为空")
    head = code.split(None, 1)[0].lower() if code.split() else ""
    allowed = (
        "graph",
        "flowchart",
        "sequencediagram",
        "classdiagram",
        "statediagram",
        "erdiagram",
        "journey",
        "gantt",
        "pie",
        "mindmap",
        "timeline",
        "gitgraph",
    )
    if not any(head.startswith(p) for p in allowed):
        raise ValueError(f"Mermaid 须以合法图类型开头，当前: {head or '(空)'}")
    return code


def build_user_prompt(
    user_request: str,
    context: Optional[str] = None,
    *,
    allow_demo: bool = False,
) -> str:
    req = (user_request or "").strip()
    ctx = (context or "").strip()
    if not req:
        if allow_demo:
            parts = [
                "【演示模式】用户未提供具体主题或材料。"
                "请生成高质量 Demo（生物医学/组会典型场景），次要细节可用 [待补充：…] 占位。"
            ]
            if ctx:
                parts.append(f"\n补充上下文：\n{ctx}")
            return "\n".join(parts)
        raise ValueError("请提供 user_request 或分析主题说明")
    parts = [f"用户需求：\n{req}"]
    if ctx:
        parts.append(f"\n补充上下文：\n{ctx}")
    return "\n".join(parts)
