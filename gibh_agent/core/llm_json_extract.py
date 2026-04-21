# -*- coding: utf-8 -*-
"""
从 LLM 原始文本中稳健提取 JSON 对象（兼容 Kimi 等模型附加 Markdown 围栏、前后缀说明）。
供 SemanticRouter、RouterAgent、Orchestrator 意图分类等共用。
"""
from __future__ import annotations

import json
import logging
import re
from typing import Any, Dict, Optional

logger = logging.getLogger(__name__)

_FENCE_PATTERN = re.compile(r"```(?:json)?\s*([\s\S]*?)\s*```", re.IGNORECASE)


def _balanced_brace_slice(s: str) -> Optional[str]:
    """从首个 `{` 起做括号深度扫描，截取最外层平衡 JSON 对象子串。"""
    start = s.find("{")
    if start < 0:
        return None
    depth = 0
    in_string = False
    escape = False
    quote = ""
    i = start
    n = len(s)
    while i < n:
        ch = s[i]
        if in_string:
            if escape:
                escape = False
            elif ch == "\\":
                escape = True
            elif ch == quote:
                in_string = False
            i += 1
            continue
        if ch in ('"', "'"):
            in_string = True
            quote = ch
            i += 1
            continue
        if ch == "{":
            depth += 1
        elif ch == "}":
            depth -= 1
            if depth == 0:
                return s[start : i + 1]
        i += 1
    return None


def _try_parse_dict(blob: str) -> Optional[Dict[str, Any]]:
    if not blob or not blob.strip():
        return None
    s = blob.strip()
    for candidate in (s, _balanced_brace_slice(s) or ""):
        if not candidate:
            continue
        try:
            obj = json.loads(candidate)
            if isinstance(obj, dict):
                return obj
        except json.JSONDecodeError:
            continue
    return None


def extract_json_object_from_llm_text(text: str) -> Optional[Dict[str, Any]]:
    """
    从模型输出中解析出单个 JSON 对象（dict）。

    处理：Markdown ```json 围栏、前后自然语言、跨行 JSON。
    """
    raw = (text or "").strip()
    if not raw:
        return None

    # 1) 优先尝试围栏内文本
    m = _FENCE_PATTERN.search(raw)
    if m:
        inner = m.group(1).strip()
        got = _try_parse_dict(inner)
        if got is not None:
            return got

    # 2) 全文 + 平衡括号切片
    got = _try_parse_dict(raw)
    if got is not None:
        return got

    # 3) 弱兜底：首 `{` 到末 `}`
    start, end = raw.find("{"), raw.rfind("}")
    if start != -1 and end > start:
        try:
            obj = json.loads(raw[start : end + 1])
            if isinstance(obj, dict):
                return obj
        except json.JSONDecodeError:
            logger.debug("llm_json_extract: fallback slice failed prefix=%r", raw[:160])

    return None
