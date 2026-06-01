# -*- coding: utf-8 -*-
"""阵列二轻量技能 · 结构化返回载荷（表格 / 指标卡 / Markdown）。"""
from __future__ import annotations

import re
from typing import Any, Dict, List, Optional, Sequence

EMPTY_DEFAULT_MESSAGE = (
    "未能根据您提供的条件检索到结果，请检查拼写或尝试其他关键词。"
)

_TIMEOUT_USER_MESSAGE = (
    "外部数据库请求超时，请稍后重试或缩小检索范围。"
)


def table_data_from_rows(
    columns: Sequence[str],
    rows: Sequence[Dict[str, Any]],
) -> Dict[str, Any]:
    return {"columns": list(columns), "rows": list(rows)}


def metrics_cards_from_pairs(pairs: Sequence[tuple[str, Any, Optional[str]]]) -> List[Dict[str, Any]]:
    """[(label, value, unit_optional), ...]"""
    out: List[Dict[str, Any]] = []
    for item in pairs:
        if len(item) < 2:
            continue
        label, value = item[0], item[1]
        unit = item[2] if len(item) > 2 else ""
        card: Dict[str, Any] = {"label": str(label), "value": value}
        if unit:
            card["unit"] = str(unit)
        out.append(card)
    return out


def _scalar_key_from_label(label: str) -> str:
    s = re.sub(r"[^\w\u4e00-\u9fff]+", "_", (label or "").strip().lower())
    s = re.sub(r"_+", "_", s).strip("_")
    return s or "metric"


def empty_result(message: Optional[str] = None, **data: Any) -> Dict[str, Any]:
    payload = dict(data)
    msg = (message or EMPTY_DEFAULT_MESSAGE).strip()
    return {
        "status": "empty",
        "message": msg,
        "markdown": f"> {msg}\n",
        "data": payload,
    }


def error_result(message: str, **data: Any) -> Dict[str, Any]:
    out: Dict[str, Any] = {"status": "error", "message": message}
    if data:
        out["data"] = dict(data)
    return out


def timeout_result(**data: Any) -> Dict[str, Any]:
    return error_result(_TIMEOUT_USER_MESSAGE, **data)


def success_payload(
    message: str,
    *,
    markdown: str = "",
    hits: Optional[List[Dict[str, Any]]] = None,
    table_data: Optional[Dict[str, Any]] = None,
    metrics_cards: Optional[List[Dict[str, Any]]] = None,
    **scalars: Any,
) -> Dict[str, Any]:
    """
    组装强类型载荷：顶层 + data 双写，供 SkillAgent 右栏 _skillVis* 解析。

    - hits：对象数组 → 科研表格
    - metrics_cards + 扁平 scalars → Metrics Board
    - table_data：columns/rows 契约（与组学工具一致）
    """
    data: Dict[str, Any] = dict(scalars)
    if hits is not None:
        data["hits"] = hits
    if table_data is not None:
        data["table_data"] = table_data
    if metrics_cards:
        data["metrics_cards"] = metrics_cards
        for card in metrics_cards:
            label = str(card.get("label") or "")
            if not label:
                continue
            key = _scalar_key_from_label(label)
            data[key] = card.get("value")

    out: Dict[str, Any] = {
        "status": "success",
        "message": message,
        "data": data,
    }
    if markdown:
        out["markdown"] = markdown
    if hits is not None:
        out["hits"] = hits
    if table_data is not None:
        out["table_data"] = table_data
    if metrics_cards:
        out["metrics_cards"] = metrics_cards
    return out
