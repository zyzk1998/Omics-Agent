# -*- coding: utf-8 -*-
"""
为 SemanticRouter / 澄清模式动态注入「工具视野」：从 WorkflowRegistry 与 ToolRegistry
拉取真实已注册的 SOP 与工具摘要，避免模型「盲人摸象」编造能力。
"""
from __future__ import annotations

import logging
from collections import defaultdict
from typing import Dict, List

logger = logging.getLogger(__name__)

_MAX_SOP_DESC_CHARS = 140
_MAX_TOOL_IDS_PER_CATEGORY = 8
_MAX_CATEGORIES = 24
_MAX_TOTAL_CHARS = 2800


def _ensure_tools_loaded() -> None:
    try:
        from gibh_agent.core.tool_registry import registry as tool_registry

        if len(tool_registry._tools) > 0:
            return
        from gibh_agent.tools import load_all_tools

        load_all_tools()
    except Exception as e:
        logger.debug("capabilities_digest: 工具懒加载跳过: %s", e)


def build_capabilities_digest_for_llm() -> str:
    """
    生成注入路由/澄清 Prompt 的短文本（中文面向模型）。

    - SOP：WorkflowRegistry.list_workflows() 域名 + 描述截断
    - 工具：按 category 聚合 tool_id 节选，避免单条过长
    """
    lines: List[str] = []

    try:
        from gibh_agent.core.workflows import WorkflowRegistry

        wr = WorkflowRegistry()
        sops: Dict[str, str] = wr.list_workflows() or {}
        lines.append("【已注册 SOP 工作流（域名 → 摘要）】")
        for name in sorted(sops.keys()):
            desc = (sops.get(name) or "").strip().replace("\n", " ")
            if len(desc) > _MAX_SOP_DESC_CHARS:
                desc = desc[: _MAX_SOP_DESC_CHARS] + "…"
            lines.append(f"- {name}：{desc if desc else '（无描述）'}")
    except Exception as e:
        logger.warning("capabilities_digest: WorkflowRegistry 失败: %s", e)
        lines.append("【已注册 SOP 工作流】（暂不可用）")

    _ensure_tools_loaded()
    try:
        from gibh_agent.core.tool_registry import registry as tool_registry

        by_cat: Dict[str, List[str]] = defaultdict(list)
        for tool_id, meta in tool_registry._tools.items():
            cat = (getattr(meta, "category", None) or "General").strip() or "General"
            by_cat[cat].append(tool_id)

        lines.append("【已注册工具（按类别节选 tool_id，勿编造未列出之名）】")
        shown_cats = 0
        for cat in sorted(by_cat.keys()):
            if shown_cats >= _MAX_CATEGORIES:
                lines.append(f"- … 另有 {len(by_cat) - shown_cats} 个类别已省略")
                break
            ids = sorted(by_cat[cat])[:_MAX_TOOL_IDS_PER_CATEGORY]
            suffix = "…" if len(by_cat[cat]) > len(ids) else ""
            lines.append(f"- {cat}：{', '.join(ids)}{suffix}")
            shown_cats += 1
    except Exception as e:
        logger.warning("capabilities_digest: ToolRegistry 失败: %s", e)
        lines.append("【已注册工具】（暂不可用）")

    text = "\n".join(lines).strip()
    if len(text) > _MAX_TOTAL_CHARS:
        text = text[:_MAX_TOTAL_CHARS] + "\n…（已截断）"
    return text
