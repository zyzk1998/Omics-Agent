# -*- coding: utf-8 -*-
"""将 ToolRegistry 中的工具转为 OpenAI Chat Completions `tools` 格式。"""
from typing import Any, Dict, List, Optional

from gibh_agent.core.tool_registry import registry


def tool_names_to_openai_tools(names: List[str]) -> List[Dict[str, Any]]:
    out: List[Dict[str, Any]] = []
    for name in names:
        meta = registry.get_metadata(name)
        if not meta:
            continue
        try:
            schema = meta.args_schema.model_json_schema()
        except Exception:
            schema = {"type": "object", "properties": {}}
        out.append(
            {
                "type": "function",
                "function": {
                    "name": meta.name,
                    "description": (meta.description or "")[:4096],
                    "parameters": schema,
                },
            }
        )
    return out
