# -*- coding: utf-8 -*-
"""
MCP tools/list 的 inputSchema → OpenAI Chat Completions tools[] 的 parameters 适配。

MCP 与 OpenAI 均使用 JSON Schema 子集；本模块做保守归一化，避免空 schema 或缺失 type 导致上游拒收。
"""
from __future__ import annotations

import re
from typing import Any, Dict, List, Set, Tuple


def registry_name_from_mcp_tool_name(mcp_tool_name: str, used: Set[str]) -> str:
    """与 gibh_agent.mcp_client._registry_name 一致，保证动态注入与 connect() 注册名对齐。"""
    base = re.sub(r"[^a-zA-Z0-9_]+", "_", mcp_tool_name).strip("_") or "tool"
    name = f"hpc_mcp_{base}"
    if name not in used:
        used.add(name)
        return name
    i = 2
    while f"{name}_{i}" in used:
        i += 1
    cand = f"{name}_{i}"
    used.add(cand)
    return cand


def normalize_openai_parameters(input_schema: Any) -> Dict[str, Any]:
    """
    将 MCP inputSchema 转为可作为 OpenAI function.parameters 的对象。
    - 已是 dict 且含 type/properties → 浅拷贝并补全缺省字段
    - 空或非法 → {"type": "object", "properties": {}}
    """
    if not isinstance(input_schema, dict):
        return {"type": "object", "properties": {}}

    # 已是完整 JSON Schema 对象
    t = input_schema.get("type")
    if t == "object" or ("properties" in input_schema):
        out = dict(input_schema)
        out.setdefault("type", "object")
        if "properties" not in out or not isinstance(out["properties"], dict):
            out["properties"] = {}
        req = out.get("required")
        if req is not None and not isinstance(req, list):
            out.pop("required", None)
        return out

    return {"type": "object", "properties": {}}


def mcp_tool_record_to_openai_function(
    name: str,
    description: str,
    input_schema: Any,
    openai_function_name: str,
) -> Dict[str, Any]:
    """单条 MCP 工具 → OpenAI tools 数组元素。"""
    params = normalize_openai_parameters(input_schema)
    desc = (description or "").strip() or f"超算 MCP 工具 `{name}`（经网关代理）。"
    if len(desc) > 4096:
        desc = desc[:4093] + "..."
    return {
        "type": "function",
        "function": {
            "name": openai_function_name,
            "description": desc,
            "parameters": params,
        },
    }


def mcp_tools_payload_to_openai_tools(
    tools: List[Dict[str, Any]],
) -> Tuple[List[Dict[str, Any]], Dict[str, str]]:
    """
    :param tools: 自网关 /tools 的 items，需含 name；可选 description、inputSchema
    :return: (openai_tools, registry_name -> 远端 MCP tool name)
    """
    used: Set[str] = set()
    openai: List[Dict[str, Any]] = []
    reg_to_mcp: Dict[str, str] = {}

    for item in tools:
        if not isinstance(item, dict):
            continue
        mcp_name = (item.get("name") or "").strip()
        if not mcp_name:
            continue
        desc = (item.get("description") or item.get("title") or "") or ""
        schema = item.get("inputSchema")
        if schema is None:
            schema = item.get("input_schema")
        reg = registry_name_from_mcp_tool_name(mcp_name, used)
        reg_to_mcp[reg] = mcp_name
        openai.append(
            mcp_tool_record_to_openai_function(
                mcp_name,
                str(desc),
                schema,
                reg,
            )
        )

    return openai, reg_to_mcp
