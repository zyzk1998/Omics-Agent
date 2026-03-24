# -*- coding: utf-8 -*-
"""
全网实时检索 MCP（与 tools/ 目录隔离）。

优先级：
1. Tavily（TAVILY_API_KEY 存在时优先；tavily-python 官方 SDK）
2. ddgs（零配置回退）
3. Serper / Bing（可选 API Key，见 .env.example）
"""
from __future__ import annotations

import json
import logging
import os
from typing import Any, Dict, List, Optional

import requests

from gibh_agent.core.tool_registry import registry

logger = logging.getLogger(__name__)


def _search_ddgs(query: str, max_results: int = 5) -> List[Dict[str, Any]]:
    """使用 ddgs 包（推荐，2025+ 维护中的 duckduckgo 搜索）。"""
    from ddgs import DDGS

    out: List[Dict[str, Any]] = []
    with DDGS() as ddgs:
        for r in ddgs.text(query, max_results=max_results):
            out.append(
                {
                    "title": (r.get("title") or "")[:500],
                    "href": (r.get("href") or r.get("url") or "")[:2000],
                    "body": (r.get("body") or "")[:1500],
                }
            )
    return out


def _search_tavily(query: str, max_results: int = 5) -> List[Dict[str, Any]]:
    """优先 Tavily 官方 SDK；未安装 tavily-python 时回退到 REST（仍使用 TAVILY_API_KEY）。"""
    key = os.getenv("TAVILY_API_KEY", "").strip()
    if not key:
        raise ValueError("TAVILY_API_KEY not set")
    data: Dict[str, Any] = {}
    try:
        from tavily import TavilyClient

        client = TavilyClient(api_key=key)
        raw = client.search(query, max_results=max_results, search_depth="basic")
        data = raw if isinstance(raw, dict) else {}
    except ImportError:
        resp = requests.post(
            "https://api.tavily.com/search",
            json={
                "api_key": key,
                "query": query,
                "max_results": max_results,
                "search_depth": "basic",
            },
            timeout=30,
        )
        resp.raise_for_status()
        data = resp.json()
    results = data.get("results") or []
    return [
        {
            "title": (x.get("title") or "")[:500],
            "href": (x.get("url") or "")[:2000],
            "body": (x.get("content") or "")[:1500],
        }
        for x in results[:max_results]
    ]


def _search_serper(query: str, max_results: int = 5) -> List[Dict[str, Any]]:
    key = os.getenv("SERPER_API_KEY", "").strip()
    if not key:
        raise ValueError("SERPER_API_KEY not set")
    resp = requests.post(
        "https://google.serper.dev/search",
        json={"q": query, "num": max_results},
        headers={"X-API-KEY": key},
        timeout=30,
    )
    resp.raise_for_status()
    data = resp.json()
    organic = data.get("organic") or []
    return [
        {
            "title": (x.get("title") or "")[:500],
            "href": (x.get("link") or "")[:2000],
            "body": (x.get("snippet") or "")[:1500],
        }
        for x in organic[:max_results]
    ]


def _search_bing(query: str, max_results: int = 5) -> List[Dict[str, Any]]:
    key = os.getenv("BING_SEARCH_API_KEY", "").strip()
    endpoint = os.getenv("BING_SEARCH_ENDPOINT", "").strip()
    if not key or not endpoint:
        raise ValueError("BING_SEARCH_API_KEY / BING_SEARCH_ENDPOINT not set")
    resp = requests.get(
        endpoint,
        headers={"Ocp-Apim-Subscription-Key": key},
        params={"q": query, "count": max_results, "mkt": "zh-CN"},
        timeout=30,
    )
    resp.raise_for_status()
    data = resp.json()
    web_pages = (data.get("webPages") or {}).get("value") or []
    return [
        {
            "title": (x.get("name") or "")[:500],
            "href": (x.get("url") or "")[:2000],
            "body": (x.get("snippet") or "")[:1500],
        }
        for x in web_pages[:max_results]
    ]


@registry.register(
    name="mcp_web_search",
    description=(
        "在互联网上搜索与查询相关的最新网页摘要与链接。"
        "当用户需要实时信息、最新新闻、近期论文或超出模型知识截止日期的内容时使用。"
        "返回若干条标题、链接与摘要，供你综合回答。"
    ),
    category="MCP",
    output_type="json",
)
def mcp_web_search(query: str, max_results: int = 5) -> Dict[str, Any]:
    """
    执行全网检索。

    Args:
        query: 搜索关键词或完整问句（建议用英文或中文关键词以提高命中率）
        max_results: 返回条数上限（默认 5）
    """
    max_results = max(1, min(int(max_results or 5), 10))
    q = (query or "").strip()
    if not q:
        return {"status": "error", "error": "query 不能为空", "results": []}

    errors: List[str] = []

    if os.getenv("TAVILY_API_KEY", "").strip():
        try:
            rows = _search_tavily(q, max_results=max_results)
            logger.info("✅ [mcp_web_search] 使用后端 tavily 返回 %s 条", len(rows))
            return {
                "status": "success",
                "backend": "tavily",
                "query": q,
                "results": rows,
            }
        except Exception as e:
            errors.append(f"tavily: {e}")
            logger.warning("⚠️ [mcp_web_search] Tavily 失败，将尝试 ddgs 等回退: %s", e)

    for label, fn in (
        ("ddgs", _search_ddgs),
        ("serper", _search_serper),
        ("bing", _search_bing),
    ):
        try:
            rows = fn(q, max_results=max_results)
            logger.info("✅ [mcp_web_search] 使用后端 %s 返回 %s 条", label, len(rows))
            return {
                "status": "success",
                "backend": label,
                "query": q,
                "results": rows,
            }
        except Exception as e:
            errors.append(f"{label}: {e}")
            logger.debug("[mcp_web_search] %s 失败: %s", label, e)

    logger.error("❌ [mcp_web_search] 全部后端失败: %s", errors)
    return {
        "status": "error",
        "error": "；".join(errors[-4:]),
        "results": [],
    }
