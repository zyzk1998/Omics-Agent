#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Skill Factory — CSV 工具清单联网校验。

对 CSV 每行工具名称与描述，通过 Tavily（或 ddgs 回退）检索后由 LLM 给出修正建议，
输出 verified_tools.csv。

用法:
  PYTHONPATH=. python scripts/validate_tool_list.py \\
    --csv docs/生物信息分析工具.csv --limit 20

环境变量:
  TAVILY_API_KEY        优先联网搜索
  SKILL_FACTORY_MODEL   LLM 校验模型，默认 deepseek-chat
"""
from __future__ import annotations

import argparse
import json
import logging
import sys
import time
from pathlib import Path
from typing import Any, Dict, List, Optional

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from scripts.skill_factory_common import (  # noqa: E402
    DEFAULT_FACTORY_MODEL,
    factory_llm_client,
    parse_csv_row,
)

logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")

VALIDATION_SYSTEM_PROMPT = """你是生物信息学与科研软件事实核查专家。

用户会提供：
1) CSV 中的一条工具记录；
2) 联网搜索摘要（可能不完整）。

请输出**仅一个 JSON 对象**（不要 Markdown 围栏），字段：
{
  "tool_name_verified": "核实后的标准名称",
  "description_verified": "修正后的中文简介（1-3句）",
  "latest_version_or_ref": "官方文档/版本信息，未知则写 unknown",
  "parameters_hint": "典型输入参数要点（路径类参数请用 file_path / input_dir / sequence_or_path 等规范名）",
  "confidence": "high|medium|low",
  "issues": "若原描述有误或工具不存在，简要说明；否则空字符串",
  "sources": ["url1", "url2"]
}
"""


def web_search_snippets(query: str, max_results: int = 5) -> List[Dict[str, str]]:
    """Tavily → ddgs 回退（与 mcp_web_search 同源逻辑）。"""
    from gibh_agent.mcp.web_search_mcp import _search_ddgs, _search_tavily
    import os

    q = (query or "").strip()
    if not q:
        return []
    if os.getenv("TAVILY_API_KEY", "").strip():
        try:
            return _search_tavily(q, max_results=max_results)
        except Exception as e:
            logger.warning("Tavily 失败，回退 ddgs: %s", e)
    try:
        return _search_ddgs(q, max_results=max_results)
    except Exception as e:
        logger.warning("ddgs 搜索失败: %s", e)
        return []


def format_search_context(hits: List[Dict[str, str]]) -> str:
    if not hits:
        return "（无联网检索结果，仅基于模型知识，confidence 应标为 low）"
    lines = []
    for i, h in enumerate(hits, 1):
        lines.append(
            f"[{i}] {h.get('title', '')}\n"
            f"    URL: {h.get('href', '')}\n"
            f"    {h.get('body', '')[:800]}"
        )
    return "\n".join(lines)


def validate_row_with_llm(
    client,
    parsed: Dict[str, str],
    search_hits: List[Dict[str, str]],
    *,
    model: str,
) -> Dict[str, Any]:
    from gibh_agent.core.llm_json_extract import extract_json_object_from_llm_text

    user_msg = (
        f"工具记录:\n"
        f"  名称: {parsed['display_name']}\n"
        f"  学科: {parsed['discipline']}\n"
        f"  子分类: {parsed['sub_category']}\n"
        f"  类型: {parsed['skills_type']}\n"
        f"  tool_chain_key: {parsed['tool_chain_key']}\n"
        f"  原简介: {parsed['description']}\n\n"
        f"联网检索摘要:\n{format_search_context(search_hits)}"
    )
    resp = client.chat(
        messages=[
            {"role": "system", "content": VALIDATION_SYSTEM_PROMPT},
            {"role": "user", "content": user_msg},
        ],
        model=model,
        temperature=0.1,
        max_tokens=2048,
    )
    text = (resp.choices[0].message.content or "").strip()
    obj = extract_json_object_from_llm_text(text)
    if not isinstance(obj, dict):
        raise ValueError(f"无法解析校验 JSON: {text[:200]}")
    return obj


def main() -> int:
    parser = argparse.ArgumentParser(description="校验生物信息工具 CSV（联网 + LLM）")
    parser.add_argument("--csv", type=Path, default=ROOT / "docs" / "生物信息分析工具.csv")
    parser.add_argument(
        "--out",
        type=Path,
        default=ROOT / "docs" / "verified_tools.csv",
    )
    parser.add_argument("--model", default=DEFAULT_FACTORY_MODEL)
    parser.add_argument("--start", type=int, default=0)
    parser.add_argument("--limit", type=int, default=0)
    parser.add_argument("--sleep", type=float, default=1.5)
    parser.add_argument("--search-max", type=int, default=5)
    args = parser.parse_args()

    if not args.csv.is_file():
        logger.error("CSV 不存在: %s", args.csv)
        return 1

    df = pd.read_csv(args.csv, encoding="utf-8-sig").fillna("")
    end = len(df) if args.limit <= 0 else min(len(df), args.start + args.limit)
    df_slice = df.iloc[args.start : end]

    client = factory_llm_client(args.model)
    out_rows: List[Dict[str, Any]] = []

    for idx, row in df_slice.iterrows():
        parsed = parse_csv_row(row.to_dict())
        name = parsed["display_name"]
        if not name:
            continue
        query = f"{name} bioinformatics tool official documentation parameters"
        if parsed["tool_chain_key"]:
            query += f" {parsed['tool_chain_key']}"

        logger.info("校验 [%s] %s …", idx, name)
        hits = web_search_snippets(query, max_results=args.search_max)

        try:
            verdict = validate_row_with_llm(client, parsed, hits, model=args.model)
        except Exception as e:
            logger.error("  LLM 校验失败: %s", e)
            verdict = {
                "tool_name_verified": name,
                "description_verified": parsed["description"],
                "latest_version_or_ref": "unknown",
                "parameters_hint": "",
                "confidence": "low",
                "issues": str(e),
                "sources": [],
            }

        out_rows.append(
            {
                **row.to_dict(),
                "verified_name": verdict.get("tool_name_verified", name),
                "verified_description": verdict.get("description_verified", ""),
                "latest_version_or_ref": verdict.get("latest_version_or_ref", ""),
                "parameters_hint": verdict.get("parameters_hint", ""),
                "validation_confidence": verdict.get("confidence", "low"),
                "validation_issues": verdict.get("issues", ""),
                "validation_sources": json.dumps(verdict.get("sources") or [], ensure_ascii=False),
            }
        )

        if args.sleep > 0:
            time.sleep(args.sleep)

    out_df = pd.DataFrame(out_rows)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    out_df.to_csv(args.out, index=False, encoding="utf-8-sig")
    logger.info("已写入 %s (%d 行)", args.out, len(out_df))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
