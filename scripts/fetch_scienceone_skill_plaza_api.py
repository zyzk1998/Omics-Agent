#!/usr/bin/env python3
"""
从磐石登录态 API 分页拉取技能广场 / 工具宇宙目录，合并进 scienceone_toolchain_tools.json。

官网 ToolChain entry bundle 仅内嵌约 120 条算子元数据；生物医药+化学 400+ 技能通常在：
  - GET/POST https://www.scienceone.cn/agent/skill/page
  - GET/POST https://www.scienceone.cn/agent/tools/tool_universe/page

需浏览器登录后复制 Cookie 到环境变量 SCIENCEONE_COOKIE（勿提交 Git）。

用法（仓库根目录）:
    export SCIENCEONE_COOKIE='session=...; other=...'
    PYTHONPATH=. python3 scripts/fetch_scienceone_skill_plaza_api.py
    PYTHONPATH=. python3 scripts/fetch_scienceone_skill_plaza_api.py --dry-run
"""
from __future__ import annotations

import argparse
import json
import os
import ssl
import sys
import urllib.error
import urllib.request
from pathlib import Path
from typing import Any, Dict, List, Optional

ROOT = Path(__file__).resolve().parents[1]
OUT_JSON = ROOT / "gibh_agent/data/scienceone_toolchain_tools.json"
BASE = "https://www.scienceone.cn/agent"

ENDPOINTS = (
    "/skill/page",
    "/tools/tool_universe/page",
)


def _request(
    path: str,
    *,
    cookie: str,
    body: Optional[dict] = None,
    timeout: int = 60,
) -> Any:
    url = BASE + path
    data = json.dumps(body or {"page": 1, "pageSize": 200}).encode("utf-8")
    headers = {
        "User-Agent": "GIBH-Agent-ScienceOneSkillPlaza/1.0",
        "Content-Type": "application/json",
        "Cookie": cookie,
    }
    req = urllib.request.Request(url, data=data, headers=headers, method="POST")
    ctx = ssl.create_default_context()
    with urllib.request.urlopen(req, timeout=timeout, context=ctx) as resp:
        raw = resp.read().decode("utf-8", errors="replace")
    payload = json.loads(raw)
    if isinstance(payload, dict) and payload.get("status") not in (None, 1, 200, "success"):
        raise RuntimeError(f"{path} 非成功响应: {payload!r:.500}")
    return payload


def _extract_rows(payload: Any) -> List[dict]:
    """兼容常见分页字段：data.records / data.list / data.items / data。"""
    if not isinstance(payload, dict):
        return []
    data = payload.get("data")
    if isinstance(data, list):
        return [r for r in data if isinstance(r, dict)]
    if not isinstance(data, dict):
        return []
    for key in ("records", "list", "items", "rows", "content"):
        block = data.get(key)
        if isinstance(block, list):
            return [r for r in block if isinstance(r, dict)]
    return []


def _row_to_tool(row: dict, *, source_api: str) -> Optional[dict]:
    key = (
        row.get("tool_chain_key")
        or row.get("toolChainKey")
        or row.get("chainKey")
        or row.get("key")
        or row.get("code")
        or row.get("id")
    )
    name = (
        row.get("display_name")
        or row.get("displayName")
        or row.get("name")
        or row.get("title")
        or row.get("skillName")
    )
    desc = (
        row.get("long_description")
        or row.get("description")
        or row.get("intro")
        or row.get("summary")
        or ""
    )
    if not key and not name:
        return None
    return {
        "tool_chain_key": str(key or name).strip(),
        "display_name": str(name or key).strip(),
        "long_description": str(desc).strip(),
        "source_block": source_api,
        "science_tool_dispatch_case": "",
        "inferred_await_function": "",
    }


def merge_tools(existing: List[dict], incoming: List[dict]) -> List[dict]:
    by_key: Dict[str, dict] = {}
    for row in existing:
        k = str(row.get("tool_chain_key") or "").strip()
        if k:
            by_key[k] = row
    for row in incoming:
        k = str(row.get("tool_chain_key") or "").strip()
        if not k:
            continue
        prev = by_key.get(k)
        if prev is None:
            by_key[k] = row
            continue
        if len(str(row.get("long_description") or "")) > len(
            str(prev.get("long_description") or "")
        ):
            prev["long_description"] = row["long_description"]
        if row.get("display_name") and not str(prev.get("display_name", "")).isascii():
            prev["display_name"] = row["display_name"]
        prev.setdefault("alias_sources", [])
        src = row.get("source_block")
        aliases = prev["alias_sources"]
        if isinstance(aliases, list) and src and src not in aliases:
            aliases.append(src)
    return sorted(by_key.values(), key=lambda x: (x.get("source_block", ""), x.get("display_name", "")))


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--dry-run", action="store_true", help="只打印各接口条数，不写 JSON")
    ap.add_argument("--cookie", default=os.environ.get("SCIENCEONE_COOKIE", ""), help="登录 Cookie")
    args = ap.parse_args()
    cookie = (args.cookie or "").strip()
    if not cookie:
        print(
            "未设置 SCIENCEONE_COOKIE：请登录 https://www.scienceone.cn 后在开发者工具复制 Cookie。",
            file=sys.stderr,
        )
        sys.exit(2)

    fetched: List[dict] = []
    for path in ENDPOINTS:
        try:
            payload = _request(path, cookie=cookie)
        except urllib.error.HTTPError as e:
            print(f"{path} HTTP {e.code}: {e.read(300)!r}", file=sys.stderr)
            continue
        except Exception as e:
            print(f"{path} 失败: {e}", file=sys.stderr)
            continue
        rows = _extract_rows(payload)
        tools = [t for r in rows if (t := _row_to_tool(r, source_api=path.lstrip("/")))]
        print(f"{path}: API 行 {len(rows)} → 可映射工具 {len(tools)}")
        fetched.extend(tools)

    if args.dry_run:
        print(f"合计可合并 {len(fetched)} 条（未写盘）")
        return

    if not OUT_JSON.is_file():
        print(f"缺少 {OUT_JSON}，请先运行 refresh_scienceone_toolchain_names.py", file=sys.stderr)
        sys.exit(1)

    payload = json.loads(OUT_JSON.read_text(encoding="utf-8"))
    existing = payload.get("tools") if isinstance(payload.get("tools"), list) else []
    merged = merge_tools(existing, fetched)
    payload["tools"] = merged
    payload["skill_plaza_api_merged"] = True
    payload["skill_plaza_api_fetched_count"] = len(fetched)
    payload["tools_parsed_total"] = len(merged)
    OUT_JSON.write_text(json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8")
    print(f"Wrote {OUT_JSON} (merged total {len(merged)}, +{len(fetched)} from API)")


if __name__ == "__main__":
    main()
