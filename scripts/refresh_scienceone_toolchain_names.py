#!/usr/bin/env python3
"""
从 ScienceOne ToolChain 前端 bundle 解析工具元数据（展示名、长描述、部分 ScienceToolType 分支的首个 await 函数名等），
写入 `gibh_agent/data/scienceone_toolchain_tools.json`。

说明：bundle 为打包后的 JavaScript，不包含各技能服务端 Python 源码；可用于对齐官方文案与定位前端 API 包装函数名。

页面：https://www.scienceone.cn/toolchain/#/science/tool
"""
from __future__ import annotations

import argparse
import json
import re
import sys
import urllib.request
from pathlib import Path
from typing import Optional

ROOT = Path(__file__).resolve().parents[1]
OUT_JSON = ROOT / "gibh_agent/data/scienceone_toolchain_tools.json"
LEGACY_JSON = ROOT / "gibh_agent/data/scienceone_toolchain_display_names.json"
BASE = "https://www.scienceone.cn"
INDEX_URL = f"{BASE}/toolchain/"

PYTHON_NOTE = (
    "ToolChain 浏览器端为打包后的 JavaScript；本解析未在 bundle 中发现随包分发的技能侧 .py 源码或内联 Python。"
    "各工具通常通过 HTTP 调用远端服务；Python 实现若在服务端则需另行对接口/文档或授权渠道获取。"
)


def _fetch(url: str, timeout: int = 60) -> str:
    req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0 GIBH-Agent-ToolChainMeta/1.0"})
    with urllib.request.urlopen(req, timeout=timeout) as resp:
        return resp.read().decode("utf-8", errors="replace")


def resolve_entry_js_url(html: str) -> str:
    m = re.search(r'src="(/toolchain/chat/assets/entry/index[a-z0-9]+\.js)"', html)
    if not m:
        raise RuntimeError("无法在 toolchain 首页解析 entry JS 路径")
    return BASE + m.group(1)


def _extract_balanced(js: str, start_brace_index: int) -> Optional[str]:
    if start_brace_index < 0 or start_brace_index >= len(js) or js[start_brace_index] != "{":
        return None
    depth = 0
    j = start_brace_index
    while j < len(js):
        if js[j] == "{":
            depth += 1
        elif js[j] == "}":
            depth -= 1
            if depth == 0:
                return js[start_brace_index + 1 : j]
        j += 1
    return None


def parse_toolnames(js: str) -> list[dict[str, str]]:
    start = js.find("toolNames:{")
    if start < 0:
        return []
    i = start + len("toolNames:{") - 1
    body = _extract_balanced(js, i)
    if body is None:
        return []
    pairs = re.findall(r'([A-Za-z_][\w]*)\s*:\s*"((?:\\.|[^"\\])*)"', body)
    return [{"key": k, "display": v.replace('\\"', '"').strip()} for k, v in pairs]


def parse_tool_descriptions(js: str) -> dict[str, str]:
    start = js.find("toolDescriptions:{")
    if start < 0:
        return {}
    i = start + len("toolDescriptions:{") - 1
    body = _extract_balanced(js, i)
    if body is None:
        return {}
    pairs = re.findall(r'([A-Za-z][\w]*)\s*:\s*"((?:\\.|[^"\\])*)"', body)
    return {k: v.replace('\\"', '"').strip() for k, v in pairs}


def parse_chemistry_tools(js: str) -> list[dict[str, str]]:
    label = "chemistryTools:{"
    start = js.find(label)
    if start < 0:
        return []
    i = start + len(label) - 1
    body = _extract_balanced(js, i)
    if body is None:
        return []
    pat = r'(\w+):\{name:"((?:\\.|[^"\\])*)",description:"((?:\\.|[^"\\])*)"'
    out = []
    for k, name, desc in re.findall(pat, body):
        out.append(
            {
                "key": k,
                "display": name.replace('\\"', '"').strip(),
                "description": desc.replace('\\"', '"').strip(),
            }
        )
    return out


def parse_science_tool_dispatch(js: str) -> dict[str, str]:
    """ScienceToolType.Case 分支内首个 await foo( 的 foo；部分分支为同步逻辑则可能为空。"""
    parts = re.split(r"case ScienceToolType\.(\w+):", js)
    out: dict[str, str] = {}
    for i in range(1, len(parts), 2):
        case = parts[i]
        chunk = parts[i + 1][:2000]
        m = re.search(r"await\s+([a-zA-Z_$][\w$]*)\s*\(", chunk)
        out[case] = m.group(1) if m else ""
    return out


def _dispatch_case_for_key(tool_key: str) -> str:
    """chemistry 多为 camelCase，与 switch 中 PascalCase 对齐。"""
    if not tool_key:
        return tool_key
    if tool_key[0].islower():
        return tool_key[0].upper() + tool_key[1:]
    return tool_key


def build_tools(js: str, entry_url: str) -> list[dict[str, object]]:
    desc_map = parse_tool_descriptions(js)
    dispatch = parse_science_tool_dispatch(js)
    tools: list[dict[str, object]] = []

    for row in parse_toolnames(js):
        k = row["key"]
        d = row["display"]
        long_desc = desc_map.get(k, "")
        case = k if k in dispatch else _dispatch_case_for_key(k)
        fn = dispatch.get(case, "") if case in dispatch else ""
        tools.append(
            {
                "tool_chain_key": k,
                "display_name": d,
                "long_description": long_desc,
                "source_block": "toolNames",
                "science_tool_dispatch_case": case if case in dispatch else "",
                "inferred_await_function": fn,
            }
        )

    for row in parse_chemistry_tools(js):
        k = row["key"]
        d = row["display"]
        long_desc = row.get("description") or ""
        case = _dispatch_case_for_key(k)
        fn = dispatch.get(case, "") if case in dispatch else ""
        tools.append(
            {
                "tool_chain_key": k,
                "display_name": d,
                "long_description": long_desc,
                "source_block": "chemistryTools",
                "science_tool_dispatch_case": case if case in dispatch else "",
                "inferred_await_function": fn,
            }
        )

    return tools


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--js-bundle",
        type=Path,
        help="离线：直接指定已下载的 ToolChain entry index*.js，跳过网络",
    )
    args = ap.parse_args()

    if args.js_bundle:
        js = args.js_bundle.read_text(encoding="utf-8", errors="replace")
        entry_url = f"file://{args.js_bundle.resolve()}"
    else:
        html = _fetch(INDEX_URL)
        entry_url = resolve_entry_js_url(html)
        js = _fetch(entry_url)

    tools = build_tools(js, entry_url)
    payload = {
        "source_page": "https://www.scienceone.cn/toolchain/#/science/tool",
        "fetched_entry_js": entry_url,
        "python_source_availability_note": PYTHON_NOTE,
        "tools": tools,
    }
    OUT_JSON.parent.mkdir(parents=True, exist_ok=True)
    with OUT_JSON.open("w", encoding="utf-8") as f:
        json.dump(payload, f, ensure_ascii=False, indent=2)

    # 兼容旧导出脚本：保留扁平列表结构
    legacy = {
        "source_page": payload["source_page"],
        "fetched_entry_js": entry_url,
        "tool_names_entries": [
            {"key": t["tool_chain_key"], "display": t["display_name"]}
            for t in tools
            if t["source_block"] == "toolNames"
        ],
        "chemistry_tool_entries": [
            {"key": t["tool_chain_key"], "display": t["display_name"]}
            for t in tools
            if t["source_block"] == "chemistryTools"
        ],
    }
    with LEGACY_JSON.open("w", encoding="utf-8") as f:
        json.dump(legacy, f, ensure_ascii=False, indent=2)

    print(f"Wrote {OUT_JSON} ({len(tools)} tools), legacy {LEGACY_JSON}")


if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        print(f"错误: {e}", file=sys.stderr)
        sys.exit(1)
