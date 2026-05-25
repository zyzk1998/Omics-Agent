#!/usr/bin/env python3
"""
从 ScienceOne ToolChain 前端 bundle 解析「能力引擎矩阵」工具元数据，写入
`gibh_agent/data/scienceone_toolchain_tools.json`。

2025+ bundle 变更：旧版 `toolNames` / `toolDescriptions` 已移除，改为：
  - `chemistryTools` / `materialTools`（i18n 内嵌 name + description）
  - `ScienceToolType` 枚举（e.AlphaFold2="AlphaFold2" 等，覆盖生信/材料/通用工具）

页面：https://www.scienceone.cn/toolchain/#/science/tool
"""
from __future__ import annotations

import argparse
import json
import re
import sys
import urllib.request
from pathlib import Path
from typing import Dict, List, Optional, Tuple

ROOT = Path(__file__).resolve().parents[1]
OUT_JSON = ROOT / "gibh_agent/data/scienceone_toolchain_tools.json"
LEGACY_JSON = ROOT / "gibh_agent/data/scienceone_toolchain_display_names.json"
BASE = "https://www.scienceone.cn"
INDEX_URL = f"{BASE}/toolchain/"

PYTHON_NOTE = (
    "ToolChain 浏览器端为打包后的 JavaScript；本解析未在 bundle 中发现随包分发的技能侧 .py 源码或内联 Python。"
    "各工具通常通过 HTTP 调用远端服务；Python 实现若在服务端则需另行对接口/文档或授权渠道获取。"
)

_NAMED_TOOL_RE = re.compile(
    r"([\w]+):\{name:\"((?:\\.|[^\"\\])*)\",description:\"((?:\\.|[^\"\\])*)\"",
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


def parse_nested_tools_block(js: str, label: str) -> List[dict[str, str]]:
    """解析 chemistryTools / materialTools 等 {key:{name,description}} 块。"""
    needle = f"{label}:{{"
    start = js.find(needle)
    if start < 0:
        return []
    i = start + len(needle) - 1
    body = _extract_balanced(js, i)
    if body is None:
        return []
    out: List[dict[str, str]] = []
    for k, name, desc in _NAMED_TOOL_RE.findall(body):
        out.append(
            {
                "key": k,
                "display": name.replace('\\"', '"').strip(),
                "description": desc.replace('\\"', '"').strip(),
            }
        )
    return out


def parse_all_named_tool_blocks(js: str) -> Dict[str, Tuple[str, str]]:
    """全文件扫描 name+description（去重，保留首条）。"""
    out: Dict[str, Tuple[str, str]] = {}
    for k, name, desc in _NAMED_TOOL_RE.findall(js):
        if k not in out:
            out[k] = (
                name.replace('\\"', '"').strip(),
                desc.replace('\\"', '"').strip(),
            )
    return out


def parse_science_tool_type_enum(js: str) -> List[Tuple[str, str]]:
    """
    解析 ScienceToolType 枚举：e.EnumName="chainKey"。
    返回 [(enum_case, tool_chain_key), ...]
    """
    anchor = js.find('e.ESMVariants="esmVariants"')
    if anchor < 0:
        anchor = js.find('e.AlphaFold2="AlphaFold2"')
    if anchor < 0:
        return []
    window = js[max(0, anchor - 300) : anchor + 14000]
    pairs = re.findall(r"e\.([A-Za-z][\w]*)=\"([a-zA-Z][\w]*)\"", window)
    # 过滤 UI/路由类非工具枚举
    skip = frozenset(
        {
            "Web",
            "three3D",
            "Code",
            "Process",
            "Bar",
            "Running",
            "Failed",
            "Completed",
            "success",
            "failed",
            "running",
            "pending",
            "Science",
            "Explore",
            "Deep",
            "Rag",
            "Database",
            "ScienceClaw",
            "display",
            "none",
            "workflow",
            "Iframe",
            "Form",
            "Self",
            "Common",
            "Collection",
        }
    )
    out: List[Tuple[str, str]] = []
    seen: set = set()
    for enum_case, chain_key in pairs:
        if enum_case in skip or chain_key in skip:
            continue
        if chain_key in seen:
            continue
        seen.add(chain_key)
        out.append((enum_case, chain_key))
    return out


def parse_science_tool_dispatch(js: str) -> Dict[str, str]:
    """ScienceToolType switch 分支内首个 await 函数名（若 bundle 仍含该结构）。"""
    parts = re.split(r"case ScienceToolType\.(\w+):", js)
    out: Dict[str, str] = {}
    for i in range(1, len(parts), 2):
        case = parts[i]
        chunk = parts[i + 1][:2000]
        m = re.search(r"await\s+([a-zA-Z_$][\w$]*)\s*\(", chunk)
        out[case] = m.group(1) if m else ""
    return out


def _dispatch_case_for_key(tool_key: str) -> str:
    if not tool_key:
        return tool_key
    if tool_key[0].islower():
        return tool_key[0].upper() + tool_key[1:]
    return tool_key


def _seed_enrichment_maps() -> Tuple[Dict[str, str], Dict[str, str], Dict[str, str]]:
    """chain_key / display_name -> description；display -> chain_key。"""
    try:
        sys.path.insert(0, str(ROOT))
        from gibh_agent.db.panshi_skill_meta import PANSHI_TOOL_CHAIN_KEY_BY_DISPLAY_NAME
        from gibh_agent.db.seed_skills import get_all_system_skills_list
    except Exception:
        return {}, {}, {}

    desc_by_key: Dict[str, str] = {}
    desc_by_display: Dict[str, str] = {}
    display_by_key: Dict[str, str] = {}
    for skill in get_all_system_skills_list():
        name = (skill.get("name") or "").strip()
        desc = (skill.get("description") or "").strip()
        if name and desc:
            desc_by_display[name] = desc
        key = PANSHI_TOOL_CHAIN_KEY_BY_DISPLAY_NAME.get(name, "")
        if key:
            display_by_key[key] = name
            if desc:
                desc_by_key[key] = desc
    return desc_by_key, desc_by_display, display_by_key


def _humanize_enum(enum_case: str) -> str:
    s = re.sub(r"([a-z])([A-Z])", r"\1 \2", enum_case)
    return s.replace("_", " ").strip()


def build_tools(js: str, entry_url: str) -> List[dict[str, object]]:
    dispatch = parse_science_tool_dispatch(js)
    named_global = parse_all_named_tool_blocks(js)
    seed_desc_by_key, seed_desc_by_display, seed_display_by_key = _seed_enrichment_maps()

    merged: Dict[str, dict[str, object]] = {}

    def upsert(
        *,
        tool_chain_key: str,
        display_name: str,
        long_description: str,
        source_block: str,
        enum_case: str = "",
    ) -> None:
        k = (tool_chain_key or "").strip()
        if not k:
            return
        case = enum_case or (_dispatch_case_for_key(k) if k in dispatch else _dispatch_case_for_key(k))
        fn = dispatch.get(case, "") if case in dispatch else ""
        row = {
            "tool_chain_key": k,
            "display_name": (display_name or k).strip(),
            "long_description": (long_description or "").strip(),
            "source_block": source_block,
            "science_tool_dispatch_case": case if case in dispatch else "",
            "inferred_await_function": fn,
        }
        prev = merged.get(k)
        if prev is None:
            merged[k] = row
            return
        # 优先保留更长描述与中文展示名
        if len(row["long_description"]) > len(str(prev.get("long_description") or "")):
            prev["long_description"] = row["long_description"]
        if row["display_name"] and not str(prev.get("display_name", "")).isascii():
            prev["display_name"] = row["display_name"]
        elif prev.get("display_name") == prev.get("tool_chain_key") and row["display_name"]:
            prev["display_name"] = row["display_name"]
        if not prev.get("science_tool_dispatch_case") and row["science_tool_dispatch_case"]:
            prev["science_tool_dispatch_case"] = row["science_tool_dispatch_case"]

    for block_label, source in (
        ("materialTools", "materialTools"),
        ("chemistryTools", "chemistryTools"),
    ):
        for row in parse_nested_tools_block(js, block_label):
            upsert(
                tool_chain_key=row["key"],
                display_name=row["display"],
                long_description=row.get("description") or "",
                source_block=source,
            )

    for enum_case, chain_key in parse_science_tool_type_enum(js):
        named = named_global.get(chain_key) or named_global.get(enum_case)
        if named:
            disp, desc = named
        else:
            # 尝试 camelCase / PascalCase 变体
            alt_keys = [
                chain_key,
                enum_case,
                chain_key[0].lower() + chain_key[1:] if chain_key else "",
            ]
            named = None
            for ak in alt_keys:
                if ak in named_global:
                    named = named_global[ak]
                    break
            disp = named[0] if named else (
                seed_display_by_key.get(chain_key)
                or seed_display_by_key.get(enum_case)
                or _humanize_enum(enum_case)
            )
            desc = (
                (named[1] if named else "")
                or seed_desc_by_key.get(chain_key, "")
                or seed_desc_by_key.get(enum_case, "")
                or seed_desc_by_display.get(disp, "")
            )
        upsert(
            tool_chain_key=chain_key,
            display_name=disp,
            long_description=desc,
            source_block="scienceToolTypeEnum",
            enum_case=enum_case,
        )

    tools = _dedupe_case_variant_keys(list(merged.values()))
    tools.sort(key=lambda x: (x.get("source_block", ""), x.get("display_name", "")))
    return tools


_SOURCE_PRIORITY = {"chemistryTools": 0, "materialTools": 1, "scienceToolTypeEnum": 2}


def _dedupe_case_variant_keys(rows: List[dict[str, object]]) -> List[dict[str, object]]:
    """合并仅大小写/Pascal 差异的 tool_chain_key（保留 chemistryTools 的 camelCase）。"""
    buckets: Dict[str, List[dict[str, object]]] = {}
    for row in rows:
        k = str(row.get("tool_chain_key") or "")
        bucket_id = re.sub(r"[^a-z0-9]", "", k.lower())
        if not bucket_id:
            continue
        buckets.setdefault(bucket_id, []).append(row)

    out: List[dict[str, object]] = []
    for variants in buckets.values():
        if len(variants) == 1:
            out.append(variants[0])
            continue
        variants.sort(
            key=lambda r: (
                _SOURCE_PRIORITY.get(str(r.get("source_block") or ""), 9),
                0 if str(r.get("tool_chain_key", ""))[:1].islower() else 1,
                -len(str(r.get("long_description") or "")),
            )
        )
        primary = dict(variants[0])
        for other in variants[1:]:
            if len(str(other.get("long_description") or "")) > len(
                str(primary.get("long_description") or "")
            ):
                primary["long_description"] = other["long_description"]
            if str(other.get("display_name", "")) and not str(primary.get("display_name", "")).isascii():
                primary["display_name"] = other["display_name"]
            alt = str(other.get("tool_chain_key") or "")
            if alt and alt != primary.get("tool_chain_key"):
                aliases = primary.setdefault("alias_keys", [])
                if isinstance(aliases, list) and alt not in aliases:
                    aliases.append(alt)
        out.append(primary)
    return out


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
        bundle_cache = ROOT / ".tmp_s1entry.js"
        if bundle_cache.is_file() and bundle_cache.stat().st_size > 100_000:
            js = bundle_cache.read_text(encoding="utf-8", errors="replace")
            entry_url = f"file://{bundle_cache.resolve()}"
        else:
            html = _fetch(INDEX_URL)
            entry_url = resolve_entry_js_url(html)
            js = _fetch(entry_url)
            bundle_cache.write_text(js, encoding="utf-8")

    tools = build_tools(js, entry_url)
    by_src: Dict[str, int] = {}
    for t in tools:
        by_src[str(t.get("source_block") or "?")] = by_src.get(str(t.get("source_block") or "?"), 0) + 1
    payload = {
        "source_page": "https://www.scienceone.cn/toolchain/#/science/tool",
        "fetched_entry_js": entry_url,
        "bundle_parser_version": "2.0_scienceToolType_enum",
        "python_source_availability_note": PYTHON_NOTE,
        "catalog_coverage_note": (
            "本 JSON 仅含 ToolChain 前端 entry bundle 内嵌的 chemistryTools/materialTools/ScienceToolType 枚举，"
            "当前约 120 条；不含需登录的 /agent/skill/page 或 /agent/tools/tool_universe/ 技能广场全量目录。"
            "官网宣传的 2000+ 工具与生物医药/化学 400+ 技能，多数在服务端分页接口，未随 SPA 静态包下发。"
        ),
        "tools_parsed_total": len(tools),
        "tools_by_source_block": by_src,
        "skill_plaza_api_hint": (
            "可选：设置环境变量 SCIENCEONE_COOKIE（浏览器登录 scienceone.cn 后的 Cookie 串），"
            "运行 scripts/fetch_scienceone_skill_plaza_api.py 合并技能广场分页数据。"
        ),
        "tools": tools,
    }
    OUT_JSON.parent.mkdir(parents=True, exist_ok=True)
    with OUT_JSON.open("w", encoding="utf-8") as f:
        json.dump(payload, f, ensure_ascii=False, indent=2)

    legacy = {
        "source_page": payload["source_page"],
        "fetched_entry_js": entry_url,
        "science_tool_type_entries": [
            {"key": t["tool_chain_key"], "display": t["display_name"]}
            for t in tools
            if t.get("source_block") == "scienceToolTypeEnum"
        ],
        "chemistry_tool_entries": [
            {"key": t["tool_chain_key"], "display": t["display_name"]}
            for t in tools
            if t.get("source_block") == "chemistryTools"
        ],
        "material_tool_entries": [
            {"key": t["tool_chain_key"], "display": t["display_name"]}
            for t in tools
            if t.get("source_block") == "materialTools"
        ],
    }
    with LEGACY_JSON.open("w", encoding="utf-8") as f:
        json.dump(legacy, f, ensure_ascii=False, indent=2)

    by_src = {}
    for t in tools:
        by_src[t.get("source_block", "?")] = by_src.get(t.get("source_block", "?"), 0) + 1
    print(f"Wrote {OUT_JSON} ({len(tools)} tools), legacy {LEGACY_JSON}")
    print("by source_block:", by_src)


if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        print(f"错误: {e}", file=sys.stderr)
        sys.exit(1)
