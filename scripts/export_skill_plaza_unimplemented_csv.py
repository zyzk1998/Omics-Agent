#!/usr/bin/env python3
"""从 seed_skills 导出技能广场中「未实现 [Skill_Route]」技能清单，按 A/B/C 分类写入 CSV。"""
from __future__ import annotations

import argparse
import csv
import json
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from gibh_agent.db.seed_skills import get_all_system_skills_list  # noqa: E402

TOOLS_JSON = ROOT / "gibh_agent/data/scienceone_toolchain_tools.json"
LEGACY_NAMES_JSON = ROOT / "gibh_agent/data/scienceone_toolchain_display_names.json"

PYTHON_FALLBACK_NOTE = (
    "未加载 ToolChain 元数据 JSON。说明：浏览器端 bundle 为 JavaScript，通常不含技能侧 Python 源码；"
    "请运行 scripts/refresh_scienceone_toolchain_names.py 生成 scienceone_toolchain_tools.json。"
)


def _normalize_display_key(name: str) -> str:
    """与 ToolChain 展示名比对时忽略大小写与空白差异。"""
    return re.sub(r"\s+", "", (name or "").strip().lower())


def load_scienceone_display_lookup() -> dict[str, str]:
    """norm_key -> ToolChain 官方展示名（canonical，保留大小写与内部空格）。"""
    if TOOLS_JSON.is_file():
        with TOOLS_JSON.open(encoding="utf-8") as f:
            data = json.load(f)
        out: dict[str, str] = {}
        for t in data.get("tools") or []:
            disp = (t.get("display_name") or "").strip()
            if not disp:
                continue
            k = _normalize_display_key(disp)
            if k not in out:
                out[k] = disp
        return out
    if LEGACY_NAMES_JSON.is_file():
        with LEGACY_NAMES_JSON.open(encoding="utf-8") as f:
            data = json.load(f)
        out = {}
        for block in ("tool_names_entries", "chemistry_tool_entries"):
            for item in data.get(block) or []:
                disp = (item.get("display") or "").strip()
                if not disp:
                    continue
                k = _normalize_display_key(disp)
                if k not in out:
                    out[k] = disp
        return out
    return {}


def load_tool_meta_maps() -> tuple[str, dict[str, dict]]:
    """(python_note, norm_display -> tool record)"""
    if not TOOLS_JSON.is_file():
        return PYTHON_FALLBACK_NOTE, {}
    with TOOLS_JSON.open(encoding="utf-8") as f:
        data = json.load(f)
    note = (data.get("python_source_availability_note") or "").strip() or PYTHON_FALLBACK_NOTE
    by_norm: dict[str, dict] = {}
    for t in data.get("tools") or []:
        disp = (t.get("display_name") or "").strip()
        if not disp:
            continue
        nk = _normalize_display_key(disp)
        if nk not in by_norm:
            by_norm[nk] = t
    return note, by_norm


def toolchain_display_for_seed_name(seed_name: str, lookup: dict[str, str]) -> tuple[str, bool]:
    """
    返回 (展示名, 是否在 ToolChain 配置中找到)。
    未找到时展示名=种子原名（多为多模态组学工作流等 ToolChain 无同名卡片）。
    """
    raw = (seed_name or "").strip()
    k = _normalize_display_key(raw)
    if k in lookup:
        return lookup[k], True
    return raw, False


def meta_row_for_display(canonical_display: str, meta_by_norm: dict[str, dict], py_note: str) -> dict[str, str]:
    t = meta_by_norm.get(_normalize_display_key(canonical_display))
    if not t:
        return {
            "ToolChain_Key": "",
            "ToolChain官方长描述": "",
            "配置来源": "",
            "ScienceToolType分支": "",
            "前端推断await函数": "",
            "Python源码可得性说明": py_note,
        }
    case = (t.get("science_tool_dispatch_case") or "").strip()
    fn = (t.get("inferred_await_function") or "").strip()
    src = (t.get("source_block") or "").strip()
    long_desc = (t.get("long_description") or "").replace("\n", " ").strip()
    return {
        "ToolChain_Key": (t.get("tool_chain_key") or "").strip(),
        "ToolChain官方长描述": long_desc,
        "配置来源": src,
        "ScienceToolType分支": case,
        "前端推断await函数": fn,
        "Python源码可得性说明": py_note,
    }


ROUTE_RE = re.compile(r"\[Skill_Route:\s*(\w+)\s*\]")

C_NAME_FRAGMENTS = (
    "Evo2",
    "ESM3",
    "ESMFold",
    "ESM-Variants",
    "AlphaFold2",
    "DiffSBDD",
    "DiffDock",
    "OligoFormer",
    "ProtGPT",
    "RFdiffusion",
    "BioGPT",
    "LigandMPNN",
    "抗体人源化",
    "抗体序列生成",
    "ADMET",
)
C_DESC_KEYWORDS = (
    "深度学习",
    "扩散模型",
    "语言模型",
    "蛋白语言模型",
    "基础模型",
    "生成模型",
)


def has_skill_route(prompt: str) -> bool:
    return bool(prompt and ROUTE_RE.search(prompt))


def classify(name: str, sub_category: str, description: str) -> tuple[str, str]:
    n = name or ""
    d = description or ""
    s = sub_category or ""

    if s == "信息检索":
        return "A", "公共数据库与轻量检索"
    if "数据库查询" in n:
        return "A", "公共数据库与轻量检索"
    if n in (
        "基因蛋白信息查询器",
        "蛋白质资料提取工具",
        "获取mRNA序列工具",
        "化学元素查询",
        "查询GSEA支持的数据库工具",
    ):
        return "A", "公共数据库与轻量检索"
    if n == "MSA search":
        return "A", "公共数据库与轻量检索"

    for frag in C_NAME_FRAGMENTS:
        if frag in n:
            return "C", "AI 深度学习大模型"
    for kw in C_DESC_KEYWORDS:
        if kw in d:
            return "C", "AI 深度学习大模型"

    return "B", "传统生信算法与统计分析"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        default=ROOT / "skill_plaza_unimplemented_ABC.csv",
        help="输出 CSV 路径（默认：仓库根目录）",
    )
    args = parser.parse_args()

    lookup = load_scienceone_display_lookup()
    py_note, meta_by_norm = load_tool_meta_maps()

    rows = []
    for skill in get_all_system_skills_list():
        name = (skill.get("name") or "").strip()
        pt = skill.get("prompt_template") or ""
        if has_skill_route(pt):
            continue
        display, matched = toolchain_display_for_seed_name(name, lookup)
        code, cat_label = classify(name, skill.get("sub_category") or "", skill.get("description") or "")
        note = "无 [Skill_Route] 直连工具（占位或通用编排）"
        if not matched:
            note += "；ToolChain 科学工具列表无同名展示项（多为组学工作流模板）"

        extra = meta_row_for_display(display, meta_by_norm, py_note)
        if not matched:
            extra = {
                "ToolChain_Key": "",
                "ToolChain官方长描述": "",
                "配置来源": "",
                "ScienceToolType分支": "",
                "前端推断await函数": "",
                "Python源码可得性说明": py_note,
            }

        rows.append(
            {
                "类别代码": code,
                "类别名称": cat_label,
                "技能名称": display,
                "种子名称": name if name != display else "",
                "主分类": skill.get("main_category") or "",
                "子分类": skill.get("sub_category") or "",
                "简介": (skill.get("description") or "").replace("\n", " ").strip(),
                "备注": note,
                **extra,
            }
        )

    order = {"A": 0, "B": 1, "C": 2}
    rows.sort(key=lambda r: (order.get(r["类别代码"], 9), r["技能名称"]))

    fieldnames = [
        "类别代码",
        "类别名称",
        "技能名称",
        "种子名称",
        "主分类",
        "子分类",
        "简介",
        "备注",
        "ToolChain_Key",
        "ToolChain官方长描述",
        "配置来源",
        "ScienceToolType分支",
        "前端推断await函数",
        "Python源码可得性说明",
    ]

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="", encoding="utf-8-sig") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        w.writerows(rows)

    print(f"Wrote {len(rows)} rows -> {args.output}")


if __name__ == "__main__":
    main()
