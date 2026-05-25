#!/usr/bin/env python3
"""
解析用户从磐石技能广场复制的清单文本 → 结构化 JSON。

输入格式（块之间空行分隔）：
  - 首行：「生物医药：技能名」或「化学：技能名」
  - 或首行「模型」+ 下一行技能名
  - 或首行即为技能名
  - 中间为简介
  - 末行可选：支持对话调用 / 接入中... / 模型

输出：gibh_agent/data/scienceone_skill_plaza_catalog.json

用法:
  PYTHONPATH=. python3 scripts/parse_scienceone_skill_plaza_paste.py
  PYTHONPATH=. python3 scripts/parse_scienceone_skill_plaza_paste.py --input docs/raw/foo.txt
"""
from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_INPUT = ROOT / "docs/raw/scienceone_skill_plaza_user_paste.txt"
OUT_JSON = ROOT / "gibh_agent/data/scienceone_skill_plaza_catalog.json"
TOOLCHAIN_JSON = ROOT / "gibh_agent/data/scienceone_toolchain_tools.json"

CATEGORY_PREFIX_RE = re.compile(r"^(生物医药|化学)[：:]\s*(.+)$")
CATEGORY_ONLY_RE = re.compile(r"^(生物医药|化学)[：:]\s*$")
STATUS_LINES = frozenset(
    {
        "支持对话调用",
        "接入中...",
        "接入中",
        "模型",
    }
)
STANDALONE_MODEL_LINE = "模型"


def _split_blocks(text: str) -> List[str]:
    blocks: List[str] = []
    cur: List[str] = []
    for line in text.replace("\r\n", "\n").split("\n"):
        if not line.strip():
            if cur:
                blocks.append("\n".join(cur).strip())
                cur = []
            continue
        cur.append(line.rstrip())
    if cur:
        blocks.append("\n".join(cur).strip())
    return [b for b in blocks if b]


def _parse_block(block: str, default_category: str) -> Tuple[Optional[Dict[str, Any]], str]:
    """返回 (entry, updated_default_category)。"""
    lines = [ln.strip() for ln in block.split("\n") if ln.strip()]
    if not lines:
        return None, default_category

    category = default_category
    name = ""
    desc_lines: List[str] = []
    skills_type = ""
    invoke_status = ""

    idx = 0
    first = lines[0]
    if CATEGORY_ONLY_RE.match(first):
        return None, first.replace("：", "").replace(":", "").strip()
    m = CATEGORY_PREFIX_RE.match(first)
    if m:
        category = m.group(1)
        name = m.group(2).strip()
        idx = 1
    elif first == STANDALONE_MODEL_LINE and len(lines) > 1:
        skills_type = "模型"
        name = lines[1].strip()
        idx = 2
    else:
        name = first
        idx = 1

    while idx < len(lines):
        ln = lines[idx]
        if ln in STATUS_LINES:
            if ln == STANDALONE_MODEL_LINE:
                skills_type = "模型"
            elif ln.startswith("接入"):
                invoke_status = ln
            else:
                invoke_status = ln
            idx += 1
            continue
        desc_lines.append(ln)
        idx += 1

    if not name or name in ("我手动复制给你",):
        return None, category

    description = "\n".join(desc_lines).strip()
    if not skills_type:
        if invoke_status == "支持对话调用" and re.search(
            r"数据库|Browser|查询|检索|PubMed|UniProt|GWAS|GEO|dbSNP|ClinVar|Reactome|InterPro|UCSC",
            name + description,
            re.I,
        ):
            skills_type = "数据库调用"
        elif re.search(
            r"Fold|Dock|GPT|Diff|AlphaFold|ESM|Evo2|BepiPred|RFdiffusion|BioGPT|GenMol|LigandMPNN|"
            r"OligoFormer|ProtGPT|基础模型|扩散模型|语言模型|预训练",
            name + description,
            re.I,
        ):
            skills_type = "模型"
        else:
            skills_type = "脚本"

    entry: Dict[str, Any] = {
        "display_name": name,
        "long_description": description,
        "main_category": category,
        "skills_type": skills_type,
        "invoke_status": invoke_status,
        "source_block": "skill_plaza_user_paste",
    }
    return entry, category


def _load_toolchain_key_index() -> Dict[str, str]:
    """display_name (norm) -> tool_chain_key。"""
    if not TOOLCHAIN_JSON.is_file():
        return {}
    try:
        payload = json.loads(TOOLCHAIN_JSON.read_text(encoding="utf-8"))
    except (json.JSONDecodeError, OSError):
        return {}
    tools = payload.get("tools") or []
    out: Dict[str, str] = {}
    for t in tools:
        if not isinstance(t, dict):
            continue
        key = (t.get("tool_chain_key") or "").strip()
        disp = (t.get("display_name") or "").strip()
        if key and disp:
            out[_norm(disp)] = key
            out[_norm(key)] = key
    try:
        sys.path.insert(0, str(ROOT))
        from gibh_agent.db.panshi_skill_meta import PANSHI_TOOL_CHAIN_KEY_BY_DISPLAY_NAME

        for disp, key in PANSHI_TOOL_CHAIN_KEY_BY_DISPLAY_NAME.items():
            if disp and key:
                out[_norm(disp)] = key
    except Exception:
        pass
    return out


def _norm(s: str) -> str:
    return re.sub(r"\s+", "", (s or "").strip().lower())


def _guess_tool_chain_key(name: str, key_index: Dict[str, str]) -> str:
    n = _norm(name)
    if n in key_index:
        return key_index[n]
    # 英文名 → camelCase 猜测
    if re.search(r"[A-Za-z]", name):
        parts = re.findall(r"[A-Za-z0-9]+", name)
        if parts:
            return parts[0].lower() + "".join(p.title() for p in parts[1:])
    return ""


def parse_catalog(text: str) -> List[Dict[str, Any]]:
    key_index = _load_toolchain_key_index()
    default_cat = "生物医药"
    entries: List[Dict[str, Any]] = []
    for block in _split_blocks(text):
        row, default_cat = _parse_block(block, default_cat)
        if not row:
            continue
        key = _guess_tool_chain_key(row["display_name"], key_index)
        row["tool_chain_key"] = key
        entries.append(row)
    return entries


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    ap.add_argument("--out", type=Path, default=OUT_JSON)
    args = ap.parse_args()
    if not args.input.is_file():
        print(f"缺少输入文件: {args.input}", file=sys.stderr)
        sys.exit(1)
    text = args.input.read_text(encoding="utf-8")
    entries = parse_catalog(text)
    by_cat: Dict[str, int] = {}
    for e in entries:
        by_cat[e.get("main_category", "?")] = by_cat.get(e.get("main_category", "?"), 0) + 1
    payload = {
        "source": "docs/raw/scienceone_skill_plaza_user_paste.txt",
        "parser_version": "1.0_skill_plaza_paste",
        "entries_total": len(entries),
        "entries_by_main_category": by_cat,
        "note": "用户从磐石技能广场复制的全量清单；tool_chain_key 已与 ToolChain bundle/种子映射对齐或留空待补。",
        "entries": entries,
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8")
    print(f"Wrote {args.out} ({len(entries)} entries)")
    print("by category:", by_cat)


if __name__ == "__main__":
    main()
