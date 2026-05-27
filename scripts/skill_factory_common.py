# -*- coding: utf-8 -*-
"""Skill Factory 脚本共用：LLM 客户端、代码提取、slug、依赖汇总与 CSV 行解析。"""
from __future__ import annotations

import ast
import os
import re
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List, Optional, Set, Tuple

from gibh_agent.core.llm_client import LLMClientFactory

DEFAULT_FACTORY_MODEL = os.getenv("SKILL_FACTORY_MODEL", "deepseek-chat")


def factory_llm_client(model: Optional[str] = None):
    return LLMClientFactory.create_for_model(
        (model or DEFAULT_FACTORY_MODEL).strip(),
        temperature=float(os.getenv("SKILL_FACTORY_TEMPERATURE", "0.2")),
        max_tokens=int(os.getenv("SKILL_FACTORY_MAX_TOKENS", "8192")),
    )


def extract_python_from_llm(text: str) -> str:
    """从 LLM 回复中提取 Python 源码（优先 ```python 围栏）。"""
    raw = (text or "").strip()
    if not raw:
        return ""
    m = re.search(r"```(?:python)?\s*\n(.*?)```", raw, re.DOTALL | re.IGNORECASE)
    if m:
        return m.group(1).strip()
    if raw.startswith("```"):
        raw = re.sub(r"^```\w*\n?", "", raw)
        raw = re.sub(r"\n?```\s*$", "", raw)
    return raw.strip()


def camel_to_snake(name: str) -> str:
    s = re.sub(r"([a-z0-9])([A-Z])", r"\1_\2", (name or "").strip())
    s = re.sub(r"[^a-zA-Z0-9_]+", "_", s)
    return s.lower().strip("_")


def skill_id_from_row(row: Dict[str, Any]) -> str:
    key = (row.get("tool_chain_key") or row.get("tool_chain") or "").strip()
    name = (row.get("技能名称") or row.get("name") or "").strip()
    if key:
        sid = camel_to_snake(key)
    else:
        sid = camel_to_snake(name)
    if not sid:
        sid = "unnamed_skill"
    if not sid[0].isalpha():
        sid = f"skill_{sid}"
    return sid[:64]


def output_filename_for_skill_id(skill_id: str) -> str:
    safe = re.sub(r"[^a-zA-Z0-9_]", "_", (skill_id or "skill").strip())[:80]
    return f"skill_{safe}.py"


def parse_csv_row(row: Dict[str, Any]) -> Dict[str, str]:
    return {
        "discipline": (row.get("学科") or "").strip(),
        "display_name": (row.get("技能名称") or row.get("name") or "").strip(),
        "sub_category": (row.get("子分类") or "").strip(),
        "skills_type": (row.get("skills类型") or row.get("skills_type") or "").strip(),
        "description": (row.get("简介") or row.get("description") or "").strip(),
        "tool_chain_key": (row.get("tool_chain_key") or "").strip(),
    }


def row_to_prompt_block(parsed: Dict[str, str], skill_id: str) -> str:
    return (
        f"skill_id: {skill_id}\n"
        f"display_name: {parsed['display_name']}\n"
        f"discipline / category: {parsed['discipline']}\n"
        f"sub_category: {parsed['sub_category']}\n"
        f"skills_type: {parsed['skills_type']}\n"
        f"tool_chain_key: {parsed['tool_chain_key']}\n"
        f"description: {parsed['description']}\n"
    )


# --- 依赖提取（AST 优先，正则兜底）---

_DEPENDENCY_LIST_RE = re.compile(
    r"__dependencies__\s*(?::\s*[^=]+)?=\s*(\[[\s\S]*?\])",
    re.MULTILINE,
)
_STRING_LITERAL_RE = re.compile(r"['\"]([^'\"]+)['\"]")


def _string_list_from_ast_node(node: ast.AST) -> List[str]:
    if isinstance(node, ast.List):
        out: List[str] = []
        for elt in node.elts:
            if isinstance(elt, ast.Constant) and isinstance(elt.value, str):
                out.append(elt.value.strip())
            elif isinstance(elt, ast.Str):  # noqa: UP038 — py<3.8 compat
                out.append(elt.s.strip())
        return out
    if isinstance(node, ast.Constant) and isinstance(node.value, list):
        return [str(x).strip() for x in node.value if str(x).strip()]
    return []


def extract_dependencies_from_source(source: str) -> List[str]:
    """从技能模块源码解析 ``__dependencies__`` 字符串列表。"""
    text = source or ""
    if not text.strip():
        return []

    try:
        tree = ast.parse(text)
    except SyntaxError:
        tree = None

    if tree is not None:
        for node in ast.walk(tree):
            if not isinstance(node, ast.ClassDef):
                continue
            for item in node.body:
                target_name: Optional[str] = None
                value_node: Optional[ast.AST] = None
                if isinstance(item, ast.Assign):
                    for t in item.targets:
                        if isinstance(t, ast.Name) and t.id == "__dependencies__":
                            target_name = "__dependencies__"
                            value_node = item.value
                            break
                elif isinstance(item, ast.AnnAssign) and isinstance(item.target, ast.Name):
                    if item.target.id == "__dependencies__":
                        target_name = "__dependencies__"
                        value_node = item.value
                if target_name and value_node is not None:
                    deps = _string_list_from_ast_node(value_node)
                    if deps:
                        return deps

    m = _DEPENDENCY_LIST_RE.search(text)
    if m:
        blob = m.group(1)
        return [s.strip() for s in _STRING_LITERAL_RE.findall(blob) if s.strip()]
    return []


def extract_dependencies_from_skill_file(py_path: Path) -> List[str]:
    try:
        return extract_dependencies_from_source(py_path.read_text(encoding="utf-8"))
    except OSError:
        return []


def collect_skill_dependency_map(skills_dir: Path) -> Dict[str, List[str]]:
    """skill_*.py → 依赖列表（仅含非空 __dependencies__）。"""
    out: Dict[str, List[str]] = {}
    if not skills_dir.is_dir():
        return out
    for py_path in sorted(skills_dir.glob("skill_*.py")):
        deps = extract_dependencies_from_skill_file(py_path)
        if deps:
            out[py_path.name] = deps
    return out


def write_skill_dependency_report(
    skills_dir: Path,
    report_path: Path,
    *,
    only_files: Optional[List[Path]] = None,
) -> Tuple[int, int]:
    """
    汇总依赖并写入 ``docs/skill_missing_dependencies.txt``。

    Returns:
        (去重依赖条数, 涉及技能文件数)
    """
    per_file: Dict[str, List[str]] = {}
    paths: List[Path]
    if only_files:
        paths = [p for p in only_files if p.is_file()]
    else:
        paths = sorted(skills_dir.glob("skill_*.py")) if skills_dir.is_dir() else []

    for py_path in paths:
        deps = extract_dependencies_from_skill_file(py_path)
        if deps:
            per_file[py_path.name] = deps

    unique: Set[str] = set()
    for deps in per_file.values():
        unique.update(deps)

    by_kind: Dict[str, List[str]] = {"pip": [], "apt": [], "conda": [], "other": []}
    for dep in sorted(unique):
        if ":" in dep:
            kind, _rest = dep.split(":", 1)
            bucket = kind.lower() if kind.lower() in by_kind else "other"
        else:
            bucket = "other"
        by_kind.setdefault(bucket, []).append(dep)

    ts = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M:%S UTC")
    lines: List[str] = [
        "# Skill Factory — 镜像待补齐依赖汇总",
        f"# 生成时间: {ts}",
        f"# 扫描目录: {skills_dir.resolve()}",
        f"# 涉及技能文件: {len(per_file)} | 去重依赖条目: {len(unique)}",
        "#",
        "# 格式: pip:包名 | apt:deb包 | conda:包名 | 其它自定义前缀",
        "# 运维: 请 DevOps 据此集中更新 services/api/Dockerfile / worker 镜像后 rebuild",
        "",
    ]

    if not unique:
        lines.append("（当前扫描范围内无 __dependencies__ 声明，或均为空列表 []）")
    else:
        lines.append("## 按安装通道分组（去重）")
        for kind in ("pip", "apt", "conda", "other"):
            items = by_kind.get(kind) or []
            if not items:
                continue
            lines.append(f"\n### {kind.upper()}")
            for item in items:
                lines.append(item)

        lines.append("\n## 按技能文件")
        for fname in sorted(per_file.keys()):
            lines.append(f"\n### {fname}")
            for d in per_file[fname]:
                lines.append(f"  - {d}")

    report_path.parent.mkdir(parents=True, exist_ok=True)
    report_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return len(unique), len(per_file)
