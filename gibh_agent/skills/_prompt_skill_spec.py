# -*- coding: utf-8 -*-
"""从 prompt_specs/SKILL.md 加载规范并驱动 LLM 执行 Prompt-Engineered 软技能。"""
from __future__ import annotations

import logging
from pathlib import Path
from typing import Any, Dict, Optional

import os

from gibh_agent.skills._prompt_skill_llm import build_user_prompt, llm_chat_text

logger = logging.getLogger(__name__)

_SPECS_PACKAGE_DIR = Path(__file__).resolve().parent / "prompt_specs"
_DOCS_PROMPT_DIR = Path(__file__).resolve().parents[2] / "docs" / "skills" / "prompt"

# 非终结态技能回退时使用（终结态技能请走 _prompt_skill_terminal.py）
_SYSTEM_SUFFIX = """
---
【系统执行附加约束】
1. 严格遵循上文 SKILL 文档中的工作流、语气与输出格式。
2. 默认交付完整 Markdown 正文，禁止编造用户未提供的量化结果、实验数据、邮件事实或文件内容。
"""


def _runtime_spec_file(spec_id: str) -> Optional[Path]:
    """查找运行时 prompt 目录中的 SKILL.md（不创建目录，避免宿主机误触 /app/uploads）。"""
    roots: list[Path] = []
    explicit = (os.environ.get("SKILLS_ASSETS_ROOT") or "").strip()
    if explicit:
        roots.append(Path(os.path.abspath(explicit)) / "prompt")
    upload = (os.environ.get("UPLOAD_DIR") or "").strip()
    if upload:
        roots.append(Path(os.path.abspath(upload)) / "skills_assets" / "prompt")
    repo_upload = Path(__file__).resolve().parents[2] / "data" / "uploads" / "skills_assets" / "prompt"
    roots.append(repo_upload)
    for root in roots:
        p = root / spec_id / "SKILL.md"
        if p.is_file():
            return p
    return None


def resolve_spec_path(spec_id: str) -> Optional[Path]:
    sid = (spec_id or "").strip()
    if not sid:
        return None
    for p in (
        _SPECS_PACKAGE_DIR / sid / "SKILL.md",
        _DOCS_PROMPT_DIR / f"{sid}.md",
    ):
        if p.is_file():
            return p
    return _runtime_spec_file(sid)


def load_prompt_spec(spec_id: str) -> str:
    path = resolve_spec_path(spec_id)
    if not path:
        return ""
    try:
        return path.read_text(encoding="utf-8").strip()
    except OSError as exc:
        logger.warning("读取 SKILL 规范失败 %s: %s", path, exc)
        return ""


def _read_text_file(file_path: str, *, max_chars: int = 80_000) -> str:
    p = Path(file_path).expanduser()
    if not p.is_file():
        return ""
    try:
        text = p.read_text(encoding="utf-8", errors="replace")
    except OSError:
        return ""
    if len(text) > max_chars:
        return text[:max_chars] + "\n\n…（附件已截断）"
    return text


def _extract_pdf_text(file_path: str, *, max_chars: int = 60_000) -> str:
    p = Path(file_path).expanduser()
    if not p.is_file():
        return ""
    try:
        from pypdf import PdfReader  # type: ignore
    except ImportError:
        try:
            import PyPDF2  # type: ignore

            reader = PyPDF2.PdfReader(str(p))
            parts = []
            for page in reader.pages:
                parts.append(page.extract_text() or "")
            text = "\n".join(parts)
        except Exception:
            logger.debug("PDF 解析依赖不可用（可 pip install pypdf）")
            return ""
    else:
        reader = PdfReader(str(p))
        parts = [page.extract_text() or "" for page in reader.pages]
        text = "\n".join(parts)
    text = (text or "").strip()
    if len(text) > max_chars:
        return text[:max_chars] + "\n\n…（PDF 文本已截断）"
    return text


def run_prompt_spec_skill(
    spec_id: str,
    *,
    user_request: str = "",
    context: str = "",
    file_path: str = "",
    extra_system: str = "",
    max_tokens: int = 8192,
    temperature: float = 0.35,
) -> Dict[str, Any]:
    """执行基于 SKILL.md 的软技能，返回标准工具结构（含 markdown）。"""
    spec = load_prompt_spec(spec_id)
    if not spec:
        return {
            "status": "error",
            "message": f"未找到技能规范文件（spec_id={spec_id}）。请确认 prompt_specs/{spec_id}/SKILL.md 已部署。",
        }

    system = spec[:100_000] + _SYSTEM_SUFFIX
    if extra_system:
        system += "\n" + extra_system.strip()

    user_parts: list[str] = []
    req = (user_request or "").strip()
    ctx = (context or "").strip()
    fp = (file_path or "").strip()

    if req:
        user_parts.append(f"用户需求：\n{req}")
    if ctx:
        user_parts.append(f"补充上下文：\n{ctx}")

    if fp:
        if spec_id == "pdf_extractor" or fp.lower().endswith(".pdf"):
            extracted = _extract_pdf_text(fp)
            kind = "PDF 提取文本"
        else:
            extracted = _read_text_file(fp)
            kind = "附件文本"
        if extracted:
            user_parts.append(f"{kind}（来自 file_path={fp}）：\n{extracted}")
        else:
            user_parts.append(
                f"用户提供的 file_path={fp}（当前环境未能自动读取；"
                "请根据路径说明处理步骤，或请用户在 user_request 中粘贴正文。）"
            )

    if not user_parts:
        try:
            user_parts.append(build_user_prompt("请按 SKILL 文档执行广场一键体验或默认演示任务。"))
        except ValueError:
            return {"status": "error", "message": "请提供 user_request、context 或 file_path。"}

    user = "\n\n".join(user_parts)
    try:
        md = llm_chat_text(system, user, max_tokens=max_tokens, temperature=temperature)
    except Exception as exc:
        logger.exception("prompt spec skill %s LLM 失败", spec_id)
        return {"status": "error", "message": f"技能执行失败：{exc}"}

    return {
        "status": "success",
        "message": "已根据技能规范生成结果，见下方 Markdown。",
        "markdown": md,
        "spec_id": spec_id,
    }
