# -*- coding: utf-8 -*-
"""Prompt 软技能 · 多轮记忆闭环与终结态（槽位填充 + 强制交付）。"""
from __future__ import annotations

import logging
import re
from typing import Any, Dict, List, Optional, Sequence

from gibh_agent.skills._prompt_skill_completeness import (
    COMPLETENESS_AND_PLACEHOLDER_RULES,
    MARKDOWN_RESUME_LAYOUT_RULES,
    MARKDOWN_WEEKLY_LAYOUT_RULES,
)
from gibh_agent.skills._prompt_skill_llm import build_user_prompt, llm_chat_json
from gibh_agent.skills._prompt_skill_schemas import PromptSkillTerminalResponse
from gibh_agent.skills._prompt_skill_spec import load_prompt_spec

logger = logging.getLogger(__name__)

_TERMINAL_SYSTEM_RULES = """
---
【系统执行 · 终结态 JSON 协议 — 覆盖 SKILL 中「待确认项」类写法】
你必须且只能输出一个 JSON 对象（禁止 Markdown 代码围栏、禁止前言后语），字段：
{
  "action": "deliver" 或 "missing_params",
  "markdown": "当 action=deliver 时为完整终稿 Markdown；当 action=missing_params 时可为空字符串",
  "missing_params": ["仅 action=missing_params 时填写，最多 2 条简短缺失项"],
  "summary": "一行中文状态说明"
}

终结法则（最高优先级）：
1. 综合「历史对话」与当前 user_request/context/file_path。完整度 ≥70% 或用户要求直接写/不用问了 → 必须 action=deliver；无实质输入 → 输出 Demo 演示终稿。
2. action=deliver 时：missing_params 必须为空数组；正文可用 `[待补充：…]` 占位符；禁止「待确认项」「请补充」「请提供」等追问式标题或编号清单；禁止以问句结尾。
3. action=missing_params **仅当完整度 <40%** 且非 Demo 场景；missing_params **最多 2 条**，每条 ≤40 字，只问最核心缺口。
4. 禁止编造用户未提供的量化结果、实验数据、邮件事实或文件内容（占位符处除外）。
""" + COMPLETENESS_AND_PLACEHOLDER_RULES

_TERMINAL_PHRASES = (
    "没有补充",
    "无补充",
    "就这些",
    "就这样",
    "直接写",
    "直接生成",
    "不用问了",
    "别问了",
    "不要追问",
    "无需补充",
    "可以了",
    "够了",
    "就这样吧",
    "no more",
    "that's all",
    "just write",
    "go ahead",
)

_CONFIRM_SECTION_RE = re.compile(
    r"(?:^|\n)\s*#{0,3}\s*(?:待确认项|待补充项|需要您补充|请补充|Please confirm)[^\n]*\n[\s\S]*?(?=\n#{1,3}\s|\Z)",
    re.IGNORECASE,
)


def _format_history_block(history: Optional[Sequence[Dict[str, Any]]], *, max_chars: int = 12_000) -> str:
    if not history:
        return ""
    lines: List[str] = []
    for msg in history:
        if not isinstance(msg, dict):
            continue
        role = (msg.get("role") or "user").strip()
        content = msg.get("content")
        if isinstance(content, list):
            text = " ".join(
                str(p.get("text") or p) for p in content if isinstance(p, dict)
            ).strip()
        else:
            text = str(content or "").strip()
        if not text or text.startswith("【系统上下文"):
            continue
        if len(text) > 2000:
            text = text[:2000] + "…"
        lines.append(f"[{role}]: {text}")
    if not lines:
        return ""
    block = "\n".join(lines)
    if len(block) > max_chars:
        block = block[-max_chars:]
    return f"【历史对话摘录（供整合，勿重复追问已覆盖信息）】\n{block}\n\n"


def user_requests_terminal(text: str) -> bool:
    t = (text or "").lower()
    return any(p in t for p in _TERMINAL_PHRASES)


def strip_confirmation_sections(md: str) -> str:
    s = (md or "").strip()
    if not s:
        return s
    return _CONFIRM_SECTION_RE.sub("\n", s).strip()


def _coverage_heuristic(spec_id: str, combined: str) -> bool:
    """≥70 分启发式：达到则强制 deliver，覆盖模型的 missing_params。"""
    t = (combined or "").lower()
    if not t:
        return False
    if len(t) >= 48:
        return True
    if spec_id == "weekly_report_writer":
        keys = ("本周", "完成", "进行", "下周", "计划", "项目", "周期", "周报", "进展", "待办", "演示")
        return sum(1 for k in keys if k in t) >= 1
    if spec_id == "tailored_resume":
        keys = ("职位", "简历", "经历", "技能", "jd", "job", "python", "sql", "年", "演示", "analyst")
        return sum(1 for k in keys if k in t) >= 1
    return len(t) >= 60


def _parse_terminal_response(raw: Dict[str, Any]) -> PromptSkillTerminalResponse:
    if "data" in raw and isinstance(raw["data"], dict):
        raw = {**raw, **raw["data"]}
    action = (raw.get("action") or raw.get("phase") or "").strip().lower()
    if action not in ("deliver", "missing_params"):
        if raw.get("markdown") and not raw.get("missing_params"):
            action = "deliver"
        else:
            action = "missing_params"
    return PromptSkillTerminalResponse(
        action=action,  # type: ignore[arg-type]
        markdown=str(raw.get("markdown") or ""),
        missing_params=list(raw.get("missing_params") or []),
        summary=str(raw.get("summary") or raw.get("message") or ""),
    )


def execute_terminal_spec_from_kwargs(
    spec_id: str,
    filled: Dict[str, Any],
    kwargs: Dict[str, Any],
    *,
    max_tokens: int = 8192,
) -> Dict[str, Any]:
    """BaseSkill.execute 通用入口：合并 conversation_history 后走终结态。"""
    history: Optional[List[Dict[str, Any]]] = None
    raw_hist = kwargs.get("conversation_history")
    if isinstance(raw_hist, list):
        history = [h for h in raw_hist if isinstance(h, dict)]
    extra = str(kwargs.get("extra_system") or "").strip()
    return run_terminal_prompt_spec_skill(
        spec_id,
        user_request=str(filled.get("user_request") or ""),
        context=str(filled.get("context") or ""),
        file_path=str(filled.get("file_path") or ""),
        conversation_history=history,
        extra_system=extra,
        max_tokens=max_tokens,
    )


def run_terminal_prompt_spec_skill(
    spec_id: str,
    *,
    user_request: str = "",
    context: str = "",
    file_path: str = "",
    conversation_history: Optional[Sequence[Dict[str, Any]]] = None,
    extra_system: str = "",
    max_tokens: int = 8192,
) -> Dict[str, Any]:
    """执行带终结态的 SKILL.md 软技能，返回 deliver 或 missing_params。"""
    from gibh_agent.skills._prompt_skill_spec import (
        _extract_pdf_text,
        _read_text_file,
        run_prompt_spec_skill,
    )

    spec = load_prompt_spec(spec_id)
    if not spec:
        return {
            "status": "error",
            "message": f"未找到技能规范（spec_id={spec_id}）。请部署 prompt_specs/{spec_id}/SKILL.md",
        }

    req = (user_request or "").strip()
    ctx = (context or "").strip()
    fp = (file_path or "").strip()
    history_block = _format_history_block(conversation_history)
    combined = "\n".join(x for x in (history_block, req, ctx) if x)

    force_deliver = user_requests_terminal(req) or user_requests_terminal(ctx) or user_requests_terminal(
        history_block
    )

    system = spec[:100_000] + _TERMINAL_SYSTEM_RULES
    if spec_id == "weekly_report_writer":
        system += "\n" + MARKDOWN_WEEKLY_LAYOUT_RULES
    elif spec_id == "tailored_resume":
        system += "\n" + MARKDOWN_RESUME_LAYOUT_RULES
    if extra_system:
        system += "\n" + extra_system.strip()
    if force_deliver:
        system += (
            "\n【强制】用户已表示无需继续补充或要求直接成稿，你必须 action=deliver，"
            "整合历史与当前全部信息输出完整 Markdown 终稿。"
        )

    user_parts: List[str] = []
    if history_block:
        user_parts.append(history_block.strip())
    if req:
        user_parts.append(f"当前用户需求：\n{req}")
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
            user_parts.append(f"{kind}（file_path={fp}）：\n{extracted}")
        else:
            user_parts.append(f"file_path={fp}（未能自动读取，请基于已有文字处理）")

    if not user_parts:
        try:
            user_parts.append(build_user_prompt("请按 SKILL 执行广场演示或默认任务。"))
        except ValueError:
            return {"status": "error", "message": "请提供 user_request、context 或 file_path。"}

    user = "\n\n".join(user_parts)

    try:
        raw = llm_chat_json(system, user, max_tokens=max_tokens)
        terminal = _parse_terminal_response(raw)
    except Exception as exc:
        logger.exception("terminal prompt skill %s 结构化失败，回退非终结执行", spec_id)
        return run_prompt_spec_skill(
            spec_id,
            user_request=req,
            context=ctx,
            file_path=fp,
            extra_system=(
                "【回退】请直接输出完整 Markdown 终稿，不要列出待确认项或追问。"
            ),
            max_tokens=max_tokens,
        )

    if force_deliver or (
        _coverage_heuristic(spec_id, combined) and terminal.action == "missing_params"
    ):
        terminal = PromptSkillTerminalResponse(
            action="deliver",
            markdown=terminal.markdown or "",
            missing_params=[],
            summary=terminal.summary or "已根据现有信息整合输出",
        )

    if terminal.action == "missing_params" and not terminal.missing_params:
        if spec_id == "tailored_resume":
            default_missing = ["目标职位或 JD 关键词（一句话即可）"]
        elif spec_id == "weekly_report_writer":
            default_missing = ["报告周期（起止日期）与本周 1–2 条核心进展"]
        else:
            default_missing = ["请补充最关键的一项信息（一句话即可）"]
        terminal = PromptSkillTerminalResponse(
            action="missing_params",
            markdown="",
            missing_params=default_missing[:2],
            summary="信息不足，需补充关键参数",
        )

    if terminal.action == "deliver":
        md = strip_confirmation_sections(terminal.markdown)
        if not md.strip():
            md = terminal.summary or "（模型未返回正文，请重试或补充需求。）"
        summary = (terminal.summary or "终稿已生成").strip()
        fe = {
            "phase": "deliver",
            "markdown": md,
            "spec_id": spec_id,
        }
        return {
            "status": "success",
            "phase": "deliver",
            "message": summary,
            "markdown": md,
            "data": fe,
            "spec_id": spec_id,
        }

    missing = terminal.missing_params[:2]
    brief = "；".join(missing)
    return {
        "status": "success",
        "phase": "missing_params",
        "message": f"尚需补充：{brief}",
        "missing_params": missing,
        "markdown": "",
        "data": {"phase": "missing_params", "missing_params": missing, "spec_id": spec_id},
        "spec_id": spec_id,
    }
