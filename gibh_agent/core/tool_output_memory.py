# -*- coding: utf-8 -*-
"""
工具结构化输出 → 会话记忆注入：从磁盘读取 chem_summary.json 等，供下一轮 LLM 引用。

与 chem_rdkit_tools._safe_load_chem_summary_json 对齐的体积上限，避免拖垮主进程。
"""
from __future__ import annotations

import json
import logging
import re
from pathlib import Path
from typing import Any, Dict, List, Optional

logger = logging.getLogger(__name__)

_DEFAULT_MAX_INJECTION_CHARS = 14_000
_DEFAULT_MAX_FILE_BYTES = 2_000_000


def _read_text_file_safe(path: Path, max_bytes: int) -> Optional[str]:
    try:
        raw = path.read_bytes()
        if len(raw) > max_bytes:
            logger.warning("tool_output_memory: 文件过大跳过 %s (%s bytes)", path, len(raw))
            return None
        return raw.decode("utf-8", errors="replace")
    except (OSError, MemoryError) as e:
        logger.warning("tool_output_memory: 读取失败 %s: %s", path, e)
        return None


def _pick_json_candidates_from_data(data: Dict[str, Any]) -> List[Path]:
    out: List[Path] = []
    od = data.get("output_dir") or data.get("output_dir_path")
    if od:
        base = Path(str(od)).expanduser().resolve()
        preferred = base / "chem_summary.json"
        if preferred.is_file():
            out.append(preferred)
    jps = data.get("json_paths")
    if isinstance(jps, list):
        for jp in jps:
            if not jp:
                continue
            p = Path(str(jp)).expanduser().resolve()
            if p.is_file() and p.suffix.lower() == ".json":
                name = p.name.lower()
                if "chem_summary" in name or "summary" in name:
                    out.append(p)
                elif not out:
                    out.append(p)
    # 去重保序
    seen = set()
    uniq: List[Path] = []
    for p in out:
        s = str(p)
        if s not in seen:
            seen.add(s)
            uniq.append(p)
    return uniq


def _fallback_compact_json(result: Dict[str, Any], max_chars: int) -> str:
    """无磁盘摘要时，将工具返回体压成可投喂 JSON（去掉 markdown 等冗长字段）。"""
    slim: Dict[str, Any] = {}
    for k, v in result.items():
        if k in ("markdown", "image_urls"):
            continue
        slim[k] = v
    try:
        s = json.dumps(slim, ensure_ascii=False, indent=2, default=str)
    except Exception:
        s = str(slim)
    if len(s) > max_chars:
        s = s[: max_chars - 80] + "\n…(truncated)…"
    return s


def build_tool_output_memory_text(
    tool_name: str,
    result: Dict[str, Any],
    *,
    max_chars: int = _DEFAULT_MAX_INJECTION_CHARS,
) -> str:
    """
    生成可写入会话历史 / 系统上下文的纯文本块（中文说明 + JSON）。

    成功路径下优先读取输出目录中的 chem_summary.json 或 json_paths 中的摘要文件；
    否则回退为精简后的整包 JSON。
    """
    if not isinstance(result, dict) or result.get("status") != "success":
        return ""

    body_parts: List[str] = []
    data = result.get("data")
    loaded_from_file = False

    if isinstance(data, dict):
        for jp in _pick_json_candidates_from_data(data):
            txt = _read_text_file_safe(jp, _DEFAULT_MAX_FILE_BYTES)
            if not txt:
                continue
            snippet = txt.strip()
            if len(snippet) > max_chars - 400:
                snippet = snippet[: max_chars - 420] + "\n…(truncated)…"
            body_parts.append(f"文件 `{jp.name}` 核心内容：\n```json\n{snippet}\n```")
            loaded_from_file = True
            break

    if not loaded_from_file:
        fb = _fallback_compact_json(result, max_chars=max_chars - 200)
        if fb:
            body_parts.append(f"工具返回体（精简 JSON）：\n```json\n{fb}\n```")

    if not body_parts:
        return ""

    header = (
        f"【系统上下文｜上一轮工具 `{tool_name}` 结构化输出】\n"
        "以下数据已由服务端从本次运行结果中读取，你可直接据此回答用户的追问（勿声称无法读取生成文件）。\n"
    )
    return header + "\n".join(body_parts)


def strip_redundant_tool_memory_from_assistant_text(text: str) -> str:
    """
    若 assistant 正文重复包含整块「系统上下文」摘要，可在展示层裁剪（当前未默认调用）。
    """
    if not text or "【系统上下文｜上一轮工具" not in text:
        return text
    return re.sub(
        r"\n【系统上下文｜上一轮工具 `[^`]+` 结构化输出】[\s\S]*$",
        "",
        text,
        count=1,
    )
