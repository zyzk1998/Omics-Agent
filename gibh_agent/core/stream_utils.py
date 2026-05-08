"""
Stream and text parsing utilities for LLM responses.
- stream_from_llm_chunks：丢弃 delta.reasoning_content 与 THINK_TAG_PAIRS 标记的思考正文，不向用户 SSE 外露；外露仅 message suggest suggestions。
- stream_and_extract_json：Planner 流式分拣——思考块不落 json_buffer；流结束后剥 fence + 括号平衡 json.loads。

工具执行期面向用户的过程日志（非 LLM token）：使用 gibh_agent.core.tool_stream_log.emit_tool_log，经 Executor/Orchestrator 转为 SSE status；本模块不重复拼装 SSE。
"""
import json
import re
import logging
from typing import Any, AsyncIterator, List, Optional, Set, Tuple, Union

logger = logging.getLogger(__name__)

SUGGEST_START = "<<<SUGGESTIONS>>>"
SUGGEST_END = "<<<END_SUGGESTIONS>>>"
# 新版追问：与前端 handleServerEvent('suggest') 对齐；流式由 stream_from_llm_chunks 拦截，不进入 message
TAG_SUGGEST_OPEN = "<suggest>"
TAG_SUGGEST_CLOSE = "</suggest>"
# 兼容多厂商：content 中出现的思考块起止标签（去重保序）
_THINK_PAIRS_ORDERED: Tuple[Tuple[str, str], ...] = (
    ('<think>', '</think>'),
    ("<redacted_reasoning>", "</redacted_reasoning>"),
    ("<reasoning>", "</reasoning>"),
    ("<thought>", "</thought>"),
    ("<thinking>", "</thinking>"),
)

def _uniq_think_pairs(pairs: Tuple[Tuple[str, str], ...]) -> Tuple[Tuple[str, str], ...]:
    seen = set()
    out: List[Tuple[str, str]] = []
    for o, c in pairs:
        k = (o, c)
        if k not in seen:
            seen.add(k)
            out.append(k)
    return tuple(out)


THINK_TAG_PAIRS: Tuple[Tuple[str, str], ...] = _uniq_think_pairs(_THINK_PAIRS_ORDERED)
# 技能填参等提示语仍用首对标签占位（与旧测试/提示词一致）
THINK_OPEN, THINK_CLOSE = THINK_TAG_PAIRS[0]
THINK_OPEN_LEN = len(THINK_OPEN)
THINK_CLOSE_LEN = len(THINK_CLOSE)
THINK_MAX_OPEN_LEN = max(len(o) for o, _ in THINK_TAG_PAIRS) if THINK_TAG_PAIRS else 0
THINK_MAX_CLOSE_LEN = max(len(c) for _, c in THINK_TAG_PAIRS) if THINK_TAG_PAIRS else 0


def _find_earliest_think_open(buf: str) -> Optional[Tuple[int, str, str]]:
    best_idx: Optional[int] = None
    best_o = ""
    best_c = ""
    for o, c in THINK_TAG_PAIRS:
        i = buf.find(o)
        if i < 0:
            continue
        if best_idx is None or i < best_idx:
            best_idx = i
            best_o, best_c = o, c
    if best_idx is None:
        return None
    return (best_idx, best_o, best_c)


def strip_think_markup_for_user(text: str) -> str:
    """剥离已知思考/推理标签整块内容，用于非流式文本与快照落库清洗。"""
    if not text or not isinstance(text, str):
        return ""
    t = text
    for open_m, close_m in THINK_TAG_PAIRS:
        pat = re.escape(open_m) + r"[\s\S]*?" + re.escape(close_m)
        t = re.sub(pat, "", t, flags=re.IGNORECASE)
    # 流式中断：无闭合结束标签时，从起始标签起删至文末（防泄漏）
    seen_open: Set[str] = set()
    for open_m, _ in THINK_TAG_PAIRS:
        if open_m in seen_open:
            continue
        seen_open.add(open_m)
        pat_tail = re.escape(open_m) + r"[\s\S]*$"
        t = re.sub(pat_tail, "", t, flags=re.IGNORECASE)
    return t.strip()


def _delta_get(delta: Any, key: str) -> Optional[str]:
    """从 delta 安全取值，兼容对象（getattr）与字典（.get）。"""
    if delta is None:
        return None
    if isinstance(delta, dict):
        val = delta.get(key)
        return str(val) if val is not None else None
    val = getattr(delta, key, None)
    return str(val) if val is not None else None


def parse_suggest_questions_json(inner: str) -> Optional[List[str]]:
    """解析 <suggest> 与 </suggest> 之间的 JSON，要求根对象含 questions 数组。"""
    if not inner or not str(inner).strip():
        return None
    s = str(inner).strip()
    try:
        obj = json.loads(s)
    except (json.JSONDecodeError, TypeError):
        balanced = _extract_balanced_json_object(s)
        if not balanced:
            return None
        try:
            obj = json.loads(balanced)
        except (json.JSONDecodeError, TypeError):
            return None
    if not isinstance(obj, dict):
        return None
    qs = obj.get("questions")
    if not isinstance(qs, list):
        return None
    out = [str(q).strip() for q in qs if q is not None and str(q).strip()]
    return out or None


def strip_suggestions_from_text(text: str) -> Tuple[str, Optional[List[str]]]:
    """
    Remove the suggestions block from a full LLM response and return cleaned text + parsed list.
    Use for non-streaming responses (diagnosis, report).

    Args:
        text: Full response string (may contain <suggest>{"questions":[...]}</suggest>
        and/or <<<SUGGESTIONS>>>["Q1","Q2"]<<<END_SUGGESTIONS>>>).

    Returns:
        (cleaned_text, suggestions_list or None). If no block found, cleaned_text is text unchanged.
    """
    if not text:
        return (text, None)
    t = text
    merged: Optional[List[str]] = None

    if TAG_SUGGEST_OPEN in t and TAG_SUGGEST_CLOSE in t:
        io = t.find(TAG_SUGGEST_OPEN)
        ic = t.find(TAG_SUGGEST_CLOSE)
        if io != -1 and ic != -1 and io < ic:
            ic_full = ic + len(TAG_SUGGEST_CLOSE)
            inner = t[io + len(TAG_SUGGEST_OPEN) : ic].strip()
            qs = parse_suggest_questions_json(inner)
            if qs:
                merged = list(qs)
            t = (t[:io] + t[ic_full:]).strip()

    if SUGGEST_START not in t or SUGGEST_END not in t:
        return (t, merged)

    idx_start = t.find(SUGGEST_START)
    idx_end = t.find(SUGGEST_END) + len(SUGGEST_END)
    before = t[:idx_start].rstrip()
    after = t[idx_end:].lstrip()
    cleaned = (before + "\n\n" + after).strip() if after else before
    try:
        json_str = t[idx_start + len(SUGGEST_START) : idx_end - len(SUGGEST_END)].strip()
        suggestions = json.loads(json_str)
        if isinstance(suggestions, list) and suggestions:
            old = [str(x).strip() for x in suggestions if str(x).strip()]
            if merged is None:
                merged = old
            else:
                merged = merged + [x for x in old if x not in merged]
            return (cleaned, merged)
    except (json.JSONDecodeError, TypeError):
        pass
    return (cleaned, merged)


async def stream_from_llm_chunks(
    chunk_iter: AsyncIterator[Any],
    model_name: Optional[str] = None,
) -> AsyncIterator[Tuple[str, object]]:
    """
    OpenAI-compat 流式状态机：
    1) delta.reasoning_content：不入下游（仅存 debug 日志），不向用户 message 透出。
    2) delta.content：对 THINK_TAG_PAIRS 多块标签做多路匹配；标签内字节一律丢弃不外发；
       外露文本仅 yield ("message", ...)；跨 chunk 用 THINK_MAX_OPEN_LEN / pending_close 尾保留防抖。
    3) <suggest>…</suggest>、<<<SUGGESTIONS>>> 行为不变。
    """
    content_buffer = ""
    in_think_tag = False
    pending_close: Optional[str] = None

    async for chunk in chunk_iter:
        choices = getattr(chunk, "choices", None) if not isinstance(chunk, dict) else chunk.get("choices")
        if not choices or len(choices) == 0:
            continue
        first = choices[0]
        delta = getattr(first, "delta", None) if not isinstance(first, dict) else first.get("delta")

        reasoning = _delta_get(delta, "reasoning_content")
        if reasoning:
            logger.debug(
                "stream_from_llm_chunks: dropped reasoning_content (~%d chars)",
                len(reasoning),
            )

        content = _delta_get(delta, "content")
        if content:
            content_buffer += content

        while True:
            if TAG_SUGGEST_OPEN in content_buffer and TAG_SUGGEST_CLOSE in content_buffer:
                io = content_buffer.find(TAG_SUGGEST_OPEN)
                ic = content_buffer.find(TAG_SUGGEST_CLOSE)
                if io != -1 and ic != -1 and io < ic:
                    ic_full = ic + len(TAG_SUGGEST_CLOSE)
                    inner = content_buffer[io + len(TAG_SUGGEST_OPEN) : ic].strip()
                    qs = parse_suggest_questions_json(inner)
                    if qs:
                        yield ("suggest", {"questions": qs})
                    content_buffer = content_buffer[:io] + content_buffer[ic_full:]
                    continue

            if SUGGEST_END in content_buffer and SUGGEST_START in content_buffer:
                idx_start = content_buffer.find(SUGGEST_START)
                idx_end = content_buffer.find(SUGGEST_END) + len(SUGGEST_END)
                if idx_start < idx_end:
                    try:
                        json_str = content_buffer[
                            idx_start + len(SUGGEST_START) : idx_end - len(SUGGEST_END)
                        ].strip()
                        suggestions = json.loads(json_str)
                        if isinstance(suggestions, list) and suggestions:
                            yield ("suggestions", suggestions)
                    except (json.JSONDecodeError, TypeError):
                        pass
                    content_buffer = content_buffer[:idx_start] + content_buffer[idx_end:]
                    continue

            if in_think_tag:
                pc = pending_close or ""
                pcl = len(pc)
                if pcl and pc in content_buffer:
                    idx = content_buffer.find(pc)
                    dropped = content_buffer[:idx]
                    if dropped:
                        logger.debug(
                            "stream_from_llm_chunks: dropped think (~%d chars)",
                            len(dropped),
                        )
                    content_buffer = content_buffer[idx + pcl :]
                    in_think_tag = False
                    pending_close = None
                    continue
                if pcl > 0 and len(content_buffer) > pcl:
                    drop_part = content_buffer[:-pcl]
                    content_buffer = content_buffer[-pcl:]
                    if drop_part:
                        logger.debug(
                            "stream_from_llm_chunks: dropped partial think (~%d chars)",
                            len(drop_part),
                        )
                    break
                break

            hit = _find_earliest_think_open(content_buffer)
            if hit is not None:
                idx, open_m, close_m = hit
                message_part = content_buffer[:idx]
                if message_part:
                    yield ("message", {"content": message_part})
                content_buffer = content_buffer[idx + len(open_m) :]
                in_think_tag = True
                pending_close = close_m
                continue

            if TAG_SUGGEST_OPEN in content_buffer:
                io = content_buffer.find(TAG_SUGGEST_OPEN)
                ic_rel = content_buffer.find(TAG_SUGGEST_CLOSE, io + len(TAG_SUGGEST_OPEN))
                if ic_rel == -1:
                    msg_part = content_buffer[:io]
                    if msg_part:
                        yield ("message", {"content": msg_part})
                    content_buffer = content_buffer[io:]
                    break
            if SUGGEST_START in content_buffer:
                safe_end = content_buffer.find(SUGGEST_START)
                to_yield = content_buffer[:safe_end]
                if to_yield:
                    yield ("message", {"content": to_yield})
                content_buffer = content_buffer[safe_end:]
                break
            if THINK_MAX_OPEN_LEN > 0 and len(content_buffer) > THINK_MAX_OPEN_LEN:
                to_yield = content_buffer[:-THINK_MAX_OPEN_LEN]
                content_buffer = content_buffer[-THINK_MAX_OPEN_LEN:]
                if to_yield:
                    yield ("message", {"content": to_yield})
                break

            break

    if content_buffer:
        if in_think_tag:
            logger.debug(
                "stream_from_llm_chunks: dropped trailing think buffer (~%d chars)",
                len(content_buffer),
            )
            return

        remaining = content_buffer
        if TAG_SUGGEST_OPEN in remaining:
            remaining = remaining.split(TAG_SUGGEST_OPEN, 1)[0]
        remaining = remaining.split(SUGGEST_START, 1)[0]
        if remaining:
            yield ("message", {"content": remaining})


async def stream_with_suggestions(
    chunk_iter: AsyncIterator[str],
) -> AsyncIterator[Tuple[str, object]]:
    """
    Consume an async generator of content chunks (e.g. from LLM stream), strip the
    <<<SUGGESTIONS>>>...<<<END_SUGGESTIONS>>> block from the stream, and yield:
    - ("message", {"content": str}) for visible text
    - ("suggestions", list) when the block is complete.

    Caller should format and yield SSE (e.g. format_sse(event_type, data)).

    Args:
        chunk_iter: Async iterator yielding content chunks (strings).

    Yields:
        ("message", {"content": str}) or ("suggestions", list).
    """
    stream_buffer = ""
    last_yielded = 0

    async for delta_content in chunk_iter:
        if not delta_content:
            continue
        stream_buffer += delta_content

        while True:
            if TAG_SUGGEST_OPEN in stream_buffer and TAG_SUGGEST_CLOSE in stream_buffer:
                io = stream_buffer.find(TAG_SUGGEST_OPEN)
                ic = stream_buffer.find(TAG_SUGGEST_CLOSE)
                if io != -1 and ic != -1 and io < ic:
                    ic_full = ic + len(TAG_SUGGEST_CLOSE)
                    to_yield = stream_buffer[last_yielded:io]
                    if to_yield:
                        yield ("message", {"content": to_yield})
                    inner = stream_buffer[io + len(TAG_SUGGEST_OPEN) : ic].strip()
                    qs = parse_suggest_questions_json(inner)
                    if qs:
                        yield ("suggest", {"questions": qs})
                    stream_buffer = stream_buffer[ic_full:]
                    last_yielded = 0
                    continue
            if SUGGEST_END in stream_buffer and SUGGEST_START in stream_buffer:
                idx_start = stream_buffer.find(SUGGEST_START)
                idx_end = stream_buffer.find(SUGGEST_END) + len(SUGGEST_END)
                to_yield = stream_buffer[last_yielded:idx_start]
                if to_yield:
                    yield ("message", {"content": to_yield})
                try:
                    json_str = stream_buffer[
                        idx_start + len(SUGGEST_START) : idx_end - len(SUGGEST_END)
                    ].strip()
                    suggestions = json.loads(json_str)
                    if isinstance(suggestions, list) and suggestions:
                        yield ("suggestions", suggestions)
                except (json.JSONDecodeError, TypeError):
                    pass
                stream_buffer = stream_buffer[idx_end:]
                last_yielded = 0
                continue
            if TAG_SUGGEST_OPEN in stream_buffer:
                io = stream_buffer.find(TAG_SUGGEST_OPEN)
                if stream_buffer.find(TAG_SUGGEST_CLOSE, io + len(TAG_SUGGEST_OPEN)) == -1:
                    to_yield = stream_buffer[last_yielded:io]
                    if to_yield:
                        yield ("message", {"content": to_yield})
                    last_yielded = io
                    break
            if SUGGEST_START in stream_buffer:
                safe_end = stream_buffer.find(SUGGEST_START)
                to_yield = stream_buffer[last_yielded:safe_end]
                if to_yield:
                    yield ("message", {"content": to_yield})
                last_yielded = safe_end
                break
            to_yield = stream_buffer[last_yielded:]
            if to_yield:
                yield ("message", {"content": to_yield})
            last_yielded = len(stream_buffer)
            break

    # Flush remaining safe content
    if last_yielded < len(stream_buffer):
        remaining = stream_buffer[last_yielded:]
        if TAG_SUGGEST_OPEN in remaining:
            remaining = remaining.split(TAG_SUGGEST_OPEN, 1)[0]
        remaining = remaining.split(SUGGEST_START, 1)[0]
        if remaining:
            yield ("message", {"content": remaining})


async def stream_with_thought_and_suggestions(
    chunk_iter: AsyncIterator[str],
) -> AsyncIterator[Tuple[str, object]]:
    """
    兼容旧接口名：对大模型输出的「字符串增量」剔除思考标签内字节（不外发），仅下发 message / suggest / suggestions。
    """
    content_buffer = ""
    in_think_tag = False
    pending_close: Optional[str] = None

    async for chunk in chunk_iter:
        if not chunk:
            continue
        content_buffer += chunk

        while True:
            if TAG_SUGGEST_OPEN in content_buffer and TAG_SUGGEST_CLOSE in content_buffer:
                io = content_buffer.find(TAG_SUGGEST_OPEN)
                ic = content_buffer.find(TAG_SUGGEST_CLOSE)
                if io != -1 and ic != -1 and io < ic:
                    ic_full = ic + len(TAG_SUGGEST_CLOSE)
                    inner = content_buffer[io + len(TAG_SUGGEST_OPEN) : ic].strip()
                    qs = parse_suggest_questions_json(inner)
                    if qs:
                        yield ("suggest", {"questions": qs})
                    content_buffer = content_buffer[:io] + content_buffer[ic_full:]
                    continue

            if SUGGEST_END in content_buffer and SUGGEST_START in content_buffer:
                idx_start = content_buffer.find(SUGGEST_START)
                idx_end = content_buffer.find(SUGGEST_END) + len(SUGGEST_END)
                if idx_start < idx_end:
                    try:
                        json_str = content_buffer[
                            idx_start + len(SUGGEST_START) : idx_end - len(SUGGEST_END)
                        ].strip()
                        suggestions = json.loads(json_str)
                        if isinstance(suggestions, list) and suggestions:
                            yield ("suggestions", suggestions)
                    except (json.JSONDecodeError, TypeError):
                        pass
                    content_buffer = content_buffer[:idx_start] + content_buffer[idx_end:]
                    continue

            if in_think_tag:
                pc = pending_close or ""
                pcl = len(pc)
                if pcl and pc in content_buffer:
                    idx = content_buffer.find(pc)
                    content_buffer = content_buffer[idx + pcl :]
                    in_think_tag = False
                    pending_close = None
                    continue
                if pcl > 0 and len(content_buffer) > pcl:
                    content_buffer = content_buffer[-pcl:]
                break

            hit = _find_earliest_think_open(content_buffer)
            if hit is not None:
                idx, open_m, close_m = hit
                message_part = content_buffer[:idx]
                if message_part:
                    yield ("message", {"content": message_part})
                content_buffer = content_buffer[idx + len(open_m) :]
                in_think_tag = True
                pending_close = close_m
                continue

            if TAG_SUGGEST_OPEN in content_buffer:
                io = content_buffer.find(TAG_SUGGEST_OPEN)
                ic_rel = content_buffer.find(TAG_SUGGEST_CLOSE, io + len(TAG_SUGGEST_OPEN))
                if ic_rel == -1:
                    msg_part = content_buffer[:io]
                    if msg_part:
                        yield ("message", {"content": msg_part})
                    content_buffer = content_buffer[io:]
                    break
            if SUGGEST_START in content_buffer:
                safe_end = content_buffer.find(SUGGEST_START)
                to_yield = content_buffer[:safe_end]
                if to_yield:
                    yield ("message", {"content": to_yield})
                content_buffer = content_buffer[safe_end:]
                break
            if THINK_MAX_OPEN_LEN > 0 and len(content_buffer) > THINK_MAX_OPEN_LEN:
                to_yield = content_buffer[:-THINK_MAX_OPEN_LEN]
                content_buffer = content_buffer[-THINK_MAX_OPEN_LEN:]
                if to_yield:
                    yield ("message", {"content": to_yield})
                break

            break

    if content_buffer:
        if in_think_tag:
            return
        remaining = content_buffer
        if TAG_SUGGEST_OPEN in remaining:
            remaining = remaining.split(TAG_SUGGEST_OPEN, 1)[0]
        remaining = remaining.split(SUGGEST_START, 1)[0]
        if remaining:
            yield ("message", {"content": remaining})


def _strip_json_code_fences(text: str) -> str:
    """去除 ``` / ```json 包裹，供流结束后 json.loads 使用。"""
    s = (text or "").strip()
    if not s:
        return s
    # 整块 ```json ... ``` 或 ``` ... ```
    m = re.match(r"^```(?:json)?\s*\r?\n?(.*)\r?\n?```\s*$", s, re.DOTALL | re.IGNORECASE)
    if m:
        return m.group(1).strip()
    if s.startswith("```"):
        lines = s.split("\n", 1)
        s = lines[1] if len(lines) > 1 else ""
        s = re.sub(r"\n?```\s*$", "", s, flags=re.DOTALL).strip()
    return s


def _extract_balanced_json_object(text: str) -> Optional[str]:
    """从文本中切出第一个花括号平衡的 JSON 对象子串（字符串内括号不计入深度）。"""
    start = text.find("{")
    if start < 0:
        return None
    depth = 0
    in_string = False
    escape = False
    for i in range(start, len(text)):
        c = text[i]
        if in_string:
            if escape:
                escape = False
            elif c == "\\":
                escape = True
            elif c == '"':
                in_string = False
            continue
        if c == '"':
            in_string = True
            continue
        if c == "{":
            depth += 1
        elif c == "}":
            depth -= 1
            if depth == 0:
                return text[start : i + 1]
    return None


def _extract_balanced_json_array(text: str) -> Optional[str]:
    """从文本中切出第一个方括号平衡的 JSON 数组子串。"""
    start = text.find("[")
    if start < 0:
        return None
    depth = 0
    in_string = False
    escape = False
    for i in range(start, len(text)):
        c = text[i]
        if in_string:
            if escape:
                escape = False
            elif c == "\\":
                escape = True
            elif c == '"':
                in_string = False
            continue
        if c == '"':
            in_string = True
            continue
        if c == "[":
            depth += 1
        elif c == "]":
            depth -= 1
            if depth == 0:
                return text[start : i + 1]
    return None


def _slice_json_by_brace_fallback(text: str) -> Optional[str]:
    """首 `{` 与末 `}` 之间的暴力切片（模型夹带自然语言时的最后手段）。"""
    t = text or ""
    start_idx = t.find("{")
    end_idx = t.rfind("}")
    if start_idx != -1 and end_idx != -1 and end_idx > start_idx:
        return t[start_idx : end_idx + 1]
    return None


def _slice_json_by_bracket_fallback(text: str) -> Optional[str]:
    """首 `[` 与末 `]` 之间的暴力切片（根为数组时）。"""
    t = text or ""
    start_idx = t.find("[")
    end_idx = t.rfind("]")
    if start_idx != -1 and end_idx != -1 and end_idx > start_idx:
        return t[start_idx : end_idx + 1]
    return None


def _prepare_json_text_for_loads(cleaned: str) -> List[str]:
    """
    生成若干候选子串供 json.loads 依次尝试。
    若首个 `{` 在首个 `[` 之前（典型「废话 + JSON 对象」），**不**加入根级 array 候选，
    避免残缺 object 内的 `\"key\": []` 被误解析为根 `[]`。
    纯数组根（如 `[{...}]`）则保留 array 候选。
    """
    t = (cleaned or "").strip()
    if not t:
        return []
    seen = set()
    out: List[str] = []

    def _add(s: Optional[str]) -> None:
        if not s:
            return
        s = s.strip()
        if s and s not in seen:
            seen.add(s)
            out.append(s)

    i_brace = t.find("{")
    i_bracket = t.find("[")
    object_first = i_brace != -1 and (i_bracket == -1 or i_brace < i_bracket)

    _add(_extract_balanced_json_object(t))
    if not object_first:
        _add(_extract_balanced_json_array(t))
    _add(_slice_json_by_brace_fallback(t))
    if not object_first:
        _add(_slice_json_by_bracket_fallback(t))
    _add(t)
    return out


JsonExtractResult = Union[dict, list]


async def stream_and_extract_json(
    chunk_iter: AsyncIterator[Any],
) -> AsyncIterator[Tuple[str, object]]:
    """
    Planner JSON 提取：剥离 delta.reasoning_content 与各类思考标签内字节（不外发）；仅向外 yield json/json_error。
    """

    is_thinking = False
    pending_close: Optional[str] = None
    json_buffer = ""
    work_buf = ""

    async for chunk in chunk_iter:
        try:
            choices = getattr(chunk, "choices", None) if not isinstance(chunk, dict) else chunk.get("choices")
            if not choices or len(choices) == 0:
                continue
            first = choices[0]
            delta = getattr(first, "delta", None) if not isinstance(first, dict) else first.get("delta")

            reasoning = _delta_get(delta, "reasoning_content")
            if reasoning:
                logger.debug(
                    "stream_and_extract_json: dropped reasoning_content (~%d chars)",
                    len(reasoning),
                )

            content = _delta_get(delta, "content")
            if not content:
                continue

            work_buf += content

            while True:
                if is_thinking:
                    pc = pending_close or ""
                    pcl = len(pc)
                    if pcl and pc in work_buf:
                        idx = work_buf.find(pc)
                        work_buf = work_buf[idx + pcl :]
                        is_thinking = False
                        pending_close = None
                        continue
                    if pcl > 0 and len(work_buf) > pcl:
                        work_buf = work_buf[-pcl:]
                    break

                hit = _find_earliest_think_open(work_buf)
                if hit is not None:
                    idx, open_m, close_m = hit
                    before = work_buf[:idx]
                    if before:
                        json_buffer += before
                    work_buf = work_buf[idx + len(open_m) :]
                    is_thinking = True
                    pending_close = close_m
                    continue

                if THINK_MAX_OPEN_LEN > 0 and len(work_buf) > THINK_MAX_OPEN_LEN:
                    safe = work_buf[:-THINK_MAX_OPEN_LEN]
                    work_buf = work_buf[-THINK_MAX_OPEN_LEN:]
                    if safe:
                        json_buffer += safe
                    break

                break
        except Exception as e:
            logger.exception("stream_and_extract_json: chunk 处理异常，已转为 json_error: %s", e)
            yield (
                "json_error",
                {
                    "error": "stream_chunk_error",
                    "message": str(e),
                    "raw_preview": (json_buffer + work_buf)[:2000],
                },
            )
            return

    # 流结束：冲刷 work_buf（思考半成品不入库）
    if work_buf:
        if not is_thinking:
            json_buffer += work_buf
        work_buf = ""

    cleaned = _strip_json_code_fences(json_buffer)
    if not cleaned:
        yield (
            "json_error",
            {
                "error": "empty_json_buffer",
                "message": "流结束后 JSON 缓冲区为空（模型未输出可解析正文或仅含思考标签）",
                "raw_preview": (json_buffer or "")[:2000],
            },
        )
        return

    candidates = _prepare_json_text_for_loads(cleaned)
    last_decode_error: Optional[json.JSONDecodeError] = None
    for cand in candidates:
        try:
            parsed: JsonExtractResult = json.loads(cand)
            if not isinstance(parsed, (dict, list)):
                yield (
                    "json_error",
                    {
                        "error": "json_type_unexpected",
                        "message": f"根类型必须为 object 或 array，实际: {type(parsed).__name__}",
                        "raw_preview": cleaned[:2000],
                    },
                )
                return
            # 残缺 object 中常含 "key": [] 片段；勿将内层 [] 误当根 JSON（否则 Planner 会误判为成功、阻断重试）
            if isinstance(parsed, list) and not cand.lstrip().startswith("["):
                continue
            yield ("json", parsed)
            return
        except json.JSONDecodeError as e:
            last_decode_error = e
            continue
        except Exception as e:
            logger.exception("stream_and_extract_json: json.loads 意外异常: %s", e)
            yield (
                "json_error",
                {
                    "error": "parse_unexpected_error",
                    "message": str(e),
                    "raw_preview": cleaned[:2000],
                },
            )
            return

    yield (
        "json_error",
        {
            "error": "json_decode_error",
            "message": str(last_decode_error) if last_decode_error else "所有候选子串均无法解析为 JSON",
            "raw_preview": cleaned[:2000],
        },
    )
