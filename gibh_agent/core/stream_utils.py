"""
Stream and text parsing utilities for LLM responses.
- 通用双轨（无模型硬编码）：reasoning_content → thought；content 内 <think>...</think> → 状态机分流，缓冲区防跨 chunk 截断；<<<SUGGESTIONS>>> 解析不变。
- stream_and_extract_json：Planner 流式升级专用分拣器——thought 与 JSON 正文严格隔离；流结束后先剥 fence，再按平衡括号或首 `{`/末 `}` 暴力切出 JSON，最后 json.loads。
"""
import json
import re
import logging
from typing import AsyncIterator, List, Optional, Tuple, Any, Union

logger = logging.getLogger(__name__)

SUGGEST_START = "<<<SUGGESTIONS>>>"
SUGGEST_END = "<<<END_SUGGESTIONS>>>"
THINK_OPEN = "<think>"
THINK_CLOSE = "</think>"
THINK_OPEN_LEN = len(THINK_OPEN)   # 7，用于跨 chunk 截断时保留尾部
THINK_CLOSE_LEN = len(THINK_CLOSE)  # 8


def _delta_get(delta: Any, key: str) -> Optional[str]:
    """从 delta 安全取值，兼容对象（getattr）与字典（.get）。"""
    if delta is None:
        return None
    if isinstance(delta, dict):
        val = delta.get(key)
        return str(val) if val is not None else None
    val = getattr(delta, key, None)
    return str(val) if val is not None else None


def strip_suggestions_from_text(text: str) -> Tuple[str, Optional[List[str]]]:
    """
    Remove the suggestions block from a full LLM response and return cleaned text + parsed list.
    Use for non-streaming responses (diagnosis, report).

    Args:
        text: Full response string (may contain <<<SUGGESTIONS>>>["Q1","Q2"]<<<END_SUGGESTIONS>>>).

    Returns:
        (cleaned_text, suggestions_list or None). If no block found, cleaned_text is text unchanged.
    """
    if not text or SUGGEST_START not in text or SUGGEST_END not in text:
        return (text, None)
    idx_start = text.find(SUGGEST_START)
    idx_end = text.find(SUGGEST_END) + len(SUGGEST_END)
    before = text[:idx_start].rstrip()
    after = text[idx_end:].lstrip()
    cleaned = (before + "\n\n" + after).strip() if after else before
    try:
        json_str = text[idx_start + len(SUGGEST_START) : idx_end - len(SUGGEST_END)].strip()
        suggestions = json.loads(json_str)
        if isinstance(suggestions, list) and suggestions:
            return (cleaned, suggestions)
    except (json.JSONDecodeError, TypeError):
        pass
    return (cleaned, None)


async def stream_from_llm_chunks(
    chunk_iter: AsyncIterator[Any],
    model_name: Optional[str] = None,
) -> AsyncIterator[Tuple[str, object]]:
    """
    通用双轨状态机（鸭子类型，不依赖 model_name）：
    1) reasoning_content 有值 → 直接 yield ("thought", ...)
    2) content 用缓冲区 + in_think_tag：<think> 内 → thought，否则 → message；
       通过保留尾部 THINK_OPEN_LEN/THINK_CLOSE_LEN 字符应对跨 chunk 截断（如 <th + ink>）。
    3) <<<SUGGESTIONS>>> 块：完整块时解析并 yield，前面未消费部分保留在 buffer 中继续走 <think> 状态机。
    """
    content_buffer = ""
    in_think_tag = False

    async for chunk in chunk_iter:
        choices = getattr(chunk, "choices", None) if not isinstance(chunk, dict) else chunk.get("choices")
        if not choices or len(choices) == 0:
            continue
        first = choices[0]
        delta = getattr(first, "delta", None) if not isinstance(first, dict) else first.get("delta")

        reasoning = _delta_get(delta, "reasoning_content")
        if reasoning:
            yield ("thought", {"content": reasoning})

        content = _delta_get(delta, "content")
        if content:
            content_buffer += content

        while True:
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
                if THINK_CLOSE in content_buffer:
                    idx = content_buffer.find(THINK_CLOSE)
                    thought_text = content_buffer[:idx]
                    if thought_text:
                        yield ("thought", {"content": thought_text})
                    content_buffer = content_buffer[idx + THINK_CLOSE_LEN:]
                    in_think_tag = False
                    continue
                if len(content_buffer) > THINK_CLOSE_LEN:
                    safe = content_buffer[:-THINK_CLOSE_LEN]
                    content_buffer = content_buffer[-THINK_CLOSE_LEN:]
                    if safe:
                        yield ("thought", {"content": safe})
                break
            else:
                if THINK_OPEN in content_buffer:
                    idx = content_buffer.find(THINK_OPEN)
                    message_part = content_buffer[:idx]
                    if message_part:
                        yield ("message", {"content": message_part})
                    content_buffer = content_buffer[idx + THINK_OPEN_LEN:]
                    in_think_tag = True
                    continue
                if SUGGEST_START in content_buffer:
                    safe_end = content_buffer.find(SUGGEST_START)
                    to_yield = content_buffer[:safe_end]
                    if to_yield:
                        yield ("message", {"content": to_yield})
                    content_buffer = content_buffer[safe_end:]
                    break
                if len(content_buffer) > THINK_OPEN_LEN:
                    to_yield = content_buffer[:-THINK_OPEN_LEN]
                    content_buffer = content_buffer[-THINK_OPEN_LEN:]
                    if to_yield:
                        yield ("message", {"content": to_yield})
                break

    if content_buffer:
        if in_think_tag:
            yield ("thought", {"content": content_buffer})
        else:
            remaining = content_buffer.split(SUGGEST_START)[0]
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
        remaining = stream_buffer[last_yielded:].split(SUGGEST_START)[0]
        if remaining:
            yield ("message", {"content": remaining})


async def stream_with_thought_and_suggestions(
    chunk_iter: AsyncIterator[str],
) -> AsyncIterator[Tuple[str, object]]:
    """
    解析 LLM 流：<think> 内内容作为 event: thought 实时推送，其余作为 message；
    同时处理 <<<SUGGESTIONS>>> 块。用于 DeepSeek-R1 等 CoT 模型。
    """
    buffer = ""
    in_think = False
    last_thought_pos = 0
    last_message_pos = 0

    async for chunk in chunk_iter:
        if not chunk:
            continue
        buffer += chunk

        while True:
            if in_think:
                if THINK_CLOSE in buffer:
                    idx = buffer.find(THINK_CLOSE)
                    thought_content = buffer[last_thought_pos:idx]
                    if thought_content:
                        yield ("thought", {"content": thought_content})
                    buffer = buffer[idx + len(THINK_CLOSE):]
                    last_thought_pos = 0
                    last_message_pos = 0
                    in_think = False
                    continue
                else:
                    new_part = buffer[last_thought_pos:]
                    if new_part:
                        yield ("thought", {"content": new_part})
                    last_thought_pos = len(buffer)
                    break
            else:
                if THINK_OPEN in buffer:
                    idx = buffer.find(THINK_OPEN)
                    before = buffer[last_message_pos:idx]
                    if before:
                        yield ("message", {"content": before})
                    buffer = buffer[idx + len(THINK_OPEN):]
                    last_message_pos = 0
                    last_thought_pos = 0
                    in_think = True
                    continue
                else:
                    if SUGGEST_END in buffer and SUGGEST_START in buffer:
                        idx_s = buffer.find(SUGGEST_START)
                        idx_e = buffer.find(SUGGEST_END) + len(SUGGEST_END)
                        if buffer[last_message_pos:idx_s]:
                            yield ("message", {"content": buffer[last_message_pos:idx_s]})
                        try:
                            json_str = buffer[idx_s + len(SUGGEST_START) : idx_e - len(SUGGEST_END)].strip()
                            suggestions = json.loads(json_str)
                            if isinstance(suggestions, list) and suggestions:
                                yield ("suggestions", suggestions)
                        except (json.JSONDecodeError, TypeError):
                            pass
                        buffer = buffer[idx_e:]
                        last_message_pos = 0
                        continue
                    if SUGGEST_START in buffer:
                        safe = buffer.find(SUGGEST_START)
                        if buffer[last_message_pos:safe]:
                            yield ("message", {"content": buffer[last_message_pos:safe]})
                        last_message_pos = safe
                    else:
                        if buffer[last_message_pos:]:
                            yield ("message", {"content": buffer[last_message_pos:]})
                        last_message_pos = len(buffer)
                    break

    if in_think and last_thought_pos < len(buffer):
        yield ("thought", {"content": buffer[last_thought_pos:]})
    elif not in_think and last_message_pos < len(buffer):
        remaining = buffer[last_message_pos:].split(SUGGEST_START)[0]
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
    流式状态机分拣器（Stateful Stream Dispatcher）：与 Planner 业务解耦，专用于「思考实时透出 + JSON 尾部结算」。

    行为摘要：
    1. 每个 chunk：若存在 delta.reasoning_content，立即 yield ("thought", {"content": ...})（原生字段优先）。
    2. delta.content：维护 is_thinking；`</think>` 与 `</think>` 之间仅 yield thought，之外仅累积 json_buffer，禁止交叉污染。
    3. 流结束：清理 fence 后，按平衡 {} / [] 或首末括号暴力切片生成候选串依次 json.loads；
       成功 yield ("json", obj)；均失败 yield ("json_error", {...})，不抛未捕获异常。

    Yields:
        ("thought", {"content": str})  — 推理片段（可多次）
        ("json", dict | list)           — 解析成功（通常一次）
        ("json_error", dict)            — 结构化错误，含 error / message / raw_preview
    """
    is_thinking = False
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
                yield ("thought", {"content": reasoning})

            content = _delta_get(delta, "content")
            if not content:
                continue

            work_buf += content

            while True:
                if is_thinking:
                    if THINK_CLOSE in work_buf:
                        idx = work_buf.find(THINK_CLOSE)
                        thought_seg = work_buf[:idx]
                        if thought_seg:
                            yield ("thought", {"content": thought_seg})
                        work_buf = work_buf[idx + THINK_CLOSE_LEN :]
                        is_thinking = False
                        continue
                    if len(work_buf) > THINK_CLOSE_LEN:
                        safe = work_buf[:-THINK_CLOSE_LEN]
                        work_buf = work_buf[-THINK_CLOSE_LEN:]
                        if safe:
                            yield ("thought", {"content": safe})
                    break
                else:
                    if THINK_OPEN in work_buf:
                        idx = work_buf.find(THINK_OPEN)
                        before = work_buf[:idx]
                        if before:
                            json_buffer += before
                        work_buf = work_buf[idx + THINK_OPEN_LEN :]
                        is_thinking = True
                        continue
                    if len(work_buf) > THINK_OPEN_LEN:
                        safe = work_buf[:-THINK_OPEN_LEN]
                        work_buf = work_buf[-THINK_OPEN_LEN:]
                        if safe:
                            json_buffer += safe
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

    # 流结束：冲刷 work_buf
    if work_buf:
        if is_thinking:
            yield ("thought", {"content": work_buf})
        else:
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
