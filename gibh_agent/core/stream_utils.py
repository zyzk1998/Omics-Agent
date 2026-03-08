"""
Stream and text parsing utilities for LLM responses.
Unified handling of <<<SUGGESTIONS>>>...<<<END_SUGGESTIONS>>> so no raw block reaches the frontend.
"""
import json
from typing import AsyncIterator, List, Optional, Tuple

SUGGEST_START = "<<<SUGGESTIONS>>>"
SUGGEST_END = "<<<END_SUGGESTIONS>>>"


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
