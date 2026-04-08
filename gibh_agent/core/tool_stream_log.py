"""
工具执行期流式日志（基建）：通过 ContextVar 将可选 sink 注入到同步/异步工具体内，
供 emit_tool_log 将中间态文案交给 Orchestrator，再格式化为 SSE status 下行。

- 未设置 sink 时 emit_tool_log 为 no-op，旧工具零改动即可运行。
- sink 签名: (content: str, state: str) -> None，须线程安全（如 queue.Queue.put）。
"""
from __future__ import annotations

import asyncio
import contextvars
import logging
from contextlib import contextmanager
from typing import Callable, Iterator, Optional

logger = logging.getLogger(__name__)

ToolLogSink = Callable[[str, str], None]

_tool_log_sink: contextvars.ContextVar[Optional[ToolLogSink]] = contextvars.ContextVar(
    "gibh_tool_log_sink", default=None
)


def emit_tool_log(content: str, *, state: str = "running") -> None:
    """
    在工具内部调用：将一行用户可见过程说明交给当前请求绑定的 sink。
    同步、非阻塞；若未绑定 sink 则忽略。
    """
    if not content or not str(content).strip():
        return
    sink = _tool_log_sink.get()
    if sink is None:
        return
    try:
        sink(str(content).strip(), str(state or "running"))
    except Exception:
        logger.exception("tool_log sink raised; continuing tool execution")


async def aemit_tool_log(content: str, *, state: str = "running") -> None:
    """异步工具内可选用：与 emit_tool_log 等价，多一次 await 以便与 async 代码风格一致。"""
    emit_tool_log(content, state=state)
    await asyncio.sleep(0)


@contextmanager
def tool_log_emitter_scope(sink: Optional[ToolLogSink]) -> Iterator[None]:
    """在 execute_step 内包裹单次工具调用，绑定/解绑 ContextVar。"""
    token = _tool_log_sink.set(sink)
    try:
        yield
    finally:
        _tool_log_sink.reset(token)


def get_bound_tool_log_sink() -> Optional[ToolLogSink]:
    """测试或诊断用：读取当前上下文是否已绑定 sink。"""
    return _tool_log_sink.get()
