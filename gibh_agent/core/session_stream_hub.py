"""
会话级 SSE 事件中枢：主连接发布增量，切回 running 会话时可重放缓冲并订阅后续事件。
"""
from __future__ import annotations

import asyncio
import logging
from collections import deque
from typing import AsyncIterator, Deque, Dict, List, Optional

logger = logging.getLogger(__name__)

_BUFFER_MAX = 800
_CLOSE = object()


class SessionStreamHub:
    def __init__(self) -> None:
        self._buffers: Dict[str, Deque[str]] = {}
        self._subscribers: Dict[str, List[asyncio.Queue]] = {}
        self._lock = asyncio.Lock()

    async def publish(self, session_id: str, event_chunk: str) -> None:
        sid = (session_id or "").strip()
        if not sid or not event_chunk:
            return
        async with self._lock:
            buf = self._buffers.setdefault(sid, deque(maxlen=_BUFFER_MAX))
            buf.append(event_chunk)
            queues = list(self._subscribers.get(sid, []))
        for q in queues:
            try:
                q.put_nowait(event_chunk)
            except asyncio.QueueFull:
                logger.warning("[SessionStreamHub] subscriber queue full session=%s", sid)

    async def close_session(self, session_id: str) -> None:
        sid = (session_id or "").strip()
        if not sid:
            return
        async with self._lock:
            queues = list(self._subscribers.pop(sid, []))
        for q in queues:
            try:
                q.put_nowait(_CLOSE)
            except asyncio.QueueFull:
                pass

    def snapshot_buffer(self, session_id: str) -> List[str]:
        sid = (session_id or "").strip()
        buf = self._buffers.get(sid)
        if not buf:
            return []
        return list(buf)

    async def subscribe(self, session_id: str) -> AsyncIterator[str]:
        sid = (session_id or "").strip()
        if not sid:
            return
        q: asyncio.Queue = asyncio.Queue(maxsize=256)
        async with self._lock:
            self._subscribers.setdefault(sid, []).append(q)
        for chunk in self.snapshot_buffer(sid):
            yield chunk
        try:
            while True:
                item = await q.get()
                if item is _CLOSE:
                    break
                yield str(item)
        finally:
            async with self._lock:
                subs = self._subscribers.get(sid, [])
                if q in subs:
                    subs.remove(q)


_hub: Optional[SessionStreamHub] = None


def get_session_stream_hub() -> SessionStreamHub:
    global _hub
    if _hub is None:
        _hub = SessionStreamHub()
    return _hub
