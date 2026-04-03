# -*- coding: utf-8 -*-
"""MCP Gateway：唯一持有官方 mcp SDK，向主 API 暴露 HTTP /tools /call /reconnect。"""
from __future__ import annotations

import asyncio
import logging
import os
from contextlib import asynccontextmanager
from typing import Any, Dict, List, Optional

from fastapi import FastAPI, HTTPException
from pydantic import BaseModel, Field

logger = logging.getLogger(__name__)

try:
    from mcp import ClientSession
    from mcp.client.sse import sse_client
    from mcp.client.streamable_http import streamable_http_client
except ImportError as e:  # pragma: no cover
    raise RuntimeError("mcp SDK is required in mcp-gateway image") from e


def _normalize_transport(transport: str) -> str:
    t = (transport or "").strip().lower().replace("_", "").replace("-", "")
    if t in ("streamablehttp", "http"):
        return "streamablehttp"
    if t == "sse":
        return "sse"
    return t or "streamablehttp"


def _call_result_to_json(result: Any) -> Dict[str, Any]:
    is_err = bool(getattr(result, "isError", None) or getattr(result, "is_error", None))
    chunks: List[Dict[str, Any]] = []
    for block in getattr(result, "content", None) or []:
        txt = getattr(block, "text", None)
        if txt is not None:
            chunks.append({"type": "text", "text": txt})
        else:
            chunks.append({"type": type(block).__name__, "data": str(block)})
    out: Dict[str, Any] = {
        "status": "error" if is_err else "success",
        "content": chunks,
    }
    if is_err:
        detail_parts: List[str] = []
        for c in chunks:
            if isinstance(c, dict) and c.get("type") == "text" and c.get("text"):
                detail_parts.append(str(c["text"]))
        detail = "\n".join(detail_parts).strip()
        if detail:
            max_len = int(os.getenv("HPC_MCP_ERROR_DETAIL_MAX", "1200"))
            clipped = detail if len(detail) <= max_len else detail[:max_len] + "…"
            out["message"] = f"远端工具返回 isError：{clipped}"
        else:
            out["message"] = "远端工具返回 isError"
    return out


class ReconnectBody(BaseModel):
    url: str = Field(default="http://192.168.32.31:8001/mcp")
    transport: str = Field(default="streamableHttp")


class CallBody(BaseModel):
    tool_name: str
    arguments: Dict[str, Any] = Field(default_factory=dict)


class GatewayState:
    def __init__(self) -> None:
        self._task: Optional[asyncio.Task] = None
        self._session: Any = None
        self._invoke_lock = asyncio.Lock()
        self._connected: bool = False
        self._last_error: Optional[str] = None
        self._tools: List[Dict[str, str]] = []
        self.url: str = (os.getenv("HPC_MCP_URL") or "http://192.168.32.31:8001/mcp").rstrip("/")
        self.transport: str = os.getenv("HPC_MCP_TRANSPORT") or "streamableHttp"

    @property
    def connected(self) -> bool:
        return self._connected and self._session is not None

    def status_dict(self) -> Dict[str, Any]:
        return {
            "connected": self.connected,
            "url": self.url,
            "transport": self.transport,
            "tools": [t["name"] for t in self._tools],
            "tool_details": list(self._tools),
            "last_error": self._last_error,
        }

    def _set_tools_from_list_response(self, listed: Any) -> None:
        raw = getattr(listed, "tools", None) or []
        out: List[Dict[str, str]] = []
        for tool in raw:
            name = getattr(tool, "name", None) or ""
            if not name:
                continue
            desc = getattr(tool, "description", None) or ""
            out.append({"name": name, "description": desc})
        self._tools = out

    async def _cancel_worker(self) -> None:
        t = self._task
        self._task = None
        if t and not t.done():
            t.cancel()
            try:
                await t
            except asyncio.CancelledError:
                pass
            except Exception as e:  # noqa: BLE001
                logger.debug("[mcp-gateway] join worker: %s", e)
        self._connected = False
        self._session = None
        self._tools = []

    async def connection_worker(self) -> None:
        self._connected = False
        self._session = None
        self._last_error = None

        nt = _normalize_transport(self.transport)
        url = self.url
        logger.info("[mcp-gateway] connecting url=%s transport=%s", url, nt)

        try:
            if nt == "streamablehttp":
                cm = streamable_http_client(url)
            elif nt == "sse":
                cm = sse_client(
                    url,
                    timeout=float(os.getenv("HPC_MCP_SSE_TIMEOUT", "30")),
                    sse_read_timeout=float(os.getenv("HPC_MCP_SSE_READ_TIMEOUT", "300")),
                )
            else:
                self._last_error = f"不支持的 transport: {self.transport}"
                logger.warning("[mcp-gateway] %s", self._last_error)
                return

            async with cm as streams:
                if nt == "streamablehttp":
                    read_stream, write_stream, _ = streams
                else:
                    read_stream, write_stream = streams
                async with ClientSession(read_stream, write_stream) as session:
                    await asyncio.wait_for(session.initialize(), timeout=30.0)
                    listed = await asyncio.wait_for(session.list_tools(), timeout=60.0)
                    self._session = session
                    self._set_tools_from_list_response(listed)
                    self._connected = True
                    self._last_error = None
                    logger.info("[mcp-gateway] connected, tools=%s", len(self._tools))
                    try:
                        await asyncio.Event().wait()
                    except asyncio.CancelledError:
                        raise
        except asyncio.CancelledError:
            logger.info("[mcp-gateway] connection task cancelled")
            raise
        except Exception as e:  # noqa: BLE001
            self._last_error = str(e)
            logger.warning("[mcp-gateway] connection failed: %s", e, exc_info=True)
        finally:
            self._connected = False
            self._session = None
            self._tools = []

    async def start(self) -> None:
        if self._task and not self._task.done():
            return
        self._task = asyncio.create_task(self.connection_worker(), name="mcp-gateway-connection")

    async def reconnect(self, url: str, transport: str) -> None:
        await self._cancel_worker()
        self.url = url.rstrip("/")
        self.transport = transport
        await self.start()

    async def disconnect(self) -> None:
        await self._cancel_worker()
        self._last_error = None

    async def call_tool(self, tool_name: str, arguments: Dict[str, Any]) -> Dict[str, Any]:
        if not self.connected or not self._session:
            raise HTTPException(status_code=503, detail="MCP session not connected")
        async with self._invoke_lock:
            try:
                res = await asyncio.wait_for(
                    self._session.call_tool(tool_name, arguments=arguments or {}),
                    timeout=float(os.getenv("HPC_MCP_CALL_TIMEOUT", "600")),
                )
                return _call_result_to_json(res)
            except asyncio.TimeoutError:
                raise HTTPException(status_code=504, detail="MCP call timeout") from None
            except HTTPException:
                raise
            except Exception as e:  # noqa: BLE001
                logger.exception("[mcp-gateway] call_tool failed: %s", tool_name)
                raise HTTPException(status_code=502, detail=str(e)) from e


state: Optional[GatewayState] = None


@asynccontextmanager
async def lifespan(app: FastAPI):
    global state
    logging.basicConfig(level=logging.INFO)
    state = GatewayState()
    app.state.gateway = state
    await state.start()
    yield
    if state:
        await state.disconnect()
        state = None


app = FastAPI(title="GIBH MCP Gateway", lifespan=lifespan)


@app.get("/health")
async def health():
    return {"status": "ok"}


@app.get("/status")
async def http_status():
    g = getattr(app.state, "gateway", None)
    if not g:
        return {"connected": False, "last_error": "gateway not initialized"}
    return g.status_dict()


@app.get("/tools")
async def http_tools():
    g = getattr(app.state, "gateway", None)
    if not g:
        return {"connected": False, "tools": []}
    return {"connected": g.connected, "tools": list(g._tools)}


@app.post("/reconnect")
async def http_reconnect(body: ReconnectBody):
    g = getattr(app.state, "gateway", None)
    if not g:
        raise HTTPException(status_code=503, detail="gateway not initialized")
    await g.reconnect(body.url.strip(), body.transport.strip())
    return {"status": "success", "message": "已调度网关重连"}


@app.post("/disconnect")
async def http_disconnect():
    g = getattr(app.state, "gateway", None)
    if not g:
        return {"status": "ok"}
    await g.disconnect()
    return {"status": "success", "message": "已断开 MCP"}


@app.post("/call")
async def http_call(body: CallBody):
    g = getattr(app.state, "gateway", None)
    if not g:
        raise HTTPException(status_code=503, detail="gateway not initialized")
    name = (body.tool_name or "").strip()
    if not name:
        raise HTTPException(status_code=400, detail="tool_name required")
    return await g.call_tool(name, body.arguments or {})
