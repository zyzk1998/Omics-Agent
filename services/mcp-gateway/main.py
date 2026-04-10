# -*- coding: utf-8 -*-
"""MCP Gateway：唯一持有官方 mcp SDK，向主 API 暴露 HTTP /tools /call /reconnect。"""
from __future__ import annotations

import asyncio
import logging
import os
import time
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


def _format_connection_exception(exc: BaseException) -> str:
    """
    将连接失败时的异常链压成一段字符串，供 503 Details 透传。
    MCP streamable_http + anyio TaskGroup 失败时抛出 ExceptionGroup：单独 str(e) 往往只有
    「unhandled errors in a TaskGroup (1 sub-exception)」，必须展开 .exceptions 才有
    httpx.ConnectError / httpcore.ConnectError。
    优先用 PEP 654 的 ``exceptions`` 元组判断，避免仅依赖 isinstance（边缘环境更稳）。
    """
    subs = getattr(exc, "exceptions", None)
    if isinstance(subs, tuple) and len(subs) > 0:
        try:
            sub_txt = [f"sub[{i}] {_format_connection_exception(x)}" for i, x in enumerate(subs)]
            return f"{type(exc).__name__}: {exc!s} || " + " || ".join(sub_txt)
        except Exception:  # noqa: BLE001
            return f"{type(exc).__name__}: {exc!s} || sub_errors_repr={subs!r}"
    if isinstance(exc, BaseExceptionGroup):
        return f"{type(exc).__name__}: {exc!s}"
    parts: List[str] = [f"{type(exc).__name__}: {exc!s}"]
    cur: Optional[BaseException] = exc.__cause__
    depth = 0
    while cur is not None and depth < 8:
        parts.append(f"caused_by {type(cur).__name__}: {cur!s}")
        cur = cur.__cause__
        depth += 1
    curx: Optional[BaseException] = exc.__context__
    if curx is not None and curx is not exc.__cause__:
        parts.append(f"context {type(curx).__name__}: {curx!s}")
    return " || ".join(parts)


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
        self._tools: List[Dict[str, Any]] = []
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

    def _tool_to_payload(self, tool: Any) -> Optional[Dict[str, Any]]:
        """MCP Tool → {name, description, inputSchema}，供主进程动态注入 LLM tools schema。"""
        schema: Dict[str, Any] = {}
        if isinstance(tool, dict):
            name = (tool.get("name") or "").strip()
            if not name:
                return None
            desc = tool.get("description") or ""
            raw_s = tool.get("inputSchema") or tool.get("input_schema")
            if isinstance(raw_s, dict):
                schema = raw_s
            return {"name": name, "description": str(desc or "").strip(), "inputSchema": schema}
        name = getattr(tool, "name", None) or ""
        if not name:
            return None
        desc = getattr(tool, "description", None) or ""
        for attr in ("inputSchema", "input_schema"):
            v = getattr(tool, attr, None)
            if isinstance(v, dict):
                schema = v
                break
        if not schema and hasattr(tool, "model_dump"):
            try:
                dumped = tool.model_dump()
                s = dumped.get("inputSchema") or dumped.get("input_schema")
                if isinstance(s, dict):
                    schema = s
            except Exception:  # noqa: BLE001
                pass
        return {"name": str(name).strip(), "description": str(desc or "").strip(), "inputSchema": schema}

    def _set_tools_from_list_response(self, listed: Any) -> None:
        raw = getattr(listed, "tools", None) or []
        out: List[Dict[str, Any]] = []
        for tool in raw:
            row = self._tool_to_payload(tool)
            if row:
                out.append(row)
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
            # 必须保留类型名 + 原文 + 异常链；仅 str(e) 在 ExceptionGroup 场景下会丢失子异常
            self._last_error = _format_connection_exception(e)
            # 不在此条附带 exc_info：否则默认 Formatter 会先打出一行与 str(e) 雷同的废话，掩盖已展开的 _last_error
            logger.warning("[mcp-gateway] connection failed: %s", self._last_error)
            logger.debug("[mcp-gateway] connection failed traceback", exc_info=True)
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

    async def wait_until_connected(self, timeout: float) -> bool:
        """在 timeout 秒内等待 connection_worker 把会话拉起来；失败或 worker 已退出则返回 False。"""
        deadline = time.monotonic() + max(0.1, float(timeout))
        while time.monotonic() < deadline:
            if self.connected and self._session is not None:
                return True
            t = self._task
            if t is not None and t.done():
                return False
            await asyncio.sleep(0.25)
        return bool(self.connected and self._session is not None)

    def describe_session_unavailable_detail(
        self,
        *,
        wait_timed_out: bool = False,
        wait_timeout_sec: Optional[float] = None,
    ) -> str:
        """
        503 返回给主 API / 大模型的唯一事实源。
        格式固定为：MCP session not connected. Details: <异常原文或等价说明> | <运维附加上下文>
        """
        err = (self._last_error or "").strip()
        if not err:
            if wait_timed_out and wait_timeout_sec is not None:
                wr = self._task is not None and not self._task.done()
                err = (
                    f"Timeout waiting for MCP session after {wait_timeout_sec:g}s "
                    f"(connection_worker_running={wr}; no Exception reached except-block yet — "
                    f"often means hung TCP/handshake to upstream)"
                )
            else:
                err = "no Exception text recorded on gateway (session never became ready)"
        t = self._task
        worker = "not_started" if t is None else ("finished_or_crashed" if t.done() else "running")
        # 主句：用户要求的固定前缀 + 真实异常链原文
        detail = f"MCP session not connected. Details: {err}"
        suffix = f" | upstream_url={self.url} | transport={self.transport} | worker_state={worker}"
        el = err.lower()
        if "refused" in el or "connection refused" in el:
            suffix += " | hint=tcp_connection_refused"
        elif "timed out" in el or "timeout" in el:
            suffix += " | hint=timeout_or_handshake_stall"
        elif "name or service not known" in el or "getaddrinfo" in el or "nodename nor servname" in el:
            suffix += " | hint=dns_resolution_failed"
        elif "certificate" in el or "ssl" in el or "tls" in el:
            suffix += " | hint=tls_or_cert"
        detail = detail + suffix
        max_len = int(os.getenv("HPC_MCP_503_DETAIL_MAX", "2500"))
        if len(detail) > max_len:
            detail = detail[: max_len - 1] + "…"
        return detail

    async def call_tool(self, tool_name: str, arguments: Dict[str, Any]) -> Dict[str, Any]:
        if self._task is None or (self._task.done() and not self.connected):
            await self.start()
        if not self.connected or not self._session:
            wait_sec = float(os.getenv("HPC_MCP_CALL_READY_TIMEOUT", "120"))
            ok = await self.wait_until_connected(wait_sec)
            if not ok or not self.connected or not self._session:
                raise HTTPException(
                    status_code=503,
                    detail=self.describe_session_unavailable_detail(
                        wait_timed_out=not ok,
                        wait_timeout_sec=wait_sec if not ok else None,
                    ),
                )
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
    startup_wait = float(os.getenv("HPC_MCP_STARTUP_WAIT_TIMEOUT", "90"))
    if startup_wait > 0:
        ok = await state.wait_until_connected(startup_wait)
        if not ok:
            logger.warning(
                "[mcp-gateway] 启动后 %.0fs 内未连上上游 MCP（请检查 HPC_MCP_URL / 网络）: %s",
                startup_wait,
                state._last_error,
            )
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
        raise HTTPException(
            status_code=503,
            detail="MCP gateway not initialized (FastAPI app.state.gateway missing).",
        )
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
        raise HTTPException(
            status_code=503,
            detail="MCP gateway not initialized (cannot POST /call: app.state.gateway missing).",
        )
    name = (body.tool_name or "").strip()
    if not name:
        raise HTTPException(status_code=400, detail="tool_name required")
    return await g.call_tool(name, body.arguments or {})
