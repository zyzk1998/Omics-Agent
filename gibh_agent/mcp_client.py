# -*- coding: utf-8 -*-
"""
超算 MCP 客户端：经 mcp-gateway HTTP 转发，主进程不依赖官方 mcp 包。

- 连接失败仅记日志，不阻塞 FastAPI。
- 需网关可达（MCP_GATEWAY_URL）；工具列表来自 GET /tools，调用经 POST /call。
"""
from __future__ import annotations

import asyncio
import json
import logging
import os
import re
import time
from typing import Any, Callable, Dict, List, Optional, Tuple

import httpx

from gibh_agent.core.mcp_schema_adapter import mcp_tools_payload_to_openai_tools
from gibh_agent.core.tool_registry import registry
from gibh_agent.core.utils import safe_tool_execution

logger = logging.getLogger(__name__)


def _gateway_call_timeout() -> httpx.Timeout:
    """
    POST /call 专用超时（与 AsyncClient 默认长读分离），防止网关/上游挂死导致 SSE 假死。
    可通过 HPC_MCP_INVOKE_TIMEOUT_SEC 调整；默认 30s。
    """
    sec = max(1.0, float(os.getenv("HPC_MCP_INVOKE_TIMEOUT_SEC", "30")))
    connect = min(15.0, sec)
    write = min(30.0, sec)
    return httpx.Timeout(connect=connect, read=sec, write=write, pool=5.0)


def _args_preview_for_log(arguments: dict, max_len: int = 400) -> str:
    try:
        s = json.dumps(arguments, ensure_ascii=False, default=str)
    except Exception:  # noqa: BLE001
        s = str(arguments)
    if len(s) > max_len:
        return s[: max_len - 3] + "..."
    return s


def _registry_name(mcp_tool_name: str, used: set) -> str:
    base = re.sub(r"[^a-zA-Z0-9_]+", "_", mcp_tool_name).strip("_") or "tool"
    name = f"hpc_mcp_{base}"
    if name not in used:
        used.add(name)
        return name
    i = 2
    while f"{name}_{i}" in used:
        i += 1
    cand = f"{name}_{i}"
    used.add(cand)
    return cand


def _gateway_error_detail_from_body(data: Dict[str, Any], response_text: str) -> str:
    """与 hpc_orchestrator._gateway_http_error_detail 对齐，避免 detail 非 str 时丢长文案。"""
    d = data.get("detail")
    if isinstance(d, str) and d.strip():
        return d.strip()
    if isinstance(d, list):
        parts: List[str] = []
        for item in d:
            if isinstance(item, dict):
                loc = item.get("loc")
                msg = item.get("msg") or item.get("message")
                if loc is not None and msg is not None:
                    parts.append(f"{loc}: {msg}")
                else:
                    parts.append(json.dumps(item, ensure_ascii=False))
            else:
                parts.append(str(item))
        if parts:
            return "; ".join(parts)
    if d is not None and not isinstance(d, (str, list)):
        try:
            return json.dumps(d, ensure_ascii=False)
        except Exception:
            return str(d)
    return (response_text or "").strip() or "(empty gateway response body)"


class HPCMCPManager:
    """超算 MCP：通过 mcp-gateway 同步工具并代理调用。"""

    def __init__(self) -> None:
        self._gateway: str = (os.getenv("MCP_GATEWAY_URL") or "http://mcp-gateway:8002").rstrip("/")
        self._url: str = (os.getenv("HPC_MCP_URL") or "http://127.0.0.1:8001/mcp").rstrip("/")
        self._transport: str = os.getenv("HPC_MCP_TRANSPORT") or "streamableHttp"
        self._client: Optional[httpx.AsyncClient] = None
        self._last_error: Optional[str] = None
        self._registered_names: List[str] = []
        self._remote_tool_names: List[str] = []
        self._mcp_name_by_registry: Dict[str, str] = {}
        # ReAct 动态 schema：按 cache_key 短期缓存 (monotonic_deadline, openai_tools, reg_to_mcp)
        self._openai_schema_cache: Dict[str, Tuple[float, List[Dict[str, Any]], Dict[str, str]]] = {}
        self._openai_schema_cache_lock: Optional[asyncio.Lock] = None

    def _connected_local(self) -> bool:
        return bool(self._remote_tool_names)

    @property
    def connected(self) -> bool:
        return self._connected_local()

    async def _get_client(self) -> httpx.AsyncClient:
        if self._client is None or self._client.is_closed:
            timeout = httpx.Timeout(
                connect=15.0,
                read=float(os.getenv("HPC_MCP_CALL_TIMEOUT", "600")),
                write=60.0,
                pool=5.0,
            )
            self._client = httpx.AsyncClient(timeout=timeout)
        return self._client

    async def _close_client(self) -> None:
        if self._client is not None and not self._client.is_closed:
            await self._client.aclose()
        self._client = None

    async def fetch_mcp_tool_schemas(self) -> List[Dict[str, Any]]:
        """
        经网关 GET /tools 拉取上游 MCP tools/list 等价数据（name / description / inputSchema）。
        失败返回空列表，不抛异常（由调用方降级）。
        """
        try:
            c = await self._get_client()
            r = await c.get(f"{self._gateway}/tools", timeout=20.0)
            if not r.is_success:
                logger.warning("[HPC-MCP] fetch_mcp_tool_schemas HTTP %s", r.status_code)
                return []
            body = r.json()
            tools = body.get("tools") or []
            if not isinstance(tools, list):
                return []
            return [t for t in tools if isinstance(t, dict)]
        except Exception as e:  # noqa: BLE001
            logger.warning("[HPC-MCP] fetch_mcp_tool_schemas 失败: %s", e, exc_info=True)
            return []

    def resolve_registry_tool_to_mcp(self, registry_tool_name: str) -> Optional[str]:
        """OpenAI / 本地注册名 hpc_mcp_* → 远端 MCP tool 名（供网关 POST /call）。"""
        return self._mcp_name_by_registry.get(registry_tool_name)

    async def get_cached_openai_tools_for_session(
        self,
        cache_key: str,
        *,
        force_refresh: bool = False,
    ) -> Tuple[List[Dict[str, Any]], Dict[str, str], bool]:
        """
        拉取并转换 MCP schema → OpenAI tools[]；带 TTL 会话级缓存，防同会话重复打网关。

        Returns
        -------
        (openai_tool_dicts, reg_name_to_mcp_name, cache_hit)
        拉取失败或列表为空时返回 ([], {}, False)，调用方应回退到 registry.metadata。
        """
        ttl = max(5.0, float(os.getenv("HPC_MCP_SCHEMA_CACHE_TTL", "90")))
        now = time.monotonic()
        if self._openai_schema_cache_lock is None:
            self._openai_schema_cache_lock = asyncio.Lock()
        async with self._openai_schema_cache_lock:
            if not force_refresh and cache_key in self._openai_schema_cache:
                exp, o_tools, rmap = self._openai_schema_cache[cache_key]
                if now < exp and o_tools:
                    for reg, mcp in rmap.items():
                        self._mcp_name_by_registry.setdefault(reg, mcp)
                    return list(o_tools), dict(rmap), True

            raw = await self.fetch_mcp_tool_schemas()
            if not raw:
                return [], {}, False
            openai_tools, reg_map = mcp_tools_payload_to_openai_tools(raw)
            for reg, mcp in reg_map.items():
                self._mcp_name_by_registry.setdefault(reg, mcp)
            self._openai_schema_cache[cache_key] = (now + ttl, list(openai_tools), dict(reg_map))
            if len(self._openai_schema_cache) > 200:
                dead = [k for k, (ex, _, _) in self._openai_schema_cache.items() if ex < now]
                for k in dead[:100]:
                    self._openai_schema_cache.pop(k, None)
            return openai_tools, reg_map, False

    async def status_dict_async(self) -> Dict[str, Any]:
        gw_status: Dict[str, Any] = {}
        try:
            c = await self._get_client()
            r = await c.get(f"{self._gateway}/status")
            if r.is_success:
                gw_status = r.json()
        except Exception as e:  # noqa: BLE001
            gw_status = {"gateway_error": str(e)}
        gw_ok = bool(gw_status.get("connected"))
        return {
            "connected": gw_ok if "connected" in gw_status else self._connected_local(),
            "url": self._url,
            "transport": self._transport,
            "gateway": self._gateway,
            "tools": list(self._remote_tool_names),
            "registered": list(self._registered_names),
            "last_error": self._last_error,
            "mcp_available": True,
            # 来自 mcp-gateway /status：上游 MCP 是否已配置鉴权头（排障用，不含密钥）
            "upstream_auth_configured": gw_status.get("upstream_auth_configured"),
            "gateway_status": gw_status,
        }

    def _clear_registered_tools(self) -> None:
        if self._registered_names:
            registry.unregister_many(self._registered_names)
        self._registered_names.clear()
        self._mcp_name_by_registry.clear()
        self._remote_tool_names.clear()

    def _register_tools_from_payload(self, tools: List[Dict[str, Any]]) -> None:
        self._clear_registered_tools()
        used_names: set = set()
        for item in tools:
            mcp_name = (item.get("name") or "").strip()
            if not mcp_name:
                continue
            desc = (item.get("description") or "").strip() or (
                f"超算 MCP 工具 `{mcp_name}`；参数请传入 arguments 字典（与远端 schema 一致）。"
            )
            reg_name = _registry_name(mcp_name, used_names)
            self._remote_tool_names.append(mcp_name)
            self._mcp_name_by_registry[reg_name] = mcp_name

            def _make_wrapper(
                _mcp: str,
                _reg: str,
                _desc: str,
                mgr: HPCMCPManager,
            ) -> Callable[..., Any]:
                @registry.register(
                    name=_reg,
                    description=_desc,
                    category="HPC-MCP",
                    output_type="json",
                )
                @safe_tool_execution
                async def _hpc_proxy(arguments: Optional[dict] = None) -> Dict[str, Any]:
                    return await mgr._invoke_remote(_mcp, arguments or {})

                return _hpc_proxy

            _make_wrapper(mcp_name, reg_name, desc, self)
            self._registered_names.append(reg_name)
            logger.info("✅ [HPC-MCP] 已注册工具 %s -> %s", mcp_name, reg_name)

    async def _invoke_remote(self, mcp_tool_name: str, arguments: dict) -> Dict[str, Any]:
        args = arguments or {}
        t_call = _gateway_call_timeout()
        read_sec = max(1.0, float(os.getenv("HPC_MCP_INVOKE_TIMEOUT_SEC", "30")))
        logger.info(
            "[MCP Execution] → POST %s/call tool=%r timeout_read≈%.1fs args_preview=%s",
            self._gateway,
            mcp_tool_name,
            read_sec,
            _args_preview_for_log(args),
        )
        try:
            c = await self._get_client()
            r = await c.post(
                f"{self._gateway}/call",
                json={"tool_name": mcp_tool_name, "arguments": args},
                timeout=t_call,
            )
            try:
                data = r.json()
            except Exception:  # noqa: BLE001
                data = {"detail": r.text}
            if r.status_code >= 400:
                msg = _gateway_error_detail_from_body(data if isinstance(data, dict) else {}, r.text)
                err_text = (
                    f"Tool Execution Failed: MCP gateway HTTP {r.status_code}. Detail: {msg}"
                )
                logger.info(
                    "[MCP Execution] ← tool=%r http=%s (error body ok) summary=%s",
                    mcp_tool_name,
                    r.status_code,
                    (msg[:120] + "…") if len(msg) > 120 else msg,
                )
                return {
                    "status": "error",
                    "message": err_text,
                }
            if not isinstance(data, dict):
                logger.info("[MCP Execution] ← tool=%r invalid JSON envelope", mcp_tool_name)
                return {
                    "status": "error",
                    "message": "Tool Execution Failed: invalid gateway response body (not a JSON object).",
                }
            if data.get("status") == "error":
                logger.info(
                    "[MCP Execution] ← tool=%r remote status=error message=%s",
                    mcp_tool_name,
                    (data.get("message") or "")[:200],
                )
                return {
                    "status": "error",
                    "message": data.get("message")
                    or "Tool Execution Failed: remote MCP tool returned error.",
                    "content": data.get("content") or [],
                }
            logger.info("[MCP Execution] ← tool=%r status=success", mcp_tool_name)
            return {
                "status": "success",
                "content": data.get("content") or [],
            }
        except httpx.TimeoutException as e:
            sec = max(1.0, float(os.getenv("HPC_MCP_INVOKE_TIMEOUT_SEC", "30")))
            msg = f"Tool Execution Failed: Gateway timeout after {sec:g}s ({type(e).__name__})."
            logger.warning("[MCP Execution] ← tool=%r %s", mcp_tool_name, msg)
            return {"status": "error", "message": msg}
        except httpx.RequestError as e:
            msg = f"Tool Execution Failed: gateway request error ({type(e).__name__}: {e})."
            logger.warning("[MCP Execution] ← tool=%r %s", mcp_tool_name, msg)
            return {"status": "error", "message": msg}
        except Exception as e:  # noqa: BLE001
            logger.exception("[HPC-MCP] gateway call 未预期异常: %s", mcp_tool_name)
            return {
                "status": "error",
                "message": f"Tool Execution Failed: unexpected error ({type(e).__name__}: {e}).",
            }

    async def call_tool(self, mcp_tool_name: str, arguments: Optional[dict] = None) -> Dict[str, Any]:
        """公开的 MCP 工具调用入口，供 API 路由按工具名转发。"""
        return await self._invoke_remote((mcp_tool_name or "").strip(), arguments or {})

    async def connect(self, url: str, transport: str = "streamableHttp") -> None:
        """通知网关重连超算 MCP，并轮询 /status + /tools 注册本地代理工具。"""
        self._clear_registered_tools()
        self._last_error = None
        self._url = (url or self._url).rstrip("/")
        self._transport = transport or self._transport

        try:
            c = await self._get_client()
            rr = await c.post(
                f"{self._gateway}/reconnect",
                json={"url": self._url, "transport": self._transport},
                timeout=30.0,
            )
            if not rr.is_success:
                self._last_error = rr.text or f"网关 reconnect HTTP {rr.status_code}"
                logger.warning("⚠️ [HPC-MCP] reconnect 失败: %s", self._last_error)
                return
        except Exception as e:  # noqa: BLE001
            self._last_error = str(e)
            logger.warning("⚠️ [HPC-MCP] reconnect 请求异常: %s", e, exc_info=True)
            return

        deadline = time.monotonic() + float(os.getenv("HPC_MCP_GATEWAY_POLL_TIMEOUT", "120"))
        poll_interval = float(os.getenv("HPC_MCP_GATEWAY_POLL_INTERVAL", "2"))
        c = await self._get_client()
        while time.monotonic() < deadline:
            try:
                st = await c.get(f"{self._gateway}/status", timeout=15.0)
                if st.is_success:
                    body = st.json()
                    if body.get("connected"):
                        tr = await c.get(f"{self._gateway}/tools", timeout=15.0)
                        tools: List[Dict[str, Any]] = []
                        if tr.is_success:
                            tj = tr.json()
                            tools = tj.get("tools") or []
                        self._register_tools_from_payload(tools if isinstance(tools, list) else [])
                        self._last_error = None
                        logger.info("✅ [HPC-MCP] 经网关已同步工具数=%s", len(self._remote_tool_names))
                        return
                    self._last_error = body.get("last_error")
            except Exception as e:  # noqa: BLE001
                self._last_error = str(e)
                logger.debug("[HPC-MCP] 轮询网关: %s", e)
            await asyncio.sleep(poll_interval)

        self._last_error = self._last_error or "网关 MCP 连接超时（轮询 /status）"
        logger.warning("⚠️ [HPC-MCP] %s", self._last_error)

    async def disconnect(self) -> None:
        try:
            c = await self._get_client()
            await c.post(f"{self._gateway}/disconnect", timeout=15.0)
        except Exception as e:  # noqa: BLE001
            logger.debug("[HPC-MCP] gateway disconnect: %s", e)
        self._clear_registered_tools()
        await self._close_client()

    async def shutdown(self) -> None:
        await self.disconnect()


hpc_mcp_manager = HPCMCPManager()
