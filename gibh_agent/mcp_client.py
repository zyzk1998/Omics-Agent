# -*- coding: utf-8 -*-
"""
超算 MCP 客户端：经 mcp-gateway HTTP 转发，主进程不依赖官方 mcp 包。

- 连接失败仅记日志，不阻塞 FastAPI。
- 需网关可达（MCP_GATEWAY_URL）；工具列表来自 GET /tools，调用经 POST /call。
"""
from __future__ import annotations

import asyncio
import logging
import os
import re
import time
from typing import Any, Callable, Dict, List, Optional

import httpx

from gibh_agent.core.tool_registry import registry
from gibh_agent.core.utils import safe_tool_execution

logger = logging.getLogger(__name__)


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


class HPCMCPManager:
    """超算 MCP：通过 mcp-gateway 同步工具并代理调用。"""

    def __init__(self) -> None:
        self._gateway: str = (os.getenv("MCP_GATEWAY_URL") or "http://mcp-gateway:8002").rstrip("/")
        self._url: str = (os.getenv("HPC_MCP_URL") or "http://192.168.32.31:8001/mcp").rstrip("/")
        self._transport: str = os.getenv("HPC_MCP_TRANSPORT") or "streamableHttp"
        self._client: Optional[httpx.AsyncClient] = None
        self._last_error: Optional[str] = None
        self._registered_names: List[str] = []
        self._remote_tool_names: List[str] = []
        self._mcp_name_by_registry: Dict[str, str] = {}

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
        try:
            c = await self._get_client()
            r = await c.post(
                f"{self._gateway}/call",
                json={"tool_name": mcp_tool_name, "arguments": arguments or {}},
            )
            try:
                data = r.json()
            except Exception:  # noqa: BLE001
                data = {"detail": r.text}
            if r.status_code >= 400:
                detail = data.get("detail") if isinstance(data, dict) else None
                msg = detail if isinstance(detail, str) else r.text
                return {
                    "status": "error",
                    "message": f"HPC Error: gateway HTTP {r.status_code} ({msg}), please inform the user.",
                }
            if not isinstance(data, dict):
                return {
                    "status": "error",
                    "message": "HPC Error: invalid gateway response, please inform the user.",
                }
            if data.get("status") == "error":
                return {
                    "status": "error",
                    "message": data.get("message") or "HPC Error: remote tool error, please inform the user.",
                    "content": data.get("content") or [],
                }
            return {
                "status": "success",
                "content": data.get("content") or [],
            }
        except httpx.TimeoutException:
            return {
                "status": "error",
                "message": "HPC Error: Connection timeout, please inform the user.",
            }
        except Exception as e:  # noqa: BLE001
            logger.exception("[HPC-MCP] gateway call 失败: %s", mcp_tool_name)
            return {
                "status": "error",
                "message": f"HPC Error: gateway call failed ({e}), please inform the user.",
            }

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
