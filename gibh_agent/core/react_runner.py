# -*- coding: utf-8 -*-
"""
DeepReActRunner：基于 astream 的 ReAct 循环，支持流式 reasoning 透出与 tool_calls 碎片聚合。

- 流式阶段仅累积 tool_calls，流结束后再执行工具（见 DeepReActRunner.run）。
- 工具执行包裹 try/except，错误以 tool 消息回注模型供自主反思。
- 虚拟工具 ask_human_for_clarification：不查注册表，发射 status 后直接结束 run（挂起）。
"""
from __future__ import annotations

import inspect
import json
import logging
import traceback
from typing import Any, Awaitable, Callable, Dict, List, Optional

from gibh_agent.core.stream_utils import _delta_get

logger = logging.getLogger(__name__)

# 与 Orchestrator 注入的 OpenAI tools 名称保持一致（虚拟工具，无 registry 实现）
ASK_HUMAN_FOR_CLARIFICATION_TOOL_NAME = "ask_human_for_clarification"

EmitCallback = Callable[[str, Dict[str, Any]], Awaitable[None]]


def _tool_args_preview(args: Dict[str, Any], max_len: int = 400) -> str:
    try:
        s = json.dumps(args, ensure_ascii=False, default=str)
    except Exception:  # noqa: BLE001
        s = str(args)
    return (s[: max_len - 3] + "...") if len(s) > max_len else s


class DeepReActCircuitBreakerError(RuntimeError):
    """ReAct 循环超过 max_steps 时抛出，用于强制熔断。"""


def _normalize_tool_call_index(raw: Any) -> int:
    if raw is None:
        return 0
    try:
        return int(raw)
    except (TypeError, ValueError):
        return 0


def _accumulate_streaming_tool_calls(
    accumulated: Dict[int, Dict[str, str]],
    chunk: Any,
) -> None:
    """
    将单个 ChatCompletionChunk 中 delta.tool_calls 的碎片按 index 聚合到 accumulated。
    每个槽位: id, name, arguments（均为字符串片段拼接）。
    """
    choices = getattr(chunk, "choices", None) if not isinstance(chunk, dict) else chunk.get("choices")
    if not choices:
        return
    first = choices[0]
    delta = getattr(first, "delta", None) if not isinstance(first, dict) else first.get("delta")
    if delta is None:
        return

    tcalls = getattr(delta, "tool_calls", None) if not isinstance(delta, dict) else delta.get("tool_calls")
    if not tcalls:
        return

    for tc in tcalls:
        idx = _normalize_tool_call_index(
            getattr(tc, "index", None) if not isinstance(tc, dict) else tc.get("index")
        )
        if idx not in accumulated:
            accumulated[idx] = {"id": "", "name": "", "arguments": ""}
        slot = accumulated[idx]

        tid = getattr(tc, "id", None) if not isinstance(tc, dict) else tc.get("id")
        if tid:
            slot["id"] = str(tid)

        fn = getattr(tc, "function", None) if not isinstance(tc, dict) else tc.get("function")
        if fn is None:
            continue
        name_piece = getattr(fn, "name", None) if not isinstance(fn, dict) else fn.get("name")
        if name_piece:
            slot["name"] += str(name_piece)
        args_piece = getattr(fn, "arguments", None) if not isinstance(fn, dict) else fn.get("arguments")
        if args_piece:
            slot["arguments"] += str(args_piece)


def _streaming_chunks_to_ordered_tool_calls(
    accumulated: Dict[int, Dict[str, str]],
) -> List[Dict[str, str]]:
    """按 index 排序，过滤无名称的占位槽位。"""
    out: List[Dict[str, str]] = []
    for idx in sorted(accumulated.keys()):
        slot = accumulated[idx]
        name = (slot.get("name") or "").strip()
        if not name:
            continue
        out.append(
            {
                "id": (slot.get("id") or "").strip(),
                "name": name,
                "arguments": slot.get("arguments") or "",
            }
        )
    return out


async def execute_tool(tool_fn: Any, args: Dict[str, Any]) -> Any:
    """
    统一执行注册工具：协程则 await，否则同步调用（与 Orchestrator 策略一致）。
    """
    if tool_fn is None:
        raise ValueError("tool_fn is None")
    if inspect.iscoroutinefunction(tool_fn):
        return await tool_fn(**args)
    return tool_fn(**args)


class DeepReActRunner:
    """
    使用 llm_client.astream 的 Deep ReAct 运行器。

    Parameters
    ----------
    llm_client:
        需实现 ``async def astream(messages, **kwargs) -> AsyncIterator[chunk]``（如 LLMClient）。
    registry:
        需实现 ``get_tool(name: str) -> Optional[Callable]``（如 ToolRegistry）。
    max_steps:
        单次 run 内允许的大模型轮数上限（含首轮），超出则抛 DeepReActCircuitBreakerError。
    """

    def __init__(
        self,
        llm_client: Any,
        registry: Any,
        max_steps: int = 10,
    ) -> None:
        self._llm = llm_client
        self._registry = registry
        self.max_steps = max(1, int(max_steps))

    async def run(
        self,
        messages: List[Dict[str, Any]],
        emit_callback: EmitCallback,
        *,
        tools: Optional[List[Dict[str, Any]]] = None,
        tool_choice: Any = "auto",
        model: Optional[str] = None,
        **llm_kwargs: Any,
    ) -> str:
        """
        驱动多轮 ReAct：就地扩展 ``messages``（调用方保留同一 list 引用即可恢复上下文）。

        emit_callback 约定：
        - ``("thought", {"content": str})``：reasoning_content 流式片段
        - ``("message", {"content": str})``：delta.content 可见正文流式片段（与 fast 模式 SSE 对齐）
        - ``("process_log", {"message": str})``：过程日志（算子执行等）
        - ``("status", {"content": str, "state": str})``：状态（如 human_input_required）

        Returns
        -------
        str
            ``"complete"``：模型最终一轮无工具调用，已追加 assistant 消息。
            ``"human_input_required"``：命中虚拟工具 ``ask_human_for_clarification``，已挂起。
        """
        step = 0
        while True:
            if step >= self.max_steps:
                raise DeepReActCircuitBreakerError(
                    f"DeepReAct 已超过 max_steps={self.max_steps}，触发强制熔断。"
                )
            step += 1

            accumulated: Dict[int, Dict[str, str]] = {}
            content_parts: List[str] = []

            stream_kwargs = dict(llm_kwargs)
            if model is not None:
                stream_kwargs["model"] = model
            if tools is not None:
                stream_kwargs["tools"] = tools
                stream_kwargs["tool_choice"] = tool_choice

            async for chunk in self._llm.astream(messages, **stream_kwargs):
                choices = getattr(chunk, "choices", None) if not isinstance(chunk, dict) else chunk.get(
                    "choices"
                )
                if not choices:
                    continue
                first = choices[0]
                delta = getattr(first, "delta", None) if not isinstance(first, dict) else first.get("delta")

                reasoning = _delta_get(delta, "reasoning_content")
                if reasoning:
                    await emit_callback("thought", {"content": reasoning})

                content = _delta_get(delta, "content")
                if content:
                    content_parts.append(content)
                    await emit_callback("message", {"content": content})

                _accumulate_streaming_tool_calls(accumulated, chunk)

            assistant_text = "".join(content_parts)
            ordered_calls = _streaming_chunks_to_ordered_tool_calls(accumulated)

            if not ordered_calls:
                messages.append(
                    {
                        "role": "assistant",
                        "content": assistant_text if assistant_text.strip() else None,
                    }
                )
                return "complete"

            openai_tool_calls: List[Dict[str, Any]] = []
            for i, tc in enumerate(ordered_calls):
                call_id = tc["id"] or f"call_{step}_{i}"
                openai_tool_calls.append(
                    {
                        "id": call_id,
                        "type": "function",
                        "function": {
                            "name": tc["name"],
                            "arguments": tc["arguments"] if tc["arguments"].strip() else "{}",
                        },
                    }
                )

            messages.append(
                {
                    "role": "assistant",
                    "content": assistant_text if assistant_text.strip() else None,
                    "tool_calls": openai_tool_calls,
                }
            )

            for i, tc in enumerate(ordered_calls):
                tool_name = tc["name"]
                args_str = tc["arguments"] if tc["arguments"].strip() else "{}"
                call_id = tc["id"] or f"call_{step}_{i}"

                if tool_name == ASK_HUMAN_FOR_CLARIFICATION_TOOL_NAME:
                    hit_options: List[str] = []
                    try:
                        _args = json.loads(args_str)
                        if isinstance(_args, dict):
                            raw_opt = _args.get("options")
                            if isinstance(raw_opt, list):
                                hit_options = [
                                    str(x).strip()
                                    for x in raw_opt
                                    if x is not None and str(x).strip()
                                ]
                    except (json.JSONDecodeError, TypeError, ValueError):
                        pass
                    await emit_callback(
                        "status",
                        {
                            "content": "需要您的确认或补充信息...",
                            "state": "human_input_required",
                        },
                    )
                    if hit_options:
                        await emit_callback("suggestions", {"suggestions": hit_options})
                    return "human_input_required"

                await emit_callback(
                    "process_log",
                    {"message": f"⏳ 正在执行算子: {tool_name}"},
                )
                try:
                    args = json.loads(args_str)
                except json.JSONDecodeError:
                    args = {}

                tool_fn = self._registry.get_tool(tool_name)
                # 动态 MCP schema 注入场景：工具可能尚未在 Registry 注册，经网关 POST /call 代理执行
                if tool_fn is None and isinstance(tool_name, str) and tool_name.startswith("hpc_mcp_"):
                    try:
                        from gibh_agent.mcp_client import hpc_mcp_manager

                        mcp_remote = hpc_mcp_manager.resolve_registry_tool_to_mcp(tool_name)
                        if mcp_remote:
                            logger.info(
                                "[MCP Execution] DeepReAct dispatch registry=%r → mcp=%r args_preview=%s",
                                tool_name,
                                mcp_remote,
                                _tool_args_preview(args),
                            )
                            try:
                                result = await hpc_mcp_manager._invoke_remote(mcp_remote, args)
                            except Exception as gw_exc:  # noqa: BLE001
                                err_line = (
                                    f"Tool Execution Failed: unexpected exception in _invoke_remote "
                                    f"({type(gw_exc).__name__}: {gw_exc})."
                                )
                                logger.exception(
                                    "[MCP Execution] DeepReAct registry=%r mcp=%r %s",
                                    tool_name,
                                    mcp_remote,
                                    err_line,
                                )
                                await emit_callback(
                                    "process_log",
                                    {
                                        "message": f"❌ {tool_name}（MCP）异常已捕获并喂回模型",
                                    },
                                )
                                messages.append(
                                    {
                                        "role": "tool",
                                        "tool_call_id": call_id,
                                        "content": err_line,
                                    }
                                )
                                continue

                            if not isinstance(result, dict):
                                err_line = "Tool Execution Failed: MCP gateway returned non-object response."
                                logger.warning(
                                    "[MCP Execution] DeepReAct registry=%r mcp=%r invalid_result_type=%s",
                                    tool_name,
                                    mcp_remote,
                                    type(result).__name__,
                                )
                                await emit_callback(
                                    "process_log",
                                    {"message": f"⚠️ {tool_name}（MCP）响应异常，已喂回模型"},
                                )
                                messages.append(
                                    {
                                        "role": "tool",
                                        "tool_call_id": call_id,
                                        "content": err_line,
                                    }
                                )
                                continue

                            payload = json.dumps(result, ensure_ascii=False, default=str)
                            if result.get("status") == "error":
                                logger.info(
                                    "[MCP Execution] DeepReAct done registry=%r mcp=%r status=error message=%s",
                                    tool_name,
                                    mcp_remote,
                                    (result.get("message") or "")[:300],
                                )
                                await emit_callback(
                                    "process_log",
                                    {
                                        "message": f"⚠️ {tool_name}（MCP）业务/网关错误，已喂回模型供反思",
                                    },
                                )
                            else:
                                logger.info(
                                    "[MCP Execution] DeepReAct done registry=%r mcp=%r status=success",
                                    tool_name,
                                    mcp_remote,
                                )
                                await emit_callback(
                                    "process_log",
                                    {"message": f"✅ {tool_name}（MCP 网关）执行成功"},
                                )
                            messages.append(
                                {
                                    "role": "tool",
                                    "tool_call_id": call_id,
                                    "content": payload,
                                }
                            )
                            continue
                    except Exception as imp_exc:  # noqa: BLE001
                        logger.warning("[DeepReAct] HPC MCP 代理分支加载失败: %s", imp_exc)

                if tool_fn is None:
                    await emit_callback(
                        "process_log",
                        {"message": f"❌ 未知算子: {tool_name}"},
                    )
                    messages.append(
                        {
                            "role": "tool",
                            "tool_call_id": call_id,
                            "content": f"Error: unknown tool {tool_name!r}",
                        }
                    )
                    continue

                try:
                    result = await execute_tool(tool_fn, args)
                    payload = json.dumps(result, ensure_ascii=False, default=str)
                    await emit_callback(
                        "process_log",
                        {"message": f"✅ {tool_name} 执行成功"},
                    )
                    messages.append(
                        {
                            "role": "tool",
                            "tool_call_id": call_id,
                            "content": payload,
                        }
                    )
                except Exception as e:
                    error_trace = traceback.format_exc()
                    await emit_callback(
                        "process_log",
                        {"message": f"❌ {tool_name} 执行异常，已触发自主反思..."},
                    )
                    messages.append(
                        {
                            "role": "tool",
                            "tool_call_id": call_id,
                            "content": f"Error: {str(e)}\nTraceback:\n{error_trace}",
                        }
                    )
