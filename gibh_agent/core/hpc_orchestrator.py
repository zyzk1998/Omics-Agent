# -*- coding: utf-8 -*-
"""
HPC 直连编排（遗留模块）：chat SSE 已不再引用本模块；统一走 AgentOrchestrator + DeepReAct + MCP 动态工具。

保留供运维脚本/单测或手工复现网关 /call 行为时参考。

- 经 MCP Gateway POST /call 调用远端工具。
- 默认工具名为 execute_hpc_command，arguments 为 {"command": "<用户输入>"}，与 curl POST /call 语义对齐。
- 若整段用户输入为 JSON 且含 tool_name + arguments，则原样透传（与 curl body 等价）。
- 环境变量 HPC_ORCH_CHAT_TOOL=hpc_mcp_chat 时仍使用旧版 NL 大包字段（需远端注册该工具）。
- 产出与主链路一致的 SSE（message / status / done / state_snapshot），供时光机入库。
"""
from __future__ import annotations

import json
import logging
import os
import re
import time
from typing import Any, AsyncIterator, Dict, List, Optional

import httpx

from gibh_agent.core.utils import sanitize_for_json

logger = logging.getLogger(__name__)

HPC_ORCH_CHAT_TOOL = (os.getenv("HPC_ORCH_CHAT_TOOL") or "execute_hpc_command").strip()

_NO_ARG_HPC_TOOLS = frozenset({"test_hpc_connection", "get_hpc_system_info"})
MCP_GATEWAY_URL = (os.getenv("MCP_GATEWAY_URL") or "http://mcp-gateway:8002").rstrip("/")

CLI_CONTRACT_PROMPT = """你是一个运行在国家级超算上的高级 CLI Agent。请遵循以下严格边界：
1. **CLI 优先**：优先使用 shell、slurm (sbatch/squeue) 等命令行工具解决问题。
2. **严禁阻塞**：对于预估耗时超过 1 分钟的重型计算任务，**绝对禁止**同步等待！必须将其提交到后台队列，并立刻向用户返回任务的 Job ID，以及查询状态的具体命令。
3. **结果解释**：执行完 CLI 命令后，必须用简明扼要的自然语言向用户解释结果。"""

_UPLOAD_PATH_RE = re.compile(r"/app/uploads/\S+")


def _format_history_messages(messages: List[Dict[str, Any]]) -> str:
    lines: List[str] = []
    for item in messages or []:
        if not isinstance(item, dict):
            continue
        role = (item.get("role") or item.get("type") or "user").strip()
        content = item.get("content")
        if content is None:
            content = item.get("message") or ""
        if isinstance(content, dict):
            content = content.get("text") or json.dumps(content, ensure_ascii=False)
        elif not isinstance(content, str):
            content = str(content)
        if content.strip():
            lines.append(f"{role}: {content.strip()}")
    return "\n".join(lines) if lines else "(无历史记录)"


def _collect_upload_paths(uploaded_files: List[Any]) -> List[str]:
    out: List[str] = []
    for f in uploaded_files or []:
        p = ""
        if isinstance(f, dict):
            p = f.get("path") or f.get("file_path") or ""
        elif isinstance(f, str):
            p = f
        p = (p or "").strip()
        if p and p not in out:
            out.append(p)
    return out


def _paths_declaration(paths: List[str]) -> str:
    if not paths:
        return "(本回合未挂载网关可读的本地文件；若需数据请让用户上传或通过 scp/rsync 同步到超算侧。)"
    joined = "\n".join(f"- {p}" for p in paths)
    return (
        f"{joined}\n（以上为网关服务器上的绝对路径；若超算节点无法直接读取，请使用 scp/rsync 拉取后再处理。）"
    )


def _try_parse_explicit_mcp_payload(text: str) -> Optional[tuple[str, Dict[str, Any]]]:
    """若用户消息整段为 JSON 且含 tool_name、arguments，则与 curl /call 的 body 等价。"""
    s = (text or "").strip()
    if len(s) < 2 or not s.startswith("{"):
        return None
    try:
        data = json.loads(s)
    except json.JSONDecodeError:
        return None
    if not isinstance(data, dict):
        return None
    tn = data.get("tool_name")
    if not isinstance(tn, str) or not tn.strip():
        return None
    args = data.get("arguments")
    if args is None:
        args = {}
    if not isinstance(args, dict):
        return None
    return (tn.strip(), args)


def _build_nl_tool_arguments(
    *,
    user_latest: str,
    history_txt: str,
    paths: List[str],
    paths_decl: str,
    session_id: str,
    user_id: str,
) -> Dict[str, Any]:
    structured_block = (
        f"[网关服务器文件路径]: {paths_decl}\n\n"
        f"[历史对话上下文]:\n{history_txt}\n\n"
        f"[用户最新指令]: {user_latest}"
    )
    full_prompt = f"{CLI_CONTRACT_PROMPT}\n\n---\n{structured_block}"
    return {
        "prompt": full_prompt,
        "user_message": user_latest,
        "history_formatted": history_txt,
        "uploaded_files_paths": paths,
        "gateway_paths_note": paths_decl,
        "session_id": session_id,
        "user_id": user_id,
        "system_directive": CLI_CONTRACT_PROMPT,
    }


def _arguments_for_named_tool(
    tool: str,
    *,
    user_latest: str,
    history_txt: str,
    paths: List[str],
    paths_decl: str,
    session_id: str,
    user_id: str,
) -> Dict[str, Any]:
    t = (tool or "").strip()
    if t == "hpc_mcp_chat":
        return _build_nl_tool_arguments(
            user_latest=user_latest,
            history_txt=history_txt,
            paths=paths,
            paths_decl=paths_decl,
            session_id=session_id,
            user_id=user_id,
        )
    if t == "execute_hpc_command":
        return {"command": user_latest}
    if t in _NO_ARG_HPC_TOOLS:
        return {}
    logger.warning(
        "[hpc_orchestrator] 未内置参数模板的工具 %s，将发送空 arguments；"
        "请改用整段 JSON 消息指定 tool_name/arguments，或设置 HPC_ORCH_CHAT_TOOL",
        t,
    )
    return {}


def _format_sse(event_type: str, data: Dict[str, Any]) -> str:
    try:
        safe = sanitize_for_json(data if isinstance(data, dict) else {"payload": data})
        json_data = json.dumps(safe, ensure_ascii=False, allow_nan=False)
    except Exception as e:  # noqa: BLE001
        logger.error("hpc_orchestrator SSE 序列化失败: %s", e, exc_info=True)
        json_data = json.dumps({"error": "sse_serialize_error", "message": str(e)}, ensure_ascii=False)
    return f"event: {event_type}\ndata: {json_data}\n\n"


def _sniff_assets_from_text(text: str) -> List[Dict[str, Any]]:
    seen = set()
    assets: List[Dict[str, Any]] = []
    for m in _UPLOAD_PATH_RE.finditer(text or ""):
        p = m.group(0).rstrip(".,;)")
        if p not in seen:
            seen.add(p)
            assets.append({"path": p, "source": "upload_sniff"})
    return assets


def _merge_state_for_emit(
    state: Dict[str, Any],
    event_type: str,
    data: Dict[str, Any],
) -> None:
    if event_type == "message":
        state["text"] = (state.get("text") or "") + (data.get("content") or "")
    elif event_type == "status" and data:
        _c = data.get("content", "")
        if not (isinstance(_c, str) and _c.startswith("[AgenticLog] ")):
            state.setdefault("process_log", []).append(
                {"content": _c, "state": data.get("state", "running")}
            )


def _gateway_http_error_detail(body: Any, response_text: str) -> str:
    """
    从网关 JSON 体提取 FastAPI/Starlette 的 detail。
    若 detail 为 list（校验错误形态）或 dict，禁止退化成整段 r.text 导致丢失结构化长文案。
    """
    if isinstance(body, dict):
        d = body.get("detail")
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


async def _call_gateway_chat(*, tool_name: str, arguments: Dict[str, Any]) -> Dict[str, Any]:
    timeout = httpx.Timeout(connect=10.0, read=300.0, write=60.0, pool=10.0)
    async with httpx.AsyncClient(timeout=timeout) as client:
        r = await client.post(
            f"{MCP_GATEWAY_URL}/call",
            json={"tool_name": tool_name, "arguments": arguments or {}},
        )
        try:
            body = r.json()
        except Exception:  # noqa: BLE001
            body = {"detail": r.text}
        if r.status_code >= 400:
            msg = _gateway_http_error_detail(body if isinstance(body, dict) else {}, r.text)
            return {"status": "error", "message": f"网关 HTTP {r.status_code}: {msg}", "content": []}
        if not isinstance(body, dict):
            return {"status": "error", "message": "网关返回非 JSON 对象", "content": []}
        return body


def _content_to_plain_text(content: Any) -> str:
    if not content:
        return ""
    parts: List[str] = []
    if isinstance(content, list):
        for block in content:
            if isinstance(block, dict):
                t = block.get("text")
                if t is not None:
                    parts.append(str(t))
                else:
                    parts.append(str(block))
            else:
                parts.append(str(block))
    elif isinstance(content, str):
        parts.append(content)
    else:
        parts.append(str(content))
    return "\n".join(parts).strip()


def _wrap_geek_terminal_markdown(text: str) -> str:
    """
    将网关/工具返回的正文包装为 Markdown 代码块。
    JSON 可解析时：若为含 output 键的 dict（execute_hpc_command 常见形态），将 output/error 以真实换行写入 bash 围栏；
    否则使用缩进 JSON 围栏。不可解析时原文写入 bash 围栏。
    """
    raw_s = text if isinstance(text, str) else str(text or "")
    trimmed = raw_s.strip()
    if not trimmed:
        return ""
    try:
        parsed = json.loads(trimmed)
        if isinstance(parsed, dict) and "output" in parsed:
            out_raw = parsed.get("output", "")
            err_raw = parsed.get("error", "")
            out_text = "" if out_raw is None else str(out_raw)
            err_text = "" if err_raw is None else str(err_raw)
            segments: List[str] = []
            if out_text.strip():
                segments.append(out_text.rstrip("\n"))
            if err_text.strip():
                if segments:
                    segments.append("")
                segments.append("[STDERR]:")
                segments.append(err_text.rstrip("\n"))
            if not segments:
                final_text = "[命令执行成功，无输出内容]"
            else:
                final_text = "\n".join(segments)
            return f"```bash\n{final_text}\n```"
        formatted = json.dumps(parsed, indent=2, ensure_ascii=False)
        return f"```json\n{formatted}\n```"
    except json.JSONDecodeError:
        body = raw_s.rstrip("\n")
        return f"```bash\n{body}\n```"


async def stream_hpc_direct_chat(
    *,
    query: str,
    history: Optional[List[Dict[str, Any]]] = None,
    messages: Optional[List[Dict[str, Any]]] = None,
    uploaded_files: Optional[List[Any]] = None,
    session_id: str = "",
    user_id: str = "",
) -> AsyncIterator[str]:
    """
    异步生成器：产出与主聊天一致的 SSE 字符串序列。
    """
    history_src = messages if (messages and len(messages) > 0) else (history or [])
    history_txt = _format_history_messages(list(history_src))
    paths = _collect_upload_paths(list(uploaded_files or []))
    paths_decl = _paths_declaration(paths)

    raw_query = (query or "").strip()
    user_latest = raw_query or "(空指令)"

    explicit = _try_parse_explicit_mcp_payload(raw_query)
    if explicit:
        use_tool, use_args = explicit
    else:
        use_tool = HPC_ORCH_CHAT_TOOL
        use_args = _arguments_for_named_tool(
            use_tool,
            user_latest=user_latest,
            history_txt=history_txt,
            paths=paths,
            paths_decl=paths_decl,
            session_id=session_id,
            user_id=user_id,
        )

    state_snapshot: Dict[str, Any] = {
        "text": "",
        "reasoning": "",
        "workflow": None,
        "steps": [],
        "process_log": [],
        "report": None,
        "duration": None,
        "_start_time": time.time(),
        "hpc_mcp_direct": True,
        "hpc_tool": use_tool,
        "hpc_explicit_json": bool(explicit),
    }

    yield _format_sse(
        "status",
        {"content": "正在通过 MCP 网关连接 HPC…", "state": "running"},
    )
    _merge_state_for_emit(state_snapshot, "status", {"content": "正在通过 MCP 网关连接 HPC…", "state": "running"})

    reply_text = ""
    gateway_plain: Optional[str] = None
    try:
        if use_tool == "execute_hpc_command":
            cmd = str((use_args or {}).get("command", "")).strip()
            if not cmd or cmd == "(空指令)":
                reply_text = (
                    "请输入要在超算执行的 shell 命令（将作为 execute_hpc_command 的 command）。\n"
                    "或与 curl 等价地发送整段 JSON，例如：\n"
                    '{"tool_name":"execute_hpc_command","arguments":{"command":"hostname"}}\n'
                    "无参工具示例：\n"
                    '{"tool_name":"test_hpc_connection","arguments":{}}'
                )
        if not reply_text:
            raw = await _call_gateway_chat(tool_name=use_tool, arguments=use_args)
            if raw.get("status") == "error":
                gateway_plain = raw.get("message") or "超算 MCP 调用失败"
            else:
                gateway_plain = _content_to_plain_text(raw.get("content"))
    except httpx.TimeoutException:
        reply_text = (
            "⏳ 超算节点正在执行耗时操作或网络延迟，连接已转入后台。请稍后主动查询进度。"
        )
        logger.warning("[hpc_orchestrator] gateway timeout after read=300s")
    except httpx.RequestError as e:
        reply_text = (
            "⏳ 超算节点正在执行耗时操作或网络延迟，连接已转入后台。请稍后主动查询进度。"
            f"（详情：{e!s}）"
        )
        logger.warning("[hpc_orchestrator] gateway network error: %s", e)
    except Exception as e:  # noqa: BLE001
        reply_text = f"HPC 网关调用异常：{e!s}"
        logger.exception("[hpc_orchestrator] unexpected error")

    if gateway_plain is not None:
        wrapped = _wrap_geek_terminal_markdown(gateway_plain)
        reply_text = wrapped if wrapped else (
            "（超算侧未返回可读文本；请确认网关 GET /tools 中是否存在该工具名，"
            "或改用 JSON 消息显式指定 tool_name/arguments。）"
        )

    if not reply_text:
        reply_text = (
            "（超算侧未返回可读文本；请确认网关 GET /tools 中是否存在该工具名，"
            "或改用 JSON 消息显式指定 tool_name/arguments。）"
        )

    yield _format_sse("message", {"content": reply_text})
    _merge_state_for_emit(state_snapshot, "message", {"content": reply_text})

    assets = _sniff_assets_from_text(reply_text)
    if assets:
        state_snapshot["assets"] = assets

    elapsed = time.time() - float(state_snapshot.get("_start_time") or time.time())
    state_snapshot["duration"] = round(elapsed, 3)
    state_snapshot.pop("_start_time", None)

    yield _format_sse("status", {"content": "回答完成", "state": "completed"})
    _merge_state_for_emit(state_snapshot, "status", {"content": "回答完成", "state": "completed"})
    yield _format_sse("done", {"status": "success"})
    yield _format_sse("state_snapshot", dict(state_snapshot))
