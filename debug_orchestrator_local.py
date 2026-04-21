#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
深水炸弹：本地复现 /api/chat 等价路径，强制打印完整 Traceback。
用法（项目根目录）:
  python debug_orchestrator_local.py
可选:
  export GIBH_ENABLE_SEMANTIC_ROUTER=0   # 若语义路由干扰，可关
"""
from __future__ import annotations

import asyncio
import json
import os
import re
import sys
import traceback
from pathlib import Path

_here = Path(__file__).resolve().parent
# 脚本拷到容器 /tmp 时不能以 /tmp 为仓库根；优先同目录含 gibh_agent，否则 GIBH_AGENT_ROOT（默认 /app）
if (_here / "gibh_agent").is_dir():
    ROOT = _here
else:
    ROOT = Path(os.environ.get("GIBH_AGENT_ROOT", "/app")).resolve()
    if not (ROOT / "gibh_agent").is_dir():
        ROOT = _here
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

_ENV_PATH = ROOT / ".env"
try:
    from dotenv import load_dotenv

    load_dotenv(_ENV_PATH)
except ImportError:
    # 无 python-dotenv 时尽力解析 KEY=VAL 行
    if _ENV_PATH.is_file():
        for line in _ENV_PATH.read_text(encoding="utf-8", errors="replace").splitlines():
            s = line.strip()
            if not s or s.startswith("#") or "=" not in s:
                continue
            k, _, v = s.partition("=")
            k, v = k.strip(), v.strip().strip('"').strip("'")
            if k and k not in os.environ:
                os.environ[k] = v


def _parse_sse_blocks(sse_text: str) -> tuple[list[str], list[str], list[str], list[dict]]:
    """返回 (message_chunks, thought_chunks, event_types_seen, error_payloads)。"""
    messages: list[str] = []
    thoughts: list[str] = []
    events: list[str] = []
    errors: list[dict] = []
    cur_event = ""
    for block in sse_text.split("\n\n"):
        if not block.strip():
            continue
        for line in block.split("\n"):
            if line.startswith("event:"):
                cur_event = line[6:].strip()
                if cur_event and cur_event not in events:
                    events.append(cur_event)
            if not line.startswith("data:"):
                continue
            raw = line[5:].strip()
            if raw in ("", "[DONE]"):
                continue
            try:
                obj = json.loads(raw)
            except json.JSONDecodeError:
                continue
            if not isinstance(obj, dict):
                continue
            if obj.get("error"):
                errors.append(obj)
            c = obj.get("content")
            if isinstance(c, str) and c:
                if cur_event == "thought":
                    thoughts.append(c)
                elif cur_event == "message":
                    messages.append(c)
                elif cur_event == "error":
                    pass
                else:
                    # 无 event 行时，有 content 的当作 message 兼容
                    if cur_event == "" and "error" not in obj:
                        messages.append(c)
    return messages, thoughts, events, errors


async def main() -> int:
    print("=" * 72)
    print("debug_orchestrator_local: 加载 .env 后直连 stream_process(chat)")
    print(f"  ROOT={ROOT}")
    print(f"  DEEPSEEK_API_KEY set: {bool((os.getenv('DEEPSEEK_API_KEY') or '').strip())}")
    print(f"  LLM_CLOUD_PROVIDER={os.getenv('LLM_CLOUD_PROVIDER', '')!r}")
    print("=" * 72)

    try:
        from gibh_agent.main import GIBHAgent
        from gibh_agent.core.orchestrator import AgentOrchestrator
    except Exception as e:
        print("❌ 导入失败:", e)
        traceback.print_exc()
        return 2

    try:
        agent = GIBHAgent(str(ROOT / "gibh_agent" / "config" / "settings.yaml"))
        orch = AgentOrchestrator(agent, upload_dir=str(ROOT / "tests" / "_debug_uploads"))
    except Exception as e:
        print("❌ 初始化 GIBHAgent / Orchestrator 失败:", e)
        traceback.print_exc()
        return 2

    payload = {
        "query": "你好，请用一句话回复。",
        "files": [],
        "history": [],
        "model_name": "deepseek-reasoner",
        "session_id": "debug-orchestrator-local",
        "user_id": "debug-user",
        "owner_id": None,
        "db": None,
        "enabled_mcps": ["web_search", "authority_db"],
        "thinking_mode": "fast",
    }

    raw_stream_parts: list[str] = []
    message_parts: list[str] = []
    thought_parts: list[str] = []
    all_events: list[str] = []
    error_events: list[dict] = []

    try:
        async for chunk in orch.stream_process(**payload):
            if not isinstance(chunk, str):
                continue
            raw_stream_parts.append(chunk)
            msgs, ths, evs, errs = _parse_sse_blocks(chunk)
            message_parts.extend(msgs)
            thought_parts.extend(ths)
            for e in evs:
                if e not in all_events:
                    all_events.append(e)
            error_events.extend(errs)
            for piece in msgs:
                sys.stdout.write(piece)
                sys.stdout.flush()
        print()
    except Exception:
        print("\n❌ stream_process 外层异常:")
        traceback.print_exc()
        return 1

    if error_events:
        print("\n❌ SSE error 事件:")
        for ev in error_events:
            print(json.dumps(ev, ensure_ascii=False, indent=2))
        return 1

    joined_msg = "".join(message_parts).strip()
    joined_th = "".join(thought_parts).strip()
    print("\n--- SSE event 类型 ---\n", ", ".join(all_events) or "(none)")
    if joined_msg:
        print("\n✅ message 文本（节选）:", joined_msg[:800])
        return 0
    if joined_th:
        print("\n✅ 仅收到 thought/reasoning（节选，Reasoner 常见）:", joined_th[:800])
        return 0

    tail = "".join(raw_stream_parts)[-2500:]
    print("\n⚠️ 未解析到 message/thought。原始流尾部:\n", tail)
    return 1


if __name__ == "__main__":
    try:
        rc = asyncio.run(main())
    except Exception:
        print("❌ asyncio.run 级异常:")
        traceback.print_exc()
        rc = 1
    sys.exit(rc)
