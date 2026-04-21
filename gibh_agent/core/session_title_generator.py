"""
首轮会话标题：流式主响应结束后异步调用轻量模型生成，不阻塞 TTFT。
"""
from __future__ import annotations

import asyncio
import logging
import os
import re
from typing import Optional

from gibh_agent.core.llm_client import LLMClientFactory
from gibh_agent.db.connection import SessionLocal
from gibh_agent.db.models import Session as SessionModel

logger = logging.getLogger(__name__)

_PROMPT = (
    "请根据以下用户的输入，提取一个精炼、专业的标题（不超过 10 个字），不要加标点符号。用户输入：\n{user_query}"
)


def _sanitize_title(raw: str) -> str:
    s = (raw or "").strip()
    s = re.sub(r"[\s\n\r\t]+", "", s)
    s = re.sub(
        r'[。，、！？;；:：""''《》<>\[\]（）()·—…\.,!\?\-_]+',
        "",
        s,
    )
    if not s:
        return ""
    if len(s) > 10:
        s = s[:10]
    return s


def generate_session_title_sync(session_id: str, user_query: str, owner_id: str) -> Optional[str]:
    """同步生成并写入 sessions.title；在线程池中调用。"""
    q = (user_query or "").strip()
    if not q or not session_id or not owner_id:
        return None
    model = (os.getenv("SESSION_TITLE_MODEL") or "deepseek-chat").strip()
    try:
        client = LLMClientFactory.create_for_model(
            model,
            temperature=0.3,
            max_tokens=64,
            timeout=45.0,
        )
        messages = [{"role": "user", "content": _PROMPT.format(user_query=q[:6000])}]
        completion = client.chat(messages, stream=False, temperature=0.3, max_tokens=64)
        text = ""
        if completion.choices:
            msg = completion.choices[0].message
            if msg is not None:
                text = (getattr(msg, "content", None) or "").strip()
        title = _sanitize_title(text)
        if not title:
            return None
        db = SessionLocal()
        try:
            row = (
                db.query(SessionModel)
                .filter(SessionModel.id == session_id, SessionModel.owner_id == owner_id)
                .first()
            )
            if not row:
                return None
            row.title = title[:512]
            db.commit()
            logger.info("✅ [SessionTitle] 已更新 session=%s title=%s", session_id, title)
            return title
        finally:
            db.close()
    except Exception as e:
        logger.warning("⚠️ [SessionTitle] 生成失败: %s", e, exc_info=True)
        return None


async def run_session_title_job(session_id: str, user_query: str, owner_id: str) -> Optional[str]:
    """异步包装：LLM 与 DB 在线程池执行，不阻塞事件循环。"""
    return await asyncio.to_thread(
        generate_session_title_sync, session_id, user_query, owner_id
    )
