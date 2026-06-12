# -*- coding: utf-8 -*-
"""会话挂载目录命名与文件索引。"""
from __future__ import annotations

from datetime import datetime

from gibh_agent.core.storage.session_paths import (
    compact_folder_timestamp,
    resolve_session_folder_name,
)


def test_session_folder_appends_compact_timestamp():
    name = resolve_session_folder_name(
        "可视化标注工具",
        "sess-abc",
        folder_timestamp="202606111628",
    )
    assert name == "可视化标注工具-202606111628"


def test_session_folder_fallback_to_id_without_timestamp():
    name = resolve_session_folder_name("", "sess-abc")
    assert name == "sess-abc"


def test_session_folder_sanitizes_invalid_chars():
    name = resolve_session_folder_name(
        'bad/name:test',
        "sess-1",
        folder_timestamp="202606111628",
    )
    assert "/" not in name
    assert ":" not in name
    assert name.endswith("-202606111628")


def test_compact_folder_timestamp_from_datetime():
    ts = compact_folder_timestamp(datetime(2026, 6, 11, 16, 28))
    assert ts == "202606111628"
