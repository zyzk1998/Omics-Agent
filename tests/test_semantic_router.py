# -*- coding: utf-8 -*-
"""
SemanticRouter 单元测试（Mock LLM，不触网）。

用例 A/B/C 对应「总指挥」关心的三条业务语义：HPC、缺文件澄清、有文件走任务。
Fail-safe 用例验证：解析失败 + 低置信 等路径最终降级为 clarify。

说明：当前环境未必安装 pytest-asyncio，故用 asyncio.run 执行协程，仍属基于 pytest 的异步逻辑验证。
"""
from __future__ import annotations

import asyncio
import json
import sys
from pathlib import Path
from typing import Any
from unittest.mock import AsyncMock, MagicMock

import pytest

ROOT = Path(__file__).resolve().parent.parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from gibh_agent.core.semantic_router import (  # noqa: E402
    RouteKind,
    RouterInput,
    SemanticRouter,
)


def _fake_completion(content: str) -> Any:
    """构造与 OpenAI ChatCompletion 兼容的最小只读结构（仅 tests 使用）。"""
    msg = MagicMock()
    msg.content = content
    ch = MagicMock()
    ch.message = msg
    comp = MagicMock()
    comp.choices = [ch]
    return comp


@pytest.fixture
def mock_llm() -> MagicMock:
    llm = MagicMock()
    llm.achat = AsyncMock()
    return llm


def test_route_hpc_pestat_compute_scheduler_on(mock_llm):
    """
    用例 A：超算查询语义 + compute_scheduler 开启 → 应路由到 hpc。

    真实环境中由模型根据 System Prompt 产出；此处 Mock 模拟「理想模型」输出，
    验证本模块能正确解析 JSON 并返回 RouterOutput。
    """

    async def _body():
        payload = {
            "route": "hpc",
            "confidence": 0.92,
            "rationale_short": "用户查询 pestat 且调度 MCP 已开启，属集群状态类问题。",
        }
        mock_llm.achat.return_value = _fake_completion(
            json.dumps(payload, ensure_ascii=False)
        )

        router = SemanticRouter(mock_llm)
        inp = RouterInput(
            query="pestat",
            file_status={},  # 无文件不妨碍 HPC 状态查询
            mcp_status={"compute_scheduler": True},
            session_flags={},
        )
        out = await router.decide_route(inp)
        assert out.route == RouteKind.hpc
        assert out.confidence >= 0.7
        mock_llm.achat.assert_awaited()

    asyncio.run(_body())


def test_route_clarify_metabolomics_qc_no_files(mock_llm):
    """
    用例 B：代谢组 QC 属执行型诉求，但无文件 → 应 clarify（缺前提）。

    Mock 返回 clarify 模拟模型遵守 Prompt 红线；同时覆盖「无 file_status 键」场景。
    """

    async def _body():
        payload = {
            "route": "clarify",
            "confidence": 0.88,
            "rationale_short": "需要数据文件才能做代谢组 QC，请上传或指定路径。",
        }
        mock_llm.achat.return_value = _fake_completion(
            json.dumps(payload, ensure_ascii=False)
        )

        router = SemanticRouter(mock_llm)
        inp = RouterInput(
            query="帮我做个代谢组QC",
            file_status={"has_files": False},
            mcp_status={},
            session_flags={},
        )
        out = await router.decide_route(inp)
        assert out.route == RouteKind.clarify

    asyncio.run(_body())


def test_route_task_diff_analysis_with_fastq(mock_llm):
    """
    用例 C：明确要跑差异分析且已有 fastq → task。

    有文件时执行意图应走任务编排通道，而非 clarify。
    """

    async def _body():
        payload = {
            "route": "task",
            "confidence": 0.91,
            "rationale_short": "用户请求差异分析且已提供 .fastq 数据，可走分析任务链。",
        }
        mock_llm.achat.return_value = _fake_completion(
            json.dumps(payload, ensure_ascii=False)
        )

        router = SemanticRouter(mock_llm)
        inp = RouterInput(
            query="跑一下这两个数据的差异分析",
            file_status={"has_files": True, "types": [".fastq"]},
            mcp_status={},
            session_flags={},
        )
        out = await router.decide_route(inp)
        assert out.route == RouteKind.task

    asyncio.run(_body())


def test_failsafe_retry_then_clarify(mock_llm):
    """
    Fail-safe：第一次坏 JSON / 第二次仍不可信 → 最终 clarify。

    第一次：不可解析；第二次：合法 JSON 但置信度 0.2（低于 0.7）。
    预期两次 achat 后降级 clarify，且不应抛异常。
    """

    async def _body():
        bad = _fake_completion("这不是 JSON")
        low_conf = json.dumps(
            {
                "route": "task",
                "confidence": 0.2,
                "rationale_short": "模型胡猜",
            },
            ensure_ascii=False,
        )
        mock_llm.achat.side_effect = [bad, _fake_completion(low_conf)]

        router = SemanticRouter(mock_llm)
        inp = RouterInput(query="任意", file_status={}, mcp_status={}, session_flags={})
        out = await router.decide_route(inp)
        assert out.route == RouteKind.clarify
        assert out.confidence == 0.0
        assert mock_llm.achat.await_count == 2

    asyncio.run(_body())


def test_failsafe_retry_recover_after_bad_json(mock_llm):
    """第一次解析失败，第二次返回高置信 hpc → 应直接成功，不再 clarify。"""

    async def _body():
        good = json.dumps(
            {
                "route": "hpc",
                "confidence": 0.95,
                "rationale_short": "第二次重试成功",
            },
            ensure_ascii=False,
        )
        mock_llm.achat.side_effect = [
            _fake_completion("```\noops\n```"),
            _fake_completion(good),
        ]

        router = SemanticRouter(mock_llm)
        inp = RouterInput(
            query="squeue",
            file_status={},
            mcp_status={"compute_scheduler": True},
            session_flags={},
        )
        out = await router.decide_route(inp)
        assert out.route == RouteKind.hpc
        assert mock_llm.achat.await_count == 2

    asyncio.run(_body())
