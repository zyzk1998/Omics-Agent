#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
BepiPred3 技能全链路实弹演习：Orchestrator 聊天模式 + LLM 工具调用 + bepipred3_prediction。

运行:
  pytest -s tests/test_bepipred3_skill.py
  python tests/test_bepipred3_skill.py

环境:
  需已配置 SILICONFLOW_API_KEY（与 gibh_agent/config/settings.yaml 一致），无需启动 FastAPI。

【BepiPred3 执行模式】
- 默认：`protein_tools.bepipred3_prediction` 在子进程中调用本地 `third_party/BepiPred-3.0`（或 BEPIPRED3_ROOT）下的 CLI，结果写入 RESULTS_DIR 并通过 /results 静态路径暴露链接。
- 无本地算子时测试会失败或走错误提示；可临时 `export BEPIPRED3_USE_REMOTE_API=1` 回退磐石 HTTP（需网络）。
"""
from __future__ import annotations

import asyncio
import json
import logging
import os
import re
import sys
import inspect
from pathlib import Path
from typing import Any, Dict, List

# 仓库根目录
ROOT = Path(__file__).resolve().parent.parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

logging.basicConfig(level=logging.INFO, format="%(levelname)s %(name)s: %(message)s")
logger = logging.getLogger("test_bepipred3_skill")

TEST_PROMPT = """[系统注入：用户调用了 BepiPred3 B细胞表位预测 技能]
请帮我预测以下蛋白质序列的B细胞表位，并按照高置信度（前20%）的标准筛选，不使用顺序平滑。序列如下：
>7lj4_B
RSTTLLALLALVLLYVSGALVFRALEQPHEQQAQRELGEVREKFLRAHPCVSDQELGLLIKEVADALGGGADPETQSTSHSAWDLGSAFFFSGTIITTIGYGNVALRTDAGRLFCIFYALVGIPLFGDILLAGVGDRLGSSLRHGIGHIEAIFLKWHVPPELVRVLSAEMLFLLIGCLLFVLTPTFVFCYMEDWSKLEAIYFVIVTLTTVGFGDYVAGADPRQDSPAYQPLVWFWILLGLPAYFASVLTTIGNWLRVVS
"""

EXPECTED_HEADER = ">7lj4_B"
EXPECTED_SEQ_PREFIX = "RSTTLLALLALVLLYVSGALVFRALEQPHEQQAQRELGEVREKFLRAHPCVSDQELGLLIKEVADALGGGADPETQSTSHSAWDLGSAFFFSGTIITTIGYGNVALRTDAGRLFCIFYALVGIPLFGDILLAGVGDRLGSSLRHGIGHIEAIFLKWHVPPELVRVLSAEMLFLLIGCLLFVLTPTFVFCYMEDWSKLEAIYFVIVTLTTVGFGDYVAGADPRQDSPAYQPLVWFWILLGLPAYFASVLTTIGNWLRVVS"


def _require_llm_env() -> None:
    if not (os.getenv("SILICONFLOW_API_KEY") or "").strip():
        try:
            import pytest

            pytest.skip("未设置 SILICONFLOW_API_KEY，跳过实弹 LLM 调用")
        except ImportError:
            print("SKIP: 未设置 SILICONFLOW_API_KEY，跳过实弹 LLM 调用", file=sys.stderr)
            sys.exit(0)


def _parse_sse_accumulated_text(raw: str) -> str:
    """从 SSE 原始串中粗略拼接 message/thought 的可见文本（用于断言 URL 是否被模型复述）。"""
    texts: List[str] = []
    for block in raw.split("\n\n"):
        if "event: message" in block or "event: thought" in block:
            for line in block.split("\n"):
                if line.startswith("data:"):
                    try:
                        payload = json.loads(line[5:].strip())
                        c = (payload.get("content") or "") if isinstance(payload, dict) else ""
                        if c:
                            texts.append(str(c))
                    except json.JSONDecodeError:
                        continue
    return "\n".join(texts)


async def _run_orchestrator_stream(
    orchestrator: Any,
    *,
    model_name: str,
) -> str:
    parts: List[str] = []
    async for chunk in orchestrator.stream_process(
        query=TEST_PROMPT,
        files=[],
        history=[],
        session_id="e2e-bepipred3-test",
        model_name=model_name,
        user_id="test",
        owner_id=None,
        db=None,
    ):
        if isinstance(chunk, str):
            parts.append(chunk)
    return "".join(parts)


async def _async_main() -> None:
    _require_llm_env()

    # 可写 Chroma 目录（沙箱 / CI 下避免 ./data/chroma_tools 不可写导致落入 Task 后失败）
    _chroma_dir = ROOT / "tests" / "_chroma_tools"
    _chroma_dir.mkdir(parents=True, exist_ok=True)
    os.environ["CHROMA_TOOLS_PERSIST_DIR"] = str(_chroma_dir)

    # 确保工具已注册
    import gibh_agent.tools.protein_tools  # noqa: F401
    from gibh_agent.core.orchestrator import AgentOrchestrator
    from gibh_agent.core.tool_registry import registry
    from gibh_agent.main import GIBHAgent

    config_path = ROOT / "gibh_agent" / "config" / "settings.yaml"
    if not config_path.is_file():
        raise RuntimeError(f"缺少配置文件: {config_path}")

    tool_calls: List[Dict[str, Any]] = []
    orig_get_tool = registry.get_tool

    def tracking_get_tool(name: str):
        fn = orig_get_tool(name)
        if name != "bepipred3_prediction" or fn is None:
            return fn

        if inspect.iscoroutinefunction(fn):

            async def wrapped(**kwargs):
                tool_calls.append(dict(kwargs))
                return await fn(**kwargs)

        else:

            def wrapped(**kwargs):
                tool_calls.append(dict(kwargs))
                return fn(**kwargs)

        return wrapped

    registry.get_tool = tracking_get_tool  # type: ignore[assignment]

    try:
        agent = GIBHAgent(str(config_path))
        orch = AgentOrchestrator(agent, upload_dir=str(ROOT / "tests" / "_asset_tmp"))
        # stream_process 默认 qwen3.5-plus 在硅基流动上常无效；与 settings.yaml / SILICONFLOW_MODEL 对齐
        model_name = (
            (os.getenv("BEP_PRED_TEST_MODEL") or "").strip()
            or (os.getenv("SILICONFLOW_MODEL") or "").strip()
            or "deepseek-ai/DeepSeek-R1"
        )

        logger.info("使用模型 model_name=%s", model_name)
        raw_sse = await _run_orchestrator_stream(orch, model_name=model_name)
    finally:
        registry.get_tool = orig_get_tool  # type: ignore[assignment]

    print("\n========== SSE 原始输出（节选前后）==========")
    head = raw_sse[:2500] if len(raw_sse) > 2500 else raw_sse
    print(head)
    if len(raw_sse) > 2500:
        print("\n... [省略中间] ...\n")
        print(raw_sse[-2500:])
    print("========== 结束 ==========\n")

    visible = _parse_sse_accumulated_text(raw_sse) + raw_sse

    # 断言 1：工具被调用
    assert tool_calls, "未记录到 bepipred3_prediction 调用（registry.get_tool 包装未触发）"
    kwargs0 = tool_calls[0]

    # 断言 2：参数
    assert kwargs0.get("top_epitope_percentage_cutoff") == "top_20", kwargs0
    assert kwargs0.get("use_sequential_smoothing") is False, kwargs0
    seq = (kwargs0.get("sequence_or_path") or "").strip()
    assert EXPECTED_HEADER in seq, f"FASTA 头缺失: keys={list(kwargs0.keys())}"
    assert EXPECTED_SEQ_PREFIX in seq.replace("\n", "").replace("\r", ""), "序列主体未传入工具"

    # 断言 3：成功链路含本地 /results 链接、远程 API 域名、或失败链路友好提示
    assert any(
        token in visible
        for token in (
            "/results/bepipred3/",
            "output_interactive_figures.html",
            "raw_output.csv",
            "bepipred3_results.zip",
            "120.220.102.26",
            "result.html",
            "result.csv",
            "技能执行未完成",
            "远程计算服务当前资源繁忙",
            "查看交互式图表",
            "下载原始数据",
        )
    ), "最终输出中未出现结果 URL 或友好错误提示，请检查工具返回与 SkillAgent 渲染"

    logger.info("全链路断言通过: tool_calls=%s", len(tool_calls))


def test_bepipred3_skill_e2e() -> None:
    """pytest 入口（同步包装 asyncio）。"""
    _require_llm_env()
    asyncio.run(_async_main())


if __name__ == "__main__":
    asyncio.run(_async_main())
    print("OK: BepiPred3 技能全链路实弹通过", file=sys.stderr)
