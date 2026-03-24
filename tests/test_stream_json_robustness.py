#!/usr/bin/env python3
"""
流式 JSON 极限压测与回归：碎片化 chunk、RNA/代谢组 generate_plan 结构、_analyze_user_intent 异步生成器契约。

运行：
  python tests/test_stream_json_robustness.py
或：
  pytest tests/test_stream_json_robustness.py -v --asyncio-mode=auto
"""
from __future__ import annotations

import asyncio
import json
import logging
import random
import sys
from io import StringIO
from pathlib import Path
from typing import Any, AsyncIterator, Dict, List, Optional, Tuple

# 项目根
_ROOT = Path(__file__).resolve().parent.parent
if str(_ROOT) not in sys.path:
    sys.path.insert(0, str(_ROOT))

from gibh_agent.core.stream_utils import THINK_CLOSE, THINK_OPEN, stream_and_extract_json
from gibh_agent.core.planner import SOPPlanner
from gibh_agent.core.tool_retriever import ToolRetriever
from gibh_agent.core.workflows.registry import WorkflowRegistry


def _chunk_dict(content: str, reasoning: Optional[str] = None) -> Dict[str, Any]:
    delta: Dict[str, Any] = {"content": content}
    if reasoning is not None:
        delta["reasoning_content"] = reasoning
    return {"choices": [{"delta": delta}]}


async def _micro_chunk_stream(
    text: str,
    rng: random.Random,
    low: int = 1,
    high: int = 2,
) -> AsyncIterator[Dict[str, Any]]:
    """模拟恶劣网络：每次仅推送 1～2 个字符（内容在 delta.content）。"""
    i = 0
    n = len(text)
    while i < n:
        take = rng.randint(low, high) if high >= low else low
        piece = text[i : i + take]
        yield _chunk_dict(piece)
        i += take


async def _async_json_fragmentation_stress() -> None:
    """体检一：前后废话 + think 标签 + 复杂 JSON，1～2 字符/chunk，必须解析成功。"""
    complex_obj: Dict[str, Any] = {
        "workflow_name": "RNA 压测工作流",
        "steps": [
            {
                "id": f"step_{k}",
                "tool_id": f"tool_{k}",
                "params": {"nested": {"arr": [1, 2, 3], "s": "x\"}"}},
            }
            for k in range(5)
        ],
        "meta": {"tags": ["单细胞", "RNA"], "ok": True},
    }
    json_core = json.dumps(complex_obj, ensure_ascii=False)
    preamble = "好的，配置如下（请忽略废话）：\n"
    thought_inner = "模型内心：先核对用户是否要完整流程…"
    suffix = "\n以上就是全部内容，谢谢。"
    full_text = f"{preamble}{THINK_OPEN}{thought_inner}{THINK_CLOSE}{json_core}{suffix}"

    rng = random.Random(42)
    stream = _micro_chunk_stream(full_text, rng, low=1, high=2)

    thoughts: List[str] = []
    parsed: Optional[Any] = None
    errors: List[Any] = []

    async for ev, data in stream_and_extract_json(stream):
        if ev == "thought":
            c = (data or {}).get("content") if isinstance(data, dict) else ""
            thoughts.append(str(c))
        elif ev == "json":
            parsed = data
        elif ev == "json_error":
            errors.append(data)

    assert not errors, f"不应产生 json_error，实际: {errors}"
    assert parsed == complex_obj, "解析结果必须与原始对象一致"
    merged_thought = "".join(thoughts)
    assert thought_inner in merged_thought or len(merged_thought) > 0, "应透出思考片段"


async def _async_rna_generate_plan_structure_offline() -> None:
    """
    体检二（离线）：不调用真实大模型；domain_name + 全量 target_steps 走 generate_plan，
    验证 RNA 无文件模板为 13 步（已剔除 CellRanger 两键），结构字段齐全。
    """
    class _DummyLLM:
        async def astream(self, *a, **kw):
            if False:
                yield _chunk_dict("")

        async def achat(self, *a, **kw):
            raise RuntimeError("本测例不应调用 achat")

    registry = WorkflowRegistry()
    rna_wf = registry.get_workflow("RNA")
    assert rna_wf is not None
    full_keys = list(rna_wf.steps_dag.keys())
    planner = SOPPlanner(ToolRetriever(), _DummyLLM())

    result = None
    async for ev, data in planner.generate_plan(
        user_query="帮我做单细胞转录组测序分析（RNA）",
        file_metadata=None,
        domain_name="RNA",
        target_steps=full_keys,
        is_template=True,
    ):
        if ev == "workflow":
            result = data

    assert result is not None
    assert result.get("type") == "workflow_config"
    assert result.get("template_mode") is True
    wd = result.get("workflow_data") or {}
    steps = wd.get("steps") or []
    # 无 FASTQ：RNA 模板剔除 rna_cellranger_count / rna_convert_cellranger_to_h5ad → 15 - 2 = 13
    assert len(steps) == 13, f"期望 13 步（非 FASTQ 全量），实际 {len(steps)}"
    for s in steps:
        assert s.get("id"), "每步需有 id"
        assert s.get("tool_id"), "每步需有 tool_id"
        assert "params" in s
    assert wd.get("workflow_name") or wd.get("name")


async def _async_metabolomics_generate_plan_structure_offline() -> None:
    """代谢组学：同结构烟测（无 LLM）。"""
    class _DummyLLM:
        async def astream(self, *a, **kw):
            if False:
                yield _chunk_dict("")

        async def achat(self, *a, **kw):
            raise RuntimeError("本测例不应调用 achat")

    registry = WorkflowRegistry()
    wf = registry.get_workflow("Metabolomics")
    assert wf is not None
    full_keys = list(wf.steps_dag.keys())
    planner = SOPPlanner(ToolRetriever(), _DummyLLM())
    result = None
    async for ev, data in planner.generate_plan(
        user_query="代谢组学分析",
        file_metadata=None,
        domain_name="Metabolomics",
        target_steps=full_keys,
        is_template=True,
    ):
        if ev == "workflow":
            result = data
    assert result and result.get("type") == "workflow_config"
    assert result.get("template_mode") is True
    assert len((result.get("workflow_data") or {}).get("steps") or []) > 0


async def _async_analyze_user_intent_yields_intent_result() -> None:
    """_analyze_user_intent 为异步生成器：至少 yield 一次 intent_result。"""
    from types import SimpleNamespace

    class _LLM:
        async def achat(self, *a, **kw):
            return SimpleNamespace(
                choices=[
                    SimpleNamespace(
                        message=SimpleNamespace(
                            content='{"target_steps": ["rna_pca"], "skip_steps": []}'
                        )
                    )
                ]
            )

        async def astream(self, *a, **kw):
            raise RuntimeError("不应走 astream")

    registry = WorkflowRegistry()
    rna_wf = registry.get_workflow("RNA")
    planner = SOPPlanner(ToolRetriever(), _LLM())

    events: List[Tuple[str, Any]] = []
    async for ev, payload in planner._analyze_user_intent("只做 PCA", rna_wf):
        events.append((ev, payload))

    intent_results = [e for e in events if e[0] == "intent_result"]
    assert len(intent_results) == 1
    ir = intent_results[0][1]
    assert isinstance(ir, dict)
    assert "rna_pca" in (ir.get("target_steps") or [])


async def _run_all() -> None:
    await _async_json_fragmentation_stress()
    await _async_rna_generate_plan_structure_offline()
    await _async_metabolomics_generate_plan_structure_offline()
    await _async_analyze_user_intent_yields_intent_result()
    print("✅ test_stream_json_robustness：全部异步测例通过")


def test_json_fragmentation_stress() -> None:
    asyncio.run(_async_json_fragmentation_stress())


def test_rna_generate_plan_structure_offline() -> None:
    asyncio.run(_async_rna_generate_plan_structure_offline())


def test_metabolomics_generate_plan_structure_offline() -> None:
    asyncio.run(_async_metabolomics_generate_plan_structure_offline())


def test_analyze_user_intent_yields_intent_result() -> None:
    asyncio.run(_async_analyze_user_intent_yields_intent_result())


def main() -> int:
    asyncio.run(_run_all())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
