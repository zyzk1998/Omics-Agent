#!/usr/bin/env python3
"""
三大组学（基因组 / 蛋白组 / 表观）闭环测试（无 HTTP）：

  FileInspector 体检
  -> SOPPlanner.generate_plan（已传 domain_name + target_steps，不触发 LLM）
  -> WorkflowExecutor.execute_workflow（真实工具注册表）

**首步质控真实计算**：`genomics_raw_qc` / `epigenomics_raw_qc_trimming` 流式读取 FASTQ(GZ)
并统计 reads、GC%、读长；`proteomics_raw_qc_conversion` 扫描 mzML 谱图数。
**下游步骤**在缺外部 CLI 时，由 `omics_derived_analysis` 基于同一输入文件做 **FASTQ 流式质量统计**
或 **mzML 元数据聚合** 的确定性代理输出（可复现、可审计，非 GATK/DIA-NN 真值）；真重型算子仍由 Worker/TaaS 承接。

运行（仓库根目录）:
  cd /home/ubuntu/GIBH-AGENT-V2 && PYTHONPATH=. python3 tests/test_omics_three_modalities_e2e.py
"""
from __future__ import annotations

import asyncio
import os
import sys
from pathlib import Path
from typing import Any, AsyncIterator, Dict, List, Optional

_ROOT = Path(__file__).resolve().parent.parent
if str(_ROOT) not in sys.path:
    sys.path.insert(0, str(_ROOT))

# 可写结果目录，避免依赖容器内 /app/results
os.environ.setdefault("RESULTS_DIR", str(_ROOT / "results" / "e2e_omics_three"))
# 宿主机无 /app/data 写权限时，避免 security_config 刷屏；亦可依赖模块内 ~/.cache 回退
os.environ.setdefault(
    "SECURITY_KEYS_PATH",
    str(_ROOT / "results" / "e2e_omics_three" / "security" / "signing_keys.json"),
)


def _find_any_markdown(steps_details: List[Dict[str, Any]]) -> Optional[str]:
    for d in steps_details or []:
        inner = d.get("step_result") or {}
        data = inner.get("data")
        if isinstance(data, dict):
            md = data.get("markdown")
            if isinstance(md, str) and md.strip():
                return md
    return None


def _qc_metrics_for_tool(steps_details: List[Dict[str, Any]], tool_id: str) -> Optional[Dict[str, Any]]:
    for d in steps_details or []:
        if d.get("tool_id") != tool_id:
            continue
        inner = d.get("step_result") or {}
        data = inner.get("data")
        if isinstance(data, dict):
            qm = data.get("qc_metrics")
            if isinstance(qm, dict) and qm:
                return qm
    return None


class _DummyLLM:
    """generate_plan 在已传 domain_name + target_steps 时不应调用 LLM；若误调则失败。"""

    async def astream(self, *a: Any, **kw: Any) -> AsyncIterator[Dict[str, Any]]:
        if False:  # pragma: no cover
            yield {}

    async def achat(self, *a: Any, **kw: Any) -> Any:
        raise RuntimeError("E2E 不应调用 achat")

    async def astream_chat(self, *a: Any, **kw: Any) -> AsyncIterator[Any]:
        if False:  # pragma: no cover
            yield None


CASES: List[Dict[str, Any]] = [
    {
        "domain": "genomics",
        "rel_path": "test_data/genomics/sample1_R1.fastq.gz",
        "prompt": "请帮我用这个 FASTQ 跑一下基因组胚系分析全流程",
        "qc_tool": "genomics_raw_qc",
        "qc_kind": "fastq",
    },
    {
        "domain": "proteomics",
        "rel_path": "test_data/proteomics/BSA1_F1.mzML",
        "prompt": "请用这份 mzML 跑蛋白质组学搜库定量全流程",
        "qc_tool": "proteomics_raw_qc_conversion",
        "qc_kind": "mzml",
    },
    {
        "domain": "epigenomics",
        "rel_path": "test_data/epigenomics/SRR1822153_1.fastq.gz",
        "prompt": "请用这个 FASTQ 跑 ATAC-seq/ChIP 主线表观组学全流程",
        "qc_tool": "epigenomics_raw_qc_trimming",
        "qc_kind": "fastq",
    },
]


async def _run_domain(domain: str, abs_file: str, prompt: str) -> Dict[str, Any]:
    from gibh_agent.tools import load_all_tools

    load_all_tools()

    from gibh_agent.core.file_inspector import FileInspector
    from gibh_agent.core.planner import SOPPlanner
    from gibh_agent.core.tool_retriever import ToolRetriever
    from gibh_agent.core.executor import WorkflowExecutor
    from gibh_agent.core.workflows.registry import WorkflowRegistry

    inspector = FileInspector(upload_dir=str(_ROOT))
    file_metadata = inspector.inspect_file(abs_file)
    assert file_metadata.get("status") == "success", file_metadata

    registry = WorkflowRegistry()
    wf = registry.get_workflow(domain)
    assert wf is not None, domain

    planner = SOPPlanner(ToolRetriever(), _DummyLLM())  # type: ignore[arg-type]
    keys = list(wf.steps_dag.keys())

    plan_result: Optional[Dict[str, Any]] = None
    async for ev, data in planner.generate_plan(
        user_query=prompt,
        file_metadata=file_metadata,
        domain_name=domain,
        target_steps=keys,
        is_template=False,
    ):
        if ev == "workflow":
            plan_result = data

    assert isinstance(plan_result, dict), plan_result
    assert plan_result.get("type") == "workflow_config", plan_result
    wd = plan_result.get("workflow_data") or {}
    assert wd.get("steps"), "workflow_data.steps 为空"

    ex = WorkflowExecutor(upload_dir=str(_ROOT))
    report = ex.execute_workflow(wd, file_paths=[abs_file])
    return report


async def main_async() -> int:
    ok = 0
    for i, case in enumerate(CASES, 1):
        fp = (_ROOT / case["rel_path"]).resolve()
        print("=" * 80)
        print(f"[{i}/3] domain={case['domain']} file={fp}")
        print("=" * 80)
        assert fp.is_file(), f"缺少测试文件: {fp}"

        report = await _run_domain(case["domain"], str(fp), case["prompt"])
        st = report.get("status")
        details = report.get("steps_details") or []

        if st != "success":
            print(f"❌ workflow_status={st!r}")
            print(report)
            raise SystemExit(1)

        md = _find_any_markdown(details)
        if not md:
            print("❌ 未在任何步骤的 step_result.data 中发现 markdown 字段")
            raise SystemExit(1)

        qtool = case["qc_tool"]
        qkind = case["qc_kind"]
        qm = _qc_metrics_for_tool(details, qtool)
        if not qm:
            print(f"❌ 未找到首步算子 {qtool} 的 qc_metrics（真实统计缺失）")
            raise SystemExit(1)
        if qkind == "fastq":
            nr = int(qm.get("n_reads") or 0)
            gc = float(qm.get("gc_percent") or 0.0)
            if nr < 1:
                print("❌ FASTQ qc_metrics.n_reads 无效")
                raise SystemExit(1)
            print(f"   📊 真实 FASTQ 统计: reads={nr}, GC%={gc}")
        elif qkind == "mzml":
            ns = int(qm.get("spectrum_count") or 0)
            if ns < 1:
                print("❌ mzML qc_metrics.spectrum_count 无效")
                raise SystemExit(1)
            print(f"   📊 真实 mzML 统计: spectra={ns}, MB={qm.get('file_size_mb')}")

        ok += 1
        print(
            f"✅ 闭环 ({ok}/3) | workflow_status={st} | steps={len(details)} | "
            f"markdown_len={len(md)} | 首步 {qtool} 含 qc_metrics"
        )

    print()
    print("=" * 80)
    print(
        f"结束：{ok}/3 组域 planner→executor 跑通，首步质控由真实 I/O 计算；"
        " 下游在缺 CLI 时输出基于输入的轻量代理指标。重型真管线请在 Worker 环境验收。"
    )
    print("=" * 80)
    return 0


def main() -> None:
    raise SystemExit(asyncio.run(main_async()))


if __name__ == "__main__":
    main()
