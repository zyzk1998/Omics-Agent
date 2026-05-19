#!/usr/bin/env python3
"""
基因组真实管线 E2E（无 Mock 贯通）：fastp → bwa → samtools → bcftools/GATK → 临床报告。

前置：
  - data/references/genomics/hg38.fa + bwa 索引（或运行 scripts/init_omics_mock_references.py）
  - test_data/genomics/sample1_R1.fastq.gz
  - 宿主 PATH 含 fastp/bwa/samtools/bcftools（API 镜像已 apt 安装）

运行:
  cd /home/ubuntu/GIBH-AGENT-V2 && GIBH_E2E_EXECUTOR_LOCAL=1 PYTHONPATH=. python3 tests/test_genomics_real_pipeline_e2e.py
"""
from __future__ import annotations

import asyncio
import json
import os
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional

_ROOT = Path(__file__).resolve().parent.parent
if str(_ROOT) not in sys.path:
    sys.path.insert(0, str(_ROOT))

os.environ.setdefault("RESULTS_DIR", str(_ROOT / "results" / "e2e_genomics_real"))
os.environ.setdefault("OMICS_REF_DIR", str(_ROOT / "data" / "references"))
os.environ.setdefault(
    "GIBH_REF_HG38", str(_ROOT / "data" / "references" / "genomics" / "hg38.fa")
)

FASTQ = _ROOT / "test_data" / "genomics" / "sample1_R1.fastq.gz"
PROMPT = "请帮我用这个 FASTQ 跑一下基因组胚系分析全流程"


class _DummyLLM:
    async def astream(self, *a: Any, **kw: Any):
        if False:
            yield {}

    async def achat(self, *a: Any, **kw: Any) -> Any:
        raise RuntimeError("E2E 不应调用 achat")

    async def astream_chat(self, *a: Any, **kw: Any):
        if False:
            yield None


async def _plan() -> Dict[str, Any]:
    from gibh_agent.tools import load_all_tools

    load_all_tools()
    from gibh_agent.core.file_inspector import FileInspector
    from gibh_agent.core.planner import SOPPlanner
    from gibh_agent.core.tool_retriever import ToolRetriever
    from gibh_agent.core.workflows.registry import WorkflowRegistry

    fp = str(FASTQ.resolve())
    inspector = FileInspector(upload_dir=str(_ROOT))
    meta = inspector.inspect_file(fp)
    assert meta.get("status") == "success", meta

    wf = WorkflowRegistry().get_workflow("genomics")
    assert wf is not None
    keys = list(wf.get_steps_dag().keys())

    planner = SOPPlanner(ToolRetriever(), _DummyLLM())  # type: ignore[arg-type]
    plan: Optional[Dict[str, Any]] = None
    async for ev, data in planner.generate_plan(
        user_query=PROMPT,
        file_metadata=meta,
        domain_name="genomics",
        target_steps=keys,
        is_template=False,
    ):
        if ev == "workflow":
            plan = data
    assert isinstance(plan, dict) and plan.get("type") == "workflow_config"
    return plan


def _execute(plan: Dict[str, Any], abs_file: str) -> Dict[str, Any]:
    from gibh_agent.tools import load_all_tools

    load_all_tools()
    from gibh_agent.core.executor import WorkflowExecutor

    ex = WorkflowExecutor(upload_dir=str(_ROOT))
    return ex.execute_workflow(plan.get("workflow_data") or plan, file_paths=[abs_file])


def _tool_payload(step: Dict[str, Any]) -> Dict[str, Any]:
    """Executor 步骤结果可能在 step_result.data 中。"""
    sr = step.get("step_result") or {}
    data = sr.get("data")
    if isinstance(data, dict) and data:
        return data
    return sr if isinstance(sr, dict) else {}


def _find_report_step(details: List[Dict[str, Any]]) -> Dict[str, Any]:
    for d in details:
        if d.get("tool_id") == "genomics_clinical_reporting":
            return _tool_payload(d)
    return {}


def main() -> int:
    from gibh_agent.tools.omics_genomics_real_io import discover_genomics_reference

    disc = discover_genomics_reference()
    print("参考基因组发现:", json.dumps(disc, indent=2, ensure_ascii=False))
    if not disc.get("bwa_index_ready"):
        print("❌ bwa 索引未就绪，请运行: python3 scripts/init_omics_mock_references.py")
        return 1
    if not FASTQ.is_file():
        print(f"❌ 缺少测试 FASTQ: {FASTQ}")
        return 1

    report = asyncio.run(_run())
    return report


async def _run() -> int:
    plan = await _plan()
    report = _execute(plan, str(FASTQ.resolve()))
    st = report.get("status")
    details = report.get("steps_details") or []

    for d in details:
        inner = d.get("step_result") or {}
        if inner.get("fatal_env_prereq"):
            print(f"❌ {d.get('step_id')}: host_env_prereq — {inner.get('message', '')[:300]}")
            return 1
        msg = str(inner.get("message") or "")
        if "基建贯通" in msg or "infra_pass" in msg:
            print(f"❌ {d.get('step_id')}: 仍含 Mock/基建贯通: {msg[:200]}")
            return 1
        if d.get("status") == "skipped":
            print(f"❌ {d.get('step_id')}: 不允许静默 skipped（下半场须真实 CLI 或 error+日志）")
            return 1
        if d.get("tool_id") in (
            "genomics_cnv_calling",
            "genomics_sv_calling",
            "genomics_variant_annotation",
            "genomics_acmg_classification",
        ):
            payload = _tool_payload(d)
            if not (
                payload.get("cli_command")
                or payload.get("cli_stderr_excerpt")
                or (payload.get("markdown") or "")[:80]
            ):
                print(f"❌ {d.get('step_id')}: 缺少 CLI 日志字段")
                return 1

    germline = next(
        (d for d in details if d.get("tool_id") == "genomics_germline_calling"), None
    )
    if germline:
        gr = _tool_payload(germline)
        vs = gr.get("variant_summary") or {}
        vcf = gr.get("vcf_path") or gr.get("file_path") or ""
        print(
            f"germline: status={germline.get('status')} variants={vs.get('variant_count')} vcf={vcf}"
        )
        if germline.get("status") != "success":
            print(f"❌ germline 未成功: {gr.get('message', '')[:400]}")
            return 1
        if not vcf or not os.path.isfile(str(vcf)) or os.path.getsize(str(vcf)) <= 0:
            print("❌ germline VCF 缺失或为空")
            return 1
        if gr.get("caller") not in ("bcftools", "gatk") and "bcftools" not in str(
            gr.get("message", "")
        ).lower():
            print(f"❌ germline 未使用真实 caller: {gr.get('caller')}")
            return 1

    rep = _find_report_step(details)
    md = rep.get("markdown") or ""
    if "变异记录数" not in md and not (rep.get("variant_summary") or {}):
        print("❌ 临床报告未含真实变异统计")
        print(md[:800])
        return 1

    if st != "success":
        print(f"❌ workflow status={st}")
        for d in details:
            if d.get("status") == "error":
                print(f"  - {d.get('step_id')}: {(d.get('step_result') or {}).get('message', '')[:200]}")
        return 1

    print("\033[32m✅ 基因组真实管线 E2E 通过\033[0m")
    print(f"   steps={len(details)} report excerpt:\n{md[:600]}...")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
