#!/usr/bin/env python3
"""
三大组学全 DAG 端到端（模拟前端 /api/execute 路径）。

前置：
  1. python3 scripts/init_omics_mock_references.py
  2. API 已启动（默认 http://127.0.0.1:8028）或设置 GIBH_E2E_EXECUTOR_LOCAL=1 仅跑 Executor

运行:
  cd /home/ubuntu/GIBH-AGENT-V2 && PYTHONPATH=. python3 tests/test_omics_full_modalities_e2e.py
"""
from __future__ import annotations

import asyncio
import json
import os
import sys
import urllib.error
import urllib.request
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

_ROOT = Path(__file__).resolve().parent.parent
if str(_ROOT) not in sys.path:
    sys.path.insert(0, str(_ROOT))

os.environ.setdefault("RESULTS_DIR", str(_ROOT / "results" / "e2e_omics_full"))
os.environ.setdefault(
    "SECURITY_KEYS_PATH",
    str(_ROOT / "results" / "e2e_omics_full" / "security" / "signing_keys.json"),
)
os.environ.setdefault("OMICS_REF_DIR", str(_ROOT / "data" / "references"))
os.environ.setdefault(
    "GIBH_REF_HG38", str(_ROOT / "data" / "references" / "genomics" / "hg38.fa")
)
os.environ.setdefault(
    "GIBH_BOWTIE2_HG38_INDEX",
    str(_ROOT / "data" / "references" / "epigenomics" / "bowtie2" / "hg38" / "hg38"),
)
os.environ.setdefault(
    "GIBH_PROTEOMICS_FASTA",
    str(_ROOT / "data" / "references" / "proteomics" / "uniprot_mini.fasta"),
)

API_BASE = os.environ.get("GIBH_API_BASE", "http://127.0.0.1:8028").rstrip("/")
USE_LOCAL_EXECUTOR = os.environ.get("GIBH_E2E_EXECUTOR_LOCAL", "").strip() in (
    "1",
    "true",
    "yes",
)

CASES: List[Dict[str, str]] = [
    {
        "domain": "genomics",
        "rel_path": "test_data/genomics/sample1_R1.fastq.gz",
        "prompt": "请帮我用这个 FASTQ 跑一下基因组胚系分析全流程",
    },
    {
        "domain": "proteomics",
        "rel_path": "test_data/proteomics/BSA1_F1.mzML",
        "prompt": "请用这份 mzML 跑蛋白质组学搜库定量全流程",
    },
    {
        "domain": "epigenomics",
        "rel_path": "test_data/epigenomics/SRR1822153_1.fastq.gz",
        "prompt": "请用这个 FASTQ 跑 ATAC-seq/ChIP 主线表观组学全流程",
    },
]


class _DummyLLM:
    async def astream(self, *a: Any, **kw: Any):
        if False:
            yield {}

    async def achat(self, *a: Any, **kw: Any) -> Any:
        raise RuntimeError("E2E 不应调用 achat")

    async def astream_chat(self, *a: Any, **kw: Any):
        if False:
            yield None


def _http_post_json(url: str, payload: Dict[str, Any], timeout: int = 7200) -> Dict[str, Any]:
    data = json.dumps(payload).encode("utf-8")
    req = urllib.request.Request(
        url,
        data=data,
        headers={"Content-Type": "application/json"},
        method="POST",
    )
    with urllib.request.urlopen(req, timeout=timeout) as resp:  # noqa: S310
        body = resp.read().decode("utf-8", errors="replace")
    return json.loads(body)


def _assert_no_env_block(details: List[Dict[str, Any]], domain: str) -> None:
    for d in details:
        inner = d.get("step_result") or {}
        if inner.get("fatal_env_prereq") or inner.get("error_category") == "host_env_prereq":
            raise AssertionError(
                f"[{domain}] 步骤 {d.get('step_id')} 触发环境阻断: {inner.get('message')}"
            )
        data = inner.get("data") if isinstance(inner.get("data"), dict) else {}
        if data.get("fatal_env_prereq"):
            raise AssertionError(f"[{domain}] 步骤 {d.get('step_id')} data 含 fatal_env_prereq")


async def _plan_workflow(domain: str, abs_file: str, prompt: str) -> Dict[str, Any]:
    from gibh_agent.tools import load_all_tools

    load_all_tools()
    from gibh_agent.core.file_inspector import FileInspector
    from gibh_agent.core.planner import SOPPlanner
    from gibh_agent.core.tool_retriever import ToolRetriever
    from gibh_agent.core.workflows.registry import WorkflowRegistry

    inspector = FileInspector(upload_dir=str(_ROOT))
    meta = inspector.inspect_file(abs_file)
    assert meta.get("status") == "success", meta

    wf = WorkflowRegistry().get_workflow(domain)
    assert wf is not None
    keys = list(wf.get_steps_dag().keys())

    planner = SOPPlanner(ToolRetriever(), _DummyLLM())  # type: ignore[arg-type]
    plan: Optional[Dict[str, Any]] = None
    async for ev, data in planner.generate_plan(
        user_query=prompt,
        file_metadata=meta,
        domain_name=domain,
        target_steps=keys,
        is_template=False,
    ):
        if ev == "workflow":
            plan = data
    assert isinstance(plan, dict) and plan.get("type") == "workflow_config"
    return plan


def _api_file_path(host_abs: str) -> str:
    """容器内 API 使用 /app 挂载路径；裸机直连则保持原路径。"""
    try:
        rel = Path(host_abs).resolve().relative_to(_ROOT.resolve())
        return f"/app/{rel.as_posix()}"
    except ValueError:
        return host_abs


async def _execute_via_api(wd: Dict[str, Any], abs_file: str) -> Dict[str, Any]:
    api_fp = _api_file_path(abs_file)
    payload = {
        "workflow_data": wd.get("workflow_data") or wd,
        "file_paths": [api_fp],
    }
    url = f"{API_BASE}/api/execute"
    return _http_post_json(url, payload)


def _execute_local(wd: Dict[str, Any], abs_file: str) -> Dict[str, Any]:
    from gibh_agent.tools import load_all_tools

    load_all_tools()
    from gibh_agent.core.executor import WorkflowExecutor

    ex = WorkflowExecutor(upload_dir=str(_ROOT))
    return ex.execute_workflow(wd.get("workflow_data") or wd, file_paths=[abs_file])


def _normalize_execute_report(report: Dict[str, Any]) -> Dict[str, Any]:
    """兼容 /api/execute 返回的 analysis_report 包装。"""
    if report.get("type") == "analysis_report" and isinstance(report.get("report_data"), dict):
        inner = report["report_data"]
        inner.setdefault("status", report.get("status"))
        return inner
    return report


async def _run_case(case: Dict[str, str]) -> Tuple[str, Dict[str, Any]]:
    domain = case["domain"]
    fp = (_ROOT / case["rel_path"]).resolve()
    assert fp.is_file(), f"缺少测试文件: {fp}"

    plan = await _plan_workflow(domain, str(fp), case["prompt"])
    wd = plan

    if USE_LOCAL_EXECUTOR:
        report = _execute_local(wd, str(fp))
    else:
        try:
            report = _normalize_execute_report(await _execute_via_api(wd, str(fp)))
            via = "api"
        except (urllib.error.URLError, urllib.error.HTTPError, TimeoutError) as exc:
            print(f"⚠️ API 不可达 ({API_BASE})，回退本地 Executor: {exc}")
            report = _execute_local(wd, str(fp))
            via = "local"
        else:
            print(f"   📡 经 API 执行完成 ({API_BASE})")

    return domain, report


async def main_async() -> int:
    manifest = _ROOT / "data" / "references" / "manifest.json"
    if not manifest.is_file():
        print("❌ 请先运行: python3 scripts/init_omics_mock_references.py")
        return 1

    ok = 0
    for i, case in enumerate(CASES, 1):
        print("=" * 80)
        print(f"[{i}/3] {case['domain']} | file={case['rel_path']}")
        print("=" * 80)
        domain, report = await _run_case(case)
        st = report.get("status")
        details = report.get("steps_details") or []
        n_ok = sum(1 for d in details if d.get("status") == "success")
        n_err = sum(1 for d in details if d.get("status") == "error")

        if st != "success":
            print(f"❌ workflow_status={st!r} success_steps={n_ok} error_steps={n_err}")
            cb = report.get("circuit_breaker")
            if cb:
                print("circuit_breaker:", cb)
            for d in details:
                if d.get("status") == "error":
                    print(
                        f"  - {d.get('step_id')}: {(d.get('step_result') or {}).get('message', '')[:200]}"
                    )
            return 1

        _assert_no_env_block(details, domain)
        ok += 1
        print(
            f"\033[32m✅ [{domain}] 全流程 success | steps={len(details)} "
            f"(ok={n_ok}, err={n_err}) | 无 host_env_prereq\033[0m"
        )

    print()
    print("\033[32m" + "=" * 80)
    print(f"🎉 三大模态全 DAG E2E 通过 ({ok}/3)")
    print("=" * 80 + "\033[0m")
    return 0


def main() -> None:
    raise SystemExit(asyncio.run(main_async()))


if __name__ == "__main__":
    main()
