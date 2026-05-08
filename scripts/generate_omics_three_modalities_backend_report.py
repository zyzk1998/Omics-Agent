#!/usr/bin/env python3
"""
三大模态（基因组 / 蛋白组 / 表观）工作流在后端单机执行，逐步汇总结果并写出 Markdown 报告。

用法（仓库根目录）:
  PYTHONPATH=. python3 scripts/generate_omics_three_modalities_backend_report.py

输出:
  OMICS_THREE_MODALITIES_BACKEND_RUN_REPORT.md（仓库根目录）
"""
from __future__ import annotations

import asyncio
import json
import os
import shutil
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List, Optional

_ROOT = Path(__file__).resolve().parent.parent
if str(_ROOT) not in sys.path:
    sys.path.insert(0, str(_ROOT))

os.environ.setdefault("RESULTS_DIR", str(_ROOT / "results" / "omics_backend_report"))
os.environ.setdefault(
    "SECURITY_KEYS_PATH",
    str(_ROOT / "results" / "omics_backend_report" / "security" / "signing_keys.json"),
)

CASES: List[Dict[str, Any]] = [
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

    async def achat(self, *a: Any, **kw: Any):
        raise RuntimeError("不应调用 achat")

    async def astream_chat(self, *a: Any, **kw: Any):
        if False:
            yield None


def _sniff_cli() -> Dict[str, Optional[str]]:
    names = [
        "bwa",
        "samtools",
        "fastp",
        "bowtie2",
        "macs2",
        "gatk",
        "bcftools",
        "diann",
        "msconvert",
        "percolator",
    ]
    return {n: shutil.which(n) for n in names}


def _truncate(s: Any, n: int = 420) -> str:
    if s is None:
        return ""
    t = str(s).strip().replace("\r\n", "\n")
    if len(t) <= n:
        return t
    return t[: n - 3] + "..."


def _summarize_data(data: Any) -> Dict[str, Any]:
    if not isinstance(data, dict):
        return {"raw_type": type(data).__name__}
    out: Dict[str, Any] = {
        "status": data.get("status"),
        "message_excerpt": _truncate(data.get("message"), 280),
    }
    if data.get("qc_metrics"):
        qm = data["qc_metrics"]
        if isinstance(qm, dict):
            keys = list(qm.keys())[:12]
            out["qc_metrics_keys"] = keys
            if "n_reads" in qm:
                out["qc_n_reads"] = qm.get("n_reads")
            if "gc_percent" in qm:
                out["qc_gc_percent"] = qm.get("gc_percent")
            if "spectrum_count" in qm:
                out["qc_spectrum_count"] = qm.get("spectrum_count")
    md = data.get("markdown")
    out["markdown_chars"] = len(md) if isinstance(md, str) else 0
    out["markdown_excerpt"] = _truncate(md, 500)
    imgs = data.get("image_urls")
    out["image_urls_count"] = len(imgs) if isinstance(imgs, list) else 0
    td = data.get("table_data")
    if isinstance(td, dict):
        out["table_data_keys"] = list(td.keys())[:20]
        cols = td.get("columns")
        rows = td.get("rows")
        if isinstance(cols, list):
            out["table_columns"] = cols[:16]
        if isinstance(rows, list):
            out["table_row_count"] = len(rows)
    return out


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
    if file_metadata.get("status") != "success":
        return {"status": "error", "error": "file_inspector", "detail": file_metadata}

    registry = WorkflowRegistry()
    wf = registry.get_workflow(domain)
    if wf is None:
        return {"status": "error", "error": "unknown_domain", "domain": domain}

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

    if not isinstance(plan_result, dict) or plan_result.get("type") != "workflow_config":
        return {"status": "error", "error": "planner", "detail": plan_result}

    wd = plan_result.get("workflow_data") or {}
    if not wd.get("steps"):
        return {"status": "error", "error": "empty_workflow", "detail": wd}

    ex = WorkflowExecutor(upload_dir=str(_ROOT))
    return ex.execute_workflow(wd, file_paths=[abs_file])


def _render_markdown(
    reports: List[Dict[str, Any]],
    sniff: Dict[str, Optional[str]],
    *,
    ok_all: bool,
) -> str:
    ts = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M:%S UTC")
    lines: List[str] = []
    lines.append("# 三大模态底层工具链 · 后端单机逐步执行报告")
    lines.append("")
    lines.append(f"- **生成时间**: {ts}")
    lines.append(f"- **仓库根**: `{_ROOT}`")
    lines.append(
        f"- **总结果**: {'✅ 三组域 workflow 均 `success`' if ok_all else '⚠️ 存在失败域，见下文'}"
    )
    lines.append("")
    lines.append("## 宿主 PATH 嗅探（可选重型 CLI）")
    lines.append("")
    lines.append("| 可执行文件 | 路径（未安装则为空） |")
    lines.append("|------------|------------------------|")
    for name in sorted(sniff.keys()):
        p = sniff[name]
        lines.append(f"| `{name}` | `{p or '—'}` |")
    lines.append("")
    lines.append(
        "> 说明：未安装重型 CLI 时，下游步骤使用基于 FASTQ/mzML **实测统计** 的轻量代理输出（确定性可复现）；"
        "检测到对应二进制与参考资源时仍可走 subprocess 真管线。"
    )
    lines.append("")

    for block in reports:
        domain = block["domain"]
        rel = block["rel_path"]
        lines.append(f"## 模态：`{domain}`")
        lines.append("")
        lines.append(f"- **输入文件**: `{rel}`")
        lines.append(f"- **绝对路径**: `{block.get('abs_path', '')}`")
        wf_status = block.get("workflow_status")
        lines.append(f"- **工作流总状态**: `{wf_status}`")
        if block.get("error"):
            lines.append(f"- **错误**: `{block.get('error')}`")
            lines.append("")
            lines.append("```json")
            lines.append(json.dumps(block.get("detail"), ensure_ascii=False, indent=2)[:8000])
            lines.append("```")
            lines.append("")
            continue

        lines.append(f"- **步骤数**: {block.get('step_count', 0)}")
        lines.append("")
        lines.append("| # | step_id | tool_id | 状态 | 耗时(s) | summary |")
        lines.append("|---|---------|---------|------|---------|---------|")
        for i, row in enumerate(block.get("step_rows", []), 1):
            lines.append(
                f"| {i} | `{row['step_id']}` | `{row['tool_id']}` | {row['status']} | "
                f"{row['duration']} | {_truncate(row['summary'], 120)} |"
            )
        lines.append("")
        lines.append("<details>")
        lines.append(f"<summary>逐步详情（{domain}）</summary>")
        lines.append("")
        for row in block.get("step_rows", []):
            lines.append(f"### `{row['tool_id']}` · {row['step_id']}")
            lines.append("")
            summ = _summarize_data(row.get("tool_out"))
            lines.append("```json")
            lines.append(json.dumps(summ, ensure_ascii=False, indent=2))
            lines.append("```")
            if summ.get("markdown_excerpt"):
                lines.append("")
                lines.append("**Markdown 摘录：**")
                lines.append("")
                lines.append("```markdown")
                lines.append(summ["markdown_excerpt"])
                lines.append("```")
            lines.append("")
        lines.append("</details>")
        lines.append("")

    lines.append("---")
    lines.append("")
    lines.append(
        "*本报告由 `scripts/generate_omics_three_modalities_backend_report.py` 自动生成。*"
    )
    return "\n".join(lines)


async def main_async() -> int:
    sniff = _sniff_cli()
    reports: List[Dict[str, Any]] = []
    ok_all = True

    for case in CASES:
        fp = (_ROOT / case["rel_path"]).resolve()
        domain = case["domain"]
        rec: Dict[str, Any] = {
            "domain": domain,
            "rel_path": case["rel_path"],
            "abs_path": str(fp),
        }
        if not fp.is_file():
            rec["workflow_status"] = "error"
            rec["error"] = "missing_file"
            rec["detail"] = str(fp)
            ok_all = False
            reports.append(rec)
            continue

        report = await _run_domain(domain, str(fp), case["prompt"])
        st = report.get("status")
        rec["workflow_status"] = st
        if st != "success":
            ok_all = False
            rec["error"] = "workflow_failed"
            rec["detail"] = report
            reports.append(rec)
            continue

        details = report.get("steps_details") or []
        rec["step_count"] = len(details)
        rows: List[Dict[str, Any]] = []
        for d in details:
            sr = d.get("step_result") or {}
            data = sr.get("data") if isinstance(sr, dict) else {}
            if not isinstance(data, dict):
                data = {}
            rows.append(
                {
                    "step_id": d.get("step_id", ""),
                    "tool_id": d.get("tool_id", ""),
                    "status": d.get("status", ""),
                    "duration": d.get("duration", 0),
                    "summary": d.get("summary", ""),
                    "tool_out": data,
                }
            )
        rec["step_rows"] = rows
        reports.append(rec)

    md = _render_markdown(reports, sniff, ok_all=ok_all)
    out_path = _ROOT / "OMICS_THREE_MODALITIES_BACKEND_RUN_REPORT.md"
    out_path.write_text(md, encoding="utf-8")
    print(f"Wrote {out_path}")
    return 0 if ok_all else 2


def main() -> None:
    raise SystemExit(asyncio.run(main_async()))


if __name__ == "__main__":
    main()
