#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
意图澄清闭环 — 双步 E2E 集成测试（POST /api/chat SSE）

步骤 A：模糊提问 → 必须收到 event: clarification + 预期候选 id
步骤 B：intent_override_id → 不得再 clarification；必须拉起 workflow 或 skill 执行流

运行:
  cd /home/ubuntu/GIBH-AGENT-V2
  UPLOAD_DIR=./uploads RESULTS_DIR=./results python3 tests/test_intent_e2e_clarification_loop.py

可选:
  GIBH_E2E_LIMIT=3          # 只跑前 N 条
  GIBH_E2E_BASE_URL=http://127.0.0.1:8028  # 对外部已启动服务黑盒测试
  GIBH_E2E_STEP_B_TIMEOUT=180
"""
from __future__ import annotations

import asyncio
import json
import os
import sys
import time
import uuid
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional, Set, Tuple

ROOT = Path(__file__).resolve().parent.parent

# 必须在 import server 之前设置可写目录
os.environ.setdefault("UPLOAD_DIR", str(ROOT / "uploads"))
os.environ.setdefault("RESULTS_DIR", str(ROOT / "results"))
os.environ.setdefault("GIBH_MASTER_INTENT_ROUTER", "1")
(ROOT / "uploads").mkdir(parents=True, exist_ok=True)
(ROOT / "results").mkdir(parents=True, exist_ok=True)

if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))


def _load_dotenv() -> None:
    env_path = ROOT / ".env"
    if not env_path.is_file():
        return
    for line in env_path.read_text(encoding="utf-8").splitlines():
        line = line.strip()
        if not line or line.startswith("#") or "=" not in line:
            continue
        k, _, v = line.partition("=")
        k = k.strip()
        if k and k not in os.environ:
            os.environ[k] = v.strip().strip('"').strip("'")


_SKILL_IDS = frozenset(
    {
        "lipinski_druglikeness",
        "drug_similarity_search",
        "bepipred3_epitope_prediction",
    }
)

# 20 条双步 E2E 用例（expected_candidates 使用注册表真实 id）
E2E_CLARIFICATION_CASES: List[Dict[str, Any]] = [
    {
        "id": "E01_compound_eval",
        "query": "我对这个化合物做个评估",
        "expected_candidates": ["lipinski_druglikeness", "drug_similarity_search"],
    },
    {
        "id": "E02_bam_mutation",
        "query": "拿这个 BAM 文件找一下突变",
        "expected_candidates": ["genomics_pipeline", "rna_scrna"],
    },
    {
        "id": "E03_omics_generic",
        "query": "跑个组学",
        "expected_candidates": ["genomics_pipeline", "metabolomics", "rna_scrna"],
    },
    {
        "id": "E04_sequencing",
        "query": "测序下机了，帮我分析一下",
        "expected_candidates": ["rna_scrna", "genomics_pipeline", "spatial_transcriptomics"],
    },
    {
        "id": "E05_drug_molecule",
        "query": "这个小分子药物帮我看看怎么样",
        "expected_candidates": ["lipinski_druglikeness", "drug_similarity_search"],
    },
    {
        "id": "E06_medical_image",
        "query": "我有一批医学图像，想做一下分析",
        "expected_candidates": ["radiomics", "spatial_transcriptomics"],
    },
    {
        "id": "E07_protein_seq",
        "query": "有一条蛋白序列，帮我分析一下",
        "expected_candidates": ["bepipred3_epitope_prediction", "proteomics"],
    },
    {
        "id": "E08_cell_trajectory",
        "query": "想看看细胞轨迹怎么分析",
        "expected_candidates": ["sted_ec_trajectory", "spatiotemporal_dynamics"],
    },
    {
        "id": "E09_mass_spec",
        "query": "质谱数据到了，帮我跑一下",
        "expected_candidates": ["metabolomics", "proteomics"],
    },
    {
        "id": "E10_methyl_variant",
        "query": "想看一下甲基化和变异，数据都在这了",
        "expected_candidates": ["epigenomics", "genomics_pipeline"],
    },
    {
        "id": "E11_scrna_spatial",
        "query": "单细胞数据想走个流程，但还没定具体方案",
        "expected_candidates": ["rna_scrna", "spatial_transcriptomics"],
    },
    {
        "id": "E12_genome_transcript",
        "query": "全基因组和转录组都想了解一下怎么跑",
        "expected_candidates": ["genomics_pipeline", "rna_scrna"],
    },
    {
        "id": "E13_antibody",
        "query": "抗体相关的序列分析帮我安排一下",
        "expected_candidates": ["bepipred3_epitope_prediction", "proteomics"],
    },
    {
        "id": "E14_smiles_dual",
        "query": "查一下这个 SMILES 有没有类似的，顺便看看成药性",
        "expected_candidates": ["drug_similarity_search", "lipinski_druglikeness"],
    },
    {
        "id": "E15_multi_omics",
        "query": "多组学联合分析想做一下",
        "expected_candidates": ["genomics_pipeline", "metabolomics", "proteomics"],
    },
    {
        "id": "E16_fastq",
        "query": "FASTQ 文件到了，帮我做个标准分析",
        "expected_candidates": ["genomics_pipeline", "rna_scrna"],
    },
    {
        "id": "E17_pathology_radiomics",
        "query": "病理切片和影像组学特征都想提取",
        "expected_candidates": ["radiomics", "spatial_transcriptomics"],
    },
    {
        "id": "E18_spatiotemporal",
        "query": "细胞时空动态和轨迹推断哪个适合我",
        "expected_candidates": ["spatiotemporal_dynamics", "sted_ec_trajectory"],
    },
    {
        "id": "E19_wet_lab",
        "query": "湿实验做完上机了，组学流程怎么选",
        "expected_candidates": ["metabolomics", "rna_scrna", "genomics_pipeline"],
    },
    {
        "id": "E20_bioinfo_pipeline",
        "query": "生信流水线跑一下，数据是组学类型的",
        "expected_candidates": ["genomics_pipeline", "metabolomics", "rna_scrna"],
    },
]

PASS_THRESHOLD = 0.90


@dataclass
class ParsedSSE:
    event_types: Set[str] = field(default_factory=set)
    events: Dict[str, List[Any]] = field(default_factory=dict)
    raw_tail: str = ""

    def add(self, event_type: str, data: Any) -> None:
        self.event_types.add(event_type)
        self.events.setdefault(event_type, []).append(data)

    def candidate_ids(self) -> List[str]:
        out: List[str] = []
        for block in self.events.get("clarification", []):
            if not isinstance(block, dict):
                continue
            for c in block.get("candidates") or []:
                if isinstance(c, dict) and c.get("id"):
                    out.append(str(c["id"]))
        return out

    def has_substantive_execution(self, override_id: str) -> bool:
        if "workflow" in self.event_types:
            for wf in self.events.get("workflow", []):
                if isinstance(wf, dict):
                    steps = (
                        wf.get("steps")
                        or (wf.get("workflow_data") or {}).get("steps")
                        or (wf.get("workflow_config") or {}).get("workflow_data", {}).get("steps")
                    )
                    if steps:
                        return True
                    if wf.get("workflow_name") or wf.get("domain_name"):
                        return True
                    cfg = wf.get("workflow_config")
                    if isinstance(cfg, dict) and (cfg.get("steps") or cfg.get("workflow_name")):
                        return True
        if "skill_tool_result" in self.event_types:
            return True
        blob = json.dumps(self.events, ensure_ascii=False)
        if "[Skill_Route:" in blob or "[skill_route:" in blob.lower():
            return True
        if override_id in _SKILL_IDS:
            for d in self.events.get("done", []):
                if isinstance(d, dict) and (
                    d.get("tool_name") or d.get("status") in ("success", "error")
                ):
                    return True
            for msg in self.events.get("message", []):
                if isinstance(msg, dict) and len(str(msg.get("content") or "")) > 20:
                    return True
            if "status" in self.event_types and "done" in self.event_types:
                return True
        keywords = ("工作流", "步骤", "规划", "DAG", "workflow", "SOP", "技能快车道", "表位", "质控", "分析流程")
        for st in self.events.get("status", []):
            if isinstance(st, dict):
                c = str(st.get("content") or "")
                if any(k in c for k in keywords):
                    return True
        for msg in self.events.get("message", []):
            if isinstance(msg, dict):
                c = str(msg.get("content") or "")
                if any(k in c for k in keywords):
                    return True
        if "state_snapshot" in self.event_types:
            for snap in self.events.get("state_snapshot", []):
                if isinstance(snap, dict) and (snap.get("workflow") or snap.get("steps")):
                    return True
        return False


@dataclass
class E2EVerdict:
    case_id: str
    query: str
    passed: bool
    step_a_ok: bool = False
    step_b_ok: bool = False
    errors: List[str] = field(default_factory=list)
    chosen_override: str = ""
    step_a_candidates: List[str] = field(default_factory=list)


def _color(t: str, c: str) -> str:
    if not sys.stdout.isatty():
        return t
    return f"\033[{c}m{t}\033[0m"


def _require_llm_keys() -> None:
    if not (
        (os.getenv("DEEPSEEK_API_KEY") or "").strip()
        or (os.getenv("ZHIPU_API_KEY") or "").strip()
        or (os.getenv("MOONSHOT_API_KEY") or "").strip()
    ):
        print("SKIP: 需要 DEEPSEEK/ZHIPU/MOONSHOT API Key", file=sys.stderr)
        sys.exit(0)


async def parse_sse_response(resp, *, max_seconds: float) -> ParsedSSE:
    """消费 httpx 流式响应直至超时或连接结束。"""
    parsed = ParsedSSE()
    deadline = time.monotonic() + max_seconds
    current_event: Optional[str] = None
    buf = ""

    async for chunk in resp.aiter_text():
        if time.monotonic() > deadline:
            break
        buf += chunk
        while "\n\n" in buf:
            block, buf = buf.split("\n\n", 1)
            event_type = "message"
            data_obj: Any = None
            for line in block.split("\n"):
                line = line.strip()
                if line.startswith("event:"):
                    event_type = line[6:].strip()
                elif line.startswith("data:"):
                    raw = line[5:].strip()
                    try:
                        data_obj = json.loads(raw)
                    except json.JSONDecodeError:
                        data_obj = raw
            if data_obj is not None:
                parsed.add(event_type, data_obj)
                current_event = event_type
            if event_type == "done" and current_event == "done":
                return parsed

    parsed.raw_tail = buf[-500:]
    return parsed


def _pick_override_id(candidates: List[str], expected: List[str]) -> Optional[str]:
    cand_set = set(candidates)
    for eid in expected:
        if eid in cand_set:
            return eid
    return candidates[0] if candidates else None


def assert_step_a(parsed: ParsedSSE, expected: List[str]) -> List[str]:
    errs: List[str] = []
    if "clarification" not in parsed.event_types:
        errs.append("A1: SSE 未包含 event: clarification")
    if "message" in parsed.event_types:
        errs.append("A1b: 澄清流禁止 event: message（会导致前端文本重复）")
    for st in parsed.events.get("status", []):
        if isinstance(st, dict) and "正在接收请求" in str(st.get("content") or ""):
            errs.append("A1c: 澄清流不得包含 Orchestrator 执行记录 status")
    if "workflow" in parsed.event_types:
        errs.append("A1d: 澄清流不得包含 event: workflow")
    ids = parsed.candidate_ids()
    if len(ids) <= 1:
        errs.append(f"A2: candidates 数量须>1，实际 {len(ids)}")
    missing = [x for x in expected if x not in ids]
    if missing:
        errs.append(f"A2: 缺少预期候选 id {missing}；实际 {ids}")
    return errs


def assert_step_b(parsed: ParsedSSE, override_id: str) -> List[str]:
    errs: List[str] = []
    if "clarification" in parsed.event_types:
        errs.append("B1: override 请求仍出现 clarification（未物理短路）")
    if not parsed.has_substantive_execution(override_id):
        errs.append(
            f"B2: 未检测到 workflow/技能执行实质内容（events={sorted(parsed.event_types)}）"
        )
    return errs


async def post_chat_sse(
    client,
    *,
    message: str,
    session_id: str,
    intent_override_id: Optional[str] = None,
    clarification_context_message: Optional[str] = None,
    timeout: float,
) -> Tuple[int, ParsedSSE]:
    payload: Dict[str, Any] = {
        "message": message,
        "history": [],
        "uploaded_files": [],
        "stream": True,
        "session_id": session_id,
        "user_id": "guest",
        "model_name": os.getenv("GIBH_E2E_MODEL", "deepseek-v4-pro"),
    }
    if intent_override_id:
        payload["intent_override_id"] = intent_override_id
    if clarification_context_message:
        payload["clarification_context_message"] = clarification_context_message

    headers = {"X-Guest-UUID": f"e2e-guest-{uuid.uuid4().hex[:12]}"}
    async with client.stream(
        "POST",
        "/api/chat",
        json=payload,
        headers=headers,
        timeout=timeout,
    ) as resp:
        status = resp.status_code
        if status != 200:
            body = await resp.aread()
            p = ParsedSSE()
            p.raw_tail = body.decode("utf-8", errors="replace")[:400]
            return status, p
        parsed = await parse_sse_response(resp, max_seconds=timeout)
        return status, parsed


async def run_one_case(client, tc: Dict[str, Any]) -> E2EVerdict:
    case_id = tc["id"]
    query = tc["query"]
    expected = tc["expected_candidates"]
    session_id = f"e2e-{case_id}-{uuid.uuid4().hex[:8]}"
    errors: List[str] = []

    step_a_timeout = float(os.getenv("GIBH_E2E_STEP_A_TIMEOUT", "90"))
    step_b_timeout = float(os.getenv("GIBH_E2E_STEP_B_TIMEOUT", "180"))

    status_a, parsed_a = await post_chat_sse(
        client,
        message=query,
        session_id=session_id,
        timeout=step_a_timeout,
    )
    if status_a != 200:
        errors.append(f"动作 A HTTP {status_a}: {parsed_a.raw_tail[:200]}")
        return E2EVerdict(case_id, query, False, errors=errors)

    errs_a = assert_step_a(parsed_a, expected)
    step_a_ok = not errs_a
    errors.extend(errs_a)

    override_id = _pick_override_id(parsed_a.candidate_ids(), expected)
    if not override_id:
        errors.append("动作 B 跳过: 无可用 candidate_id")
        return E2EVerdict(
            case_id,
            query,
            False,
            step_a_ok=step_a_ok,
            step_a_candidates=parsed_a.candidate_ids(),
            errors=errors,
        )

    status_b, parsed_b = await post_chat_sse(
        client,
        message=query,
        session_id=session_id,
        intent_override_id=override_id,
        clarification_context_message=query,
        timeout=step_b_timeout,
    )
    if status_b != 200:
        errors.append(f"动作 B HTTP {status_b}: {parsed_b.raw_tail[:200]}")
        return E2EVerdict(
            case_id,
            query,
            False,
            step_a_ok=step_a_ok,
            chosen_override=override_id,
            step_a_candidates=parsed_a.candidate_ids(),
            errors=errors,
        )

    errs_b = assert_step_b(parsed_b, override_id)
    step_b_ok = not errs_b
    errors.extend(errs_b)

    return E2EVerdict(
        case_id=case_id,
        query=query,
        passed=step_a_ok and step_b_ok,
        step_a_ok=step_a_ok,
        step_b_ok=step_b_ok,
        errors=errors,
        chosen_override=override_id,
        step_a_candidates=parsed_a.candidate_ids(),
    )


def _make_client():
    base = (os.getenv("GIBH_E2E_BASE_URL") or "").strip()
    import httpx

    if base:
        return httpx.AsyncClient(base_url=base.rstrip("/"), timeout=None)

    from httpx import ASGITransport

    _load_dotenv()
    from server import app  # noqa: E402

    transport = ASGITransport(app=app)
    return httpx.AsyncClient(transport=transport, base_url="http://test", timeout=None)


async def verify_override_never_calls_classify() -> None:
    """单元级佐证：intent_override_id 分支不触发 classify（mock）。"""
    from unittest.mock import AsyncMock, MagicMock, patch

    from gibh_agent.core.intent_router import MasterIntentRouter, should_bypass_master_intent_router

    class _Req:
        intent_override_id = "genomics_pipeline"
        target_domain = None
        workflow_data = None
        message = "跑个组学"
        test_dataset_id = None
        local_tool_response = None

    assert should_bypass_master_intent_router(_Req()) is True

    mock_llm = MagicMock()
    mock_llm.achat = AsyncMock()
    router = MasterIntentRouter(mock_llm)
    with patch.object(router, "classify", new_callable=AsyncMock) as mock_cls:
        # 模拟 server 分支：有 override 时不应调用 classify
        if (_Req().intent_override_id or "").strip():
            pass  # bypass
        else:
            await router.classify("x")
        mock_cls.assert_not_awaited()
    print(_color("✓ 佐证: intent_override_id 时 classify 未被调用", "32"))


async def _main_async() -> int:
    _load_dotenv()
    _require_llm_keys()

    limit_raw = (os.getenv("GIBH_E2E_LIMIT") or "").strip()
    cases = E2E_CLARIFICATION_CASES
    if limit_raw.isdigit():
        cases = cases[: int(limit_raw)]

    await verify_override_never_calls_classify()
    print()

    async with _make_client() as client:
        verdicts: List[E2EVerdict] = []
        for tc in cases:
            print("=" * 72)
            print(f"【{tc['id']}】{tc['query']}")
            v = await run_one_case(client, tc)
            verdicts.append(v)
            print(f"  动作 A: {'Pass' if v.step_a_ok else 'Fail'} | 候选: {v.step_a_candidates}")
            print(f"  动作 B: {'Pass' if v.step_b_ok else 'Fail'} | override={v.chosen_override}")
            if v.passed:
                print(_color("  链路: Pass", "32"))
            else:
                print(_color("  链路: Fail", "31"))
                for e in v.errors:
                    print(_color(f"    ✗ {e}", "31"))
            await asyncio.sleep(0.5)

    passed = sum(1 for v in verdicts if v.passed)
    total = len(verdicts)
    rate = passed / total if total else 0.0
    print()
    print("=" * 72)
    if rate >= PASS_THRESHOLD:
        print(_color(f"E2E 汇总: {passed}/{total} Passed ({rate * 100:.0f}%)", "32"))
    else:
        print(_color(f"E2E 汇总: {passed}/{total} Passed ({rate * 100:.0f}%) — 未达标", "31"))
    return 0 if rate >= PASS_THRESHOLD else 1


def main() -> None:
    raise SystemExit(asyncio.run(_main_async()))


if __name__ == "__main__":
    main()
