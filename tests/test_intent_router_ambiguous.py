#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
MasterIntentRouter 模糊意图（AMBIGUOUS）实弹评测 — 20 条边缘用例。

运行:
  python tests/test_intent_router_ambiguous.py
  pytest -s tests/test_intent_router_ambiguous.py

环境:
  需 DEEPSEEK_API_KEY / ZHIPU_API_KEY / MOONSHOT_API_KEY 之一；会真实调用 LLM。
  可选 GIBH_INTENT_ROUTER_TEST_MODEL=deepseek-v4-pro
"""
from __future__ import annotations

import asyncio
import os
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Set

ROOT = Path(__file__).resolve().parent.parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from gibh_agent.core.intent_router import (  # noqa: E402
    MasterIntentRouter,
    build_global_intent_registry,
)
from gibh_agent.core.workflows import WorkflowRegistry  # noqa: E402
from gibh_agent.core.tool_registry import registry as tool_registry  # noqa: E402

# ---------------------------------------------------------------------------
# 20 条硬核模糊用例：required_ids 中的每一项都必须出现在 candidates 里
# ---------------------------------------------------------------------------
AMBIGUOUS_TEST_CASES: List[Dict[str, Any]] = [
    {
        "id": "TC01_compound_eval",
        "message": "我对这个化合物做个评估",
        "required_ids": ["lipinski_druglikeness", "drug_similarity_search"],
        "has_files": False,
    },
    {
        "id": "TC02_bam_mutation",
        "message": "拿这个 BAM 文件找一下突变",
        "required_ids": ["genomics_pipeline", "rna_scrna"],
        "has_files": True,
    },
    {
        "id": "TC03_omics_generic",
        "message": "跑个组学",
        "required_ids": ["genomics_pipeline", "metabolomics", "rna_scrna"],
        "has_files": False,
    },
    {
        "id": "TC04_sequencing_generic",
        "message": "测序下机了，帮我分析一下",
        "required_ids": ["rna_scrna", "genomics_pipeline", "spatial_transcriptomics"],
        "has_files": True,
    },
    {
        "id": "TC05_drug_molecule",
        "message": "这个小分子药物帮我看看怎么样",
        "required_ids": ["lipinski_druglikeness", "drug_similarity_search"],
        "has_files": False,
    },
    {
        "id": "TC06_medical_image",
        "message": "我有一批医学图像，想做一下分析",
        "required_ids": ["radiomics", "spatial_transcriptomics"],
        "has_files": True,
    },
    {
        "id": "TC07_protein_sequence",
        "message": "有一条蛋白序列，帮我分析一下",
        "required_ids": ["bepipred3_epitope_prediction", "proteomics"],
        "has_files": False,
    },
    {
        "id": "TC08_cell_trajectory",
        "message": "想看看细胞轨迹怎么分析",
        "required_ids": ["sted_ec_trajectory", "spatiotemporal_dynamics"],
        "has_files": False,
    },
    {
        "id": "TC09_mass_spec",
        "message": "质谱数据到了，帮我跑一下",
        "required_ids": ["metabolomics", "proteomics"],
        "has_files": True,
    },
    {
        "id": "TC10_epigenetic_or_genomic",
        "message": "想看一下甲基化和变异，数据都在这了",
        "required_ids": ["epigenomics", "genomics_pipeline"],
        "has_files": True,
    },
    {
        "id": "TC11_scrna_or_spatial",
        "message": "单细胞数据想走个流程，但还没定具体方案",
        "required_ids": ["rna_scrna", "spatial_transcriptomics"],
        "has_files": True,
    },
    {
        "id": "TC12_genome_or_transcriptome",
        "message": "全基因组和转录组都想了解一下怎么跑",
        "required_ids": ["genomics_pipeline", "rna_scrna"],
        "has_files": False,
    },
    {
        "id": "TC13_antibody_or_epitope",
        "message": "抗体相关的序列分析帮我安排一下",
        "required_ids": ["bepipred3_epitope_prediction", "proteomics"],
        "has_files": False,
    },
    {
        "id": "TC14_similarity_or_druglike",
        "message": "查一下这个 SMILES 在数据库里有没有类似的，顺便看看成药性",
        "required_ids": ["drug_similarity_search", "lipinski_druglikeness"],
        "has_files": False,
    },
    {
        "id": "TC15_multi_omics_vague",
        "message": "多组学联合分析想做一下",
        "required_ids": ["genomics_pipeline", "metabolomics", "proteomics"],
        "has_files": False,
    },
    {
        "id": "TC16_fastq_unknown",
        "message": "FASTQ 文件到了，帮我做个标准分析",
        "required_ids": ["rna_scrna", "genomics_pipeline"],
        "has_files": True,
    },
    {
        "id": "TC17_pathology_or_radiomics",
        "message": "病理切片和影像组学特征都想提取",
        "required_ids": ["radiomics", "spatial_transcriptomics"],
        "has_files": True,
    },
    {
        "id": "TC18_temporal_cell",
        "message": "细胞时空动态和轨迹推断哪个适合我",
        "required_ids": ["spatiotemporal_dynamics", "sted_ec_trajectory"],
        "has_files": False,
    },
    {
        "id": "TC19_wet_lab_omics",
        "message": "湿实验做完上机了，组学流程怎么选",
        "required_ids": ["metabolomics", "rna_scrna", "genomics_pipeline"],
        "has_files": True,
    },
    {
        "id": "TC20_bioinfo_pipeline",
        "message": "生信流水线跑一下，数据是组学类型的",
        "required_ids": ["genomics_pipeline", "metabolomics", "rna_scrna"],
        "has_files": True,
    },
]

PASS_THRESHOLD = 0.90


@dataclass
class CaseVerdict:
    case_id: str
    message: str
    passed: bool
    errors: List[str] = field(default_factory=list)
    status: str = ""
    question: str = ""
    candidate_ids: List[str] = field(default_factory=list)
    rationale: str = ""


def _color(text: str, code: str) -> str:
    if not sys.stdout.isatty():
        return text
    return f"\033[{code}m{text}\033[0m"


def _green(t: str) -> str:
    return _color(t, "32")


def _red(t: str) -> str:
    return _color(t, "31")


def _cyan(t: str) -> str:
    return _color(t, "36")


def _yellow(t: str) -> str:
    return _color(t, "33")


def _require_llm_env() -> None:
    if not (
        (os.getenv("DEEPSEEK_API_KEY") or "").strip()
        or (os.getenv("ZHIPU_API_KEY") or "").strip()
        or (os.getenv("MOONSHOT_API_KEY") or "").strip()
    ):
        print(_yellow("SKIP: 未设置 DEEPSEEK/ZHIPU/MOONSHOT API Key"), file=sys.stderr)
        sys.exit(0)


def _load_dotenv_if_present() -> None:
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


def _question_is_friendly_interrogative(q: str) -> bool:
    q = (q or "").strip()
    if len(q) < 6:
        return False
    if "？" in q or "?" in q:
        return True
    markers = ("请问", "是否", "哪一种", "哪种", "哪类", "需要", "想", "还是", "或者", "请选择", "确认")
    return any(m in q for m in markers)


def evaluate_ambiguous_case(
    result: Any,
    *,
    required_ids: Sequence[str],
) -> List[str]:
    """三重严格断言，返回错误列表（空=通过）。"""
    errors: List[str] = []

    status = (getattr(result, "status", None) or "").strip().lower()
    if status != "clarify":
        errors.append(f'Status 断言失败: 期望 "clarify"，实际 "{status}"')

    question = (getattr(result, "question", None) or "").strip()
    if not question:
        errors.append("Question 断言失败: question 为空")
    elif not _question_is_friendly_interrogative(question):
        errors.append(f"Question 断言失败: 不像友好反问句 → {question[:80]}")

    candidates = getattr(result, "candidates", None) or []
    cand_ids: List[str] = [
        str(c.get("id", "")).strip()
        for c in candidates
        if isinstance(c, dict) and c.get("id")
    ]
    if len(cand_ids) <= 1:
        errors.append(f"Candidates 断言失败: 候选数量 {len(cand_ids)} 须 > 1")

    required_set: Set[str] = {str(x).strip() for x in required_ids if x}
    cand_set = set(cand_ids)
    missing = required_set - cand_set
    if missing:
        errors.append(
            f"Candidates 断言失败: 缺少必需 id {sorted(missing)}；"
            f"模型给出 {cand_ids}"
        )

    if status == "exact":
        errors.append("禁止 EXACT_MATCH: 模糊用例不得判 exact")
    if status == "chitchat":
        errors.append("禁止 CHITCHAT: 模糊用例有业务执行诉求，不得判 chitchat")

    return errors


def _print_registry_banner() -> None:
    wf = WorkflowRegistry()
    print(_cyan("── WorkflowRegistry ──"))
    for name, desc in (wf.list_workflows() or {}).items():
        print(f"  • {name}: {(desc or '')[:60]}")
    print(_cyan(f"── ToolRegistry 已注册工具数: {len(tool_registry._tools)} ──"))
    reg = build_global_intent_registry()
    print(_cyan(f"── MasterIntentRouter 注册表条目: {len(reg)} ──"))
    for e in reg:
        print(f"  • {e.id} ({e.route_kind.value})")


async def run_ambiguous_suite(router: MasterIntentRouter) -> List[CaseVerdict]:
    verdicts: List[CaseVerdict] = []
    for tc in AMBIGUOUS_TEST_CASES:
        case_id = tc["id"]
        message = tc["message"]
        required_ids = tc["required_ids"]
        has_files = bool(tc.get("has_files"))

        result = await router.classify(message, has_files=has_files)
        errors = evaluate_ambiguous_case(result, required_ids=required_ids)
        cand_ids = [
            str(c.get("id", ""))
            for c in (result.candidates or [])
            if isinstance(c, dict)
        ]
        verdict = CaseVerdict(
            case_id=case_id,
            message=message,
            passed=not errors,
            errors=errors,
            status=result.status,
            question=result.question or "",
            candidate_ids=cand_ids,
            rationale=result.rationale_short or "",
        )
        verdicts.append(verdict)

        print("─" * 72)
        print(f"【{case_id}】输入: {message}")
        print(f"  判定状态: {verdict.status}")
        print(f"  澄清问句: {verdict.question[:120]}")
        print(f"  候选选项: {verdict.candidate_ids}")
        if verdict.rationale:
            print(f"  理由: {verdict.rationale}")
        if verdict.passed:
            print(_green("  结果: Pass"))
        else:
            print(_red("  结果: Fail"))
            for err in errors:
                print(_red(f"    ✗ {err}"))

        await asyncio.sleep(0.35)

    return verdicts


def _create_router() -> MasterIntentRouter:
    from gibh_agent.core.llm_client import LLMClientFactory
    from gibh_agent.core.llm_cloud_providers import (
        get_default_chat_model,
        validate_and_resolve_model_name,
    )

    model_raw = (os.getenv("GIBH_INTENT_ROUTER_TEST_MODEL") or "").strip()
    model_name = validate_and_resolve_model_name(model_raw or None)
    if not model_raw:
        model_name = validate_and_resolve_model_name(get_default_chat_model())
    print(_cyan(f"使用模型: {model_name}"))
    llm = LLMClientFactory.create_for_model(model_name)
    return MasterIntentRouter(llm)


async def _main_async() -> int:
    _load_dotenv_if_present()
    _require_llm_env()
    _print_registry_banner()
    print()
    router = _create_router()
    verdicts = await run_ambiguous_suite(router)

    passed = sum(1 for v in verdicts if v.passed)
    total = len(verdicts)
    rate = passed / total if total else 0.0

    print()
    print("=" * 72)
    if rate >= PASS_THRESHOLD:
        print(_green(f"汇总: {passed}/{total} Passed ({rate * 100:.0f}%) — 达标 (≥{PASS_THRESHOLD * 100:.0f}%)"))
    else:
        print(_red(f"汇总: {passed}/{total} Passed ({rate * 100:.0f}%) — 未达标 (需 ≥{PASS_THRESHOLD * 100:.0f}%)"))
        failed = [v.case_id for v in verdicts if not v.passed]
        print(_yellow(f"失败用例: {', '.join(failed)}"))

    return 0 if rate >= PASS_THRESHOLD else 1


def main() -> None:
    raise SystemExit(asyncio.run(_main_async()))


# pytest 入口：逐条用例（共享同一 router 需 session fixture 时可再拆）
def test_ambiguous_cases_live_llm() -> None:
    """pytest -s 时跑完整 20 条（触网）。"""
    raise SystemExit(asyncio.run(_main_async()))


if __name__ == "__main__":
    main()
