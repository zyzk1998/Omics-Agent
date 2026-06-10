# -*- coding: utf-8 -*-
"""
语料库就绪 (Corpus-Ready) · 一键入库守门人。

在 POST /api/ingestion/trigger 打包前盘点快照资产；若尚无标准语料归档，
则按模态（VLM / NLP）代偿生成 corpus_archive/ 后再交付 Data Packager。
"""
from __future__ import annotations

import logging
import os
from pathlib import Path
from typing import Any, Dict, List, Optional

from gibh_agent.core.hitl_resume import _extract_steps_details, _resolve_draft_expert_report
from gibh_agent.tools.hitl_tools import parse_image_path_inputs

logger = logging.getLogger(__name__)

CORPUS_MODALITY_VLM = "vlm"
CORPUS_MODALITY_NLP = "nlp"


def _execution_snapshot(snap: Dict[str, Any]) -> Dict[str, Any]:
    ex = snap.get("execution_snapshot")
    return ex if isinstance(ex, dict) else {}


def default_corpus_archive_dir(session_id: str) -> Path:
    sid = str(session_id or "unknown").strip()
    return Path(os.getenv("RESULTS_DIR", "/app/results")).expanduser() / "corpus_archive" / sid


def corpus_archive_is_ready(archive_dir: Path) -> bool:
    return archive_dir.is_dir() and (archive_dir / "dataset.json").is_file()


def resolve_existing_corpus_archive(
    snap: Dict[str, Any],
    session_id: str,
) -> Optional[Path]:
    """从快照或默认工作区定位已生成的 corpus_archive。"""
    ex = _execution_snapshot(snap)
    for key in ("corpus_archive_path", "corpus_archive_dir"):
        raw = str(ex.get(key) or "").strip()
        if raw:
            p = Path(raw).expanduser()
            if corpus_archive_is_ready(p):
                return p.resolve()
    default_dir = default_corpus_archive_dir(session_id)
    if corpus_archive_is_ready(default_dir):
        return default_dir.resolve()
    return None


def sniff_corpus_modality(
    *,
    steps_details: Optional[List[Dict[str, Any]]] = None,
    steps_results: Optional[List[Dict[str, Any]]] = None,
    hitl_annotations: Any = None,
    hitl_image_paths: Optional[List[str]] = None,
) -> str:
    """
    嗅探任务模态：存在可归档图像 → VLM；否则 → 纯文本 NLP SFT。
    """
    from gibh_agent.agents.corpus_processing_agent import CorpusProcessingAgent

    paths = list(hitl_image_paths or [])
    paths.extend(
        CorpusProcessingAgent.collect_resolvable_image_paths(
            steps_details,
            steps_results,
        )
    )
    if paths:
        return CORPUS_MODALITY_VLM

    if hitl_annotations is not None:
        if CorpusProcessingAgent.annotations_contain_images(hitl_annotations):
            return CORPUS_MODALITY_VLM

    return CORPUS_MODALITY_NLP


def _get_corpus_agent(agent_root: Any) -> Any:
    agents = getattr(agent_root, "agents", None) or {}
    agent = agents.get("corpus_processing_agent")
    if agent:
        return agent
    llm_clients = getattr(agent_root, "llm_clients", None) or {}
    logic_llm = llm_clients.get("logic") if isinstance(llm_clients, dict) else None
    if not logic_llm:
        logic_llm = getattr(agent_root, "llm_client", None)
    if logic_llm:
        from gibh_agent.agents.corpus_processing_agent import CorpusProcessingAgent
        from gibh_agent.core.prompt_manager import create_default_prompt_manager

        return CorpusProcessingAgent(
            llm_client=logic_llm,
            prompt_manager=create_default_prompt_manager(),
        )
    return None


def _build_ingestion_summary_context(
    *,
    snap: Dict[str, Any],
    expert_report_md: str,
    hitl_annotations: Any,
    hitl_meta: Optional[Dict[str, Any]],
    skip_hitl: bool,
) -> Dict[str, Any]:
    ex = _execution_snapshot(snap)
    hitl = snap.get("hitl") if isinstance(snap.get("hitl"), dict) else {}
    reg = hitl_meta or {}
    draft = _resolve_draft_expert_report(snap) or expert_report_md
    ann_json = ""
    if hitl_annotations is not None:
        import json

        try:
            ann_json = json.dumps(hitl_annotations, ensure_ascii=False, default=str)[:12000]
        except (TypeError, ValueError):
            ann_json = ""

    stored_ann = ex.get("hitl_annotations")
    ctx: Dict[str, Any] = {
        "expert_report_draft": draft,
        "draft_expert_report": draft,
        "expert_report_markdown": expert_report_md or draft,
        "hitl_skipped": bool(skip_hitl or ex.get("hitl_skipped")),
        "hitl_resume": bool(ex.get("hitl_final")),
        "workflow_status": "ingestion_gatekeeper",
    }
    if hitl_annotations is not None:
        ctx["hitl_annotations_raw"] = hitl_annotations
        if ann_json:
            ctx["hitl_annotations_json"] = ann_json
    elif stored_ann is not None:
        ctx["hitl_annotations_raw"] = stored_ann
    hitl_images = parse_image_path_inputs(
        hitl.get("image_paths") or hitl.get("image_path") or reg.get("image_path") or ""
    )
    if hitl_images:
        ctx["hitl_image_paths"] = hitl_images
    return ctx


async def ensure_corpus_ready_for_ingestion(
    *,
    session_id: str,
    snap: Dict[str, Any],
    agent_root: Any,
    expert_report_md: str = "",
    steps_details: Optional[List[Dict[str, Any]]] = None,
    hitl_annotations: Any = None,
    hitl_meta: Optional[Dict[str, Any]] = None,
    output_dir: Optional[str] = None,
    skip_hitl: bool = False,
) -> Dict[str, Any]:
    """
    入库前卡口：若 corpus_archive 未就绪则生成；已存在则复用。
    返回 bundle 供 Data Packager 纳入 tar.gz。
    """
    sid = str(session_id or "").strip()
    steps_details = steps_details if steps_details is not None else _extract_steps_details(snap)
    steps_results: List[Dict[str, Any]] = []
    for sd in steps_details:
        sr = sd.get("step_result") if isinstance(sd.get("step_result"), dict) else sd
        if isinstance(sr, dict):
            steps_results.append(sr)

    existing = resolve_existing_corpus_archive(snap, sid)
    if existing:
        dataset = existing / "dataset.json"
        modality = CORPUS_MODALITY_VLM if (existing / "images").is_dir() else CORPUS_MODALITY_NLP
        try:
            import json

            records = json.loads(dataset.read_text(encoding="utf-8"))
            count = len(records) if isinstance(records, list) else 0
        except (OSError, json.JSONDecodeError):
            count = 0
        logger.info("[CorpusGatekeeper] reuse archive session=%s path=%s", sid, existing)
        return {
            "status": "success",
            "archive_dir": str(existing),
            "modality": modality,
            "count": count,
            "reused": True,
        }

    corpus_agent = _get_corpus_agent(agent_root)
    if not corpus_agent or not hasattr(corpus_agent, "generate_corpus_archive_bundle"):
        return {
            "status": "error",
            "message": "CorpusProcessingAgent 不可用，无法生成语料归档",
            "archive_dir": None,
            "modality": None,
            "count": 0,
            "reused": False,
        }

    summary_context = _build_ingestion_summary_context(
        snap=snap,
        expert_report_md=expert_report_md,
        hitl_annotations=hitl_annotations,
        hitl_meta=hitl_meta,
        skip_hitl=skip_hitl,
    )
    summary_context["steps_details"] = steps_details

    archive_dir = default_corpus_archive_dir(sid)
    bundle = await corpus_agent.generate_corpus_archive_bundle(
        summary_context=summary_context,
        steps_results=steps_results,
        steps_details=steps_details,
        archive_dir=archive_dir,
        expert_report_md=expert_report_md or summary_context.get("expert_report_draft") or "",
        output_dir=output_dir,
    )
    bundle["reused"] = False
    return bundle
