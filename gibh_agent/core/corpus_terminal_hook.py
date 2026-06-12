# -*- coding: utf-8 -*-
"""
Task / Skill 执行生命周期终点站：统一生成 SFT 语料并推送 corpus_asset SSE。

- 有 HITL 标注 → Gold Standard（VLM / NLP）
- 无标注（跳过 HITL 或纯自动化技能）→ Silver / Machine-Generated

落盘：corpus_archive/{session_id}/ 与 corpus_hitl/{session_id}/sft_corpus_*.json
"""
from __future__ import annotations

import json
import logging
from datetime import datetime
from pathlib import Path
from typing import Any, AsyncIterator, Callable, Dict, List, Optional

from gibh_agent.core.corpus_gatekeeper import (
    CORPUS_MODALITY_NLP,
    default_corpus_archive_dir,
    resolve_existing_corpus_archive,
)
from gibh_agent.core.execution_snapshot import update_corpus_asset_in_state
from gibh_agent.core.hitl_resume import _extract_steps_details, _resolve_draft_expert_report
from gibh_agent.core.utils import sanitize_for_json
from gibh_agent.tools.hitl_tools import parse_image_path_inputs

logger = logging.getLogger(__name__)

CORPUS_STANDARD_GOLD = "gold"
CORPUS_STANDARD_SILVER = "silver"

CORPUS_TERMINAL_STEP_ID = "downstream_corpus_sft"
CORPUS_TERMINAL_STEP_NAME = "科学语料数据加工"
CORPUS_STATUS_MSG_START = "📥 正在执行下游流程：自动化科学语料数据加工..."
CORPUS_STATUS_MSG_DONE = "✅ 语料 JSON 生成完毕"


def _execution_snapshot(snap: Dict[str, Any]) -> Dict[str, Any]:
    ex = snap.get("execution_snapshot")
    return ex if isinstance(ex, dict) else {}


def _resolve_expert_report_md(snap: Dict[str, Any], override: str = "") -> str:
    md = str(override or "").strip()
    if md:
        return md
    ex = _execution_snapshot(snap)
    md = str(ex.get("expert_report_markdown") or snap.get("expert_report_markdown") or "").strip()
    if md:
        return md
    return _resolve_draft_expert_report(snap)


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


def _build_summary_context(
    *,
    snap: Dict[str, Any],
    expert_report_md: str,
    hitl_annotations: Any,
    hitl_meta: Optional[Dict[str, Any]],
    skip_hitl: bool,
    steps_details: Optional[List[Dict[str, Any]]],
    standalone_media_paths: Optional[List[str]] = None,
) -> Dict[str, Any]:
    ex = _execution_snapshot(snap)
    hitl = snap.get("hitl") if isinstance(snap.get("hitl"), dict) else {}
    reg = hitl_meta or {}
    draft = _resolve_draft_expert_report(snap) or expert_report_md
    ann_json = ""
    if hitl_annotations is not None:
        try:
            ann_json = json.dumps(hitl_annotations, ensure_ascii=False, default=str)[:12000]
        except (TypeError, ValueError):
            ann_json = ""

    stored_ann = ex.get("hitl_annotations")
    ctx: Dict[str, Any] = {
        "expert_report_draft": draft,
        "draft_expert_report": draft,
        "expert_report_markdown": expert_report_md or draft,
        "hitl_skipped": bool(skip_hitl or ex.get("hitl_skipped") or snap.get("hitl_skipped")),
        "hitl_resume": bool(ex.get("hitl_final") or snap.get("hitl_resumed")),
        "workflow_status": "corpus_terminal",
        "steps_details": steps_details or [],
    }
    if hitl_annotations is not None:
        ctx["hitl_annotations_raw"] = hitl_annotations
        if ann_json:
            ctx["hitl_annotations_json"] = ann_json
    elif stored_ann is not None:
        ctx["hitl_annotations_raw"] = stored_ann
        try:
            ctx["hitl_annotations_json"] = json.dumps(stored_ann, ensure_ascii=False, default=str)[:12000]
        except (TypeError, ValueError):
            pass
    hitl_images = parse_image_path_inputs(
        hitl.get("image_paths") or hitl.get("image_path") or reg.get("image_path") or ""
    )
    if hitl_images:
        ctx["hitl_image_paths"] = hitl_images
    if standalone_media_paths:
        merged = list(ctx.get("hitl_image_paths") or [])
        for p in standalone_media_paths:
            ps = str(p or "").strip()
            if ps and ps not in merged:
                merged.append(ps)
        if merged:
            ctx["hitl_image_paths"] = merged
        ctx["standalone_corpus_skill"] = True
        ctx["image_paths"] = merged
    return ctx


def _build_corpus_terminal_step(*, status: str, message: str) -> Dict[str, Any]:
    st = str(status or "running").lower()
    return {
        "step_id": CORPUS_TERMINAL_STEP_ID,
        "step_name": CORPUS_TERMINAL_STEP_NAME,
        "name": CORPUS_TERMINAL_STEP_NAME,
        "tool_name": "corpus_sft_terminal",
        "status": st,
        "message": message,
        "step_result": {
            "status": st,
            "message": message,
            "data": {"phase": "corpus_terminal", "downstream": True},
        },
    }


def _upsert_corpus_terminal_step(
    steps_details: Optional[List[Dict[str, Any]]],
    step: Dict[str, Any],
) -> List[Dict[str, Any]]:
    steps = list(steps_details or [])
    for idx, item in enumerate(steps):
        if not isinstance(item, dict):
            continue
        if (
            item.get("step_id") == CORPUS_TERMINAL_STEP_ID
            or item.get("tool_name") == "corpus_sft_terminal"
        ):
            steps[idx] = {**item, **step}
            return steps
    steps.append(step)
    return steps


def _emit_corpus_steps_sse(
    *,
    emit_sse: Callable[[str, Dict[str, Any]], str],
    state_snapshot: Dict[str, Any],
    steps_details: Optional[List[Dict[str, Any]]],
    step: Dict[str, Any],
    workflow_name: str = "",
) -> List[Dict[str, Any]]:
    from gibh_agent.core.execution_snapshot import apply_execution_snapshot_to_state

    updated = _upsert_corpus_terminal_step(steps_details, step)
    apply_execution_snapshot_to_state(
        state_snapshot,
        updated,
        workflow_name=workflow_name or None,
    )
    emit_sse(
        "step_result",
        {
            "report_data": {
                "steps_details": updated,
                "workflow_name": workflow_name or "",
            }
        },
    )
    return updated


def _write_corpus_hitl_snapshot(
    *,
    session_id: str,
    records: List[Any],
    standard: str,
    modality: str,
) -> Optional[Path]:
    import os

    sid = str(session_id or "unknown").strip()
    hitl_dir = Path(os.getenv("RESULTS_DIR", "/app/results")).expanduser() / "corpus_hitl" / sid
    hitl_dir.mkdir(parents=True, exist_ok=True)
    ts = datetime.utcnow().strftime("%Y%m%d_%H%M%S")
    path = hitl_dir / f"sft_corpus_{standard}_{ts}.json"
    try:
        path.write_text(
            json.dumps(sanitize_for_json(records), ensure_ascii=False, indent=2),
            encoding="utf-8",
        )
        return path.resolve()
    except OSError as exc:
        logger.warning("[CorpusTerminal] write corpus_hitl failed: %s", exc)
        return None


def _load_existing_corpus_bundle(snap: Dict[str, Any], session_id: str) -> Optional[Dict[str, Any]]:
    existing = resolve_existing_corpus_archive(snap, session_id)
    if not existing:
        return None
    dataset_path = existing / "dataset.json"
    try:
        records = json.loads(dataset_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return None
    if not isinstance(records, list) or not records:
        return None
    modality = CORPUS_MODALITY_NLP
    if (existing / "images").is_dir():
        from gibh_agent.core.corpus_gatekeeper import CORPUS_MODALITY_VLM

        modality = CORPUS_MODALITY_VLM
    ex = _execution_snapshot(snap)
    has_ann = ex.get("hitl_annotations") is not None or snap.get("hitl_annotations_export") is not None
    standard = CORPUS_STANDARD_GOLD if has_ann or ex.get("hitl_final") else CORPUS_STANDARD_SILVER
    return {
        "status": "success",
        "records": records,
        "archive_dir": str(existing),
        "modality": modality,
        "count": len(records),
        "standard": standard,
        "reused": True,
        "ai_only": standard == CORPUS_STANDARD_SILVER,
    }


async def finalize_session_corpus_bundle(
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
    standalone_media_paths: Optional[List[str]] = None,
) -> Dict[str, Any]:
    """生成或复用语料归档，返回 bundle（含 records / standard / modality）。"""
    sid = str(session_id or "").strip() or "unknown"
    steps_details = steps_details if steps_details is not None else _extract_steps_details(snap)
    expert_report_md = _resolve_expert_report_md(snap, expert_report_md)

    # 重新标注 / Resume 带 LS 导出：必须覆写归档，不可复用首期银标 bundle
    force_regenerate = hitl_annotations is not None
    existing_bundle = None if force_regenerate else _load_existing_corpus_bundle(snap, sid)
    if existing_bundle:
        hitl_path = _write_corpus_hitl_snapshot(
            session_id=sid,
            records=existing_bundle["records"],
            standard=existing_bundle["standard"],
            modality=existing_bundle["modality"],
        )
        if hitl_path:
            existing_bundle["corpus_hitl_path"] = str(hitl_path)
        return existing_bundle

    corpus_agent = _get_corpus_agent(agent_root)
    if not corpus_agent or not hasattr(corpus_agent, "generate_corpus_archive_bundle"):
        return {
            "status": "error",
            "message": "CorpusProcessingAgent 不可用，无法生成语料",
            "records": [],
            "standard": None,
            "modality": None,
            "count": 0,
        }

    steps_results: List[Dict[str, Any]] = []
    for sd in steps_details:
        sr = sd.get("step_result") if isinstance(sd.get("step_result"), dict) else sd
        if isinstance(sr, dict):
            steps_results.append(sr)

    summary_context = _build_summary_context(
        snap=snap,
        expert_report_md=expert_report_md,
        hitl_annotations=hitl_annotations,
        hitl_meta=hitl_meta,
        skip_hitl=skip_hitl,
        steps_details=steps_details,
        standalone_media_paths=standalone_media_paths,
    )

    has_annotations = (
        hitl_annotations is not None
        or summary_context.get("hitl_annotations_raw") is not None
    )
    standard = CORPUS_STANDARD_GOLD if has_annotations else CORPUS_STANDARD_SILVER

    archive_dir = default_corpus_archive_dir(sid)
    bundle = await corpus_agent.generate_corpus_archive_bundle(
        summary_context=summary_context,
        steps_results=steps_results,
        steps_details=steps_details,
        archive_dir=archive_dir,
        expert_report_md=expert_report_md,
        output_dir=output_dir,
    )
    if bundle.get("status") != "success":
        bundle.setdefault("records", [])
        bundle.setdefault("standard", standard)
        return bundle

    records = bundle.get("records") or []
    if not records and bundle.get("dataset_path"):
        try:
            records = json.loads(Path(bundle["dataset_path"]).read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError):
            records = []

    bundle["records"] = records
    bundle["standard"] = standard
    bundle["reused"] = False
    hitl_path = _write_corpus_hitl_snapshot(
        session_id=sid,
        records=records,
        standard=standard,
        modality=str(bundle.get("modality") or CORPUS_MODALITY_NLP),
    )
    if hitl_path:
        bundle["corpus_hitl_path"] = str(hitl_path)
    return bundle


def _corpus_asset_payload(bundle: Dict[str, Any]) -> Dict[str, Any]:
    records = bundle.get("records") or []
    corpus_json = json.dumps(sanitize_for_json(records), ensure_ascii=False, indent=2)
    return sanitize_for_json(
        {
            "status": bundle.get("status"),
            "standard": bundle.get("standard"),
            "modality": bundle.get("modality"),
            "count": bundle.get("count") or len(records),
            "archive_dir": bundle.get("archive_dir"),
            "corpus_hitl_path": bundle.get("corpus_hitl_path"),
            "corpus_json": corpus_json,
            "records": records,
            "reused": bool(bundle.get("reused")),
            "ai_only": bool(bundle.get("ai_only", bundle.get("standard") == CORPUS_STANDARD_SILVER)),
            "message": bundle.get("message") or "",
        }
    )


def _resolve_session_title_from_db(db: Any, session_id: str, owner_id: str = "") -> Optional[str]:
    meta = _resolve_session_meta_from_db(db, session_id, owner_id)
    return meta.get("title") if meta else None


def _resolve_session_meta_from_db(db: Any, session_id: str, owner_id: str = "") -> Optional[Dict[str, Any]]:
    if not db or not session_id:
        return None
    try:
        from gibh_agent.db.models import Session as SessionModel

        row = db.query(SessionModel).filter(SessionModel.id == session_id).first()
        if not row:
            return None
        if owner_id and getattr(row, "owner_id", None) and row.owner_id != owner_id:
            return None
        title = str(getattr(row, "title", "") or "").strip()
        return {
            "title": title or None,
            "created_at": getattr(row, "created_at", None),
        }
    except Exception as exc:
        logger.debug("[CorpusTerminal] resolve session meta failed: %s", exc)
        return None


async def stream_finalize_session_corpus(
    *,
    emit_sse: Callable[[str, Dict[str, Any]], str],
    state_snapshot: Dict[str, Any],
    session_id: Optional[str],
    agent_root: Any,
    steps_details: Optional[List[Dict[str, Any]]] = None,
    expert_report_md: str = "",
    output_dir: Optional[str] = None,
    hitl_annotations: Any = None,
    hitl_meta: Optional[Dict[str, Any]] = None,
    skip_hitl: bool = False,
    skip_if_pending_hitl: bool = False,
    standalone_media_paths: Optional[List[str]] = None,
    db: Any = None,
    owner_id: str = "",
    workflow_name: str = "",
    workspace_context: Optional[Dict[str, Any]] = None,
    local_workspace_mounted: bool = False,
    session_title: Optional[str] = None,
) -> AsyncIterator[str]:
    """在 completed 之前生成语料、写入快照、推送 corpus_asset SSE，并可选落库。"""
    if skip_if_pending_hitl and state_snapshot.get("hitl_pending"):
        return

    sid = str(session_id or "").strip()
    if not sid:
        logger.info("[CorpusTerminal] skip: missing session_id")
        return

    steps_details = steps_details if steps_details is not None else _extract_steps_details(state_snapshot)
    wf_name = workflow_name or str(
        (state_snapshot.get("workflow") or {}).get("workflow_name") or ""
    ).strip()

    running_step = _build_corpus_terminal_step(
        status="running",
        message=CORPUS_STATUS_MSG_START,
    )
    steps_details = _emit_corpus_steps_sse(
        emit_sse=emit_sse,
        state_snapshot=state_snapshot,
        steps_details=steps_details,
        step=running_step,
        workflow_name=wf_name,
    )

    yield emit_sse(
        "status",
        {"content": CORPUS_STATUS_MSG_START, "state": "generating_corpus"},
    )

    bundle = await finalize_session_corpus_bundle(
        session_id=sid,
        snap=state_snapshot,
        agent_root=agent_root,
        expert_report_md=expert_report_md,
        steps_details=steps_details,
        hitl_annotations=hitl_annotations,
        hitl_meta=hitl_meta,
        output_dir=output_dir,
        skip_hitl=skip_hitl,
        standalone_media_paths=standalone_media_paths,
    )

    if bundle.get("status") != "success":
        logger.warning(
            "[CorpusTerminal] partial session=%s: %s",
            sid,
            bundle.get("message"),
        )
        fail_step = _build_corpus_terminal_step(
            status="error",
            message=bundle.get("message") or "语料生成未完成",
        )
        _emit_corpus_steps_sse(
            emit_sse=emit_sse,
            state_snapshot=state_snapshot,
            steps_details=steps_details,
            step=fail_step,
            workflow_name=wf_name,
        )
        yield emit_sse(
            "status",
            {
                "content": bundle.get("message") or "语料生成未完成（可稍后重试）",
                "state": "error",
            },
        )
        if db and owner_id:
            from gibh_agent.core.hitl_resume import persist_session_snapshot_to_db

            persist_session_snapshot_to_db(
                db,
                sid,
                owner_id,
                state_snapshot,
                steps_details=steps_details,
                workflow_name=wf_name,
            )
        return

    asset_payload = _corpus_asset_payload(bundle)
    update_corpus_asset_in_state(state_snapshot, asset_payload)

    session_meta = _resolve_session_meta_from_db(db, sid, owner_id) or {}
    resolved_title = session_title or session_meta.get("title")
    session_created_at = session_meta.get("created_at")
    mounted_flag = bool(
        local_workspace_mounted or state_snapshot.get("local_workspace_mounted")
    )
    try:
        from gibh_agent.core.silent_local_corpus_deploy import (
            silent_deploy_corpus_bundle_to_local_mount,
        )
        from gibh_agent.core.storage.dual_path import notify_session_files_changed

        silent_delivery = silent_deploy_corpus_bundle_to_local_mount(
            session_id=sid,
            session_title=resolved_title,
            bundle=bundle,
            workspace_context=workspace_context,
            local_workspace_mounted=mounted_flag,
            owner_id=owner_id,
            session_created_at=session_created_at,
        )
        if silent_delivery:
            state_snapshot["local_silent_deploy"] = silent_delivery
            files_notify = notify_session_files_changed(
                session_id=sid,
                changed_paths=silent_delivery.get("changed_paths") or [],
                mount_tree=silent_delivery.get("mount_tree"),
            )
            if silent_delivery.get("client_deploy_package"):
                files_notify["client_deploy_package"] = silent_delivery["client_deploy_package"]
            asset_payload["silent_deploy"] = silent_delivery
            asset_payload["files_notification"] = files_notify
            relay_msg = (
                "语料已打包，等待本机 Sidecar 落盘"
                if silent_delivery.get("client_deploy_package")
                else "已静默落盘至本地项目目录"
            )
            yield emit_sse(
                "status",
                {
                    "content": relay_msg,
                    "state": "running",
                    "files_notification": files_notify,
                },
            )
    except Exception as exc:
        logger.warning("[CorpusTerminal] silent local deploy skipped: %s", exc)

    yield emit_sse("corpus_asset", asset_payload)

    done_step = _build_corpus_terminal_step(
        status="success",
        message=CORPUS_STATUS_MSG_DONE,
    )
    steps_details = _emit_corpus_steps_sse(
        emit_sse=emit_sse,
        state_snapshot=state_snapshot,
        steps_details=steps_details,
        step=done_step,
        workflow_name=wf_name,
    )

    label = "金标准" if bundle.get("standard") == CORPUS_STANDARD_GOLD else "银标准"
    yield emit_sse(
        "status",
        {"content": CORPUS_STATUS_MSG_DONE, "state": "completed"},
    )
    yield emit_sse(
        "status",
        {
            "content": f"SFT 语料已就绪（{label}，{bundle.get('count', 0)} 条）",
            "state": "running",
        },
    )

    if db and owner_id:
        from gibh_agent.core.hitl_resume import persist_session_snapshot_to_db

        persist_session_snapshot_to_db(
            db,
            sid,
            owner_id,
            state_snapshot,
            steps_details=steps_details,
            workflow_name=wf_name,
        )
