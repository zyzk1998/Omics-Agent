# -*- coding: utf-8 -*-
"""
HITL 唤醒：从 waiting_for_hitl 恢复，拉取 LS 标注并生成《专家分析报告（最终版）》。
"""
from __future__ import annotations

import copy
import json
import logging
import asyncio
import time
from typing import Any, AsyncIterator, Dict, List, Optional, Tuple

from sqlalchemy.orm import Session as OrmSession
from sqlalchemy.orm.attributes import flag_modified

from gibh_agent.core.execution_snapshot import (
    apply_execution_snapshot_to_state,
    mirror_execution_snapshot_to_message_content,
    update_expert_report_in_state,
)
from gibh_agent.core.hitl_session_registry import get_hitl_session, register_hitl_session
from gibh_agent.core.session_runtime import SESSION_COMPLETED, SESSION_RUNNING, set_session_status
from gibh_agent.core.utils import sanitize_for_json, scrub_markdown_poison_urls
from gibh_agent.db.models import Message as MessageModel, Session as SessionModel
from gibh_agent.utils.ls_client import LabelStudioClient, LabelStudioClientError

logger = logging.getLogger(__name__)


def _load_latest_agent_snapshot(
    db: OrmSession, session_id: str, owner_id: str
) -> Tuple[Optional[MessageModel], Dict[str, Any]]:
    from gibh_agent.db.message_queries import load_latest_agent_snapshot

    return load_latest_agent_snapshot(db, session_id, owner_id)


def _enrich_hitl_registry_from_snapshot(
    session_id: str,
    snap: Dict[str, Any],
    reg: Optional[Dict[str, Any]] = None,
    hitl: Optional[Dict[str, Any]] = None,
) -> Dict[str, Any]:
    """completed 会话重入：从 execution_snapshot / workflow / hitl 补全 registry 上下文。"""
    merged: Dict[str, Any] = dict(reg or {})
    hitl = hitl if isinstance(hitl, dict) else {}
    ex = snap.get("execution_snapshot") if isinstance(snap.get("execution_snapshot"), dict) else {}
    wf = snap.get("workflow") if isinstance(snap.get("workflow"), dict) else {}
    skill = snap.get("skill_last_tool") if isinstance(snap.get("skill_last_tool"), dict) else {}

    def _pick(key: str, *sources: Any) -> None:
        if merged.get(key) is not None:
            return
        for src in sources:
            if isinstance(src, dict) and src.get(key) is not None:
                merged[key] = src[key]
                return

    for key in (
        "project_id",
        "ls_project_id",
        "ls_url",
        "ls_project_url",
        "scenario_type",
        "workflow_name",
        "output_dir",
        "agent_key",
        "skill_id",
        "image_path",
        "file_path",
    ):
        _pick(key, hitl, ex, wf, reg, skill.get("tool_result") if isinstance(skill.get("tool_result"), dict) else {})

    for list_key in ("image_paths", "uploaded_file_paths"):
        if merged.get(list_key):
            continue
        for src in (hitl, ex, reg):
            if isinstance(src, dict) and src.get(list_key):
                merged[list_key] = list(src.get(list_key) or [])
                break

    pid = merged.get("project_id") or merged.get("ls_project_id")
    if pid is not None:
        try:
            merged["project_id"] = int(pid)
        except (TypeError, ValueError):
            pass

    sid = str(session_id or "").strip()
    if sid:
        register_hitl_session(sid, merged)
    return merged


def _extract_steps_details(snap: Dict[str, Any]) -> List[Dict[str, Any]]:
    ex = snap.get("execution_snapshot") if isinstance(snap.get("execution_snapshot"), dict) else {}
    steps = ex.get("steps_details")
    if isinstance(steps, list) and steps:
        return steps
    rep = snap.get("report") if isinstance(snap.get("report"), dict) else {}
    rd = rep.get("report_data") if isinstance(rep.get("report_data"), dict) else {}
    steps2 = rd.get("steps_details")
    if isinstance(steps2, list):
        return steps2
    if isinstance(snap.get("steps"), list):
        return snap["steps"]
    return []


def _fetch_ls_annotations(project_id: int) -> Tuple[Any, str]:
    client = LabelStudioClient()
    summary_lines: List[str] = []
    export_data: Any = None
    try:
        export_data = client.export_annotations(project_id)
    except LabelStudioClientError as exc:
        logger.warning("[HITL Resume] export_annotations failed: %s", exc)
        try:
            tasks = client.list_tasks(project_id)
            export_data = tasks
        except LabelStudioClientError as exc2:
            return None, f"（未能从 Label Studio 拉取标注：{exc2}）"

    if isinstance(export_data, list):
        summary_lines.append(f"- 导出任务/标注条目数: {len(export_data)}")
        for i, item in enumerate(export_data[:5]):
            if isinstance(item, dict):
                ann = item.get("annotations") or item.get("completions")
                summary_lines.append(f"- 条目 #{i + 1}: annotations={len(ann) if isinstance(ann, list) else 0}")
    elif isinstance(export_data, dict):
        summary_lines.append(f"- 导出结构键: {', '.join(list(export_data.keys())[:8])}")
    else:
        summary_lines.append("- 导出格式: 非结构化")

    return export_data, "\n".join(summary_lines)


def _seed_state_snapshot_from_prior(prior: Dict[str, Any]) -> Dict[str, Any]:
    """HITL 唤醒/跳过：在既有快照上增量更新，禁止用空快照覆盖时光机数据。"""
    base = copy.deepcopy(prior) if isinstance(prior, dict) and prior else {}
    defaults: Dict[str, Any] = {
        "text": "",
        "reasoning": "",
        "workflow": None,
        "steps": [],
        "process_log": [],
        "report": None,
    }
    for key, val in defaults.items():
        if key not in base or base[key] is None:
            base[key] = val if not isinstance(val, (list, dict)) else copy.deepcopy(val)
    base["_start_time"] = time.time()
    return base


def _persist_agent_message_snapshot(
    msg: MessageModel,
    state_snapshot: Dict[str, Any],
    *,
    steps_details: Optional[List[Dict[str, Any]]] = None,
    workflow_name: str = "",
) -> None:
    """合并 execution_snapshot 并双写 content 根级字段。"""
    if steps_details:
        apply_execution_snapshot_to_state(
            state_snapshot,
            steps_details,
            workflow_name=workflow_name or None,
        )
    content = msg.content if isinstance(msg.content, dict) else {}
    content["state_snapshot"] = sanitize_for_json(state_snapshot)
    mirror_execution_snapshot_to_message_content(content, state_snapshot)
    ex = state_snapshot.get("execution_snapshot")
    if isinstance(ex, dict):
        for key in (
            "corpus_sft_json",
            "corpus_sft_records",
            "corpus_standard",
            "corpus_modality",
            "corpus_count",
        ):
            if ex.get(key) is not None:
                content[key] = sanitize_for_json(ex.get(key))
    msg.content = content
    flag_modified(msg, "content")


def persist_session_snapshot_to_db(
    db: OrmSession,
    session_id: str,
    owner_id: str,
    state_snapshot: Dict[str, Any],
    *,
    steps_details: Optional[List[Dict[str, Any]]] = None,
    workflow_name: str = "",
) -> bool:
    """将含 corpus 在内的 state_snapshot 写入最新 agent 消息并 commit。"""
    sid = str(session_id or "").strip()
    if not sid:
        return False
    msg, _snap = _load_latest_agent_snapshot(db, sid, owner_id)
    if not msg:
        return False
    _persist_agent_message_snapshot(
        msg,
        state_snapshot,
        steps_details=steps_details,
        workflow_name=workflow_name,
    )
    db.commit()
    return True


_CORPUS_SKILL_ID = "skill_corpus_data_processing"
_CORPUS_WORKFLOW_NAME = "科学语料数据加工"


def _is_standalone_corpus_skill(
    snap: Dict[str, Any],
    reg: Optional[Dict[str, Any]] = None,
    hitl: Optional[Dict[str, Any]] = None,
) -> bool:
    reg = reg or {}
    hitl = hitl or {}
    wf = snap.get("workflow") if isinstance(snap.get("workflow"), dict) else {}
    if reg.get("agent_key") == "corpus_processing_agent":
        return True
    if str(reg.get("skill_id") or hitl.get("skill_id") or wf.get("skill_id") or "").strip() == _CORPUS_SKILL_ID:
        return True
    for name in (
        reg.get("workflow_name"),
        hitl.get("workflow_name"),
        wf.get("workflow_name"),
    ):
        if str(name or "").strip() == _CORPUS_WORKFLOW_NAME:
            return True
    st = snap.get("skill_last_tool")
    if isinstance(st, dict) and st.get("tool_name") == _CORPUS_SKILL_ID:
        return True
    return False


def _extract_standalone_corpus_media_paths(
    snap: Dict[str, Any],
    reg: Optional[Dict[str, Any]] = None,
    hitl: Optional[Dict[str, Any]] = None,
) -> List[str]:
    """独立语料技能 Skip：从会话快照 / HITL 注册表提取用户上传的图像或文件路径。"""
    reg = reg or {}
    hitl = hitl or {}
    paths: List[str] = []
    seen: set[str] = set()

    def _add(raw: Any) -> None:
        p = str(raw or "").strip()
        if not p or p in seen:
            return
        seen.add(p)
        paths.append(p)

    for src in (reg, hitl):
        if not isinstance(src, dict):
            continue
        for item in src.get("image_paths") or []:
            _add(item)
        _add(src.get("image_path"))
        _add(src.get("file_path"))
        for item in src.get("uploaded_file_paths") or []:
            _add(item)

    ex = snap.get("execution_snapshot") if isinstance(snap.get("execution_snapshot"), dict) else {}
    for att in ex.get("input_draft_attachments") or snap.get("input_draft_attachments") or []:
        if not isinstance(att, dict):
            continue
        _add(att.get("path") or att.get("file_path") or att.get("absolute_path") or att.get("name"))

    wf = snap.get("workflow") if isinstance(snap.get("workflow"), dict) else {}
    for item in wf.get("file_paths") or wf.get("files") or []:
        if isinstance(item, dict):
            _add(item.get("path") or item.get("file_path") or item.get("name"))
        else:
            _add(item)

    st = snap.get("skill_last_tool")
    if isinstance(st, dict):
        tr = st.get("tool_result")
        if isinstance(tr, dict):
            for item in tr.get("image_paths") or []:
                _add(item)
            _add(tr.get("image_path"))
            _add(tr.get("file_path"))

    from gibh_agent.tools.hitl_tools import parse_image_path_inputs

    expanded: List[str] = []
    for p in paths:
        expanded.extend(parse_image_path_inputs(p))
    out: List[str] = []
    seen2: set[str] = set()
    for p in expanded:
        if p and p not in seen2:
            seen2.add(p)
            out.append(p)
    return out


def _resolve_draft_expert_report(snap: Dict[str, Any]) -> str:
    """提取《专家分析报告（初稿）》或当前会话中的报告正文，供 Skip 代偿语料使用。"""
    ex = snap.get("execution_snapshot") if isinstance(snap.get("execution_snapshot"), dict) else {}
    md = str(ex.get("expert_report_markdown") or snap.get("expert_report_markdown") or "").strip()
    if md:
        return md
    rep = snap.get("report") if isinstance(snap.get("report"), dict) else {}
    rd = rep.get("report_data") if isinstance(rep.get("report_data"), dict) else {}
    return str(rd.get("report") or rep.get("report") or "").strip()


def _pick_agent_for_workflow(
    agent_root: Any,
    workflow_name: str,
    steps: List[Dict],
    *,
    scenario_type: str = "",
    agent_key: str = "",
) -> Any:
    agents = getattr(agent_root, "agents", None) or {}
    if agent_key and agent_key in agents:
        return agents[agent_key]
    scenario = str(scenario_type or "").strip()
    if scenario == "generic_corpus_processing" or "语料" in (workflow_name or ""):
        if "corpus_processing_agent" in agents:
            return agents["corpus_processing_agent"]
    wnl = (workflow_name or "").lower()
    domain = "Metabolomics"
    if "rna" in wnl or "转录" in workflow_name or "单细胞" in workflow_name:
        domain = "RNA"
    elif "spatial" in wnl or "空间" in workflow_name:
        domain = "Spatial"
    elif "radiomics" in wnl or "影像" in workflow_name:
        domain = "Radiomics"
    elif "genomics" in wnl or "基因组" in workflow_name:
        domain = "genomics"
    mapping = {
        "RNA": "rna_agent",
        "Metabolomics": "metabolomics_agent",
        "Spatial": "spatial_agent",
        "Radiomics": "radiomics_agent",
        "genomics": "dna_agent",
    }
    key = mapping.get(domain)
    if key and key in agents:
        return agents[key]
    return next(iter(agents.values()), None) if agents else None


async def stream_resume_from_hitl(
    *,
    orchestrator: Any,
    session_id: str,
    owner_id: str,
    db: OrmSession,
    project_id: Optional[int] = None,
    trigger: str = "api",
) -> AsyncIterator[str]:
    """
    异步 SSE 流：唤醒 HITL → 生成最终专家报告 → 覆写语料 JSON → 更新 DB → done。
    支持 completed 会话二次重入（重新标注覆盖金标准）。
    """
    sid = str(session_id or "").strip()
    seed_snap: Dict[str, Any] = {"text": "", "process_log": []}

    try:
        session = db.query(SessionModel).filter(SessionModel.id == sid).first()
        if not session or session.owner_id != owner_id:
            yield orchestrator._emit_sse(
                seed_snap,
                "error",
                {"message": "会话不存在或无权访问"},
            )
            yield orchestrator._emit_sse(seed_snap, "done", {"status": "error"})
            return

        msg, snap = _load_latest_agent_snapshot(db, sid, owner_id)
        if not msg:
            yield orchestrator._emit_sse(
                seed_snap,
                "error",
                {"message": "未找到可恢复的 Agent 消息"},
            )
            yield orchestrator._emit_sse(seed_snap, "done", {"status": "error"})
            return

        state_snapshot: Dict[str, Any] = _seed_state_snapshot_from_prior(snap)

        def emit(event_type: str, data: Dict[str, Any]) -> str:
            return orchestrator._emit_sse(state_snapshot, event_type, data)

        yield emit("status", {"content": "正在唤醒 HITL 流程…", "state": "running"})

        hitl = snap.get("hitl") if isinstance(snap.get("hitl"), dict) else {}
        reg = _enrich_hitl_registry_from_snapshot(sid, snap, get_hitl_session(sid) or {}, hitl)
        pid = project_id or hitl.get("project_id") or hitl.get("ls_project_id") or reg.get("project_id")
        if pid is None:
            yield emit("error", {"message": "缺少 Label Studio project_id，无法拉取专家标注"})
            yield emit("done", {"status": "error"})
            return
        try:
            pid_int = int(pid)
        except (TypeError, ValueError):
            yield emit("error", {"message": f"无效的 Label Studio project_id: {pid!r}"})
            yield emit("done", {"status": "error"})
            return

        set_session_status(db, sid, SESSION_RUNNING, owner_id=owner_id)

        steps_details = _extract_steps_details(snap)
        workflow_name = (
            hitl.get("workflow_name")
            or reg.get("workflow_name")
            or (snap.get("workflow") or {}).get("workflow_name")
            or "工作流"
        )
        output_dir = reg.get("output_dir") or hitl.get("output_dir")

        yield emit(
            "status",
            {"content": f"正在从 Label Studio 拉取标注 (project {pid_int})…", "state": "running"},
        )

        annotations, ann_summary = _fetch_ls_annotations(pid_int)

        register_hitl_session(
            sid,
            {
                **reg,
                "project_id": pid_int,
                "workflow_name": workflow_name,
                "output_dir": output_dir,
                "last_resume_trigger": trigger,
            },
        )

        steps_results = []
        for sd in steps_details:
            sr = sd.get("step_result") if isinstance(sd.get("step_result"), dict) else sd
            if isinstance(sr, dict):
                steps_results.append(sr)

        summary_context: Dict[str, Any] = {
            "has_failures": any(s.get("status") == "error" for s in steps_details),
            "failed_steps": [s for s in steps_details if s.get("status") == "error"],
            "successful_steps": [s for s in steps_details if s.get("status") in ("success", "hitl_required")],
            "workflow_status": "hitl_resumed",
            "hitl_annotations_summary": ann_summary,
            "hitl_resume": True,
            "hitl_scenario": hitl.get("scenario_type") or reg.get("scenario_type"),
        }
        if annotations is not None:
            summary_context["hitl_annotations_raw"] = annotations
            try:
                summary_context["hitl_annotations_json"] = json.dumps(
                    annotations, ensure_ascii=False, default=str
                )[:12000]
            except (TypeError, ValueError):
                pass

        target_agent = _pick_agent_for_workflow(
            orchestrator.agent,
            workflow_name,
            steps_details,
            scenario_type=str(hitl.get("scenario_type") or reg.get("scenario_type") or ""),
            agent_key=str(reg.get("agent_key") or ""),
        )
        summary: Optional[str] = None

        yield emit("status", {"content": "正在生成《专家分析报告（最终版）》…", "state": "generating_report"})

        if target_agent and hasattr(target_agent, "_generate_analysis_summary"):
            omics_type = workflow_name
            try:
                summary = await target_agent._generate_analysis_summary(
                    steps_results,
                    omics_type=omics_type,
                    workflow_name=workflow_name,
                    summary_context=summary_context,
                    output_dir=output_dir,
                )
            except Exception as exc:
                logger.exception("[HITL Resume] LLM summary failed: %s", exc)

        if not summary or not str(summary).strip():
            summary = (
                "## 专家分析报告（最终版）\n\n"
                "本报告在 **Label Studio 专家标注完成** 后自动生成。\n\n"
                "### 专家复核摘要\n\n"
                f"{ann_summary}\n\n"
                "### 说明\n\n"
                "已合并专家标注导出；如需完整生物学解读，请确认 LLM 服务可用后再次提交标注。"
            )
        else:
            if "最终版" not in summary and "专家分析" not in summary:
                summary = "## 专家分析报告（最终版）\n\n" + summary.lstrip()

        summary = scrub_markdown_poison_urls(summary.strip())
        update_expert_report_in_state(state_snapshot, summary, is_final=True)
        state_snapshot["hitl_pending"] = False
        state_snapshot["hitl_resumed"] = True
        state_snapshot["hitl_reannotation_enabled"] = True
        ls_url = (
            hitl.get("ls_url")
            or hitl.get("ls_project_url")
            or reg.get("ls_url")
            or reg.get("ls_project_url")
            or ""
        )
        if not ls_url:
            ls_url = f"/label-studio/projects/{pid_int}/data"
        state_snapshot["hitl"] = sanitize_for_json(
            {
                "status": "hitl_completed",
                "project_id": pid_int,
                "ls_project_id": pid_int,
                "ls_url": ls_url,
                "ls_project_url": ls_url,
                "scenario_type": hitl.get("scenario_type") or reg.get("scenario_type"),
                "workflow_name": workflow_name,
                "output_dir": output_dir,
                "reannotation_enabled": True,
            }
        )
        if annotations is not None:
            state_snapshot["hitl_annotations_export"] = sanitize_for_json(annotations)
            ex_ann = state_snapshot.get("execution_snapshot")
            if not isinstance(ex_ann, dict):
                ex_ann = {}
                state_snapshot["execution_snapshot"] = ex_ann
            ex_ann["hitl_annotations"] = sanitize_for_json(annotations)
            ex_ann["hitl_final"] = True

        diagnosis_response = {
            "report_data": {
                "report": summary,
                "workflow_name": workflow_name,
                "hitl_final": True,
            }
        }
        yield emit("diagnosis", diagnosis_response)
        yield emit(
            "result",
            {
                "report_data": {
                    "report": summary,
                    "steps_details": steps_details,
                    "workflow_name": workflow_name,
                }
            },
        )
        async for _corpus_sse in orchestrator._stream_session_corpus_terminal(
            state_snapshot,
            session_id=sid,
            steps_details=steps_details,
            expert_report_md=summary or "",
            output_dir=output_dir,
            hitl_annotations=annotations,
            hitl_meta=reg,
            skip_if_pending_hitl=False,
            db=db,
            owner_id=owner_id,
            workflow_name=workflow_name,
        ):
            yield _corpus_sse
            await asyncio.sleep(0.01)
        _ex_after = state_snapshot.get("execution_snapshot")
        if isinstance(_ex_after, dict) and _ex_after.get("steps_details"):
            steps_details = _ex_after["steps_details"]
        yield emit("status", {"content": "专家报告（最终版）已就绪", "state": "completed"})
        yield emit("done", {"status": "success", "hitl_resumed": True})

        msg_fresh, _ = _load_latest_agent_snapshot(db, sid, owner_id)
        if msg_fresh:
            msg = msg_fresh
        _persist_agent_message_snapshot(
            msg,
            state_snapshot,
            steps_details=steps_details,
            workflow_name=workflow_name,
        )
        set_session_status(db, sid, SESSION_COMPLETED, owner_id=owner_id)
        db.commit()

        yield orchestrator._format_sse("state_snapshot", state_snapshot)
    except Exception as exc:
        logger.exception("[HITL Resume] stream_resume_from_hitl failed session=%s: %s", sid, exc)
        try:
            db.rollback()
        except Exception:
            pass
        yield orchestrator._emit_sse(seed_snap, "error", {"message": str(exc)})
        yield orchestrator._emit_sse(seed_snap, "done", {"status": "error"})


def finalize_hitl_skip(
    db: OrmSession,
    session_id: str,
    owner_id: str,
    *,
    trigger: str = "skip",
) -> Dict[str, Any]:
    """
    用户跳过 HITL：关闭复核入口、保留初稿并标记会话完成。
    SFT 语料由 stream_skip_from_hitl / 编排器终点站生成（不再延迟至入库）。
    """
    sid = str(session_id or "").strip()
    session = db.query(SessionModel).filter(SessionModel.id == sid).first()
    if not session or session.owner_id != owner_id:
        return {"status": "error", "message": "会话不存在或无权访问"}

    msg, snap = _load_latest_agent_snapshot(db, sid, owner_id)
    if not msg:
        return {"status": "error", "message": "未找到可恢复的 Agent 消息"}

    reg = get_hitl_session(sid) or {}
    hitl = snap.get("hitl") if isinstance(snap.get("hitl"), dict) else {}
    draft_md = _resolve_draft_expert_report(snap)

    register_hitl_session(
        sid,
        {**reg, "last_resume_trigger": trigger, "hitl_skipped": True},
    )

    state_snapshot = _seed_state_snapshot_from_prior(snap)
    state_snapshot["hitl_pending"] = False
    state_snapshot["hitl_skipped"] = True
    state_snapshot["hitl_resumed"] = False
    state_snapshot.pop("hitl", None)
    if draft_md:
        state_snapshot["text"] = draft_md
        update_expert_report_in_state(state_snapshot, draft_md, is_final=False)
        ex_skip = state_snapshot.get("execution_snapshot")
        if isinstance(ex_skip, dict):
            ex_skip["hitl_skipped"] = True
            ex_skip["hitl_final"] = False

    steps_details = _extract_steps_details(snap)
    workflow_name = (
        hitl.get("workflow_name")
        or reg.get("workflow_name")
        or (snap.get("workflow") or {}).get("workflow_name")
        or ""
    )
    _persist_agent_message_snapshot(
        msg,
        state_snapshot,
        steps_details=steps_details,
        workflow_name=workflow_name,
    )
    set_session_status(db, sid, SESSION_COMPLETED, owner_id=owner_id)
    db.commit()

    return {
        "status": "skipped",
        "message": "已跳过专家复核，保留初稿报告。SFT 语料将在会话收尾自动生成。",
        "hitl_skipped": True,
    }


async def stream_skip_from_hitl(
    *,
    orchestrator: Any,
    session_id: str,
    owner_id: str,
    db: OrmSession,
    trigger: str = "skip",
) -> AsyncIterator[str]:
    """轻量 SSE：同步 Skip 状态并在收尾生成银标准语料（含独立语料技能 AI 代偿）。"""
    sid = str(session_id or "").strip()

    def emit(state_snapshot: Dict[str, Any], event_type: str, data: Dict[str, Any]) -> str:
        return orchestrator._emit_sse(state_snapshot, event_type, data)

    seed: Dict[str, Any] = {
        "text": "",
        "reasoning": "",
        "workflow": None,
        "steps": [],
        "process_log": [],
        "report": None,
        "_start_time": time.time(),
    }
    yield emit(seed, "status", {"content": "已跳过专家复核，准备生成语料…", "state": "running"})
    result = finalize_hitl_skip(db, sid, owner_id, trigger=trigger)
    if result.get("status") == "error":
        yield emit(seed, "error", {"message": result.get("message")})
        yield emit(seed, "done", {"status": "error"})
        return

    msg, snap_after = _load_latest_agent_snapshot(db, sid, owner_id)
    snap_after = snap_after or {}
    state_snapshot = _seed_state_snapshot_from_prior(snap_after)
    state_snapshot["hitl_pending"] = False
    state_snapshot["hitl_skipped"] = True
    state_snapshot.pop("hitl", None)

    reg = get_hitl_session(sid) or {}
    hitl = snap_after.get("hitl") if isinstance(snap_after.get("hitl"), dict) else {}
    is_corpus_skill = _is_standalone_corpus_skill(snap_after, reg, hitl)
    media_paths = _extract_standalone_corpus_media_paths(snap_after, reg, hitl) if is_corpus_skill else []
    steps_details = _extract_steps_details(snap_after)
    workflow_name = (
        hitl.get("workflow_name")
        or reg.get("workflow_name")
        or (snap_after.get("workflow") or {}).get("workflow_name")
        or (_CORPUS_WORKFLOW_NAME if is_corpus_skill else "")
    )
    output_dir = reg.get("output_dir") or hitl.get("output_dir")
    draft_md = _resolve_draft_expert_report(snap_after)

    if draft_md:
        state_snapshot["text"] = draft_md
        yield emit(
            state_snapshot,
            "result",
            {
                "report_data": {
                    "report": draft_md,
                    "hitl_skipped": True,
                }
            },
        )
    elif is_corpus_skill:
        skip_msg = (
            "已跳过 Label Studio 手动标注，将基于您上传的图像/文件自动生成 SFT 语料 JSON。"
        )
        state_snapshot["text"] = skip_msg
        yield emit(state_snapshot, "message", {"content": skip_msg})

    hitl_meta = {**reg, **hitl}
    if media_paths:
        hitl_meta["image_paths"] = media_paths
        hitl_meta["image_path"] = media_paths[0]
        hitl_meta["uploaded_file_paths"] = media_paths

    async for _corpus_sse in orchestrator._stream_session_corpus_terminal(
        state_snapshot,
        session_id=sid,
        steps_details=steps_details,
        expert_report_md=draft_md,
        output_dir=output_dir,
        hitl_meta=hitl_meta,
        skip_hitl=True,
        skip_if_pending_hitl=False,
        standalone_media_paths=media_paths if is_corpus_skill else None,
        db=db,
        owner_id=owner_id,
        workflow_name=workflow_name,
    ):
        yield _corpus_sse
        await asyncio.sleep(0.01)

    yield emit(
        state_snapshot,
        "status",
        {
            "content": "已跳过专家复核" if not is_corpus_skill else "已自动生成语料 JSON",
            "state": "completed",
        },
    )
    yield emit(state_snapshot, "done", {"status": "success", "hitl_skipped": True})
    yield orchestrator._format_sse("state_snapshot", state_snapshot)
