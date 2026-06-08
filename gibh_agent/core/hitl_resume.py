# -*- coding: utf-8 -*-
"""
HITL 唤醒：从 waiting_for_hitl 恢复，拉取 LS 标注并生成《专家分析报告（最终版）》。
"""
from __future__ import annotations

import json
import logging
import time
from typing import Any, AsyncIterator, Dict, List, Optional, Tuple

from sqlalchemy.orm import Session as OrmSession
from sqlalchemy.orm.attributes import flag_modified

from gibh_agent.core.execution_snapshot import update_expert_report_in_state
from gibh_agent.core.hitl_session_registry import get_hitl_session, register_hitl_session
from gibh_agent.core.session_runtime import SESSION_COMPLETED, SESSION_RUNNING, set_session_status
from gibh_agent.core.utils import sanitize_for_json
from gibh_agent.db.models import Message as MessageModel, Session as SessionModel
from gibh_agent.utils.ls_client import LabelStudioClient, LabelStudioClientError

logger = logging.getLogger(__name__)


def _load_latest_agent_snapshot(
    db: OrmSession, session_id: str, owner_id: str
) -> Tuple[Optional[MessageModel], Dict[str, Any]]:
    session = db.query(SessionModel).filter(SessionModel.id == session_id).first()
    if not session or session.owner_id != owner_id:
        return None, {}
    msg = (
        db.query(MessageModel)
        .filter(MessageModel.session_id == session_id, MessageModel.role == "agent")
        .order_by(MessageModel.created_at.desc(), MessageModel.id.desc())
        .first()
    )
    if not msg:
        return None, {}
    content = msg.content if isinstance(msg.content, dict) else {}
    snap = content.get("state_snapshot") if isinstance(content.get("state_snapshot"), dict) else content
    return msg, snap if isinstance(snap, dict) else {}


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
    异步 SSE 流：唤醒 HITL → 生成最终专家报告 → 更新 DB → done。
    """
    sid = str(session_id or "").strip()
    state_snapshot: Dict[str, Any] = {
        "text": "",
        "reasoning": "",
        "workflow": None,
        "steps": [],
        "process_log": [],
        "report": None,
        "_start_time": time.time(),
    }

    def emit(event_type: str, data: Dict[str, Any]) -> str:
        return orchestrator._emit_sse(state_snapshot, event_type, data)

    yield emit("status", {"content": "正在唤醒 HITL 流程…", "state": "running"})

    session = db.query(SessionModel).filter(SessionModel.id == sid).first()
    if not session or session.owner_id != owner_id:
        yield emit("error", {"message": "会话不存在或无权访问"})
        yield emit("done", {"status": "error"})
        return

    msg, snap = _load_latest_agent_snapshot(db, sid, owner_id)
    if not msg:
        yield emit("error", {"message": "未找到可恢复的 Agent 消息"})
        yield emit("done", {"status": "error"})
        return

    reg = get_hitl_session(sid) or {}
    hitl = snap.get("hitl") if isinstance(snap.get("hitl"), dict) else {}
    pid = project_id or hitl.get("project_id") or reg.get("project_id")
    if pid is None:
        yield emit("error", {"message": "缺少 Label Studio project_id，无法拉取专家标注"})
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

    yield emit("status", {"content": f"正在从 Label Studio 拉取标注 (project {pid})…", "state": "running"})

    annotations, ann_summary = _fetch_ls_annotations(int(pid))

    register_hitl_session(
        sid,
        {
            **reg,
            "project_id": int(pid),
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
        try:
            summary_context["hitl_annotations_json"] = json.dumps(annotations, ensure_ascii=False, default=str)[:12000]
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
            "自动化 pipeline 已在 HITL 环节挂起；当前版本已合并专家标注导出。"
            "若需完整生物学解读，请确认 LLM 服务可用后再次点击「继续生成报告」。"
        )
    else:
        if "最终版" not in summary and "专家分析" not in summary:
            summary = "## 专家分析报告（最终版）\n\n" + summary.lstrip()

    update_expert_report_in_state(state_snapshot, summary.strip())
    state_snapshot["hitl_pending"] = False
    state_snapshot["hitl_resumed"] = True
    state_snapshot.pop("hitl", None)
    if annotations is not None:
        state_snapshot["hitl_annotations_export"] = sanitize_for_json(annotations)

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
    yield emit("status", {"content": "专家报告（最终版）已就绪", "state": "completed"})
    yield emit("done", {"status": "success", "hitl_resumed": True})

    content = msg.content if isinstance(msg.content, dict) else {}
    content["state_snapshot"] = sanitize_for_json(state_snapshot)
    if isinstance(content.get("execution_snapshot"), dict):
        content["execution_snapshot"]["expert_report_markdown"] = summary
        content["execution_snapshot"]["hitl_final"] = True
        if annotations is not None:
            content["execution_snapshot"]["hitl_annotations"] = sanitize_for_json(annotations)
    msg.content = content
    flag_modified(msg, "content")
    set_session_status(db, sid, SESSION_COMPLETED, owner_id=owner_id)
    db.commit()

    yield orchestrator._format_sse("state_snapshot", state_snapshot)
