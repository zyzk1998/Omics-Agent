"""
工作流执行「时光机」快照：三大法定键物理隔离，steps_details 为 Live 全量镜像（1:1 入库）。

入库载体为 messages.content（JSON），字段 execution_snapshot 独立存储：
  - data_diagnosis_html：数据诊断/校验富文本（Markdown 或 HTML 字符串）
  - expert_report_markdown：专家结题报告 Markdown（有 HITL 最终版时为最终稿）
  - expert_report_draft_markdown：专家分析报告初稿（HITL 最终版生成后保留用于时光机对比）
  - steps_details：与 Live 工作台一致的完整步骤数组（含诊断节点、管线步骤、日志与图表等）
  - hitl_final / hitl_skipped：HITL 终态标记
  - report_timestamp：报告生成时刻（YYYY-MM-DD HH:MM:SS，入库瞬间冻结，不可变）
  - input_draft_text：主页面输入框未发送自然语言草稿
  - input_draft_attachments：输入框上方待发送/已挂载文件药丸元数据列表
  - corpus_sft_json / corpus_sft_records / corpus_standard / corpus_modality：Task/Skill 终点站生成的 SFT 语料

单会话多次分析：state_snapshot.execution_snapshots = {
  "snapshots": { "<task_id>": { ...完整镜像... }, ... },
  "order": ["<task_id>", ...],
  "active_snapshot_id": "<task_id>"
}
execution_snapshot 始终指向 active_snapshot_id 对应条目（向后兼容）。
"""
from __future__ import annotations

import copy
import json
import logging
import time
from datetime import datetime, timezone
from typing import Any, Dict, List, Optional

from gibh_agent.core.file_asset_meta import normalize_attachment_list
from gibh_agent.core.utils import sanitize_for_json, scrub_markdown_poison_urls

logger = logging.getLogger(__name__)


def _utc_iso_now() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat().replace("+00:00", "Z")


def _frozen_report_timestamp_str() -> str:
    """审计栏展示用本地时刻字符串（入库瞬间冻结，禁止前端 new Date 覆盖）。"""
    return datetime.now().strftime("%Y-%m-%d %H:%M:%S")


def _freeze_report_timestamp(ex: Dict[str, Any]) -> None:
    """报告生成时间仅写入一次；report_timestamp 为法定展示字段。"""
    if not ex.get("report_timestamp"):
        ex["report_timestamp"] = _frozen_report_timestamp_str()
    if not ex.get("report_generated_at"):
        ex["report_generated_at"] = _utc_iso_now()


def ensure_report_timestamp_frozen(state_snapshot: Dict[str, Any]) -> None:
    """流水线结束时确保 execution_snapshot 已写入冻结的 report_timestamp。"""
    if not isinstance(state_snapshot, dict):
        return
    ex = state_snapshot.get("execution_snapshot")
    if not isinstance(ex, dict):
        return
    has_content = bool(
        (ex.get("steps_details") or [])
        or (ex.get("expert_report_markdown") or "").strip()
        or (ex.get("data_diagnosis_html") or "").strip()
    )
    if has_content:
        _freeze_report_timestamp(ex)
        state_snapshot["execution_snapshot"] = ex


def _diagnosis_to_text(diagnosis: Any) -> str:
    """从诊断 payload 提取可持久化正文（Markdown/纯文本）。"""
    if diagnosis is None:
        return ""
    if isinstance(diagnosis, str):
        return diagnosis.strip()
    if isinstance(diagnosis, dict):
        for key in ("message", "diagnosis_report", "report", "summary", "text"):
            val = diagnosis.get(key)
            if isinstance(val, str) and val.strip():
                return val.strip()
    return str(diagnosis).strip()


def normalize_input_draft_attachments(attachments: Any) -> List[Dict[str, Any]]:
    """规范化 Composer 附件队列元数据（含 source: hpc_asset | data_asset）。"""
    return normalize_attachment_list(attachments)


def deep_copy_steps_details(steps_details: Any) -> List[Dict[str, Any]]:
    """深拷贝步骤列表（全量镜像，不做筛选或重组）。"""
    if not isinstance(steps_details, list):
        return []
    out: List[Dict[str, Any]] = []
    for item in steps_details:
        if not isinstance(item, dict):
            continue
        try:
            out.append(copy.deepcopy(item))
        except Exception as exc:  # noqa: BLE001
            logger.warning("execution_snapshot deepcopy step failed: %s", exc)
            out.append(dict(item))
    return out


def _resolve_diagnosis_from_state(state_snapshot: Dict[str, Any]) -> str:
    if not isinstance(state_snapshot, dict):
        return ""
    ex = state_snapshot.get("execution_snapshot")
    if isinstance(ex, dict):
        stored = ex.get("data_diagnosis_html")
        if isinstance(stored, str) and stored.strip():
            return stored.strip()
    report = state_snapshot.get("report")
    if not isinstance(report, dict):
        return ""
    for key in ("diagnosis", "diagnosis_report"):
        text = _diagnosis_to_text(report.get(key))
        if text:
            return text
    rd = report.get("report_data")
    if isinstance(rd, dict):
        text = _diagnosis_to_text(rd.get("diagnosis"))
        if text:
            return text
    return ""


def _resolve_expert_markdown_from_state(state_snapshot: Dict[str, Any]) -> str:
    if not isinstance(state_snapshot, dict):
        return ""
    top = state_snapshot.get("expert_report_markdown")
    if isinstance(top, str) and top.strip():
        return top.strip()
    ex = state_snapshot.get("execution_snapshot")
    if isinstance(ex, dict):
        stored = ex.get("expert_report_markdown")
        if isinstance(stored, str) and stored.strip():
            return stored.strip()
    report = state_snapshot.get("report")
    if not isinstance(report, dict):
        return ""
    rd = report.get("report_data")
    if isinstance(rd, dict):
        for key in ("report", "summary"):
            val = rd.get(key)
            if isinstance(val, str) and val.strip():
                return val.strip()
    for key in ("report", "summary"):
        val = report.get(key)
        if isinstance(val, str) and val.strip():
            return val.strip()
    return ""


def _is_final_expert_report_markdown(markdown: str, *, is_final: bool = False) -> bool:
    if is_final:
        return True
    return "最终版" in (markdown or "")


def build_execution_snapshot(
    steps_details: List[Dict[str, Any]],
    *,
    workflow_name: Optional[str] = None,
    data_diagnosis_html: str = "",
    expert_report_markdown: str = "",
    expert_report_draft_markdown: str = "",
    input_draft_text: str = "",
    input_draft_attachments: Optional[List[Dict[str, Any]]] = None,
    hitl_final: bool = False,
    hitl_skipped: bool = False,
) -> Dict[str, Any]:
    """
    构造执行快照（写入 state_snapshot.execution_snapshot）。
    steps_details 与 Live 侧数组 1:1 对应，禁止入库前过滤或改名。
    """
    steps = deep_copy_steps_details(steps_details)
    diag = (data_diagnosis_html or "").strip()
    expert = (expert_report_markdown or "").strip()
    expert_draft = (expert_report_draft_markdown or "").strip()
    draft_text = (input_draft_text or "").strip() if isinstance(input_draft_text, str) else ""
    draft_att = normalize_input_draft_attachments(input_draft_attachments or [])
    payload: Dict[str, Any] = {
        "version": 5,
        "workflow_name": workflow_name or "",
        "data_diagnosis_html": diag,
        "expert_report_markdown": expert,
        "expert_report_draft_markdown": expert_draft,
        "steps_details": steps,
        "input_draft_text": draft_text,
        "input_draft_attachments": draft_att,
    }
    if hitl_final:
        payload["hitl_final"] = True
    if hitl_skipped:
        payload["hitl_skipped"] = True
    if workflow_name:
        payload["report_data"] = {
            "workflow_name": workflow_name,
            "steps_details": steps,
        }
    return payload


def _ensure_snapshots_collection(state_snapshot: Dict[str, Any]) -> Dict[str, Any]:
    col = state_snapshot.get("execution_snapshots")
    if not isinstance(col, dict):
        col = {"snapshots": {}, "order": [], "active_snapshot_id": None}
        state_snapshot["execution_snapshots"] = col
    if not isinstance(col.get("snapshots"), dict):
        col["snapshots"] = {}
    if not isinstance(col.get("order"), list):
        col["order"] = []
    return col


def resolve_snapshot_id(
    state_snapshot: Dict[str, Any],
    ex: Optional[Dict[str, Any]] = None,
    *,
    explicit_id: Optional[str] = None,
) -> str:
    """解析本次分析快照 ID（execution_id / task_* / 显式传入）。"""
    if explicit_id and str(explicit_id).strip():
        return str(explicit_id).strip()
    if isinstance(state_snapshot, dict):
        cur = state_snapshot.get("_current_execution_id")
        if cur and str(cur).strip():
            return str(cur).strip()
    if isinstance(ex, dict):
        sid = ex.get("snapshot_id") or ex.get("execution_id")
        if sid and str(sid).strip():
            return str(sid).strip()
    if isinstance(state_snapshot, dict):
        wf = state_snapshot.get("workflow")
        if isinstance(wf, dict):
            for key in ("execution_id", "executionId"):
                val = wf.get(key)
                if val and str(val).strip():
                    return str(val).strip()
    return f"task_{int(time.time() * 1000)}"


def register_execution_snapshot_in_collection(
    state_snapshot: Dict[str, Any],
    ex: Dict[str, Any],
    *,
    snapshot_id: Optional[str] = None,
) -> str:
    """将完整 execution_snapshot 写入会话级多快照集（append/覆盖同 ID）。"""
    if not isinstance(state_snapshot, dict) or not isinstance(ex, dict):
        return snapshot_id or ""
    sid = resolve_snapshot_id(state_snapshot, ex, explicit_id=snapshot_id)
    try:
        stored = copy.deepcopy(ex)
    except Exception as exc:  # noqa: BLE001
        logger.warning("execution_snapshots deepcopy failed: %s", exc)
        stored = dict(ex)
    stored["snapshot_id"] = sid
    col = _ensure_snapshots_collection(state_snapshot)
    col["snapshots"][sid] = stored
    if sid not in col["order"]:
        col["order"].append(sid)
    col["active_snapshot_id"] = sid
    state_snapshot["execution_snapshots"] = col
    return sid


def get_execution_snapshot_from_collection(
    state_snapshot: Dict[str, Any],
    snapshot_id: Optional[str] = None,
) -> Optional[Dict[str, Any]]:
    if not isinstance(state_snapshot, dict):
        return None
    col = state_snapshot.get("execution_snapshots")
    if not isinstance(col, dict):
        ex = state_snapshot.get("execution_snapshot")
        return ex if isinstance(ex, dict) else None
    snaps = col.get("snapshots")
    if not isinstance(snaps, dict) or not snaps:
        ex = state_snapshot.get("execution_snapshot")
        return ex if isinstance(ex, dict) else None
    sid = snapshot_id or col.get("active_snapshot_id")
    if sid and sid in snaps and isinstance(snaps[sid], dict):
        return snaps[sid]
    order = col.get("order")
    if isinstance(order, list) and order:
        last = order[-1]
        if last in snaps and isinstance(snaps[last], dict):
            return snaps[last]
    ex = state_snapshot.get("execution_snapshot")
    return ex if isinstance(ex, dict) else None


def _merge_execution_snapshot(
    state_snapshot: Dict[str, Any],
    ex: Dict[str, Any],
    *,
    snapshot_id: Optional[str] = None,
) -> None:
    sid = register_execution_snapshot_in_collection(state_snapshot, ex, snapshot_id=snapshot_id)
    state_snapshot["execution_snapshot"] = ex
    if sid:
        ex["snapshot_id"] = sid
    snap_steps = ex.get("steps_details") or []
    state_snapshot["steps"] = snap_steps

    report = state_snapshot.get("report")
    if not isinstance(report, dict):
        report = {}
    report_data = report.get("report_data")
    if not isinstance(report_data, dict):
        report_data = {}
    legacy_rd = ex.get("report_data")
    if isinstance(legacy_rd, dict):
        report_data = {**report_data, **legacy_rd}
    report_data["steps_details"] = snap_steps
    if ex.get("workflow_name"):
        report_data["workflow_name"] = ex["workflow_name"]
    report["report_data"] = report_data
    state_snapshot["report"] = report


def update_diagnosis_in_state(state_snapshot: Dict[str, Any], diagnosis: Any) -> None:
    """仅更新 data_diagnosis_html，不触碰 steps_details。"""
    if not isinstance(state_snapshot, dict):
        return
    text = _diagnosis_to_text(diagnosis)
    if not text:
        return
    ex = state_snapshot.get("execution_snapshot")
    if not isinstance(ex, dict):
        ex = build_execution_snapshot([], data_diagnosis_html=text)
    else:
        ex = {**ex, "data_diagnosis_html": text}
    _freeze_report_timestamp(ex)
    _merge_execution_snapshot(state_snapshot, ex)

    report = state_snapshot.get("report") or {}
    if not isinstance(report, dict):
        report = {}
    report["diagnosis"] = diagnosis if not isinstance(diagnosis, str) else text
    report["diagnosis_report"] = text
    rd = report.get("report_data") or {}
    if not isinstance(rd, dict):
        rd = {}
    rd["diagnosis"] = text
    report["report_data"] = rd
    state_snapshot["report"] = report


def update_expert_report_in_state(
    state_snapshot: Dict[str, Any],
    markdown: str,
    *,
    is_final: bool = False,
) -> None:
    """更新专家报告 Markdown；最终版时保留初稿至 expert_report_draft_markdown。"""
    if not isinstance(state_snapshot, dict):
        return
    md = scrub_markdown_poison_urls((markdown or "").strip())
    if not md:
        return
    final_report = _is_final_expert_report_markdown(md, is_final=is_final)
    ex = state_snapshot.get("execution_snapshot")
    if not isinstance(ex, dict):
        ex = {}
    current = str(ex.get("expert_report_markdown") or state_snapshot.get("expert_report_markdown") or "").strip()
    draft_stored = str(ex.get("expert_report_draft_markdown") or "").strip()

    if final_report:
        if not draft_stored and current and not _is_final_expert_report_markdown(current):
            ex["expert_report_draft_markdown"] = scrub_markdown_poison_urls(current)
        ex["expert_report_markdown"] = md
        ex["hitl_final"] = True
        state_snapshot["hitl_final"] = True
        state_snapshot["hitl_resumed"] = True
    else:
        ex["expert_report_markdown"] = md
        if not ex.get("hitl_final"):
            ex["expert_report_draft_markdown"] = md

    state_snapshot["expert_report_markdown"] = md
    _freeze_report_timestamp(ex)
    _merge_execution_snapshot(state_snapshot, ex)

    report = state_snapshot.get("report") or {}
    if not isinstance(report, dict):
        report = {}
    rd = report.get("report_data") or {}
    if not isinstance(rd, dict):
        rd = {}
    rd["report"] = md
    if final_report:
        rd["hitl_final"] = True
    report["report_data"] = rd
    report["report"] = md
    state_snapshot["report"] = report


def update_corpus_asset_in_state(
    state_snapshot: Dict[str, Any],
    payload: Dict[str, Any],
) -> None:
    """将 corpus_asset SSE 载荷写入 execution_snapshot（时光机 / MySQL 持久化）。"""
    if not isinstance(state_snapshot, dict) or not isinstance(payload, dict):
        return
    records = payload.get("records")
    if records is None and payload.get("corpus_sft_records") is not None:
        records = payload.get("corpus_sft_records")
    if not isinstance(records, list):
        records = []

    corpus_json = payload.get("corpus_json") or payload.get("corpus_sft_json")
    if not isinstance(corpus_json, str) or not corpus_json.strip():
        try:
            corpus_json = json.dumps(sanitize_for_json(records), ensure_ascii=False, indent=2)
        except (TypeError, ValueError):
            corpus_json = "[]"

    ex = state_snapshot.get("execution_snapshot")
    if not isinstance(ex, dict):
        ex = {}
    ex["corpus_sft_records"] = sanitize_for_json(records)
    ex["corpus_sft_json"] = corpus_json
    ex["corpus_standard"] = payload.get("standard") or payload.get("corpus_standard") or ""
    ex["corpus_modality"] = payload.get("modality") or payload.get("corpus_modality") or ""
    ex["corpus_count"] = int(payload.get("count") or len(records) or 0)
    archive_dir = payload.get("archive_dir") or payload.get("corpus_archive_path")
    if archive_dir:
        ex["corpus_archive_path"] = str(archive_dir)
    hitl_path = payload.get("corpus_hitl_path")
    if hitl_path:
        ex["corpus_hitl_path"] = str(hitl_path)

    _merge_execution_snapshot(state_snapshot, ex)
    state_snapshot["corpus_sft_json"] = corpus_json
    state_snapshot["corpus_sft_records"] = ex["corpus_sft_records"]
    state_snapshot["corpus_standard"] = ex.get("corpus_standard")
    state_snapshot["corpus_modality"] = ex.get("corpus_modality")


def mirror_execution_snapshot_to_message_content(
    content: Dict[str, Any],
    state_snapshot: Dict[str, Any],
) -> None:
    """将 state_snapshot 内 execution_snapshot 双写至 messages.content 根级（时光机入库）。"""
    if not isinstance(content, dict) or not isinstance(state_snapshot, dict):
        return
    ex = state_snapshot.get("execution_snapshot")
    if isinstance(ex, dict):
        content["execution_snapshot"] = sanitize_for_json(ex)
    col = state_snapshot.get("execution_snapshots")
    if isinstance(col, dict) and col.get("snapshots"):
        content["execution_snapshots"] = sanitize_for_json(col)
        active = col.get("active_snapshot_id")
        if active:
            content["active_snapshot_id"] = active


def _resolve_composer_draft_from_state(state_snapshot: Dict[str, Any]) -> tuple:
    if not isinstance(state_snapshot, dict):
        return "", []
    ex = state_snapshot.get("execution_snapshot")
    if isinstance(ex, dict):
        text = ex.get("input_draft_text")
        att = ex.get("input_draft_attachments")
        if isinstance(text, str) or att:
            return (text or "").strip() if isinstance(text, str) else "", normalize_input_draft_attachments(att)
    text_top = state_snapshot.get("input_draft_text")
    att_top = state_snapshot.get("input_draft_attachments")
    if isinstance(text_top, str) or att_top:
        return (
            (text_top or "").strip() if isinstance(text_top, str) else "",
            normalize_input_draft_attachments(att_top),
        )
    return "", []


def update_composer_draft_in_state(
    state_snapshot: Dict[str, Any],
    *,
    input_draft_text: str = "",
    input_draft_attachments: Optional[List[Dict[str, Any]]] = None,
) -> None:
    """仅更新主页面 Composer 草稿字段，不触碰 steps_details / 诊断 / 专家报告。"""
    if not isinstance(state_snapshot, dict):
        return
    draft_text = (input_draft_text or "").strip() if isinstance(input_draft_text, str) else ""
    draft_att = normalize_input_draft_attachments(input_draft_attachments or [])
    ex = state_snapshot.get("execution_snapshot")
    if not isinstance(ex, dict):
        ex = build_execution_snapshot([], input_draft_text=draft_text, input_draft_attachments=draft_att)
    else:
        ex = {
            **ex,
            "version": max(int(ex.get("version") or 2), 3),
            "input_draft_text": draft_text,
            "input_draft_attachments": draft_att,
        }
    _merge_execution_snapshot(state_snapshot, ex)
    state_snapshot["input_draft_text"] = draft_text
    state_snapshot["input_draft_attachments"] = draft_att


def apply_execution_snapshot_to_state(
    state_snapshot: Dict[str, Any],
    steps_details: Any,
    *,
    workflow_name: Optional[str] = None,
) -> None:
    """将 Live steps_details 全量 1:1 合并进 state_snapshot（不做步骤过滤）。"""
    if not isinstance(state_snapshot, dict):
        return
    mirrored_steps = deep_copy_steps_details(steps_details)
    diag_html = _resolve_diagnosis_from_state(state_snapshot)
    expert_md = _resolve_expert_markdown_from_state(state_snapshot)
    draft_text, draft_att = _resolve_composer_draft_from_state(state_snapshot)
    prev = state_snapshot.get("execution_snapshot")
    prev_draft = ""
    prev_hitl_final = False
    prev_hitl_skipped = False
    if isinstance(prev, dict):
        prev_draft = str(prev.get("expert_report_draft_markdown") or "").strip()
        prev_hitl_final = bool(prev.get("hitl_final"))
        prev_hitl_skipped = bool(prev.get("hitl_skipped"))
    ex = build_execution_snapshot(
        mirrored_steps,
        workflow_name=workflow_name,
        data_diagnosis_html=diag_html,
        expert_report_markdown=expert_md,
        expert_report_draft_markdown=prev_draft,
        input_draft_text=draft_text,
        input_draft_attachments=draft_att,
        hitl_final=prev_hitl_final,
        hitl_skipped=prev_hitl_skipped,
    )
    if mirrored_steps or expert_md or diag_html:
        if isinstance(prev, dict):
            if prev.get("report_timestamp"):
                ex["report_timestamp"] = prev["report_timestamp"]
            if prev.get("report_generated_at"):
                ex["report_generated_at"] = prev["report_generated_at"]
            if not draft_text and prev.get("input_draft_text"):
                ex["input_draft_text"] = str(prev.get("input_draft_text") or "").strip()
            if not draft_att and prev.get("input_draft_attachments"):
                ex["input_draft_attachments"] = normalize_input_draft_attachments(
                    prev.get("input_draft_attachments")
                )
            if prev.get("expert_report_draft_markdown"):
                ex["expert_report_draft_markdown"] = prev.get("expert_report_draft_markdown")
            if prev.get("hitl_final"):
                ex["hitl_final"] = True
            if prev.get("hitl_skipped"):
                ex["hitl_skipped"] = True
            if prev.get("hitl_annotations") is not None:
                ex["hitl_annotations"] = prev.get("hitl_annotations")
            for _corpus_key in (
                "corpus_sft_json",
                "corpus_sft_records",
                "corpus_standard",
                "corpus_modality",
                "corpus_count",
                "corpus_archive_path",
                "corpus_hitl_path",
            ):
                if prev.get(_corpus_key) is not None:
                    ex[_corpus_key] = prev.get(_corpus_key)
        if not ex.get("report_timestamp"):
            _freeze_report_timestamp(ex)
    _merge_execution_snapshot(state_snapshot, ex)


def build_session_composer_draft_payload(
    *,
    input_draft_text: str = "",
    input_draft_attachments: Optional[List[Dict[str, Any]]] = None,
) -> Dict[str, Any]:
    """会话级 Composer 草稿（与 execution_snapshot 同键名，供 PATCH /api/sessions/.../composer_draft）。"""
    return {
        "input_draft_text": (input_draft_text or "").strip() if isinstance(input_draft_text, str) else "",
        "input_draft_attachments": normalize_input_draft_attachments(input_draft_attachments or []),
    }


def extract_steps_details_from_agent_content(content: Any) -> List[Dict[str, Any]]:
    """从已入库 agent 消息 content 中解析 steps_details（全量，不筛选）。"""
    if not isinstance(content, dict):
        return []

    def _from_ex(ex: Any) -> List[Dict[str, Any]]:
        if not isinstance(ex, dict):
            return []
        sd = ex.get("steps_details")
        if isinstance(sd, list) and sd:
            return deep_copy_steps_details(sd)
        return []

    ex = content.get("execution_snapshot")
    steps = _from_ex(ex)
    if steps:
        return steps

    col = content.get("execution_snapshots")
    if isinstance(col, dict):
        snaps = col.get("snapshots")
        if isinstance(snaps, dict) and snaps:
            order = col.get("order")
            if isinstance(order, list) and order:
                for sid in reversed(order):
                    steps = _from_ex(snaps.get(sid))
                    if steps:
                        return steps
            for _sid, snap_ex in snaps.items():
                steps = _from_ex(snap_ex)
                if steps:
                    return steps

    snap = content.get("state_snapshot")
    if isinstance(snap, dict):
        steps = _from_ex(snap.get("execution_snapshot"))
        if steps:
            return steps
        col2 = snap.get("execution_snapshots")
        if isinstance(col2, dict):
            snaps2 = col2.get("snapshots")
            if isinstance(snaps2, dict) and snaps2:
                order2 = col2.get("order")
                if isinstance(order2, list) and order2:
                    for sid2 in reversed(order2):
                        steps = _from_ex(snaps2.get(sid2))
                        if steps:
                            return steps
        if isinstance(snap.get("steps"), list) and snap["steps"]:
            return deep_copy_steps_details(snap["steps"])
        rep = snap.get("report")
        if isinstance(rep, dict):
            rd = rep.get("report_data")
            if isinstance(rd, dict) and isinstance(rd.get("steps_details"), list):
                return deep_copy_steps_details(rd["steps_details"])
    return []
