# -*- coding: utf-8 -*-
"""
科学语料数据加工专用 Agent：程序化触发 Label Studio HITL，唤醒后导出 SFT JSON 语料。
"""
from __future__ import annotations

import json
import logging
import os
import re
import time
import uuid
from pathlib import Path
from typing import Any, AsyncIterator, Callable, Dict, List, Optional

from gibh_agent.agents.base_agent import BaseAgent
from gibh_agent.core.prompt_manager import PromptManager
from gibh_agent.core.llm_client import LLMClient
from gibh_agent.core.utils import sanitize_for_json
from gibh_agent.tools.hitl_tools import (
    Trigger_Expert_Annotation,
    normalize_frontend_media_path,
    normalize_hitl_payload_for_frontend,
)

logger = logging.getLogger(__name__)

_IMAGE_SUFFIXES = frozenset({
    ".png", ".jpg", ".jpeg", ".jpe", ".jfif", ".gif", ".bmp", ".webp",
    ".tif", ".tiff", ".svg", ".ico", ".heic", ".heif", ".avif",
    ".ppm", ".pgm", ".pbm",
})
_TEXT_SUFFIXES = frozenset({".txt", ".md", ".csv", ".tsv"})
_JSON_SUFFIXES = frozenset({".json"})
_SCENARIO = "generic_corpus_processing"
_SKILL_ID = "skill_corpus_data_processing"
_WORKFLOW_NAME = "科学语料数据加工"

_SFT_SYSTEM_PROMPT = """你是一个科研语料整理专家。
请将用户在 Label Studio 中标注的原始 JSON 导出结果，清洗并转换为标准的 SFT 微调语料格式。

输出要求：
1. 仅输出一个 JSON 数组，不要 Markdown 围栏或解释文字。
2. 每条记录必须包含字段：instruction、input、output（均为字符串，可为空字符串但键必须存在）。
3. 从矩形框/多边形/Choices/TextArea 等标注中归纳语义：instruction 为任务指令，input 为上下文或问题，output 为期望回答或标签文本。
4. 若无法从某条任务推断完整三元组，可合理留空 input 或将 instruction 设为通用描述。
5. 合并同一任务内多个标注为一条语料，避免重复条目。"""


def _extract_json_block(text: str) -> Dict[str, Any]:
    if not (text or "").strip():
        return {}
    m = re.search(r"```(?:json)?\s*\n(\{.*?\})\s*```", text, re.DOTALL | re.IGNORECASE)
    if not m:
        return {}
    try:
        obj = json.loads(m.group(1))
        return obj if isinstance(obj, dict) else {}
    except json.JSONDecodeError:
        return {}


class CorpusProcessingAgent(BaseAgent):
    """科学语料数据加工：硬 HITL 挂起 + 唤醒后 SFT JSON 导出。"""

    def __init__(
        self,
        llm_client: Optional[LLMClient],
        prompt_manager: PromptManager,
    ):
        super().__init__(llm_client, prompt_manager, "corpus_processing_expert")

    async def process_query(
        self,
        query: str,
        history: List[Dict[str, str]] = None,
        uploaded_files: List[Dict[str, str]] = None,
        **kwargs,
    ) -> Dict[str, Any]:
        return {
            "type": "skill",
            "skill_id": _SKILL_ID,
            "message": "请通过技能广场「科学语料数据加工」上传文件以启动 Label Studio 标注流程。",
        }

    def resolve_upload_media(
        self,
        file_paths: List[str],
        user_query: str = "",
    ) -> Dict[str, Any]:
        """从上传列表与 prompt JSON 块解析 image_path / file_path / 内联文本。"""
        tpl = _extract_json_block(user_query)
        image_path = str(tpl.get("image_path") or "").strip()
        file_path = str(tpl.get("file_path") or "").strip()
        inline_text = ""

        for raw in file_paths or []:
            p = str(raw or "").strip()
            if not p:
                continue
            rel = normalize_frontend_media_path(p) or p
            suffix = Path(rel.split("?")[0]).suffix.lower()
            if suffix in _IMAGE_SUFFIXES and not image_path:
                image_path = rel
            elif suffix in _JSON_SUFFIXES and not file_path:
                file_path = rel
            elif suffix in _TEXT_SUFFIXES and not file_path:
                file_path = rel
                try:
                    physical = self._physical_path(rel)
                    if physical.is_file():
                        inline_text = physical.read_text(encoding="utf-8", errors="replace")[:8000]
                except OSError:
                    pass
            elif not image_path and not file_path:
                file_path = rel

        return {
            "image_path": image_path,
            "file_path": file_path,
            "inline_text": inline_text,
            "project_title": str(tpl.get("project_title") or "科学语料专家标注").strip(),
        }

    @staticmethod
    def _physical_path(rel_or_abs: str) -> Path:
        raw = str(rel_or_abs or "").strip()
        if raw.startswith("/results/"):
            base = Path(os.getenv("RESULTS_DIR", "/app/results"))
            return base / raw[len("/results/") :].lstrip("/")
        if raw.startswith("/uploads/"):
            base = Path(os.getenv("UPLOAD_DIR", "/app/uploads"))
            return base / raw[len("/uploads/") :].lstrip("/")
        if raw.startswith("/app/results/"):
            return Path(raw)
        if raw.startswith("/app/uploads/"):
            return Path(raw)
        p = Path(raw).expanduser()
        if p.is_file():
            return p
        for root in (
            Path(os.getenv("UPLOAD_DIR", "/app/uploads")),
            Path(os.getenv("RESULTS_DIR", "/app/results")),
        ):
            cand = root / raw.lstrip("/")
            if cand.is_file():
                return cand
        return p

    def trigger_corpus_hitl(
        self,
        *,
        image_path: str = "",
        file_path: str = "",
        inline_text: str = "",
        project_title: str = "",
    ) -> Dict[str, Any]:
        """程序化调用 Trigger_Expert_Annotation，无需 LLM 决策。"""
        title = (project_title or "科学语料专家标注").strip()
        tasks_json = ""
        if inline_text and not image_path and not file_path:
            tasks_json = json.dumps([{"data": {"text": inline_text}}], ensure_ascii=False)

        result = Trigger_Expert_Annotation(
            scenario_type=_SCENARIO,
            project_title=title,
            image_path=image_path or "",
            file_path=file_path or "",
            tasks_json=tasks_json,
        )
        if isinstance(result, dict):
            return normalize_hitl_payload_for_frontend(result)
        return {"status": "error", "message": "Trigger_Expert_Annotation 返回异常"}

    async def stream_skill_flow(
        self,
        *,
        state_snapshot: Dict[str, Any],
        emit_sse: Callable[[str, Dict[str, Any]], str],
        file_paths: List[str],
        user_query: str,
        session_id: str,
        db: Any,
        owner_id: str,
        model_name: str,
    ) -> AsyncIterator[str]:
        """技能快车道：收到有效文件后立即挂起进入 Label Studio（硬 HITL）。"""
        self.set_request_model_name(model_name)
        yield emit_sse("status", {"content": "正在准备语料标注任务…", "state": "running"})

        media = self.resolve_upload_media(file_paths, user_query)
        if not any([media.get("image_path"), media.get("file_path"), media.get("inline_text")]):
            yield emit_sse(
                "message",
                {
                    "content": (
                        "请先上传待标注的**图像**（png/jpg 等）或**文本/JSON 任务文件**，"
                        "再点击「使用」启动科学语料数据加工。"
                    ),
                },
            )
            yield emit_sse("done", {"status": "error"})
            return

        yield emit_sse(
            "status",
            {"content": "正在创建 Label Studio 标注项目（无需大模型推理）…", "state": "running"},
        )
        hitl_result = self.trigger_corpus_hitl(
            image_path=media.get("image_path") or "",
            file_path=media.get("file_path") or "",
            inline_text=media.get("inline_text") or "",
            project_title=media.get("project_title") or "",
        )

        if hitl_result.get("status") != "hitl_required":
            err = hitl_result.get("message") or "创建 Label Studio 项目失败"
            yield emit_sse("message", {"content": f"**语料标注启动失败**：{err}"})
            yield emit_sse("done", {"status": "error"})
            return

        sid = str(session_id or "").strip()
        if db and sid:
            from gibh_agent.core.session_runtime import SESSION_WAITING_FOR_HITL, set_session_status

            set_session_status(db, sid, SESSION_WAITING_FOR_HITL, owner_id=owner_id)

        _hitl_event = {
            "ls_url": hitl_result.get("ls_project_url") or hitl_result.get("ls_url") or "",
            "project_id": hitl_result.get("ls_project_id") or hitl_result.get("project_id"),
            "scenario_type": _SCENARIO,
            "message": hitl_result.get("message") or "请在 Label Studio 中完成语料标注",
            "workflow_name": _WORKFLOW_NAME,
        }
        state_snapshot["hitl_pending"] = True
        state_snapshot["hitl"] = _hitl_event
        state_snapshot["workflow"] = {"workflow_name": _WORKFLOW_NAME, "skill_id": _SKILL_ID}

        if sid:
            from gibh_agent.core.hitl_session_registry import register_hitl_session

            register_hitl_session(
                sid,
                {
                    "project_id": _hitl_event.get("project_id"),
                    "ls_url": _hitl_event.get("ls_url"),
                    "scenario_type": _SCENARIO,
                    "workflow_name": _WORKFLOW_NAME,
                    "agent_key": "corpus_processing_agent",
                },
            )

        yield emit_sse("message", {"content": "已创建 Label Studio 语料标注项目，请在右侧进入专家工作台完成打标。"})
        yield emit_sse("hitl_action", _hitl_event)
        yield emit_sse(
            "status",
            {"content": "⏸ 流程已挂起，等待您在 Label Studio 中完成语料标注", "state": "waiting"},
        )
        yield emit_sse("done", {"status": "waiting_for_hitl", "tool_name": _SKILL_ID, "tool_result": hitl_result})

    async def _generate_analysis_summary(
        self,
        steps_results: List[Dict[str, Any]],
        omics_type: str = "Corpus",
        workflow_name: str = _WORKFLOW_NAME,
        summary_context: Optional[Dict[str, Any]] = None,
        output_dir: Optional[str] = None,
    ) -> Optional[str]:
        """HITL 唤醒：将 LS 标注导出清洗为 SFT JSON 并写入工作区。"""
        ctx = summary_context or {}
        ann_json = str(ctx.get("hitl_annotations_json") or "").strip()
        if not ann_json:
            return (
                "## 科学语料导出（最终版）\n\n"
                "未能获取 Label Studio 标注导出，请确认已在 LS 中提交标注后重试唤醒。"
            )

        user_content = (
            "【Label Studio 原始标注 JSON】\n"
            f"{ann_json[:12000]}\n\n"
            "请输出清洗后的 SFT 语料 JSON 数组。"
        )
        messages = [
            {"role": "system", "content": _SFT_SYSTEM_PROMPT},
            {"role": "user", "content": user_content},
        ]
        sft_records: List[Dict[str, str]] = []
        raw_llm = ""
        try:
            completion = await self.llm_for_request().achat(messages, temperature=0.1, max_tokens=4096)
            if completion and getattr(completion, "choices", None):
                _llm = self.llm_for_request()
                if hasattr(_llm, "extract_think_and_content"):
                    _, raw_llm = _llm.extract_think_and_content(completion)
                else:
                    block = completion.choices[0]
                    msg = getattr(block, "message", None)
                    raw_llm = getattr(msg, "content", None) if msg else getattr(block, "content", "")
            sft_records = self._parse_sft_json(raw_llm)
        except Exception as exc:
            logger.exception("CorpusProcessingAgent SFT 清洗 LLM 失败: %s", exc)

        if not sft_records:
            sft_records = self._fallback_sft_from_ls_export(ann_json)

        out_dir = self._resolve_corpus_output_dir(output_dir)
        out_dir.mkdir(parents=True, exist_ok=True)
        ts = time.strftime("%Y%m%d_%H%M%S")
        fname = f"sft_corpus_{ts}_{uuid.uuid4().hex[:8]}.json"
        out_path = out_dir / fname
        payload = sanitize_for_json(sft_records)
        out_path.write_text(json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8")

        try:
            rel = out_path.resolve().relative_to(
                Path(os.getenv("RESULTS_DIR", "/app/results")).expanduser().resolve()
            )
            json_url = f"/results/{rel.as_posix()}"
        except ValueError:
            json_url = f"/results/corpus_hitl/{fname}"

        preview = json.dumps(payload[:3], ensure_ascii=False, indent=2)
        if len(payload) > 3:
            preview += f"\n... （共 {len(payload)} 条，完整内容见下载链接）"

        return (
            "## 科学语料导出（最终版）\n\n"
            f"已将 Label Studio 专家标注清洗为 **{len(payload)}** 条 SFT 微调语料。\n\n"
            f"### 下载\n\n"
            f"- [SFT 语料 JSON]({json_url})（`instruction` / `input` / `output`）\n\n"
            f"### 预览（前 3 条）\n\n"
            f"```json\n{preview}\n```\n\n"
            "### 说明\n\n"
            "语料已保存至工作区，可直接用于大模型 SFT/RLHF 训练流水线。"
        )

    @staticmethod
    def _parse_sft_json(text: str) -> List[Dict[str, str]]:
        raw = (text or "").strip()
        if not raw:
            return []
        if raw.startswith("```"):
            raw = re.sub(r"^```(?:json)?\s*", "", raw, flags=re.IGNORECASE)
            raw = re.sub(r"\s*```$", "", raw)
        try:
            parsed = json.loads(raw)
        except json.JSONDecodeError:
            m = re.search(r"\[[\s\S]*\]", raw)
            if not m:
                return []
            try:
                parsed = json.loads(m.group(0))
            except json.JSONDecodeError:
                return []
        if isinstance(parsed, dict):
            parsed = [parsed]
        if not isinstance(parsed, list):
            return []
        out: List[Dict[str, str]] = []
        for item in parsed:
            if not isinstance(item, dict):
                continue
            out.append(
                {
                    "instruction": str(item.get("instruction") or "").strip(),
                    "input": str(item.get("input") or "").strip(),
                    "output": str(item.get("output") or "").strip(),
                }
            )
        return [x for x in out if any(x.values())]

    @staticmethod
    def _fallback_sft_from_ls_export(ann_json: str) -> List[Dict[str, str]]:
        """LLM 不可用时的最小结构化兜底。"""
        try:
            data = json.loads(ann_json)
        except json.JSONDecodeError:
            return []
        tasks = data if isinstance(data, list) else [data]
        out: List[Dict[str, str]] = []
        for task in tasks:
            if not isinstance(task, dict):
                continue
            data_block = task.get("data") if isinstance(task.get("data"), dict) else {}
            anns = task.get("annotations") or task.get("completions") or []
            labels: List[str] = []
            text_out = ""
            if isinstance(anns, list):
                for ann in anns:
                    if not isinstance(ann, dict):
                        continue
                    for res in ann.get("result") or []:
                        if not isinstance(res, dict):
                            continue
                        val = res.get("value") or {}
                        if isinstance(val, dict):
                            if val.get("text"):
                                text_out = str(val["text"][0] if isinstance(val["text"], list) else val["text"])
                            choices = val.get("choices") or []
                            if choices:
                                labels.extend(str(c) for c in choices)
            instruction = "请根据科研语料完成下列标注任务"
            inp = str(data_block.get("text") or data_block.get("image") or "")[:2000]
            output = text_out or ", ".join(labels) or json.dumps(data_block, ensure_ascii=False)[:500]
            out.append({"instruction": instruction, "input": inp, "output": output})
        return out

    @staticmethod
    def _resolve_corpus_output_dir(output_dir: Optional[str]) -> Path:
        if output_dir:
            p = Path(output_dir).expanduser()
            if p.is_dir():
                return p
        base = Path(os.getenv("RESULTS_DIR", "/app/results")).expanduser() / "corpus_hitl"
        base.mkdir(parents=True, exist_ok=True)
        return base
