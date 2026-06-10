# -*- coding: utf-8 -*-
"""
科学语料数据加工专用 Agent：程序化触发 Label Studio HITL，唤醒后导出 SFT JSON 语料。
"""
from __future__ import annotations

import base64
import json
import logging
import os
import re
import shutil
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
    parse_image_path_inputs,
    resolve_ls_import_image_payload,
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

_VLM_SFT_SYSTEM_PROMPT = """你是一名多模态数据科学家（Vision-Language Model 语料工程师），专精 LLaVA / Qwen-VL 等视觉指令微调数据构建。

你将收到 Label Studio 专家标注的结构化摘要（每条任务含 task_id、regions 坐标与标签、可选文本字段）。
请将其炼制为**多模态 VLM 指令微调语料** JSON 数组。

## 输出格式（严格遵守，仅输出 JSON 数组，禁止 Markdown 围栏与解释文字）

每条记录结构：
{
  "id": "corpus_task_01",
  "image": "<保留占位符 __IMAGE_PLACEHOLDER__，不要自行编造 Base64>",
  "conversations": [
    {
      "from": "human",
      "value": "<image>\\n请详细描述图像中被标注的各个细胞群/感兴趣区域的位置及其类别。"
    },
    {
      "from": "gpt",
      "value": "在这张图像中，我检测到了以下关键区域：\\n1. 在图像左上方（坐标：[ymin, xmin, ymax, xmax]），专家标注为 [T-Cell] 细胞群。\\n2. ..."
    }
  ]
}

## 坐标与语义铁律
1. 输入中 `bbox_norm_1000` 已是 **[ymin, xmin, ymax, xmax]**，取值范围 0–1000；输出 gpt 回答中必须原样引用或基于此描述空间位置（如「左上方」「中央偏右」+ 坐标）。
2. 将生物学/病理/空间组学术语写入 gpt 回答：细胞群、微环境、ROI、组织区域等；禁止只罗列干巴巴标签名。
3. human 的 value **必须以** `<image>\\n` 开头，提出与标注场景匹配的视觉理解问题（中文）。
4. 同一 LS 任务内多个框选合并为**一条**语料；gpt 回答用有序列表逐条描述各区域。
5. `id` 与输入 task_id 一致；`image` 字段固定写 `__IMAGE_PLACEHOLDER__`（系统会注入真实 Base64）。
6. 纯文本任务（无 regions 仅有 text_fields）可省略坐标，但仍须 conversations 双轮结构。
7. 禁止输出 instruction/input/output 旧格式；禁止空 conversations。"""

_AI_ONLY_SFT_SYSTEM_PROMPT = """你是一名多模态数据科学家。当前任务**没有经过人类专家手动打标**（用户已跳过 Label Studio 复核），
请基于系统提供的 **AI 初步分析报告** 与图像信息，构建 LLaVA / Qwen-VL 格式的多模态指令微调语料。

## 输出格式（严格遵守，仅输出 JSON 数组，禁止 Markdown 围栏与解释文字）

每条记录结构：
{
  "id": "corpus_task_01",
  "image": "<保留占位符 __IMAGE_PLACEHOLDER__，不要自行编造 Base64>",
  "conversations": [
    {
      "from": "human",
      "value": "<image>\\n请对这张科研图像给出全面的生物学解读，包括主要细胞群、空间分布与关键生物学意义。"
    },
    {
      "from": "gpt",
      "value": "基于 AI 初步分析报告（未经专家手工框选），对图像的全局描述如下：\\n..."
    }
  ]
}

## 代偿生成铁律（无 Bounding Box 时）
1. 输入中 `regions` 通常为空；**禁止编造虚假坐标框**。请使用 **全局图像描述（Global Description）** 构建 gpt 回答。
2. 将 `text_fields.analysis_draft` 中的聚类解读、标志基因、生物学结论融入 gpt 回答，使其像专家级视觉解读。
3. human 的 value **必须以** `<image>\\n` 开头；gpt 回答须完整、可独立用于 SFT，禁止「见上文」「同上」。
4. 每张输入图像对应**一条**语料；`id` 与输入 task_id 一致；`image` 固定写 `__IMAGE_PLACEHOLDER__`。
5. 禁止输出 instruction/input/output 旧格式；禁止空 conversations。"""

_NLP_SFT_SYSTEM_PROMPT = """你是一个科研数据整理专家。请将以下任务输入和输出转化为标准对话指令微调语料。

## 输出格式（严格遵守，仅输出 JSON 数组，禁止 Markdown 围栏与解释文字）

[
  {
    "instruction": "请根据以下科研分析上下文，给出专业、可复用的结论摘要。",
    "input": "<任务输入或上下文>",
    "output": "<高质量输出，融合专家报告要点>"
  }
]

## 铁律
1. 本任务**无图像资产**，禁止出现任何关于图像、坐标、Bounding Box、<image> 的描述。
2. instruction / input / output 均使用中文科研表述，output 须完整可独立用于 SFT。
3. 可拆分为多条语料，但每条须自洽；禁止空字段。"""

_IMAGE_PATH_PLACEHOLDER = "__IMAGE_PATH__"

_PLOT_IMAGE_KEYS = (
    "plot",
    "plot_path",
    "image_path",
    "umap_plot",
    "figure_path",
    "output_image",
    "image",
    "heatmap_path",
    "cluster_plot",
    "visualization_path",
)


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
        tpl_file = str(tpl.get("file_path") or "").strip()
        tpl_images = parse_image_path_inputs(
            tpl.get("image_paths") or tpl.get("image_path") or ""
        )
        image_paths: List[str] = []
        file_path = ""
        inline_text = ""

        for raw in file_paths or []:
            p = str(raw or "").strip()
            if not p:
                continue
            rel = normalize_frontend_media_path(p) or p
            suffix = Path(rel.split("?")[0]).suffix.lower()
            if suffix in _IMAGE_SUFFIXES:
                image_paths.append(rel)
            elif suffix in _JSON_SUFFIXES:
                file_path = rel
            elif suffix in _TEXT_SUFFIXES:
                file_path = rel
                try:
                    physical = self._physical_path(rel)
                    if physical.is_file():
                        inline_text = physical.read_text(encoding="utf-8", errors="replace")[:8000]
                except OSError:
                    pass
            elif not image_paths and not file_path:
                file_path = rel

        if not image_paths:
            image_paths = tpl_images
        if not file_path:
            file_path = tpl_file

        return {
            "image_paths": image_paths,
            "image_path": image_paths[0] if image_paths else "",
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
        image_paths: Optional[List[str]] = None,
        file_path: str = "",
        inline_text: str = "",
        project_title: str = "",
    ) -> Dict[str, Any]:
        """程序化调用 Trigger_Expert_Annotation，无需 LLM 决策。"""
        title = (project_title or "科学语料专家标注").strip()
        tasks_json = ""
        paths = list(image_paths or []) or parse_image_path_inputs(image_path)
        if inline_text and not paths and not file_path:
            tasks_json = json.dumps([{"data": {"text": inline_text}}], ensure_ascii=False)

        result = Trigger_Expert_Annotation(
            scenario_type=_SCENARIO,
            project_title=title,
            image_path=paths if len(paths) > 1 else (paths[0] if paths else ""),
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
        if not any([media.get("image_paths"), media.get("file_path"), media.get("inline_text")]):
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
            image_paths=media.get("image_paths") or [],
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
        """HITL 唤醒：记录专家复核完成；标准语料在【一键入库】时由 Corpus Gatekeeper 生成。"""
        ctx = summary_context or {}
        has_annotations = ctx.get("hitl_annotations_raw") is not None or bool(
            str(ctx.get("hitl_annotations_json") or "").strip()
        )
        if has_annotations:
            return (
                "## 科学语料复核（已完成）\n\n"
                "Label Studio 专家标注已合并至会话快照。\n\n"
                "标准 **VLM / NLP 语料归档**（`corpus_archive/dataset.json` + 解耦图像路径）"
                "将在您点击 **【一键入库 / 沉淀业务资产】** 时由入库守门人自动生成并打包。"
            )
        return (
            "## 科学语料导出\n\n"
            "当前无专家标注导出。标准语料将在 **一键入库** 时基于 AI 初稿按模态（有图→VLM，无图→NLP）自动生成。"
        )

    async def generate_corpus_archive_bundle(
        self,
        *,
        summary_context: Optional[Dict[str, Any]] = None,
        steps_results: Optional[List[Dict[str, Any]]] = None,
        steps_details: Optional[List[Dict[str, Any]]] = None,
        archive_dir: Path,
        expert_report_md: str = "",
        output_dir: Optional[str] = None,
    ) -> Dict[str, Any]:
        """
        Corpus-Ready 主入口：模态路由 + 资产解耦归档。

        产出目录结构::
            corpus_archive/
              dataset.json
              images/      (VLM)
              reports/expert_report.md
        """
        from gibh_agent.core.corpus_gatekeeper import (
            CORPUS_MODALITY_NLP,
            CORPUS_MODALITY_VLM,
            sniff_corpus_modality,
        )

        ctx = summary_context or {}
        ann_json = str(ctx.get("hitl_annotations_json") or "").strip()
        has_annotations = ctx.get("hitl_annotations_raw") is not None or bool(ann_json)
        hitl_images = list(ctx.get("hitl_image_paths") or [])

        modality = sniff_corpus_modality(
            steps_details=steps_details,
            steps_results=steps_results,
            hitl_annotations=ctx.get("hitl_annotations_raw") if has_annotations else None,
            hitl_image_paths=hitl_images,
        )

        archive_dir = Path(archive_dir).expanduser()
        archive_dir.mkdir(parents=True, exist_ok=True)
        images_dir = archive_dir / "images"
        reports_dir = archive_dir / "reports"
        reports_dir.mkdir(parents=True, exist_ok=True)

        report_md = (expert_report_md or ctx.get("expert_report_markdown") or ctx.get("expert_report_draft") or "").strip()
        if report_md:
            (reports_dir / "expert_report.md").write_text(report_md, encoding="utf-8")

        if modality == CORPUS_MODALITY_VLM:
            records, meta = await self._generate_vlm_slim_records(
                ctx=ctx,
                steps_results=steps_results,
                steps_details=steps_details,
                images_dir=images_dir,
                has_annotations=has_annotations,
                ann_json=ann_json,
            )
        else:
            records, meta = await self._generate_nlp_slim_records(
                ctx=ctx,
                steps_details=steps_details,
                steps_results=steps_results,
            )

        if not records:
            return {
                "status": "error",
                "message": meta.get("message") or "语料记录为空，无法归档",
                "archive_dir": str(archive_dir.resolve()),
                "modality": modality,
                "count": 0,
                "dataset_path": None,
            }

        dataset_path = archive_dir / "dataset.json"
        payload = sanitize_for_json(records)
        dataset_path.write_text(json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8")

        bundle = {
            "status": "success",
            "records": payload,
            "archive_dir": str(archive_dir.resolve()),
            "dataset_path": str(dataset_path.resolve()),
            "modality": modality,
            "count": len(payload),
            "ai_only": meta.get("ai_only", not has_annotations),
            "message": meta.get("message", ""),
        }
        self._last_corpus_archive_bundle = bundle
        return bundle

    async def generate_vlm_corpus_bundle(
        self,
        *,
        summary_context: Optional[Dict[str, Any]] = None,
        steps_results: Optional[List[Dict[str, Any]]] = None,
        output_dir: Optional[str] = None,
        steps_details: Optional[List[Dict[str, Any]]] = None,
        session_id: str = "legacy",
    ) -> Dict[str, Any]:
        """兼容旧调用方：委托 generate_corpus_archive_bundle（路径引用、无 Base64）。"""
        archive_dir = (
            Path(os.getenv("RESULTS_DIR", "/app/results")).expanduser()
            / "corpus_archive"
            / str(session_id or "legacy")
        )
        bundle = await self.generate_corpus_archive_bundle(
            summary_context=summary_context,
            steps_results=steps_results,
            steps_details=steps_details,
            archive_dir=archive_dir,
            expert_report_md=str((summary_context or {}).get("expert_report_markdown") or ""),
            output_dir=output_dir,
        )
        if bundle.get("status") == "success" and bundle.get("dataset_path"):
            try:
                rel = Path(bundle["dataset_path"]).resolve().relative_to(
                    Path(os.getenv("RESULTS_DIR", "/app/results")).expanduser().resolve()
                )
                json_url = f"/results/{rel.as_posix()}"
            except ValueError:
                json_url = None
            bundle["path"] = bundle["dataset_path"]
            bundle["json_url"] = json_url
        self._last_vlm_corpus_bundle = bundle
        return bundle

    @staticmethod
    def format_vlm_corpus_markdown(bundle: Dict[str, Any]) -> str:
        if bundle.get("status") != "success":
            title = "AI 代偿版" if bundle.get("ai_only") else "最终版"
            return f"## 科学语料导出（{title}）\n\n{bundle.get('message') or '语料生成失败。'}"
        payload = bundle.get("records") or []
        json_url = bundle.get("json_url") or ""
        preview = json.dumps(payload[:3], ensure_ascii=False, indent=2)
        if len(payload) > 3:
            preview += f"\n... （共 {len(payload)} 条，完整内容见下载链接）"
        title = "AI 代偿版" if bundle.get("ai_only") else "最终版"
        source = (
            "已将 AI 初步分析报告与步骤图像代偿清洗为"
            if bundle.get("ai_only")
            else "已将 Label Studio 专家标注清洗为"
        )
        return (
            f"## 科学语料导出（{title}）\n\n"
            f"{source} **{len(payload)}** 条 VLM 多模态微调语料。\n\n"
            f"### 下载\n\n"
            f"- [VLM 语料 JSON]({json_url})（`id` / `image` / `conversations` · LLaVA 风格）\n\n"
            f"### 预览（前 3 条）\n\n"
            f"```json\n{preview}\n```\n\n"
            "### 说明\n\n"
            "语料已保存至工作区 `corpus_archive/`，`dataset.json` 内 image 为相对路径引用。"
        )

    @staticmethod
    def format_corpus_archive_markdown(bundle: Dict[str, Any]) -> str:
        if bundle.get("status") != "success":
            return f"## 语料归档\n\n{bundle.get('message') or '语料归档失败。'}"
        modality = bundle.get("modality") or "unknown"
        label = "VLM 多模态" if modality == "vlm" else "NLP 纯文本"
        return (
            f"## 语料归档（{label}）\n\n"
            f"已生成 **{bundle.get('count', 0)}** 条语料，归档于 `corpus_archive/dataset.json`。\n\n"
            f"- 模态：`{modality}`\n"
            f"- 路径：`{bundle.get('archive_dir', '')}`\n"
            f"- 图像引用：相对路径 `images/...`（无内嵌 Base64）"
        )

    async def _generate_vlm_slim_records(
        self,
        *,
        ctx: Dict[str, Any],
        steps_results: Optional[List[Dict[str, Any]]],
        steps_details: Optional[List[Dict[str, Any]]],
        images_dir: Path,
        has_annotations: bool,
        ann_json: str,
    ) -> tuple[List[Dict[str, Any]], Dict[str, Any]]:
        if has_annotations:
            structured_tasks = self._extract_structured_ls_tasks(
                ctx.get("hitl_annotations_raw"),
                ann_json,
            )
            system_prompt = _VLM_SFT_SYSTEM_PROMPT
            user_payload = self._build_vlm_llm_user_payload(structured_tasks)
            fallback_fn = self._fallback_vlm_corpus_from_structured
            ai_only = False
        else:
            draft = str(
                ctx.get("expert_report_draft")
                or ctx.get("draft_expert_report")
                or ctx.get("expert_report_markdown")
                or ""
            ).strip()
            image_paths = self.collect_resolvable_image_paths(steps_details, steps_results)
            for raw_path in ctx.get("hitl_image_paths") or ctx.get("image_paths") or []:
                for p in parse_image_path_inputs(raw_path):
                    if p not in image_paths:
                        image_paths.append(p)
            structured_tasks = self._build_ai_only_structured_tasks(draft, image_paths, steps_details)
            system_prompt = _AI_ONLY_SFT_SYSTEM_PROMPT
            user_payload = self._build_ai_only_llm_user_payload(structured_tasks)
            fallback_fn = self._fallback_ai_only_vlm_from_structured
            ai_only = True

        if not structured_tasks:
            err = (
                "未能从 AI 初稿与步骤结果中解析可用图像，无法生成 VLM 语料。"
                if ai_only
                else "未能获取 Label Studio 标注导出。"
            )
            return [], {"message": err, "ai_only": ai_only}

        self._materialize_structured_task_images(structured_tasks, images_dir)

        messages = [
            {"role": "system", "content": system_prompt},
            {"role": "user", "content": user_payload},
        ]
        vlm_records: List[Dict[str, Any]] = []
        if self.llm_client is not None:
            raw_llm = ""
            try:
                completion = await self.llm_for_request().achat(messages, temperature=0.1, max_tokens=8192)
                if completion and getattr(completion, "choices", None):
                    _llm = self.llm_for_request()
                    if hasattr(_llm, "extract_think_and_content"):
                        _, raw_llm = _llm.extract_think_and_content(completion)
                    else:
                        block = completion.choices[0]
                        msg = getattr(block, "message", None)
                        raw_llm = getattr(msg, "content", None) if msg else getattr(block, "content", "")
                vlm_records = self._parse_vlm_corpus_json(raw_llm)
                vlm_records = self._inject_image_paths_into_vlm_records(vlm_records, structured_tasks)
            except Exception as exc:
                logger.exception("VLM slim corpus LLM failed: %s", exc)

        if not vlm_records:
            vlm_records = self._inject_image_paths_into_vlm_records(
                fallback_fn(structured_tasks),
                structured_tasks,
            )
        return vlm_records, {"ai_only": ai_only}

    async def _generate_nlp_slim_records(
        self,
        *,
        ctx: Dict[str, Any],
        steps_details: Optional[List[Dict[str, Any]]],
        steps_results: Optional[List[Dict[str, Any]]],
    ) -> tuple[List[Dict[str, Any]], Dict[str, Any]]:
        draft = str(
            ctx.get("expert_report_draft")
            or ctx.get("draft_expert_report")
            or ctx.get("expert_report_markdown")
            or ""
        ).strip()
        steps_lines: List[str] = []
        for sd in steps_details or []:
            if not isinstance(sd, dict):
                continue
            name = sd.get("step_name") or sd.get("name") or sd.get("tool_name") or "步骤"
            status = sd.get("status") or "unknown"
            steps_lines.append(f"- {name}: {status}")
        user_payload = (
            "【科研任务上下文 · 纯文本模态】\n"
            f"### 分析报告\n{draft[:8000] or '（无报告正文）'}\n\n"
            f"### 步骤摘要\n{chr(10).join(steps_lines) or '（无步骤）'}\n\n"
            "请输出标准 instruction/input/output JSON 数组；禁止任何图像相关描述。"
        )
        messages = [
            {"role": "system", "content": _NLP_SFT_SYSTEM_PROMPT},
            {"role": "user", "content": user_payload},
        ]
        records: List[Dict[str, Any]] = []
        if self.llm_client is not None:
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
                records = self._parse_nlp_corpus_json(raw_llm)
            except Exception as exc:
                logger.exception("NLP slim corpus LLM failed: %s", exc)

        if not records:
            records = self._fallback_nlp_corpus_from_context(draft, steps_lines)
        return records, {"ai_only": True, "message": ""}

    @classmethod
    def collect_resolvable_image_paths(
        cls,
        steps_details: Optional[List[Dict[str, Any]]],
        steps_results: Optional[List[Dict[str, Any]]] = None,
    ) -> List[str]:
        """收集步骤中可解析为文件的图像路径（不含纯文本任务）。"""
        seen: set = set()
        paths: List[str] = []

        def _add(raw: Any) -> None:
            p = str(raw or "").strip()
            if not p or p in seen:
                return
            if p.startswith("data:image") or cls._resolve_image_fs_path(p) is not None:
                seen.add(p)
                paths.append(p)

        for sd in steps_details or []:
            if not isinstance(sd, dict):
                continue
            for key in _PLOT_IMAGE_KEYS:
                _add(sd.get(key))
            sr = sd.get("step_result") if isinstance(sd.get("step_result"), dict) else sd
            if isinstance(sr, dict):
                for key in _PLOT_IMAGE_KEYS:
                    _add(sr.get(key))
        for sr in steps_results or []:
            if isinstance(sr, dict):
                for key in _PLOT_IMAGE_KEYS:
                    _add(sr.get(key))
        return paths

    @classmethod
    def annotations_contain_images(cls, annotations_raw: Any) -> bool:
        tasks = annotations_raw if isinstance(annotations_raw, list) else [annotations_raw]
        for task in tasks:
            if not isinstance(task, dict):
                continue
            data = task.get("data") if isinstance(task.get("data"), dict) else {}
            image = str(data.get("image") or "").strip()
            if image:
                return True
        return False

    @classmethod
    def _resolve_image_fs_path(cls, raw_path: str) -> Optional[Path]:
        raw = str(raw_path or "").strip()
        if not raw or raw.startswith("data:"):
            return None
        rel = normalize_frontend_media_path(raw)
        bases = [
            Path(os.getenv("UPLOAD_DIR", "/app/uploads")).expanduser(),
            Path(os.getenv("RESULTS_DIR", "/app/results")).expanduser(),
            Path("/app/assets").expanduser(),
        ]
        for base in bases:
            for candidate in (base / rel.lstrip("/"), base / raw.lstrip("/")):
                if candidate.is_file():
                    return candidate.resolve()
        p = Path(raw).expanduser()
        return p.resolve() if p.is_file() else None

    @classmethod
    def _materialize_image_asset(
        cls,
        src_path: str,
        images_dir: Path,
        stem: str,
    ) -> Optional[str]:
        images_dir.mkdir(parents=True, exist_ok=True)
        raw = str(src_path or "").strip()
        safe_stem = re.sub(r"[^\w.-]+", "_", stem) or "image"

        if raw.startswith("data:image"):
            m = re.match(r"data:image/([\w+.-]+);base64,(.+)", raw, re.DOTALL)
            if not m:
                return None
            ext = m.group(1).replace("jpeg", "jpg").split("+")[0] or "png"
            b64_raw = m.group(2).strip()
            try:
                data = base64.b64decode(b64_raw, validate=False)
            except (ValueError, TypeError):
                try:
                    pad = b64_raw + "=" * (-len(b64_raw) % 4)
                    data = base64.b64decode(pad)
                except (ValueError, TypeError):
                    return None
            fname = f"{safe_stem}.{ext}"
            dest = images_dir / fname
            dest.write_bytes(data)
            return f"images/{fname}"

        fs_path = cls._resolve_image_fs_path(raw)
        if fs_path and fs_path.is_file():
            ext = fs_path.suffix or ".png"
            fname = f"{safe_stem}{ext}"
            dest = images_dir / fname
            if not dest.exists():
                shutil.copy2(fs_path, dest)
            return f"images/{fname}"
        return None

    @classmethod
    def _materialize_structured_task_images(
        cls,
        structured_tasks: List[Dict[str, Any]],
        images_dir: Path,
    ) -> None:
        for task in structured_tasks:
            task_id = str(task.get("task_id") or "corpus_task_01")
            src = str(task.get("image_source_path") or task.get("image") or "").strip()
            if not src:
                continue
            rel = cls._materialize_image_asset(src, images_dir, task_id)
            if rel:
                task["image_rel_path"] = rel

    @classmethod
    def _inject_image_paths_into_vlm_records(
        cls,
        records: List[Dict[str, Any]],
        structured_tasks: List[Dict[str, Any]],
    ) -> List[Dict[str, Any]]:
        by_id = {str(t.get("task_id")): t for t in structured_tasks}
        out: List[Dict[str, Any]] = []
        for idx, rec in enumerate(records):
            rid = str(rec.get("id") or "").strip()
            src = by_id.get(rid) or (structured_tasks[idx] if idx < len(structured_tasks) else None)
            image_rel = ""
            if src:
                image_rel = str(src.get("image_rel_path") or "").strip()
                if not image_rel:
                    raw_img = str(rec.get("image") or "").strip()
                    if raw_img and not raw_img.startswith("data:") and raw_img != "__IMAGE_PLACEHOLDER__":
                        image_rel = raw_img
            if not image_rel:
                continue
            out.append(
                {
                    "id": rid or (src.get("task_id") if src else f"corpus_task_{idx + 1:02d}"),
                    "image": image_rel,
                    "conversations": rec.get("conversations") or [],
                }
            )
        return out

    @staticmethod
    def _parse_nlp_corpus_json(text: str) -> List[Dict[str, str]]:
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
            row = {
                "instruction": str(item.get("instruction") or "").strip(),
                "input": str(item.get("input") or "").strip(),
                "output": str(item.get("output") or "").strip(),
            }
            if any(row.values()):
                out.append(row)
        return out

    @staticmethod
    def _fallback_nlp_corpus_from_context(
        draft: str,
        steps_lines: List[str],
    ) -> List[Dict[str, str]]:
        inp = draft[:3000] if draft else "\n".join(steps_lines[:20])
        out_text = draft[3000:6000] if len(draft) > 3000 else (draft or "（基于步骤摘要生成的科研结论占位）")
        return [
            {
                "instruction": "请根据以下科研分析上下文，给出专业、可复用的结论摘要。",
                "input": inp or "（无输入上下文）",
                "output": out_text,
            }
        ]

    @staticmethod
    def _xywh_percent_to_norm_1000(x: float, y: float, w: float, h: float) -> List[int]:
        """LS 百分比框 (x,y,width,height) → [ymin,xmin,ymax,xmax] 0–1000。"""
        ymin = int(round(y * 10))
        xmin = int(round(x * 10))
        ymax = int(round((y + h) * 10))
        xmax = int(round((x + w) * 10))
        return [min(1000, max(0, v)) for v in (ymin, xmin, ymax, xmax)]

    @staticmethod
    def _short_image_ref(image: str) -> str:
        raw = str(image or "").strip()
        if not raw:
            return ""
        if raw.startswith("data:") and len(raw) > 80:
            return raw[:48] + f"...<base64 len={len(raw)}>"
        return raw[:200]

    @classmethod
    def _extract_region_from_ls_result(cls, res: Dict[str, Any]) -> Optional[Dict[str, Any]]:
        if not isinstance(res, dict):
            return None
        val = res.get("value") or {}
        if not isinstance(val, dict):
            return None
        labels = (
            val.get("rectanglelabels")
            or val.get("polygonlabels")
            or val.get("brushlabels")
            or val.get("labels")
            or []
        )
        label = str(labels[0]) if labels else "Unknown"
        x, y, w, h = val.get("x"), val.get("y"), val.get("width"), val.get("height")
        if all(v is not None for v in (x, y, w, h)):
            xf, yf, wf, hf = float(x), float(y), float(w), float(h)
            return {
                "label": label,
                "bbox_norm_1000": cls._xywh_percent_to_norm_1000(xf, yf, wf, hf),
                "xywh_percent": [xf, yf, wf, hf],
                "type": res.get("type") or "rectanglelabels",
            }
        if val.get("text"):
            text_val = val["text"][0] if isinstance(val["text"], list) else val["text"]
            return {"label": label, "text": str(text_val), "type": res.get("type") or "textarea"}
        choices = val.get("choices") or []
        if choices:
            return {"label": label, "choices": [str(c) for c in choices], "type": res.get("type") or "choices"}
        return {"label": label, "type": res.get("type") or "unknown"}

    @classmethod
    def _extract_structured_ls_tasks(
        cls,
        annotations_raw: Any,
        ann_json: str,
    ) -> List[Dict[str, Any]]:
        """从 LS 导出（对象或 JSON 字符串）提取任务级结构化摘要。"""
        tasks_data: Any = None
        if annotations_raw is not None:
            tasks_data = annotations_raw
        elif ann_json:
            try:
                tasks_data = json.loads(ann_json)
            except json.JSONDecodeError:
                return []
        if tasks_data is None:
            return []
        tasks = tasks_data if isinstance(tasks_data, list) else [tasks_data]
        structured: List[Dict[str, Any]] = []
        for idx, task in enumerate(tasks):
            if not isinstance(task, dict):
                continue
            data = task.get("data") if isinstance(task.get("data"), dict) else {}
            image_source = str(data.get("image") or "").strip()
            task_id = str(task.get("id") or f"corpus_task_{idx + 1:02d}")
            regions: List[Dict[str, Any]] = []
            text_fields: Dict[str, str] = {}
            for ann in task.get("annotations") or task.get("completions") or []:
                if not isinstance(ann, dict):
                    continue
                for res in ann.get("result") or []:
                    region = cls._extract_region_from_ls_result(res)
                    if region and region.get("bbox_norm_1000"):
                        regions.append(region)
                    elif region and region.get("text"):
                        text_fields.setdefault("annotation_text", str(region["text"]))
                    elif region and region.get("choices"):
                        text_fields.setdefault("choices", ", ".join(region["choices"]))
            if data.get("text"):
                text_fields.setdefault("task_text", str(data["text"]))
            if not image_source:
                continue
            structured.append(
                {
                    "task_id": task_id,
                    "image_source_path": image_source,
                    "image_ref": cls._short_image_ref(image_source),
                    "regions": regions,
                    "text_fields": text_fields,
                }
            )
        return structured

    @staticmethod
    def _build_vlm_llm_user_payload(structured_tasks: List[Dict[str, Any]]) -> str:
        slim = []
        for t in structured_tasks:
            slim.append(
                {
                    "task_id": t.get("task_id"),
                    "image_ref": t.get("image_ref"),
                    "regions": t.get("regions") or [],
                    "text_fields": t.get("text_fields") or {},
                }
            )
        return (
            "【Label Studio 结构化标注摘要】\n"
            f"{json.dumps(slim, ensure_ascii=False, indent=2)}\n\n"
            "请输出符合规范的 VLM 语料 JSON 数组；image 字段一律写 __IMAGE_PLACEHOLDER__。"
        )

    @classmethod
    def _collect_plot_image_paths_from_steps(
        cls,
        steps_details: Optional[List[Dict[str, Any]]],
        steps_results: Optional[List[Dict[str, Any]]] = None,
    ) -> List[str]:
        """从步骤详情/结果中提取可视化图像路径（供 Skip 代偿语料使用）。"""
        seen: set = set()
        paths: List[str] = []

        def _add(raw: Any) -> None:
            p = str(raw or "").strip()
            if not p or p in seen:
                return
            low = p.split("?")[0].lower()
            if p.startswith("data:image") or any(low.endswith(suf) for suf in _IMAGE_SUFFIXES):
                seen.add(p)
                paths.append(p)

        for sd in steps_details or []:
            if not isinstance(sd, dict):
                continue
            for key in _PLOT_IMAGE_KEYS:
                _add(sd.get(key))
            sr = sd.get("step_result") if isinstance(sd.get("step_result"), dict) else sd
            if isinstance(sr, dict):
                for key in _PLOT_IMAGE_KEYS:
                    _add(sr.get(key))
                for val in sr.values():
                    if isinstance(val, str):
                        _add(val)
        for sr in steps_results or []:
            if not isinstance(sr, dict):
                continue
            for key in _PLOT_IMAGE_KEYS:
                _add(sr.get(key))
        return paths

    @classmethod
    def _build_ai_only_structured_tasks(
        cls,
        draft_report: str,
        image_paths: List[str],
        steps_details: Optional[List[Dict[str, Any]]] = None,
    ) -> List[Dict[str, Any]]:
        draft = str(draft_report or "").strip()
        if not draft and steps_details:
            lines: List[str] = []
            for sd in steps_details[:12]:
                if not isinstance(sd, dict):
                    continue
                name = sd.get("step_name") or sd.get("name") or sd.get("tool_name") or "步骤"
                status = sd.get("status") or "unknown"
                lines.append(f"- {name}: {status}")
            if lines:
                draft = "AI Pipeline 执行摘要：\n" + "\n".join(lines)

        structured: List[Dict[str, Any]] = []
        for idx, img_path in enumerate(image_paths or []):
            p = str(img_path or "").strip()
            if not p:
                continue
            if not p.startswith("data:image") and cls._resolve_image_fs_path(p) is None:
                continue
            structured.append(
                {
                    "task_id": f"corpus_task_{idx + 1:02d}",
                    "image_source_path": p,
                    "image_ref": cls._short_image_ref(p),
                    "regions": [],
                    "text_fields": {
                        "analysis_draft": draft[:8000],
                        "task_text": draft[:2000],
                    },
                    "ai_only": True,
                }
            )
        return structured

    @staticmethod
    def _build_ai_only_llm_user_payload(structured_tasks: List[Dict[str, Any]]) -> str:
        slim = []
        for t in structured_tasks:
            slim.append(
                {
                    "task_id": t.get("task_id"),
                    "image_ref": t.get("image_ref"),
                    "regions": t.get("regions") or [],
                    "text_fields": t.get("text_fields") or {},
                    "ai_only": True,
                }
            )
        return (
            "【AI 初稿分析报告与图像上下文（无专家手工框选）】\n"
            f"{json.dumps(slim, ensure_ascii=False, indent=2)}\n\n"
            "请基于 AI 初步分析报告对每张图像构建全局描述型 VLM 语料；"
            "无 bbox 时在 gpt 回答中使用全局图像描述（Global Description）。"
            "image 字段一律写 __IMAGE_PLACEHOLDER__。"
        )

    @classmethod
    def _fallback_ai_only_vlm_from_structured(
        cls,
        structured_tasks: List[Dict[str, Any]],
    ) -> List[Dict[str, Any]]:
        """Skip / 无标注时 LLM 不可用：用 AI 初稿构建全局描述型 VLM 语料。"""
        out: List[Dict[str, Any]] = []
        for task in structured_tasks:
            task_id = str(task.get("task_id") or "corpus_task_01")
            image = str(task.get("image_rel_path") or "").strip()
            if not image:
                continue
            draft = str((task.get("text_fields") or {}).get("analysis_draft") or "").strip()
            human_q = (
                "<image>\n请对这张科研图像给出全面的生物学解读，"
                "包括主要细胞群分布、空间结构特征与关键生物学意义。"
            )
            if draft:
                gpt_body = (
                    "基于 AI 初步分析报告（未经专家手工框选），对图像的全局描述如下：\n"
                    f"{draft[:3500]}"
                )
            else:
                gpt_body = (
                    "基于 AI 初步分析报告（未经专家手工框选），对图像的全局描述如下：\n"
                    "图像展示了本次组学分析的关键可视化结果，请结合 pipeline 上下文理解主要细胞群与空间分布特征。"
                )
            out.append(
                {
                    "id": task_id,
                    "image": image,
                    "conversations": [
                        {"from": "human", "value": human_q},
                        {"from": "gpt", "value": gpt_body},
                    ],
                }
            )
        return out

    @staticmethod
    def _parse_vlm_corpus_json(text: str) -> List[Dict[str, Any]]:
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
        out: List[Dict[str, Any]] = []
        for item in parsed:
            if not isinstance(item, dict):
                continue
            conv = item.get("conversations")
            if not isinstance(conv, list) or len(conv) < 2:
                continue
            out.append(
                {
                    "id": str(item.get("id") or "").strip() or None,
                    "image": str(item.get("image") or "").strip(),
                    "conversations": conv,
                }
            )
        return [x for x in out if x.get("conversations")]

    @staticmethod
    def _inject_images_into_vlm_records(
        records: List[Dict[str, Any]],
        structured_tasks: List[Dict[str, Any]],
    ) -> List[Dict[str, Any]]:
        by_id = {str(t.get("task_id")): t for t in structured_tasks}
        out: List[Dict[str, Any]] = []
        for idx, rec in enumerate(records):
            rid = str(rec.get("id") or "").strip()
            src = by_id.get(rid) or (structured_tasks[idx] if idx < len(structured_tasks) else None)
            image = ""
            if src:
                image = str(src.get("image") or "").strip()
            elif rec.get("image") and str(rec["image"]) != "__IMAGE_PLACEHOLDER__":
                image = str(rec["image"])
            if not image:
                continue
            out.append(
                {
                    "id": rid or (src.get("task_id") if src else f"corpus_task_{idx + 1:02d}"),
                    "image": image,
                    "conversations": rec.get("conversations") or [],
                }
            )
        return out

    @classmethod
    def _fallback_vlm_corpus_from_structured(
        cls,
        structured_tasks: List[Dict[str, Any]],
    ) -> List[Dict[str, Any]]:
        """LLM 不可用时的确定性 VLM 语料兜底（含规范化坐标与 Base64 图像）。"""
        out: List[Dict[str, Any]] = []
        for task in structured_tasks:
            task_id = str(task.get("task_id") or "corpus_task_01")
            image = str(task.get("image_rel_path") or "").strip()
            regions = task.get("regions") or []
            lines: List[str] = []
            for i, reg in enumerate(regions, start=1):
                bbox = reg.get("bbox_norm_1000") or []
                label = reg.get("label") or "Unknown"
                if len(bbox) == 4:
                    pos = cls._describe_bbox_position(bbox)
                    lines.append(
                        f"{i}. 在图像{pos}（坐标：[{bbox[0]}, {bbox[1]}, {bbox[2]}, {bbox[3]}]），"
                        f"专家标注为 [{label}] 感兴趣区域。"
                    )
                else:
                    lines.append(f"{i}. 专家标注为 [{label}]。")
            gpt_body = (
                "在这张图像中，我检测到了以下关键区域：\n" + "\n".join(lines)
                if lines
                else "图像中暂无矩形框标注，请结合任务文本字段理解专家意图。"
            )
            text_hint = (task.get("text_fields") or {}).get("task_text", "")
            human_q = (
                "<image>\n请详细描述图像中被标注的各个细胞群/感兴趣区域的位置及其类别。"
            )
            if text_hint:
                human_q += f"\n（任务上下文：{text_hint[:500]}）"
            if not image:
                continue
            out.append(
                {
                    "id": task_id,
                    "image": image,
                    "conversations": [
                        {"from": "human", "value": human_q},
                        {"from": "gpt", "value": gpt_body},
                    ],
                }
            )
        return out

    @staticmethod
    def _describe_bbox_position(bbox: List[int]) -> str:
        if len(bbox) != 4:
            return "中"
        ymin, xmin, ymax, xmax = bbox
        cy, cx = (ymin + ymax) / 2, (xmin + xmax) / 2
        v = "上方" if cy < 333 else ("下方" if cy > 666 else "中部")
        h = "左侧" if cx < 333 else ("右侧" if cx > 666 else "中央")
        return f"{v}{h}"

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
