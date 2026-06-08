# -*- coding: utf-8 -*-
"""
Human-in-the-loop (HITL) 工具：Label Studio 专家标注白名单网关。

仅允许在预定义生信场景下触发 Label Studio；其它场景一律拒绝。
"""
from __future__ import annotations

import json
import logging
import os
from pathlib import Path
from typing import Any, Dict, List, Optional, Union

from gibh_agent.core.tool_registry import registry
from gibh_agent.core.utils import safe_tool_execution
from gibh_agent.utils.ls_client import (
    DEFAULT_IMAGE_LABEL_CONFIG,
    LabelStudioClient,
    LabelStudioClientError,
)

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# 白名单：仅下列 scenario_type 允许触发 Label Studio
# ---------------------------------------------------------------------------
HITL_SCENARIO_WHITELIST: frozenset[str] = frozenset(
    {
        "scrna_cell_type_annotation",       # 转录组 / 单细胞细胞类型注释
        "spatial_microenvironment",         # 空间组微环境 / 生态位标注
        "spatial_cell_segmentation_qc",   # 空间组细胞分割与 QC 复核
        "histopathology_region_review",   # 病理 / H&E 组织区域专家复核
        "radiomics_roi_validation",       # 影像组学 ROI 边界验证
        "generic_corpus_processing",        # 科学语料数据加工 · 通用图像/文本打标
    }
)

HITL_WHITELIST_DESCRIPTION = (
    "允许场景: "
    + ", ".join(sorted(HITL_SCENARIO_WHITELIST))
)

# 白名单场景 → Label Studio 最简标注界面 XML（创建项目时默认注入，可被 label_config 覆盖）
SCENARIO_TO_LS_XML: Dict[str, str] = {
    "scrna_cell_type_annotation": """<View>
  <Image name="image" value="$image"/>
  <RectangleLabels name="cell_type" toName="image">
    <Label value="T_cell" background="#e6194B"/>
    <Label value="B_cell" background="#3cb44b"/>
    <Label value="Macrophage" background="#ffe119"/>
    <Label value="Unknown" background="#911eb4"/>
  </RectangleLabels>
</View>""",
    "spatial_microenvironment": """<View>
  <Image name="image" value="$image"/>
  <PolygonLabels name="niche" toName="image">
    <Label value="Tumor" background="#e6194B"/>
    <Label value="Stroma" background="#3cb44b"/>
    <Label value="Immune_infiltrate" background="#4363d8"/>
    <Label value="Niche_other" background="#f58231"/>
  </PolygonLabels>
</View>""",
    "spatial_cell_segmentation_qc": """<View>
  <Image name="image" value="$image"/>
  <BrushLabels name="segmentation" toName="image">
    <Label value="Cell" background="#46f0f0"/>
    <Label value="Artifact" background="#f032e6"/>
    <Label value="Background" background="#808080"/>
  </BrushLabels>
</View>""",
    "histopathology_region_review": """<View>
  <Image name="image" value="$image"/>
  <RectangleLabels name="region" toName="image">
    <Label value="Tumor" background="#e6194B"/>
    <Label value="Normal" background="#3cb44b"/>
    <Label value="Necrosis" background="#000000"/>
    <Label value="Margin" background="#ffe119"/>
  </RectangleLabels>
</View>""",
    "radiomics_roi_validation": """<View>
  <Text name="finding" value="$text"/>
  <Labels name="qc" toName="finding">
    <Label value="Accept" background="#3cb44b"/>
    <Label value="Revise" background="#e6194B"/>
    <Label value="Need_expert" background="#ffe119"/>
  </Labels>
</View>""",
    "generic_corpus_processing": """<View>
  <Header value="科学语料标注：框选图像区域并填写 SFT 语料字段"/>
  <Image name="image" value="$image"/>
  <RectangleLabels name="region" toName="image">
    <Label value="Target" background="#e6194B"/>
    <Label value="Context" background="#3cb44b"/>
    <Label value="Ignore" background="#808080"/>
  </RectangleLabels>
  <Text name="text" value="$text"/>
  <Choices name="corpus_role" toName="text" choice="single">
    <Choice value="instruction"/>
    <Choice value="input"/>
    <Choice value="output"/>
    <Choice value="qa_pair"/>
  </Choices>
  <TextArea name="sft_output" toName="text" placeholder="填写规范化 output 文本（SFT 训练标签）" rows="4"/>
</View>""",
}

# 场景 label_config 中 $variable 占位符 → import_task 时 task.data 必填键
SCENARIO_REQUIRED_TASK_DATA_KEYS: Dict[str, tuple[str, ...]] = {
    "generic_corpus_processing": ("image", "text"),
    "radiomics_roi_validation": ("text",),
}

_CORPUS_IMAGE_ONLY_TEXT_PLACEHOLDER = (
    "（图像语料：请用上方工具框选 T-Cell / B-Cell 等区域，并在下方 Choices / TextArea 填写 SFT 语料字段）"
)


def _normalize_tasks_for_scenario(
    scenario: str,
    tasks: List[Dict[str, Any]],
) -> List[Dict[str, Any]]:
    """按场景补齐 LS task.data 必填键，避免 label_config 与导入数据不一致导致 HTTP 400。"""
    required = SCENARIO_REQUIRED_TASK_DATA_KEYS.get(str(scenario or "").strip())
    if not required:
        return tasks
    out: List[Dict[str, Any]] = []
    for raw in tasks:
        if not isinstance(raw, dict):
            continue
        task = dict(raw)
        data = task.get("data") if isinstance(task.get("data"), dict) else None
        if data is None:
            data = {
                k: v
                for k, v in task.items()
                if k not in ("id", "annotations", "completions", "predictions", "meta")
            }
            task = {"data": data}
        else:
            data = dict(data)
        for key in required:
            if key not in data or data[key] is None:
                data[key] = ""
        if scenario == "generic_corpus_processing":
            has_image = bool(str(data.get("image") or "").strip())
            has_text = bool(str(data.get("text") or "").strip())
            if has_image and not has_text:
                data["text"] = _CORPUS_IMAGE_ONLY_TEXT_PLACEHOLDER
            elif has_text and not has_image:
                data["image"] = ""
        task["data"] = data
        out.append(task)
    return out


def resolve_ls_browser_project_url(project_id: int) -> str:
    """供浏览器 iframe / 新窗口打开的 LS 项目 URL（默认同源 /label-studio/ 反代）。"""
    from gibh_agent.utils.ls_public_url import resolve_ls_browser_project_url as _resolve

    return _resolve(project_id)


def normalize_hitl_payload_for_frontend(payload: Dict[str, Any]) -> Dict[str, Any]:
    """统一 hitl_required 载荷字段，确保前端总能拿到 ls_url / project_id / 真实错误信息。"""
    if not isinstance(payload, dict):
        return {}
    out = dict(payload)
    pid = out.get("ls_project_id") or out.get("project_id")
    url = str(out.get("ls_project_url") or out.get("ls_url") or "").strip()
    err_msg = str(out.get("error_detail") or out.get("message") or "").strip()
    if pid is not None and not url:
        url = resolve_ls_browser_project_url(int(pid))
    elif pid is not None and out.get("ls_unavailable"):
        # 保留后端显式不可用标记，但仍提供浏览器 LS 链接供排查
        url = url or resolve_ls_browser_project_url(int(pid))
    elif pid is not None:
        url = resolve_ls_browser_project_url(int(pid))
    out["ls_project_id"] = pid
    out["project_id"] = pid
    out["ls_project_url"] = url
    out["ls_url"] = url
    out.setdefault("status", "hitl_required")
    if err_msg:
        out["error_detail"] = err_msg
        out["message"] = err_msg
    if out.get("ls_unavailable"):
        out["ls_unavailable"] = True
    elif url and pid is not None:
        out.pop("ls_unavailable", None)
    elif out.get("status") == "hitl_required" and not url:
        out["ls_unavailable"] = True
        if not err_msg:
            out["message"] = "Label Studio 项目链接未生成"
            out["error_detail"] = out["message"]
    return out


def normalize_frontend_media_path(path: str) -> str:
    """
    前端 / SSE / state_snapshot 专用：统一为相对路径 /results/... 或 /uploads/...。
    剥离任何绝对 URL（含 https on 8028、192.168.x.x 等毒药地址）。
    """
    raw = str(path or "").strip()
    if not raw:
        return ""
    if raw.startswith("http://") or raw.startswith("https://"):
        from urllib.parse import urlparse

        parsed = urlparse(raw)
        p = (parsed.path or "").strip()
        if p.startswith("/results/") or p.startswith("/uploads/"):
            return p + (parsed.query if parsed.query else "")
        if p.startswith("/app/results/"):
            return "/results/" + p[len("/app/results/") :].lstrip("/")
        if p.startswith("/app/uploads/"):
            return "/uploads/" + p[len("/app/uploads/") :].lstrip("/")
        logger.warning("normalize_frontend_media_path: 无法从绝对 URL 提取相对路径: %s", raw[:200])
        return p or raw
    if raw.startswith("/app/results/"):
        return "/results/" + raw[len("/app/results/") :].lstrip("/")
    if raw.startswith("/app/uploads/"):
        return "/uploads/" + raw[len("/app/uploads/") :].lstrip("/")
    if raw.startswith("results/"):
        return "/" + raw
    if raw.startswith("uploads/"):
        return "/" + raw
    return raw


def _sanitize_internal_http_base(base: str) -> str:
    """LS 容器拉取图片用基址：强制 http://，禁止 https。"""
    b = str(base or "").strip().rstrip("/")
    if not b:
        return ""
    if b.startswith("https://"):
        b = "http://" + b[len("https://") :]
        logger.warning("LS 图片基址已从 https 降级为 http: %s", b)
    if not b.startswith("http://"):
        b = "http://" + b.lstrip("/")
    return b


def _resolve_ls_image_fetch_base() -> str:
    """
    Label Studio 服务端拉取任务图片时使用的 HTTP 基址（内网专用）。
    禁止：外部 IP（192.168.x.x）、https、浏览器侧 OMICS_AGENT_WEB_URL。
    默认：Docker 内 http://nginx 或 http://api-server:8028。
    """
    explicit = os.getenv("LS_IMAGE_FETCH_BASE_URL", "").strip()
    if explicit:
        return _sanitize_internal_http_base(explicit)

    if Path("/.dockerenv").is_file():
        # 同 gibh-network：nginx 静态挂载 /results；api-server 亦挂载 StaticFiles
        for candidate in (
            "http://nginx",
            "http://api-server:8028",
            "http://host.docker.internal:8018",
        ):
            return _sanitize_internal_http_base(candidate)

    return "http://127.0.0.1:8018"


def resolve_ls_import_image_url(path: str) -> str:
    """
    仅供 Label Studio import_task 使用：内网可达的绝对 HTTP URL。
    本系统 /results/ 路径 → http://nginx/results/...；外部公网 URL 原样透传。
    """
    raw = str(path or "").strip()
    if not raw:
        return ""
    if raw.startswith("http://") or raw.startswith("https://"):
        rel = normalize_frontend_media_path(raw)
        if rel.startswith("/results/") or rel.startswith("/uploads/"):
            base = _resolve_ls_image_fetch_base()
            url = base + rel
            logger.info("HITL LS 内网图片 URL(由绝对 URL 规范化): %s → %s", raw[:120], url[:220])
            return url
        logger.info("HITL LS 外部图片 URL 透传: %s", raw[:200])
        return raw

    rel = normalize_frontend_media_path(raw)

    if rel.startswith("/results/") or rel.startswith("/uploads/"):
        base = _resolve_ls_image_fetch_base()
        url = base + rel
        logger.info("HITL LS 内网图片 URL: %s → %s", str(path)[:120], url[:220])
        return url

    p = Path(str(path or rel)).expanduser()
    if p.is_file():
        try:
            rdir = Path(os.getenv("RESULTS_DIR", "/app/results")).expanduser().resolve()
            rel_path = p.resolve().relative_to(rdir)
            base = _resolve_ls_image_fetch_base()
            url = f"{base}/results/{rel_path.as_posix()}"
            logger.info("HITL LS 内网图片 URL(物理路径): %s → %s", str(path)[:120], url[:220])
            return url
        except ValueError:
            pass
        try:
            udir = Path(os.getenv("UPLOAD_DIR", "/app/uploads")).expanduser().resolve()
            rel_path = p.resolve().relative_to(udir)
            base = _resolve_ls_image_fetch_base()
            url = f"{base}/uploads/{rel_path.as_posix()}"
            logger.info("HITL LS 内网图片 URL(上传路径): %s → %s", str(path)[:120], url[:220])
            return url
        except ValueError:
            pass

    logger.warning("resolve_ls_import_image_url: 未能映射 %s", str(path)[:200])
    return ""


def resolve_ls_accessible_image_url(path: str) -> str:
    """向后兼容别名：等同 resolve_ls_import_image_url（仅 LS 导入，非前端）。"""
    return resolve_ls_import_image_url(path)


@registry.register(
    name="Trigger_Expert_Annotation",
    description=(
        "在需要人类视觉判断或领域专家复核时，启动 Label Studio 人机协同标注。"
        "Label Studio 是开源的数据标注与复核平台，用于细胞类型注释、空间微环境划分、"
        "ROI 边界确认等需人工介入的任务。"
        f"仅可在以下白名单场景调用: {HITL_WHITELIST_DESCRIPTION}。"
        "若自动化结果置信度不足、存在歧义簇/边界，或用户明确要求专家复核，方可调用本工具。"
    ),
    category="Human-in-the-loop",
    output_type="json",
)
@safe_tool_execution
def Trigger_Expert_Annotation(
    scenario_type: str,
    project_title: str,
    file_path: str = "",
    image_path: str = "",
    label_config: str = "",
    tasks_json: str = "",
) -> Dict[str, Any]:
    """
    触发 Label Studio 专家标注流程（Human-in-the-loop）。

    Label Studio 是什么？
        Label Studio 是开源的数据标注与复核（Human-in-the-loop）平台，支持图像、
        文本、时序等多模态任务。本系统将其作为「专家介入」出口：当生信 pipeline
        的自动注释/分割/分类结果需要人类视觉或专业判断确认时，由调度器调用本工具
        创建临时标注项目并导入待复核数据。

    何时调用？
        1. scenario_type 必须属于系统白名单（见 HITL_SCENARIO_WHITELIST）。
        2. 当前步骤的结果无法仅靠算法闭环（如低置信度细胞群、空间 niche 边界争议、
           病理区域划分存疑、影像 ROI 需放射科确认）。
        3. 用户或诊断模块明确要求「专家复核 / 人工标注」。

    何时禁止调用？
        - scenario_type 不在白名单（如常规 QC 报表、纯数值统计、无需视觉判断的步骤）。
        - 尚无待标注数据（既无 file_path / image_path，也无 tasks_json）。

    参数:
        scenario_type: 白名单场景标识（如 scrna_cell_type_annotation）。
        project_title: Label Studio 项目名称。
        file_path: 可选，指向 JSON 任务列表文件（每项含 data 或扁平字段）。
        image_path: 可选，单张待标注图像的本地路径或 URL。
        label_config: 可选，Label Studio 标注界面 XML 配置。
        tasks_json: 可选，内联 JSON 字符串（任务数组），与 file_path 二选一。

    返回:
        白名单通过且 LS 项目创建成功时:
        {"status": "hitl_required", "ls_project_url": "...", "ls_project_id": ..., ...}
        白名单拒绝: {"status": "error", "message": "..."}
    """
    scenario = str(scenario_type or "").strip()
    if scenario not in HITL_SCENARIO_WHITELIST:
        return {
            "status": "error",
            "message": (
                f"场景 '{scenario or '(空)'}' 不在 Label Studio 白名单内，已拒绝触发。"
                f" 允许: {HITL_WHITELIST_DESCRIPTION}"
            ),
            "allowed_scenarios": sorted(HITL_SCENARIO_WHITELIST),
        }

    title = str(project_title or "").strip()
    if not title:
        return {"status": "error", "message": "project_title 不能为空"}

    resolved_image = resolve_ls_import_image_url(image_path) if image_path else ""
    if image_path and not resolved_image:
        return {
            "status": "error",
            "message": f"无法将图片路径映射为 Label Studio 内网 HTTP URL: {image_path}",
        }
    tasks = _normalize_tasks_for_scenario(
        scenario,
        _resolve_tasks(
            file_path=file_path,
            image_path=resolved_image or image_path,
            tasks_json=tasks_json,
        ),
    )
    if not tasks:
        return {
            "status": "error",
            "message": "需要提供待标注数据：file_path、image_path 或 tasks_json 至少一项",
        }

    resolved_label_config = (
        label_config.strip()
        or SCENARIO_TO_LS_XML.get(scenario)
        or DEFAULT_IMAGE_LABEL_CONFIG
    )

    client = LabelStudioClient()
    try:
        project = client.create_project(
            title=title,
            label_config=resolved_label_config,
            description=f"GIBH HITL · scenario={scenario}",
        )
        project_id = int(project["id"])
        import_result = client.import_task(project_id, tasks)
        task_count = (
            import_result.get("task_count")
            or import_result.get("import")
            or len(tasks)
        )
        logger.info(
            "Trigger_Expert_Annotation 建项成功 project_id=%s tasks=%s image_sample=%s",
            project_id,
            task_count,
            (tasks[0].get("data") or tasks[0]) if tasks else None,
        )
        ls_url = resolve_ls_browser_project_url(project_id)
    except LabelStudioClientError as exc:
        logger.warning("Trigger_Expert_Annotation LS 失败 scenario=%s: %s", scenario, exc)
        return {"status": "error", "message": str(exc)}

    return normalize_hitl_payload_for_frontend(
        {
            "status": "hitl_required",
            "message": f"已创建 Label Studio 标注项目，等待专家复核（场景: {scenario}）",
            "scenario_type": scenario,
            "ls_project_id": project_id,
            "ls_project_url": ls_url,
            "task_count": len(tasks),
        }
    )


def _resolve_tasks(
    *,
    file_path: str,
    image_path: str,
    tasks_json: str,
) -> List[Dict[str, Any]]:
    """从 file_path / image_path / tasks_json 解析任务列表。"""
    if tasks_json and str(tasks_json).strip():
        try:
            parsed = json.loads(tasks_json)
        except json.JSONDecodeError as exc:
            raise ValueError(f"tasks_json 不是合法 JSON: {exc}") from exc
        if isinstance(parsed, list):
            return [x for x in parsed if isinstance(x, dict)]
        if isinstance(parsed, dict):
            return [parsed]
        raise ValueError("tasks_json 须为对象或数组")

    fp = str(file_path or "").strip()
    if fp:
        path = Path(fp).expanduser()
        if not path.is_file():
            raise ValueError(f"file_path 不存在或不是文件: {path}")
        raw = json.loads(path.read_text(encoding="utf-8"))
        if isinstance(raw, list):
            return [x for x in raw if isinstance(x, dict)]
        if isinstance(raw, dict):
            return [raw]
        raise ValueError("file_path 内 JSON 须为对象或数组")

    img = str(image_path or "").strip()
    if img:
        return [{"data": {"image": img}}]

    return []


def is_hitl_scenario_allowed(scenario_type: str) -> bool:
    """供编排器/诊断模块查询白名单（无需实例化工具）。"""
    return str(scenario_type or "").strip() in HITL_SCENARIO_WHITELIST
