# -*- coding: utf-8 -*-
"""
Human-in-the-loop (HITL) 工具：Label Studio 专家标注白名单网关。

仅允许在预定义生信场景下触发 Label Studio；其它场景一律拒绝。
"""
from __future__ import annotations

import base64
import json
import logging
import mimetypes
import os
from pathlib import Path
import re
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


def get_image_base64_data_uri(file_path: str) -> str:
    """将本地图像文件编码为 Data URI，供 Label Studio import_task 内嵌渲染（绕过 CSP / 内网 DNS）。"""
    path = Path(str(file_path or "").strip()).expanduser()
    if not path.is_file():
        raise FileNotFoundError(f"找不到待标注图片: {file_path}")
    mime_type, _ = mimetypes.guess_type(str(path))
    mime_type = mime_type or "image/png"
    encoded = base64.b64encode(path.read_bytes()).decode("utf-8")
    return f"data:{mime_type};base64,{encoded}"


def _results_dir() -> Path:
    return Path(os.getenv("RESULTS_DIR", "/app/results")).expanduser().resolve()


def _uploads_dir() -> Path:
    return Path(os.getenv("UPLOAD_DIR", "/app/uploads")).expanduser().resolve()


def _static_assets_dir() -> Path:
    """Nginx 静态资源根（含技能广场演示图）。"""
    explicit = os.getenv("NGINX_HTML_ROOT", "").strip()
    if explicit:
        return Path(explicit).expanduser().resolve()
    return Path(__file__).resolve().parents[2] / "services" / "nginx" / "html"


def _find_uploaded_image_by_basename(basename: str) -> Optional[Path]:
    """在 uploads 卷内按文件名查找（用户会话路径常为 /uploads/{owner}/{ts}/file.png）。"""
    name = str(basename or "").strip()
    if not name:
        return None
    uploads = _uploads_dir()
    if not uploads.is_dir():
        return None
    try:
        for candidate in uploads.rglob(name):
            if candidate.is_file():
                return candidate.resolve()
    except OSError:
        return None
    return None


def resolve_local_image_file(path: str) -> Optional[Path]:
    """将 /uploads/、/results/、物理路径或内网 HTTP URL 解析为容器内可读图像文件。"""
    raw = str(path or "").strip()
    if not raw:
        return None
    if raw.startswith("data:"):
        return None

    if raw.startswith("http://") or raw.startswith("https://"):
        rel = normalize_frontend_media_path(raw)
        if not (rel.startswith("/results/") or rel.startswith("/uploads/")):
            return None
        raw = rel

    rel = normalize_frontend_media_path(raw)
    if rel.startswith("/results/"):
        candidate = _results_dir() / rel[len("/results/") :].lstrip("/")
        return candidate if candidate.is_file() else None
    if rel.startswith("/assets/"):
        candidate = _static_assets_dir() / rel.lstrip("/")
        return candidate if candidate.is_file() else None
    if rel.startswith("/uploads/"):
        candidate = _uploads_dir() / rel[len("/uploads/") :].lstrip("/")
        if candidate.is_file():
            return candidate
        basename = Path(rel).name
        demo = _static_assets_dir() / "images" / "demos" / "corpus" / basename
        if demo.is_file():
            return demo.resolve()
        return _find_uploaded_image_by_basename(basename)

    p = Path(raw).expanduser()
    if p.is_file():
        return p.resolve()
    if p.is_absolute():
        return p if p.is_file() else None

    for base in (_uploads_dir(), _results_dir(), Path.cwd()):
        candidate = (base / raw).resolve()
        if candidate.is_file():
            return candidate
    return None


def resolve_ls_import_image_payload(path: str) -> str:
    """
    Label Studio import_task 图像字段：本地文件 → Base64 Data URI；
    已是 data: URI 原样返回；无法解析本地文件时返回空串。
    """
    raw = str(path or "").strip()
    if not raw:
        return ""
    if raw.startswith("data:"):
        return raw

    local_file = resolve_local_image_file(raw)
    if local_file is not None:
        data_uri = get_image_base64_data_uri(str(local_file))
        logger.info(
            "HITL LS 图片 Base64 内嵌: %s → data:%s;base64,<len=%s>",
            str(path)[:120],
            data_uri[5 : data_uri.find(";")],
            len(data_uri),
        )
        return data_uri

    if raw.startswith("http://") or raw.startswith("https://"):
        rel = normalize_frontend_media_path(raw)
        if rel.startswith("/results/") or rel.startswith("/uploads/"):
            logger.warning(
                "resolve_ls_import_image_payload: 内网/本系统 URL 无本地文件，拒绝 HTTP 透传 %s",
                str(path)[:200],
            )
            return ""
        logger.info("HITL LS 外部图片 URL 透传(非本系统路径): %s", raw[:200])
        return raw

    logger.warning("resolve_ls_import_image_payload: 未能解析本地图片 %s", str(path)[:200])
    return ""


def _is_internal_ls_image_ref(value: str) -> bool:
    raw = str(value or "").strip()
    if not raw or raw.startswith("data:"):
        return False
    if raw.startswith("/uploads/") or raw.startswith("/results/") or raw.startswith("/assets/"):
        return True
    if raw.startswith("http://") or raw.startswith("https://"):
        rel = normalize_frontend_media_path(raw)
        return rel.startswith("/uploads/") or rel.startswith("/results/") or rel.startswith("/assets/")
    return False


def _embed_task_image_fields_as_base64(tasks: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
    """将任务 data.image 中的路径 / 内网 URL 转为 Base64 Data URI。"""
    out: List[Dict[str, Any]] = []
    for task in tasks:
        item = dict(task)
        data = item.get("data")
        if isinstance(data, dict) and data.get("image"):
            raw_img = str(data["image"])
            embedded = resolve_ls_import_image_payload(raw_img)
            data = dict(data)
            if embedded:
                data["image"] = embedded
            elif _is_internal_ls_image_ref(raw_img):
                data["image"] = ""
            item["data"] = data
        elif isinstance(item, dict) and item.get("image") and "data" not in item:
            raw_img = str(item["image"])
            embedded = resolve_ls_import_image_payload(raw_img)
            item = dict(item)
            if embedded:
                item["image"] = embedded
            elif _is_internal_ls_image_ref(raw_img):
                item["image"] = ""
        out.append(item)
    return out


def _validate_tasks_have_embedded_images(tasks: List[Dict[str, Any]]) -> Optional[str]:
    """导入 LS 前校验：凡含 image 字段的任务须已为 data: URI 或留空（纯文本语料）。"""
    for idx, task in enumerate(tasks):
        if not isinstance(task, dict):
            continue
        data = task.get("data") if isinstance(task.get("data"), dict) else task
        if not isinstance(data, dict):
            continue
        img = str(data.get("image") or "").strip()
        if not img:
            continue
        if img.startswith("data:"):
            continue
        if _is_internal_ls_image_ref(img) or "nginx" in img:
            return (
                f"任务 #{idx + 1} 图像未能内嵌为 Base64（拒绝 http://nginx 等浏览器不可达 URL）: "
                f"{img[:120]}"
            )
    return None


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
    Label Studio import_task 图像字段（向后兼容函数名）。
    本地 /uploads/、/results/ 及可解析物理路径 → Base64 Data URI。
    """
    return resolve_ls_import_image_payload(path)


def resolve_ls_accessible_image_url(path: str) -> str:
    """向后兼容别名：等同 resolve_ls_import_image_payload（LS 导入专用）。"""
    return resolve_ls_import_image_payload(path)


def parse_image_path_inputs(image_path: Union[str, List[str], None]) -> List[str]:
    """
    将 image_path 解析为有序路径列表。
    支持：List[str]、JSON 数组字符串、逗号/换行分隔的多路径字符串、单路径、完整 data: URI。
    注意：禁止按分号切分（会破坏 data:image/png;base64,...）。
    """
    if image_path is None:
        return []
    if isinstance(image_path, list):
        out: List[str] = []
        for item in image_path:
            out.extend(parse_image_path_inputs(item))
        return [p for p in out if p]
    raw = str(image_path).strip()
    if not raw:
        return []
    if raw.startswith("data:"):
        return [raw]
    if raw.startswith("[") and raw.endswith("]"):
        try:
            parsed = json.loads(raw)
            if isinstance(parsed, list):
                return parse_image_path_inputs(parsed)
        except json.JSONDecodeError:
            pass
    parts = re.split(r"[,\n]+", raw)
    return [p.strip() for p in parts if p and p.strip()]


def resolve_image_paths_to_base64_payloads(image_paths: List[str]) -> tuple[List[str], Optional[str]]:
    """批量将路径解析为 Base64 Data URI；失败时返回 (已解析列表, 错误信息)。"""
    resolved: List[str] = []
    for path in image_paths:
        payload = resolve_ls_import_image_payload(path)
        if not payload:
            return resolved, f"无法读取待标注图片并编码为 Base64: {path}"
        resolved.append(payload)
    return resolved, None


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
    image_path: Union[str, List[str]] = "",
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
        image_path: 可选，单张或多张待标注图像路径（List[str]、逗号分隔字符串或 JSON 数组字符串）。
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

    raw_image_paths = parse_image_path_inputs(image_path)
    resolved_images, resolve_err = resolve_image_paths_to_base64_payloads(raw_image_paths)
    if resolve_err:
        return {"status": "error", "message": resolve_err}
    tasks = _embed_task_image_fields_as_base64(
        _normalize_tasks_for_scenario(
            scenario,
            _resolve_tasks(
                file_path=file_path,
                image_paths=resolved_images,
                tasks_json=tasks_json,
            ),
        )
    )
    embed_err = _validate_tasks_have_embedded_images(tasks)
    if embed_err:
        return {"status": "error", "message": embed_err}
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
        sample = (tasks[0].get("data") or tasks[0]) if tasks else None
        if isinstance(sample, dict) and isinstance(sample.get("image"), str):
            img_val = sample["image"]
            if img_val.startswith("data:") and len(img_val) > 80:
                sample = {**sample, "image": img_val[:48] + f"...<base64 len={len(img_val)}>"}
        logger.info(
            "Trigger_Expert_Annotation 建项成功 project_id=%s tasks=%s image_sample=%s",
            project_id,
            task_count,
            sample,
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
    image_paths: Optional[List[str]] = None,
    tasks_json: str,
) -> List[Dict[str, Any]]:
    """从 file_path / image_paths / tasks_json 解析任务列表（多图 → 多 Task）。"""
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

    paths = [str(p).strip() for p in (image_paths or []) if str(p).strip()]
    if paths:
        return [{"data": {"image": img}} for img in paths]

    return []


def is_hitl_scenario_allowed(scenario_type: str) -> bool:
    """供编排器/诊断模块查询白名单（无需实例化工具）。"""
    return str(scenario_type or "").strip() in HITL_SCENARIO_WHITELIST
