"""
工作流执行前加固：将本次会话解析到的上传绝对路径写入首启步骤 params，
避免 Planner 遗留占位路径或跨会话字符串导致算子读不到真实 FASTQ/mzML。

首步除 ``file_path`` 外，对含 ``data_path`` / ``image_path`` / ``mask_path`` 的模板（如影像组学）
同步绑定，避免工具签名不认 ``file_path`` 而空转。
"""
from __future__ import annotations

import logging
import os
from pathlib import Path
from typing import MutableMapping, Optional, Sequence

logger = logging.getLogger(__name__)

_PLACEHOLDER_PATHS = frozenset(
    {
        "",
        "<PENDING_UPLOAD>",
        "<待上传数据>",
        "<user_input>",
        "<USER_INPUT>",
    }
)


def _path_slot_empty(val: object) -> bool:
    if val is None:
        return True
    s = str(val).strip()
    if not s:
        return True
    if s in _PLACEHOLDER_PATHS:
        return True
    if s.startswith("<") and s.endswith(">") and len(s) >= 3:
        inner = s[1:-1].strip().lower()
        if inner in ("pending_upload", "待上传数据", "user_input"):
            return True
    return False


def resolve_all_existing_upload_paths(
    file_paths: Sequence[str],
    upload_dir: Optional[str] = None,
) -> list[str]:
    """
    依次解析 file_paths 中每一项，返回容器内真实存在且为文件的绝对路径列表（保序去重）。
    用于影像组学等多附件：首步 data_path 可安全拼接为逗号串供 resolve_omics_paths 消费。
    """
    out: list[str] = []
    seen: set[str] = set()
    for raw in file_paths:
        s = str(raw).strip() if raw is not None else ""
        if not s:
            continue
        one = resolve_primary_upload_path([s], upload_dir=upload_dir)
        if one and one not in seen:
            seen.add(one)
            out.append(one)
    return out


def resolve_primary_upload_path(
    file_paths: Sequence[str],
    upload_dir: Optional[str] = None,
) -> Optional[str]:
    """
    在 API 容器内解析「首个真实存在的文件」绝对路径。
    依次尝试：原始字符串绝对路径、basename 落在 UPLOAD_DIR、相对 uploads 的子路径。
    """
    paths = [str(p).strip() for p in file_paths if str(p).strip()]
    if not paths:
        return None
    raw = paths[0]
    ud: Optional[Path] = None
    if upload_dir:
        try:
            ud = Path(upload_dir).resolve()
        except OSError:
            ud = None

    trials: list[Path] = []
    p = Path(raw)
    trials.append(p)
    if not p.is_absolute() and ud is not None:
        trials.append(ud / raw)
        trials.append(ud / p.name)
    if p.is_absolute():
        trials.append(p)

    seen: set[str] = set()
    for cand in trials:
        try:
            rp = cand.resolve()
        except OSError:
            continue
        key = str(rp)
        if key in seen:
            continue
        seen.add(key)
        if rp.is_file():
            logger.info("[Preflight] 解析到容器内可读文件: %s", key)
            return key

    logger.error(
        "[Preflight] 文件注入失败：容器内找不到路径 raw=%r upload_dir=%r ",
        raw,
        str(ud) if ud else None,
    )
    return None


def inject_primary_file_path_into_workflow(
    workflow_data: MutableMapping[str, object],
    file_paths: Sequence[str],
    upload_dir: Optional[str] = None,
) -> bool:
    """
    将解析到的真实上传绝对路径绑定到第一个启用步骤：
    - 始终写入 params.file_path（与 Executor 反射 / 占位符解析对齐）；
    - 若 params 含 data_path 且仍为空或占位符（如 radiomics_data_validation），写入主路径或多文件逗号拼接；
    - 若含 image_path / mask_path 且槽位为空：image_path 用主路径；mask_path 在有第二路径时写入（常见图+掩膜顺序）。

    若无法在磁盘上定位任一文件，返回 False（调用方应中止执行并向前端报错）。
    """
    all_existing = resolve_all_existing_upload_paths(file_paths, upload_dir=upload_dir)
    if not all_existing:
        return False
    primary = all_existing[0]
    if not os.path.isfile(primary):
        logger.error("[Preflight] 注入中止：路径非文件 %s", primary)
        return False

    steps = workflow_data.get("steps")
    if not isinstance(steps, list):
        logger.warning("[Preflight] workflow_data.steps 非列表，跳过注入")
        return False

    injected = False
    for st in steps:
        if not isinstance(st, dict):
            continue
        if st.get("enabled") is False or st.get("selected") is False:
            continue
        params = st.get("params")
        if not isinstance(params, dict):
            params = {}
            st["params"] = params
        params["file_path"] = primary
        multi = ",".join(all_existing) if len(all_existing) > 1 else primary
        if "data_path" in params and _path_slot_empty(params.get("data_path")):
            params["data_path"] = multi
            logger.info("[Preflight] 已绑定首步 data_path（多附件=%s）", len(all_existing) > 1)
        if "image_path" in params and _path_slot_empty(params.get("image_path")):
            params["image_path"] = primary
        if (
            "mask_path" in params
            and _path_slot_empty(params.get("mask_path"))
            and len(all_existing) > 1
        ):
            params["mask_path"] = all_existing[1]
            logger.info("[Preflight] 已绑定首步 mask_path 为第二附件")
        injected = True
        logger.info("[Preflight] 已强制首步 file_path=%s（覆盖 Planner 遗留路径）", primary)
        break

    if not injected:
        logger.warning("[Preflight] 无启用步骤，未能注入 file_path")
        return False
    return True
