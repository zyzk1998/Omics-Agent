"""
工作流执行前加固：将本次会话解析到的上传绝对路径写入首步 params，
避免 Planner 遗留占位路径或跨会话字符串导致算子读不到真实 FASTQ/mzML。
"""
from __future__ import annotations

import logging
import os
from pathlib import Path
from typing import MutableMapping, Optional, Sequence

logger = logging.getLogger(__name__)


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
    将解析到的首个真实文件绝对路径写入第一个启用步骤的 params.file_path。
    若无法在磁盘上定位文件，返回 False（调用方应中止执行并向前端报错）。
    """
    primary = resolve_primary_upload_path(file_paths, upload_dir=upload_dir)
    if not primary:
        return False
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
        injected = True
        logger.info("[Preflight] 已强制首步 file_path=%s（覆盖 Planner 遗留路径）", primary)
        break

    if not injected:
        logger.warning("[Preflight] 无启用步骤，未能注入 file_path")
        return False
    return True
