"""
插件上传 API：POST /api/plugins/upload

用户提交 → 审核暂存区（skills_review_staging）+ DB status=pending；不在广场展示、不挂载路由。
"""
from __future__ import annotations

import json
import logging
import re
import uuid as _uuid
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, Optional, Tuple

from fastapi import APIRouter, Depends, File, Form, HTTPException, UploadFile, status
from sqlalchemy.orm import Session

from gibh_agent.core.deps import get_current_user
from gibh_agent.core.skills_assets_layout import get_skills_review_staging_dir
from gibh_agent.core.user_notifications import notify_all_admins
from gibh_agent.db.connection import get_db_session
from gibh_agent.db.models import User

from gibh_agent.plugin_system.parser import parse_skill_md
from gibh_agent.plugin_system.registry import persist_plugin

logger = logging.getLogger(__name__)

router = APIRouter(prefix="/api", tags=["PluginSystem"])

_ALLOWED_DRIVER_TYPES = frozenset({"prompt", "script", "model"})
_DRIVER_TO_SKILL_TYPE = {
    "prompt": "prompt",
    "script": "script",
    "model": "model",
}


def _normalize_driver_type(raw: str) -> str:
    key = (raw or "").strip().lower()
    aliases = {
        "纯 prompt 驱动": "prompt",
        "纯prompt驱动": "prompt",
        "prompt": "prompt",
        "单脚本驱动": "script",
        "script": "script",
        "模型驱动": "model",
        "model": "model",
    }
    key = aliases.get(key, key)
    if key not in _ALLOWED_DRIVER_TYPES:
        raise HTTPException(
            status_code=status.HTTP_400_BAD_REQUEST,
            detail="driver_type 必须是 prompt、script 或 model",
        )
    return key


def _safe_user_segment(user_id: str) -> str:
    s = re.sub(r"[^a-zA-Z0-9_@.\-]+", "_", (user_id or "").strip()).strip("_")
    return (s[:128] if s else "anonymous")


def _stage_submission_files(
    *,
    user_id: str,
    driver_type: str,
    prompt_text: Optional[str],
    file_bytes: Optional[bytes],
    original_filename: Optional[str],
) -> Tuple[str, str, Dict[str, Any]]:
    """
    将提交物写入 ${UPLOAD_DIR}/skills_review_staging/{user_id}/{uuid}/。
    返回 (submission_uid, extract_dir, manifest_dict)。
    """
    submission_uid = str(_uuid.uuid4())
    dest = (get_skills_review_staging_dir() / _safe_user_segment(user_id) / submission_uid).resolve()
    dest.mkdir(parents=True, exist_ok=True)

    saved_files: list[str] = []
    pt = (prompt_text or "").strip()
    if pt:
        prompt_path = dest / "prompt.md"
        prompt_path.write_text(pt, encoding="utf-8")
        saved_files.append("prompt.md")

    if file_bytes:
        safe_name = Path(original_filename or "upload.bin").name
        if not safe_name or safe_name in (".", ".."):
            safe_name = "upload.bin"
        target = dest / safe_name
        target.write_bytes(file_bytes)
        saved_files.append(safe_name)

    manifest: Dict[str, Any] = {
        "submission_uid": submission_uid,
        "author_id": user_id,
        "driver_type": driver_type,
        "saved_files": saved_files,
        "submitted_at": datetime.now(timezone.utc).isoformat(),
        "review_status": "pending",
    }
    (dest / "manifest.json").write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2),
        encoding="utf-8",
    )
    return submission_uid, str(dest), manifest


def _infer_metadata_from_submission(
    *,
    submission_uid: str,
    driver_type: str,
    prompt_text: Optional[str],
    file_bytes: Optional[bytes],
    original_filename: Optional[str],
) -> Dict[str, Any]:
    """轻量解析 Name/Description/Parameters，仅用于审核队列展示字段。"""
    md_source = (prompt_text or "").strip()
    if not md_source and file_bytes and (original_filename or "").lower().endswith(".md"):
        md_source = file_bytes.decode("utf-8", errors="replace")

    if md_source:
        meta = parse_skill_md(md_source)
    else:
        stem = Path(original_filename or "").stem or f"submission_{submission_uid[:8]}"
        meta = {
            "name": re.sub(r"[^a-zA-Z0-9_]+", "_", stem.lower()).strip("_") or f"submission_{submission_uid[:8]}",
            "display_name": stem[:255],
            "description": f"用户提交的{driver_type}类技能（待管理员审核）",
            "parameters_schema": {"type": "object", "properties": {}, "required": []},
        }
    return meta


@router.post("/plugins/upload")
async def upload_plugin(
    driver_type: str = Form(...),
    prompt_text: Optional[str] = Form(None),
    file: Optional[UploadFile] = File(None),
    current_user: User = Depends(get_current_user),
    db: Session = Depends(get_db_session),
):
    """
    用户提交技能至管理员审核队列：落盘 skills_review_staging，DB status=pending。
    不执行脚本挂载、不注入 ToolRegistry、不在技能广场展示。
    """
    dt = _normalize_driver_type(driver_type)
    pt = (prompt_text or "").strip()
    file_bytes: Optional[bytes] = None
    filename: Optional[str] = None

    if file is not None and file.filename:
        file_bytes = await file.read()
        filename = file.filename
        if not file_bytes:
            raise HTTPException(status_code=400, detail="空文件")

    if dt == "prompt":
        if not pt and not file_bytes:
            raise HTTPException(status_code=400, detail="纯 Prompt 驱动须填写提示词或上传 .md/.zip 文件")
    else:
        if not file_bytes:
            raise HTTPException(status_code=400, detail="单脚本/模型驱动须上传 .zip 技能包")
        if not (filename or "").lower().endswith(".zip"):
            raise HTTPException(status_code=400, detail="单脚本/模型驱动仅支持 .zip 压缩包")

    if file_bytes and filename:
        lower = filename.lower()
        if dt == "prompt" and not (lower.endswith(".zip") or lower.endswith(".md")):
            raise HTTPException(status_code=400, detail="纯 Prompt 驱动上传文件仅支持 .zip 或 .md")
        if dt != "prompt" and not lower.endswith(".zip"):
            raise HTTPException(status_code=400, detail="单脚本/模型驱动仅支持 .zip")

    try:
        submission_uid, extract_dir, manifest = _stage_submission_files(
            user_id=current_user.username,
            driver_type=dt,
            prompt_text=pt or None,
            file_bytes=file_bytes,
            original_filename=filename,
        )
        meta = _infer_metadata_from_submission(
            submission_uid=submission_uid,
            driver_type=dt,
            prompt_text=pt or None,
            file_bytes=file_bytes,
            original_filename=filename,
        )
        params_schema = dict(meta.get("parameters_schema") or {"type": "object", "properties": {}})
        params_schema["_review_meta"] = {
            "driver_type": dt,
            "submission_uid": submission_uid,
            "manifest": manifest,
        }

        row = persist_plugin(
            db,
            name=meta["name"],
            display_name=meta.get("display_name"),
            description=meta.get("description"),
            parameters_schema=params_schema,
            skill_type=_DRIVER_TO_SKILL_TYPE[dt],
            script_path=None,
            extract_dir=extract_dir,
            author_id=current_user.username,
            status="pending",
        )
        display = (row.display_name or row.name or "未命名技能").strip()
        is_admin = (current_user.role or "").strip().lower() == "admin"
        if not is_admin:
            try:
                notify_all_admins(
                    db,
                    ntype="admin_plugin_submitted",
                    title="新安装技能待审核",
                    content=(
                        f"用户「{current_user.username}」提交了技能「{display}」"
                        f"（{dt}），请在管理员控制台「技能审核」中处理。（plugin_id={row.id}）"
                    ),
                    commit=False,
                )
            except Exception as ne:
                logger.warning("通知管理员 plugin 提交失败 id=%s: %s", row.id, ne)
        db.commit()
        db.refresh(row)
    except ValueError as ve:
        db.rollback()
        raise HTTPException(status_code=409, detail=str(ve)) from ve
    except HTTPException:
        raise
    except Exception as e:
        logger.exception("plugins/upload 审核提交失败: %s", e)
        raise HTTPException(status_code=400, detail=str(e)) from e

    logger.info(
        "plugin_system: 审核提交 user=%s driver=%s uid=%s plugin_id=%s path=%s",
        current_user.username,
        dt,
        submission_uid,
        row.id,
        extract_dir,
    )
    return {
        "status": "success",
        "message": "技能已提交，等待管理员审核",
        "review_status": "pending",
        "plugin_id": row.id,
        "submission_uid": submission_uid,
        "name": row.name,
        "display_name": row.display_name,
        "driver_type": dt,
        "extract_dir": extract_dir,
    }
