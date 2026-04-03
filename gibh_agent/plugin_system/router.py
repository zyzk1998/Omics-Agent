"""
插件上传 API：POST /api/plugins/upload
"""
from __future__ import annotations

import logging
import os
import tempfile
import uuid as _uuid

from fastapi import APIRouter, Depends, File, HTTPException, UploadFile, status
from sqlalchemy.orm import Session

from gibh_agent.core.deps import get_current_user
from gibh_agent.db.connection import get_db_session
from gibh_agent.db.models import User

from gibh_agent.plugin_system.parser import get_dynamic_skills_root, parse_skill_markdown_upload, parse_skill_zip
from gibh_agent.plugin_system.registry import persist_plugin

logger = logging.getLogger(__name__)

router = APIRouter(prefix="/api", tags=["PluginSystem"])


def _require_zip_or_md(filename: str) -> str:
    lower = (filename or "").lower()
    if lower.endswith(".zip"):
        return "zip"
    if lower.endswith(".md"):
        return "md"
    raise HTTPException(
        status_code=status.HTTP_400_BAD_REQUEST,
        detail="仅支持 .zip 或 .md",
    )


@router.post("/plugins/upload")
async def upload_plugin(
    file: UploadFile = File(...),
    current_user: User = Depends(get_current_user),
    db: Session = Depends(get_db_session),
):
    """ZIP：含 main.py → script；仅文档 → prompt。纯 .md → prompt，不生成假 main.py。"""
    kind = _require_zip_or_md(file.filename or "")
    suffix = ".zip" if kind == "zip" else ".md"

    data = await file.read()
    if not data:
        raise HTTPException(status_code=400, detail="空文件")

    with tempfile.NamedTemporaryFile(delete=False, suffix=suffix) as tmp:
        tmp.write(data)
        tmp_path = tmp.name

    try:
        if kind == "zip":
            parsed = parse_skill_zip(tmp_path)
        else:
            text = data.decode("utf-8", errors="replace")
            root = get_dynamic_skills_root()
            uid = str(_uuid.uuid4())
            dest = (root / uid).resolve()
            parsed = parse_skill_markdown_upload(text, dest, uid)

        try:
            row = persist_plugin(
                db,
                name=parsed["name"],
                display_name=parsed.get("display_name"),
                description=parsed.get("description"),
                parameters_schema=parsed.get("parameters_schema") or {"type": "object", "properties": {}},
                skill_type=parsed.get("skill_type") or "prompt",
                script_path=parsed.get("main_py_path"),
                extract_dir=parsed.get("extract_dir"),
                author_id=current_user.username,
            )
        except ValueError as ve:
            db.rollback()
            raise HTTPException(status_code=409, detail=str(ve)) from ve

        return {
            "status": "success",
            "plugin_id": row.id,
            "name": row.name,
            "display_name": row.display_name,
            "skill_type": row.skill_type,
            "script_path": row.script_path,
            "worker_route": row.worker_route,
            "extract_dir": row.extract_dir,
            "parameters_schema": row.parameters_schema,
        }
    except HTTPException:
        raise
    except Exception as e:
        logger.exception("plugins/upload 失败: %s", e)
        raise HTTPException(status_code=400, detail=str(e)) from e
    finally:
        try:
            os.unlink(tmp_path)
        except OSError:
            pass
