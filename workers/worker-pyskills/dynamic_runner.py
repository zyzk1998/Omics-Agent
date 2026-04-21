"""
动态技能沙盒：仅允许加载 dynamic_skills 下的 main.py 并调用 run(**kwargs)。
"""
from __future__ import annotations

import importlib.util
import logging
import os
import uuid
from pathlib import Path
from typing import Any, Dict

from fastapi import APIRouter, HTTPException
from pydantic import BaseModel, Field

logger = logging.getLogger(__name__)

router = APIRouter(tags=["DynamicSkill"])

def _allowed_dynamic_root() -> Path:
    """
    与 gibh_agent/core/skills_assets_layout.get_dynamic_skills_dir 保持一致（本镜像不含 gibh_agent 包）。
    优先级：DYNAMIC_SKILLS_DIR → ${UPLOAD_DIR}/skills_assets/script/runtime
    """
    override = os.environ.get("DYNAMIC_SKILLS_DIR")
    if override and str(override).strip():
        root = Path(os.path.abspath(override.strip()))
        root.mkdir(parents=True, exist_ok=True)
        return root
    upload = os.environ.get("UPLOAD_DIR", "/app/uploads").strip()
    root = Path(os.path.abspath(upload)) / "skills_assets" / "script" / "runtime"
    root.mkdir(parents=True, exist_ok=True)
    return root


ALLOWED_ROOT = _allowed_dynamic_root()


class DynamicRunBody(BaseModel):
    script_path: str = Field(..., description="main.py 绝对路径，必须在 dynamic_skills 下")
    kwargs: Dict[str, Any] = Field(default_factory=dict)


def _assert_script_in_sandbox(script_path: str) -> Path:
    raw = (script_path or "").strip()
    if not raw:
        raise HTTPException(status_code=400, detail="script_path 为空")
    p = Path(os.path.abspath(raw))
    root = Path(os.path.abspath(str(ALLOWED_ROOT)))
    try:
        p.relative_to(root)
    except ValueError as e:
        logger.warning("dynamic/run 路径越界: %s (root=%s)", p, root)
        raise HTTPException(status_code=403, detail="script_path 不在允许的 dynamic_skills 目录下") from e
    if not p.is_file():
        raise HTTPException(status_code=400, detail=f"文件不存在: {p}")
    if p.name != "main.py":
        raise HTTPException(status_code=400, detail="仅允许加载名为 main.py 的入口文件")
    return p


@router.post("/api/dynamic/run")
def dynamic_run(body: DynamicRunBody) -> Dict[str, Any]:
    path = _assert_script_in_sandbox(body.script_path)
    mod_name = f"_dyn_skill_{uuid.uuid4().hex}"
    spec = importlib.util.spec_from_file_location(mod_name, str(path))
    if spec is None or spec.loader is None:
        raise HTTPException(status_code=500, detail="无法创建模块 spec")
    module = importlib.util.module_from_spec(spec)
    try:
        spec.loader.exec_module(module)
    except Exception as e:
        logger.exception("exec_module failed: %s", path)
        raise HTTPException(status_code=500, detail=f"加载模块失败: {e}") from e
    run_fn = getattr(module, "run", None)
    if not callable(run_fn):
        raise HTTPException(status_code=500, detail="main.py 中未找到可调用对象 run")
    try:
        out = run_fn(**(body.kwargs or {}))
    except Exception as e:
        logger.exception("run(**kwargs) failed: %s", path)
        raise HTTPException(status_code=500, detail=f"run 执行失败: {e}") from e
    if not isinstance(out, dict):
        return {"status": "success", "message": "run 返回非 dict，已包装", "data": out}
    return out
