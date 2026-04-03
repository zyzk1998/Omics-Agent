"""
动态技能执行代理：prompt 类返回编排提示文本；script 类转发 worker-pyskills。
Executor 反射：显式形参仅 skill_name + tool_params_json（额外脚本参数以 JSON 传入）。
"""
from __future__ import annotations

import json
import logging
import os
from typing import Any, Dict

import requests

from gibh_agent.core.tool_registry import registry
from gibh_agent.core.utils import safe_tool_execution
from gibh_agent.db.connection import SessionLocal
from gibh_agent.plugin_system.registry import DynamicSkillPlugin, default_worker_base_url

logger = logging.getLogger(__name__)


@registry.register(
    name="execute_dynamic_skill",
    description=(
        "执行用户上架的动态技能：skill_name 为插件英文标识（与 SKILL 中 Name 对应 slug）；"
        "script 类将 kwargs 发往 PySkills Worker；prompt 类返回提示词与 Schema 供模型遵循。"
    ),
    category="PluginSystem",
    output_type="json",
)
@safe_tool_execution
def execute_dynamic_skill(
    skill_name: str,
    tool_params_json: str = "{}",
) -> Dict[str, Any]:
    """
    按技能类型分流执行动态插件。

    Args:
        skill_name: 动态技能库中的 ``name``（小写 slug，与上传解析一致）。
        tool_params_json: JSON 对象字符串，键值对作为 ``run(**kwargs)`` 传入 worker（仅 ``skill_type=script`` 时使用）。

    Returns:
        统一字典：至少含 ``status`` (``success``|``error``)、``message``；script 分支另含 ``worker_result``。
    """
    key = (skill_name or "").strip()
    if not key:
        return {"status": "error", "message": "skill_name 不能为空"}

    try:
        extra: Dict[str, Any] = json.loads(tool_params_json or "{}")
        if not isinstance(extra, dict):
            return {"status": "error", "message": "tool_params_json 必须是 JSON 对象"}
    except json.JSONDecodeError as e:
        return {"status": "error", "message": f"tool_params_json 解析失败: {e}"}

    db = SessionLocal()
    try:
        row = db.query(DynamicSkillPlugin).filter(DynamicSkillPlugin.name == key).first()
        if not row or (row.status or "") != "approved":
            return {"status": "error", "message": f"未找到已审核的动态技能: {key!r}"}

        st = (row.skill_type or "script").strip().lower()
        if st == "prompt":
            schema_str = json.dumps(row.parameters_schema or {}, ensure_ascii=False, indent=2)
            text = (
                "这是一个提示词技能（无远程脚本执行）。请根据以下描述、参数 Schema 与用户输入进行专业回答，"
                "不要虚构已执行代码或数值结果。\n\n"
                f"【技能说明】\n{row.description or '（无）'}\n\n"
                f"【参数 Schema（JSON）】\n{schema_str}\n"
            )
            return {
                "status": "success",
                "message": text,
                "skill_type": "prompt",
                "skill_name": row.name,
                "display_name": row.display_name,
            }

        if st != "script":
            return {"status": "error", "message": f"未知 skill_type: {st!r}"}

        sp = row.script_path
        if not sp or not str(sp).strip():
            return {"status": "error", "message": "脚本技能缺少 script_path"}

        base = default_worker_base_url()
        route = (row.worker_route or "/api/dynamic/run").strip()
        if not route.startswith("/"):
            route = "/" + route
        url = f"{base}{route}"
        timeout_s = float(os.environ.get("DYNAMIC_SKILL_HTTP_TIMEOUT", "600"))
        logger.info("execute_dynamic_skill POST %s script_path=%s", url, sp)
        resp = requests.post(
            url,
            json={"script_path": sp, "kwargs": extra},
            timeout=timeout_s,
        )
        try:
            payload = resp.json()
        except Exception:
            payload = {"raw": resp.text[:8000]}
        if resp.status_code >= 400:
            return {
                "status": "error",
                "message": payload.get("detail") if isinstance(payload, dict) else str(payload),
                "http_status": resp.status_code,
                "worker_result": payload,
            }
        return {
            "status": "success",
            "message": "worker 执行完成",
            "skill_type": "script",
            "skill_name": row.name,
            "worker_result": payload,
        }
    finally:
        db.close()
