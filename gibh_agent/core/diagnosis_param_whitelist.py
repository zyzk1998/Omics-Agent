# -*- coding: utf-8 -*-
"""
数据诊断 LLM 用的「前端可绑定参数名」白名单。

前端表单控件名为 name="param_{步骤索引}_{参数名}"，其中「参数名」与 workflow 步骤
params 字典的键一致（通常为 snake_case）。一键应用使用这些键做 querySelector。
"""
from __future__ import annotations

import logging
from typing import Any, Dict, List, Optional, Tuple

logger = logging.getLogger(__name__)


def humanize_param_key_as_ui_label(key: str) -> str:
    """与前端 index.html 中 displayName 规则对齐：下划线转空格 + 首字母大写。"""
    if not key or not isinstance(key, str):
        return str(key)
    parts = [p for p in key.replace("-", "_").split("_") if p]
    if not parts:
        return key
    return " ".join(p[:1].upper() + p[1:] if p else "" for p in parts)


def build_diagnosis_whitelist_prompt(
    workflow: Any,
    target_step_ids: Optional[List[str]],
    file_metadata: Optional[Dict[str, Any]],
) -> Tuple[List[str], str]:
    """
    从 generate_template 产物按「步骤」列出可调 Key；返回 (扁平 key 列表供解析过滤, 注入 LLM 的多行约束文案)。

    若模板为空则回退到 DAG + get_step_metadata 的 default_params（按步骤聚合）。
    """
    flat_ordered: List[str] = []
    seen: set = set()
    step_lines: List[str] = []

    tmpl: Optional[Dict[str, Any]] = None
    if workflow is not None:
        try:
            ts = target_step_ids if target_step_ids else None
            tmpl = workflow.generate_template(ts, file_metadata)
        except Exception as e:
            logger.debug("build_diagnosis_whitelist_prompt: generate_template 失败: %s", e)

    steps: List[Any] = []
    if tmpl and isinstance(tmpl.get("workflow_data"), dict):
        steps = tmpl["workflow_data"].get("steps") or []

    if steps:
        for st in steps:
            if not isinstance(st, dict):
                continue
            sid = st.get("step_id") or st.get("id") or ""
            sname = st.get("name") or st.get("step_name") or sid or "步骤"
            params = st.get("params")
            if not isinstance(params, dict) or not params:
                continue
            keys = [k for k in params if isinstance(k, str)]
            if not keys:
                continue
            keys_csv = ", ".join(keys)
            step_lines.append(
                f"- Step: {sname} (step_id={sid}) -> Valid Keys: [{keys_csv}]"
            )
            for k in keys:
                if k not in seen:
                    seen.add(k)
                    flat_ordered.append(k)
                lbl = humanize_param_key_as_ui_label(k)
                step_lines.append(
                    f"  - Key: {k}  |  UI 标签（仅作阅读，表格第一列禁止写标签）: {lbl}"
                )
    else:
        # 回退：无 steps 时按 DAG 元数据逐 step 列出
        try:
            dag = getattr(workflow, "steps_dag", None) or {}
            if target_step_ids:
                valid = [s for s in target_step_ids if s in dag]
                if not valid:
                    valid = list(dag.keys())
                resolved = workflow.resolve_dependencies(valid)
            else:
                resolved = list(dag.keys())
        except Exception:
            resolved = []
        for step_id in resolved:
            try:
                meta = workflow.get_step_metadata(step_id)
            except Exception:
                meta = {}
            sname = (meta or {}).get("name", step_id)
            params = (meta or {}).get("default_params") or {}
            if not isinstance(params, dict) or not params:
                continue
            keys = [k for k in params if isinstance(k, str)]
            if not keys:
                continue
            keys_csv = ", ".join(keys)
            step_lines.append(
                f"- Step: {sname} (step_id={step_id}) -> Valid Keys: [{keys_csv}]"
            )
            for k in keys:
                if k not in seen:
                    seen.add(k)
                    flat_ordered.append(k)
                lbl = humanize_param_key_as_ui_label(k)
                step_lines.append(
                    f"  - Key: {k}  |  UI 标签（仅作阅读，表格第一列禁止写标签）: {lbl}"
                )

    body = "\n".join(step_lines) if step_lines else "(当前未能解析出任何可调参数；请勿输出参数推荐表)"
    prompt = f"""

【致命约束：参数推荐规则】
如果你决定推荐参数，你必须且只能从下列「各步骤 Valid Keys」中出现过的 **Key** 里挑选（这些 Key 与前端工作流卡片输入框绑定字段一致；一键应用使用这些字符串匹配 name="param_{{步骤索引}}_{{Key}}"）：

{body}

注意：
1. 表格的第一列「参数名」必须严格使用上述某一 Key，**一字不差**（区分大小写，通常为 snake_case）。**绝不允许**捏造列表外不存在的参数（例如 log_transform 若未出现在上表中则禁止写入）！
2. **绝不允许**用 UI 标签、中文名、Title Case 或简写代替 Key；不满足则**宁可省略整张参数表**，仅在正文定性描述。
3. 表格第一列**禁止**用 Markdown 加粗（**）或反引号（`）包裹 Key。
"""
    return flat_ordered, prompt


def collect_param_names_from_plan_result(plan: Optional[Dict[str, Any]]) -> List[str]:
    """从 planner / generate_template 返回的 workflow 结构中收集所有 params 键（保序去重）。"""
    if not plan or not isinstance(plan, dict):
        return []
    wd = plan.get("workflow_data")
    inner: Dict[str, Any] = wd if isinstance(wd, dict) else plan
    steps = inner.get("steps") or []
    if not isinstance(steps, list):
        return []
    seen = set()
    ordered: List[str] = []
    for st in steps:
        if not isinstance(st, dict):
            continue
        params = st.get("params")
        if isinstance(params, dict):
            for k in params.keys():
                if isinstance(k, str) and k not in seen:
                    seen.add(k)
                    ordered.append(k)
        elif isinstance(params, list):
            for item in params:
                if isinstance(item, dict):
                    n = item.get("name") or item.get("key")
                    if n is not None:
                        s = str(n)
                        if s and s not in seen:
                            seen.add(s)
                            ordered.append(s)
    return ordered


def collect_param_whitelist_for_diagnosis(
    workflow: Any,
    target_step_ids: Optional[List[str]],
    file_metadata: Optional[Dict[str, Any]],
) -> List[str]:
    """
    优先用 workflow.generate_template（与前端最终 steps 结构一致），失败则回退到
    steps_dag + get_step_metadata 的 default_params 键。
    """
    if workflow is None:
        return []
    try:
        ts = target_step_ids if target_step_ids else None
        tmpl = workflow.generate_template(ts, file_metadata)
        names = collect_param_names_from_plan_result(tmpl)
        if names:
            return names
    except Exception as e:
        logger.debug("collect_param_whitelist_for_diagnosis: generate_template 失败，回退 DAG: %s", e)

    try:
        dag = getattr(workflow, "steps_dag", None) or {}
        if target_step_ids:
            valid = [s for s in target_step_ids if s in dag]
            if not valid:
                valid = list(dag.keys())
            resolved = workflow.resolve_dependencies(valid)
        else:
            resolved = list(dag.keys())
    except Exception as e:
        logger.debug("collect_param_whitelist_for_diagnosis: DAG 解析失败: %s", e)
        resolved = []

    seen = set()
    ordered: List[str] = []
    for step_id in resolved:
        try:
            meta = workflow.get_step_metadata(step_id)
        except Exception:
            meta = {}
        params = (meta or {}).get("default_params") or {}
        if isinstance(params, dict):
            for k in params.keys():
                if isinstance(k, str) and k not in seen:
                    seen.add(k)
                    ordered.append(k)
    return ordered
