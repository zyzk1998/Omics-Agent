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

# ---------------------------------------------------------------------------
# 组学三模态：LLM 参数推荐「硬兜底」——与 Genomics/Proteomics/Epigenomics 工作流
# get_step_metadata.default_params 及工具签名对齐；在 generate_template 偶发无参时仍保证白名单非空。
# ---------------------------------------------------------------------------
OMICS_GENOMICS_ALGO_PARAM_KEYS: frozenset = frozenset(
    {
        "min_read_length",
        "quality_cutoff",
        "threads",
        "mismatch_penalty",
        "gap_open_penalty",
        "optical_duplicate_distance",
        "bqsr_max_cycles",
        "min_base_quality",
        "min_mapping_quality",
        "stand_call_conf",
        "cnv_bin_width",
        "min_sv_len",
        "tranche_sensitivity",
        "pick_allele",
        "pp2_ba1_threshold",
        "report_style",
    }
)
OMICS_PROTEOMICS_ALGO_PARAM_KEYS: frozenset = frozenset(
    {
        "mz_tolerance_ppm",
        "snr_threshold",
        "fragment_tol_da",
        "missed_cleavages",
        "peptide_fdr",
        "target_fdr",
        "razor_min_peptides",
        "lfq_ratio_type",
        "knn_neighbors",
        "norm_method",
        "group_column",
        "log2fc_cutoff",
        "n_top_features",
        "enrich_padj",
        "string_score_cutoff",
        "report_depth",
    }
)
OMICS_EPIGENOMICS_ALGO_PARAM_KEYS: frozenset = frozenset(
    {
        "trim_quality_threshold",
        "threads",
        "mismatch_penalty",
        "min_mapq",
        "shift_correction_bp",
        "qvalue_threshold",
        "broad_peak",
        "idr_threshold",
        "merge_distance_bp",
        "promoter_window_bp",
        "padj_cutoff",
        "motif_e_value",
        "nuc_resolution_bp",
        "max_distance_bp",
        "min_correlation",
    }
)

OMICS_DOMAIN_ALGO_KEYS: Dict[str, frozenset] = {
    "genomics": OMICS_GENOMICS_ALGO_PARAM_KEYS,
    "proteomics": OMICS_PROTEOMICS_ALGO_PARAM_KEYS,
    "epigenomics": OMICS_EPIGENOMICS_ALGO_PARAM_KEYS,
}


def _workflow_domain_name(workflow: Any) -> str:
    try:
        n = workflow.get_name()
        return n if isinstance(n, str) else ""
    except Exception:  # noqa: BLE001
        return ""


def merge_omics_static_whitelist_keys(workflow: Any, keys: List[str]) -> List[str]:
    """将三模态硬编码算法键并入列表（去重保序）。"""
    dom = _workflow_domain_name(workflow)
    extra = OMICS_DOMAIN_ALGO_KEYS.get(dom) or frozenset()
    if not extra:
        return list(keys)
    seen = set(keys)
    out = list(keys)
    for k in sorted(extra):
        if k in seen:
            continue
        if not _is_algorithmic_param_key(k):
            continue
        seen.add(k)
        out.append(k)
    return out


# 仅路径/资产流转类参数：不应出现在「💡 参数推荐」表格中凑数（无可调算法语义时整步静默）
PIPELINE_PLUMBING_PARAM_KEYS = frozenset(
    {
        "file_path",
        "data_path",
        "input_dir",
        "matrix_dir",
        "adata_path",
        "h5ad_path",
        "image_path",
        "mask_path",
        "sequence_or_path",
        "counts_file",
        "raw_counts_path",
        "reference_path",
        "reference_id",
        "vcf_path",
        "bam_path",
        "ingress_file_path",
    }
)


def _is_algorithmic_param_key(key: str) -> bool:
    if not key or not isinstance(key, str):
        return False
    return key not in PIPELINE_PLUMBING_PARAM_KEYS


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
            keys = [k for k in params if isinstance(k, str) and _is_algorithmic_param_key(k)]
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
            keys = [k for k in params if isinstance(k, str) and _is_algorithmic_param_key(k)]
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

    dom = _workflow_domain_name(workflow)
    extra_dom = OMICS_DOMAIN_ALGO_KEYS.get(dom)
    if extra_dom:
        step_lines.append(
            f"- **{dom} 模态·诊断兜底 Key 列表（与 workflow 工具 default_params / 签名对齐；仅当该步骤实际暴露此键时可推荐）**: "
            + ", ".join(sorted(extra_dom))
        )

    flat_ordered = merge_omics_static_whitelist_keys(workflow, flat_ordered)

    body = "\n".join(step_lines) if step_lines else "(当前未能解析出任何可调参数；请勿输出参数推荐表)"
    prompt = f"""

【致命约束：参数推荐规则】
如果你决定推荐参数，你必须且只能从下列「各步骤 Valid Keys」中出现过的 **Key** 里挑选（这些 Key 与前端工作流卡片输入框绑定字段一致；一键应用使用这些字符串匹配 name="param_{{步骤索引}}_{{Key}}"）：

{body}

注意：
1. 表格的第一列「参数名」必须严格使用上述某一 Key，**一字不差**（区分大小写，通常为 snake_case）。**绝不允许**捏造列表外不存在的参数（例如 log_transform 若未出现在上表中则禁止写入）！
2. **绝不允许**用 UI 标签、中文名、Title Case 或简写代替 Key；不满足则**宁可省略整张参数表**，仅在正文定性描述。
3. 表格第一列**禁止**用 Markdown 加粗（**）或反引号（`）包裹 Key。
4. **静默规则**：若某步骤仅有 file_path / input_dir 等路径流转、且无可调算法阈值，则**禁止**为该步骤单独占一行凑数；白名单中已剔除此类 Key。
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
                if (
                    isinstance(k, str)
                    and _is_algorithmic_param_key(k)
                    and k not in seen
                ):
                    seen.add(k)
                    ordered.append(k)
        elif isinstance(params, list):
            for item in params:
                if isinstance(item, dict):
                    n = item.get("name") or item.get("key")
                    if n is None:
                        continue
                    s = str(n).strip()
                    # list 形态必须与 dict 键一致：剔除管线 plumbing（file_path 等）
                    if not s or not _is_algorithmic_param_key(s):
                        continue
                    if s not in seen:
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
            return merge_omics_static_whitelist_keys(workflow, names)
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
                if (
                    isinstance(k, str)
                    and _is_algorithmic_param_key(k)
                    and k not in seen
                ):
                    seen.add(k)
                    ordered.append(k)
    return merge_omics_static_whitelist_keys(workflow, ordered)
