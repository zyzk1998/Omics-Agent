"""
工具执行前：参数整形与路径安检门。

在 Executor / ToolRegistry 真正调用底层算子前，对路径类参数做：
1. 从 uploaded_files / step_context 智能回退补齐；
2. file_path ↔ adata_path 等同义键映射；
3. os.path.exists 物理校验，失败时返回规范 JSON 错误结构。
"""
from __future__ import annotations

import inspect
import logging
import os
from contextvars import ContextVar
from pathlib import Path
from typing import Any, Callable, Dict, List, MutableMapping, Optional, Sequence, Set, Tuple

logger = logging.getLogger(__name__)

# 供 ToolRegistry wrapper 在无显式 context 时读取（Executor 会在 execute_step 内 set）
_tool_input_context: ContextVar[Optional[Dict[str, Any]]] = ContextVar(
    "tool_input_context", default=None
)

_PLACEHOLDER_PATHS = frozenset(
    {
        "",
        "<PENDING_UPLOAD>",
        "<待上传数据>",
        "<user_input>",
        "<USER_INPUT>",
    }
)

# 架构宪法允许的路径类参数名 + 组学域扩展
_PATH_PARAM_NAMES: frozenset[str] = frozenset(
    {
        "file_path",
        "data_path",
        "adata_path",
        "h5ad_path",
        "matrix_dir",
        "input_dir",
        "image_path",
        "mask_path",
        "sequence_or_path",
        "fastqs_path",
        "fastq_path",
        "reference_path",
        "cellranger_matrix_dir",
        "transcriptome_path",
        "trajectory_data_path",
        "features_csv_path",
        "rad_score_csv_path",
    }
)

# 允许目录存在的参数（其余默认期望文件或可解析为文件的路径）
_DIR_ALLOWED_PARAMS: frozenset[str] = frozenset(
    {
        "matrix_dir",
        "input_dir",
        "fastqs_path",
        "fastq_path",
        "cellranger_matrix_dir",
        "adata_path",  # 10x 目录
        "file_path",
        "data_path",
        "transcriptome_path",
    }
)

# 输出路径：存在性校验跳过（工具会创建）
_OUTPUT_PATH_PARAMS: frozenset[str] = frozenset(
    {
        "output_path",
        "output_file",
        "output_dir",
        "output_h5ad",
        "output_h5ad_path",
    }
)


def set_tool_input_context(context: Optional[Dict[str, Any]]) -> Any:
    """设置当前协程/线程的工具输入上下文 token（Executor 调用）。"""
    return _tool_input_context.set(context)


def reset_tool_input_context(token: Any) -> None:
    _tool_input_context.reset(token)


def get_tool_input_context() -> Optional[Dict[str, Any]]:
    return _tool_input_context.get()


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


def _resolve_upload_dir(context: Optional[Dict[str, Any]]) -> Optional[Path]:
    if not context:
        return None
    for key in ("upload_dir", "UPLOAD_DIR"):
        raw = context.get(key)
        if raw:
            try:
                return Path(str(raw)).resolve()
            except OSError:
                pass
    return None


def _absolutize_path(raw: str, upload_dir: Optional[Path]) -> str:
    s = str(raw).strip()
    if not s:
        return s
    p = Path(s)
    if p.is_absolute():
        return str(p)
    if upload_dir is not None:
        return str((upload_dir / s).resolve())
    return s


def _collect_uploaded_file_entries(context: Optional[Dict[str, Any]]) -> List[Dict[str, Any]]:
    if not context:
        return []
    out: List[Dict[str, Any]] = []
    for key in ("uploaded_files", "files"):
        items = context.get(key)
        if isinstance(items, list):
            for item in items:
                if isinstance(item, dict):
                    out.append(item)
    return out


def _collect_candidate_paths(context: Optional[Dict[str, Any]]) -> List[str]:
    """从上下文收集候选绝对/相对路径（保序去重）。"""
    seen: Set[str] = set()
    candidates: List[str] = []
    upload_dir = _resolve_upload_dir(context)

    def _add(raw: object) -> None:
        if raw is None:
            return
        s = str(raw).strip()
        if not s or _path_slot_empty(s):
            return
        abs_s = _absolutize_path(s, upload_dir)
        key = abs_s.rstrip("/")
        if key in seen:
            return
        seen.add(key)
        candidates.append(abs_s)

    if context:
        omics = context.get("omics_resolved") or {}
        if isinstance(omics, dict):
            for bucket in ("10x_mtx", "h5ad", "tables", "images", "masks"):
                for p in omics.get(bucket) or []:
                    _add(p)
        for key in (
            "resolved_input_paths",
            "file_paths",
            "current_file_path",
        ):
            val = context.get(key)
            if isinstance(val, list):
                for p in val:
                    _add(p)
            elif val:
                _add(val)
        fm = context.get("file_metadata") or {}
        if isinstance(fm, dict):
            _add(fm.get("real_data_path"))
            _add(fm.get("file_path"))

    for entry in _collect_uploaded_file_entries(context):
        if entry.get("is_10x") and entry.get("group_dir"):
            _add(entry["group_dir"])
        _add(entry.get("path") or entry.get("file_path"))

    # 逗号拼接串展开
    expanded: List[str] = []
    for c in candidates:
        if "," in c or ";" in c:
            try:
                from gibh_agent.core.path_resolvers import resolve_omics_paths

                resolved = resolve_omics_paths(c)
                for bucket in ("10x_mtx", "h5ad", "tables"):
                    for p in resolved.get(bucket) or []:
                        if p not in seen:
                            seen.add(p.rstrip("/"))
                            expanded.append(p)
            except Exception:
                for part in c.replace(";", ",").split(","):
                    part = part.strip()
                    if part and part.rstrip("/") not in seen:
                        seen.add(part.rstrip("/"))
                        expanded.append(part)
        else:
            expanded.append(c)
    return expanded


def _pick_primary_path(candidates: Sequence[str], tool_id: str = "", context: Optional[Dict[str, Any]] = None) -> Optional[str]:
    """优先链式 h5ad → 10x 目录 → 其它存在路径。"""
    if context:
        chain = context.get("current_file_path")
        if isinstance(chain, str) and chain.strip().lower().endswith(".h5ad") and os.path.isfile(chain.strip()):
            if tool_id not in ("rna_data_validation", "rna_qc_filter"):
                return chain.strip()
    if not candidates:
        return None
    upload_hints = ("10x_data", "10x_mtx", "matrix.mtx")
    h5ad_files = [c for c in candidates if c.lower().endswith((".h5ad", ".h5")) and os.path.isfile(c)]
    if h5ad_files and tool_id not in ("rna_data_validation",):
        return h5ad_files[0]
    dir_candidates = [c for c in candidates if os.path.isdir(c)]
    for c in dir_candidates:
        low = c.lower()
        if any(h in low for h in upload_hints) or "matrix.mtx" in low:
            return c
    for c in candidates:
        if c.lower().endswith((".h5ad", ".h5")):
            return c
    for c in candidates:
        if os.path.exists(c):
            return c
    return candidates[0]


def _path_looks_like_upload_bundle(path_val: object) -> bool:
    if not isinstance(path_val, str) or not path_val.strip():
        return False
    p = Path(path_val.strip())
    if not p.is_dir():
        return False
    low = str(p).lower()
    if "10x_data" in low:
        return True
    try:
        names = {x.lower() for x in os.listdir(p)}
    except OSError:
        return False
    return any(
        n in names
        for n in (
            "matrix.mtx",
            "matrix.mtx.gz",
            "barcodes.tsv",
            "barcodes.tsv.gz",
            "features.tsv",
            "features.tsv.gz",
        )
    )


def _tool_path_param_names(tool_func: Callable) -> Set[str]:
    try:
        sig = inspect.signature(tool_func)
    except (TypeError, ValueError):
        return set()
    names: Set[str] = set()
    for pname in sig.parameters:
        if pname in _PATH_PARAM_NAMES:
            names.add(pname)
    return names


def _normalize_synonym_keys(
    params: MutableMapping[str, Any],
    accepted: Set[str],
) -> None:
    """file_path ↔ adata_path 等同义键互补。"""
    if "adata_path" in accepted and _path_slot_empty(params.get("adata_path")):
        for src in ("file_path", "data_path", "h5ad_path"):
            if src in params and not _path_slot_empty(params.get(src)):
                params["adata_path"] = params[src]
                logger.info("[InputValidator] 映射 %s -> adata_path", src)
                break
    if "file_path" in accepted and _path_slot_empty(params.get("file_path")):
        for src in ("adata_path", "data_path"):
            if src in params and not _path_slot_empty(params.get(src)):
                params["file_path"] = params[src]
                logger.info("[InputValidator] 映射 %s -> file_path", src)
                break
    if "fastqs_path" in accepted and _path_slot_empty(params.get("fastqs_path")):
        for src in ("file_path", "data_path", "input_dir"):
            if src in params and not _path_slot_empty(params.get(src)):
                params["fastqs_path"] = params[src]
                break
    if "cellranger_matrix_dir" in accepted and _path_slot_empty(params.get("cellranger_matrix_dir")):
        for src in ("file_path", "adata_path", "matrix_dir"):
            if src in params and not _path_slot_empty(params.get(src)):
                params["cellranger_matrix_dir"] = params[src]
                break


def _apply_path_fallbacks(
    params: MutableMapping[str, Any],
    accepted: Set[str],
    primary: Optional[str],
    tool_id: str = "",
) -> None:
    if not primary:
        return
    for pname in accepted:
        if pname in _OUTPUT_PATH_PARAMS:
            continue
        if _path_slot_empty(params.get(pname)):
            params[pname] = primary
            logger.info(
                "[InputValidator] 回退补齐 %s=%s（来自 uploaded_files 上下文）",
                pname,
                primary,
            )


def _prefer_chain_h5ad_over_upload_dir(
    params: MutableMapping[str, Any],
    accepted: Set[str],
    context: Optional[Dict[str, Any]],
    tool_id: str,
) -> None:
    """若 adata_path 仍为上传 10x 目录，但上下文已有前序 h5ad，则替换为链式路径。"""
    if tool_id in ("rna_data_validation", "rna_qc_filter"):
        return
    chain: Optional[str] = None
    if context:
        cfp = context.get("current_file_path")
        if isinstance(cfp, str) and cfp.strip().lower().endswith(".h5ad") and os.path.isfile(cfp.strip()):
            chain = cfp.strip()
    if not chain:
        return
    for pname in ("adata_path", "h5ad_path", "file_path"):
        if pname not in accepted:
            continue
        cur = params.get(pname)
        if _path_slot_empty(cur):
            continue
        if _path_looks_like_upload_bundle(cur):
            params[pname] = chain
            logger.info(
                "[InputValidator] 链式 h5ad 替换上传目录 %s: %s -> %s",
                pname,
                cur,
                chain,
            )


def _resolve_param_on_disk(
    raw: str,
    upload_dir: Optional[Path],
    param_name: str,
    context: Optional[Dict[str, Any]] = None,
) -> Tuple[Optional[str], List[str]]:
    """尝试解析为容器内真实路径；返回 (resolved_or_none, attempted)."""
    attempted: List[str] = []
    if not raw or _path_slot_empty(raw):
        return None, attempted

    ctx = context or {}
    allow_dir = param_name in _DIR_ALLOWED_PARAMS
    context_paths = _collect_candidate_paths(ctx)
    uploaded_entries = _collect_uploaded_file_entries(ctx)

    try:
        from gibh_agent.core.asset_locator import resolve_asset_path

        resolved = resolve_asset_path(
            raw,
            upload_dir=str(upload_dir) if upload_dir else None,
            owner_id=ctx.get("owner_id"),
            db=ctx.get("db"),
            context_paths=context_paths,
            uploaded_entries=uploaded_entries,
            must_exist=True,
        )
        attempted.append(resolved)
        p = Path(resolved)
        if p.is_file() or (allow_dir and p.is_dir()):
            return resolved, attempted
    except Exception as exc:
        logger.debug("[InputValidator] asset_locator 解析失败，回退本地 trials: %s", exc)

    trials: List[Path] = []
    p = Path(str(raw).strip())
    trials.append(p)
    if upload_dir is not None:
        if not p.is_absolute():
            trials.append(upload_dir / raw)
            trials.append(upload_dir / p.name)
        else:
            trials.append(upload_dir / p.name)

    seen: Set[str] = set()
    for cand in trials:
        try:
            rp = cand.resolve()
        except OSError:
            continue
        key = str(rp)
        if key in seen:
            continue
        seen.add(key)
        attempted.append(key)
        if rp.exists():
            if rp.is_file() or (allow_dir and rp.is_dir()):
                return key, attempted

    # Basename 降维模糊自愈（上下文 + owner 目录）
    from gibh_agent.core.omics_io_registry import fuzzy_heal_path_by_basename

    healed, heal_attempts = fuzzy_heal_path_by_basename(
        raw,
        upload_dir=upload_dir,
        owner_id=ctx.get("owner_id"),
        context_paths=context_paths,
        uploaded_entries=uploaded_entries,
        allow_dir=allow_dir,
    )
    attempted.extend(heal_attempts)
    if healed:
        logger.info("[InputValidator] 模糊自愈 %s: %s -> %s", param_name, raw, healed)
        return healed, attempted

    return None, attempted


def _fuzzy_heal_path_like_param_values(
    params: MutableMapping[str, Any],
    context: Optional[Dict[str, Any]],
) -> None:
    """对 params 内所有路径型字符串做 basename 自愈（不限于签名内路径参数名）。"""
    from gibh_agent.core.omics_io_registry import fuzzy_heal_path_by_basename, path_looks_like_data_file

    ctx = context or {}
    upload_dir = _resolve_upload_dir(ctx)
    context_paths = _collect_candidate_paths(ctx)
    uploaded_entries = _collect_uploaded_file_entries(ctx)
    owner_id = ctx.get("owner_id")

    for key, val in list(params.items()):
        if key in _OUTPUT_PATH_PARAMS:
            continue
        if not isinstance(val, str) or _path_slot_empty(val):
            continue
        if not path_looks_like_data_file(val):
            continue
        try:
            if Path(val.strip()).exists():
                continue
        except OSError:
            pass
        healed, _ = fuzzy_heal_path_by_basename(
            val,
            upload_dir=upload_dir,
            owner_id=owner_id,
            context_paths=context_paths,
            uploaded_entries=uploaded_entries,
            allow_dir=True,
        )
        if healed and healed != val:
            params[key] = healed
            logger.info("[InputValidator] 全局模糊自愈 %s: %s -> %s", key, val, healed)


def _validate_path_params_exist(
    params: MutableMapping[str, Any],
    accepted: Set[str],
    upload_dir: Optional[Path],
    context: Optional[Dict[str, Any]] = None,
) -> Optional[Dict[str, Any]]:
    """物理防呆：输入路径必须存在。"""
    ctx = context or {}
    missing: List[Dict[str, Any]] = []
    for pname in sorted(accepted):
        if pname in _OUTPUT_PATH_PARAMS:
            continue
        raw = params.get(pname)
        if raw is None or _path_slot_empty(raw):
            continue
        if not isinstance(raw, str):
            continue
        resolved, attempted = _resolve_param_on_disk(raw, upload_dir, pname, ctx)
        if resolved:
            if resolved != raw:
                params[pname] = resolved
            continue
        missing.append(
            {
                "param": pname,
                "value": raw,
                "attempted": attempted,
            }
        )

    if not missing:
        return None

    first = missing[0]
    checked = first.get("attempted") or [str(first.get("value", ""))]
    detail = (
        f"参数 `{first['param']}` 指向的路径在容器内不存在: {first['value']!r}。"
        f" 已尝试: {checked}"
    )
    user_msg = (
        "无法在容器内定位上传文件，请检查 Docker 挂载与路径。"
        f"（工具参数 {first['param']}={first['value']}）"
    )
    return {
        "status": "error",
        "error_category": "input_path_not_found",
        "error": detail,
        "message": user_msg,
        "user_message": user_msg,
        "path_validation": missing,
    }


def validate_and_normalize_inputs(
    tool_id: str,
    params: MutableMapping[str, Any],
    tool_func: Optional[Callable] = None,
    context: Optional[Dict[str, Any]] = None,
) -> Tuple[MutableMapping[str, Any], Optional[Dict[str, Any]]]:
    """
    工具执行前统一安检门。

    Returns:
        (normalized_params, error_dict)
        error_dict 非 None 时调用方应直接返回该结构，勿再调用底层工具。
    """
    ctx = context if context is not None else get_tool_input_context()
    if tool_func is None:
        from gibh_agent.core.tool_registry import registry

        tool_func = registry.get_tool(tool_id)

    accepted = _tool_path_param_names(tool_func) if tool_func else set(
        k for k in params if k in _PATH_PARAM_NAMES
    )
    if not accepted:
        return params, None

    _normalize_synonym_keys(params, accepted)

    candidates = _collect_candidate_paths(ctx)
    primary = _pick_primary_path(candidates, tool_id=tool_id, context=ctx)
    _apply_path_fallbacks(params, accepted, primary, tool_id=tool_id)
    _prefer_chain_h5ad_over_upload_dir(params, accepted, ctx, tool_id)

    _fuzzy_heal_path_like_param_values(params, ctx)

    upload_dir = _resolve_upload_dir(ctx)
    err = _validate_path_params_exist(params, accepted, upload_dir, ctx)
    if err:
        err["tool_id"] = tool_id
        from gibh_agent.core.omics_io_registry import is_elastic_skip_tool

        if is_elastic_skip_tool(tool_id):
            err["can_skip"] = True
        logger.error("[InputValidator] %s 路径校验失败: %s", tool_id, err.get("error"))
    return params, err


# 对外别名（任务书命名）
_validate_and_normalize_inputs = validate_and_normalize_inputs
