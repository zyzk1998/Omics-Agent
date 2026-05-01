"""
通用工作流执行器 - The Hands

动态执行工作流，不依赖硬编码的工具逻辑。
使用 ToolRegistry 查找和执行工具。
支持动态参数智能推荐与注入：按函数签名过滤并做类型转换，防止幻觉参数与 TypeError。
"""
import asyncio
import inspect
import json
import os
import logging
import time
import traceback
from typing import Dict, Any, List, Optional, Callable, get_origin, get_args
from pathlib import Path
from datetime import datetime

from .tool_registry import registry
from .tool_stream_log import tool_log_emitter_scope
from .utils import sanitize_for_json
from .file_inspector import (
    path_looks_like_medical_imaging,
    DataAsset,
    OmicsAssetManager,
    to_legacy_resolved_dict,
)

# 隐患 3：类型转换失败时的哨兵，用于丢弃非法推荐参数而非透传导致底层 TypeError
_COERCE_FAILED = object()


class SecurityException(Exception):
    """Raised when file signature verification fails before execution."""

    def __init__(self, message: str, path: Optional[str] = None):
        self.path = path
        super().__init__(message)


# 🔥 CRITICAL: 确保工具模块被加载（触发工具注册）
try:
    from ..tools import load_all_tools
    # 如果工具还未加载，立即加载
    if len(registry._tools) == 0:
        logger.info("🔍 [Executor] 工具注册表为空，正在加载工具...")
        load_all_tools()
        logger.info(f"✅ [Executor] 工具加载完成，已注册 {len(registry._tools)} 个工具")
    try:
        from ..mcp import load_mcp_plugins
        load_mcp_plugins()
    except Exception as _mcp_e:
        logger.warning("⚠️ [Executor] MCP 插件加载失败: %s", _mcp_e)
except Exception as e:
    logger.warning(f"⚠️ [Executor] 工具加载失败: {e}")

logger = logging.getLogger(__name__)


class WorkflowExecutor:
    """
    工作流执行器
    
    职责：
    1. 从 ToolRegistry 查找工具
    2. 验证参数
    3. 执行工具
    4. 处理步骤间的数据流
    5. 生成符合前端格式的执行报告
    """
    
    def __init__(self, output_dir: Optional[str] = None, upload_dir: Optional[str] = None):
        """
        初始化工作流执行器
        
        Args:
            output_dir: 输出目录（如果为 None，将在执行时创建）
            upload_dir: 上传目录（用于解析相对路径）
        """
        self.output_dir = output_dir
        self.upload_dir = upload_dir or os.getenv("UPLOAD_DIR", "/app/uploads")
        self.upload_dir = upload_dir or os.getenv("UPLOAD_DIR", "/app/uploads")
        self.results_dir = os.getenv("RESULTS_DIR", "/app/results")
        self.step_results: Dict[str, Any] = {}  # 存储步骤结果，用于数据流传递
    
    def _resolve_file_path(self, file_path: str) -> str:
        """
        🔥 SYSTEM-WIDE REFACTOR: Smart Path Resolver
        
        按照严格的逻辑顺序解析文件路径，确保不会错误地修改已存在的绝对路径。
        
        检查顺序（严格按照此顺序，不要改变）：
        1. 绝对路径且存在 -> 直接返回（不修改）
        2. Results 目录 -> 尝试在 RESULTS_DIR 中查找
        3. Uploads 目录 -> 尝试在 UPLOAD_DIR 中查找
        4. 当前工作目录 -> 尝试相对路径
        5. 失败 -> 抛出 FileNotFoundError
        
        Args:
            file_path: 原始文件路径（可能是绝对路径、相对路径或文件名）
            
        Returns:
            解析后的绝对路径
            
        Raises:
            FileNotFoundError: 如果所有检查都失败
        """
        if not file_path or file_path in ["<待上传数据>", "<PENDING_UPLOAD>", ""]:
            return file_path

        from ..utils.path_resolver import (
            is_hpc_style_path,
            is_windows_abs_path,
            normalize_duplicate_tail_filename,
        )

        file_path = normalize_duplicate_tail_filename(str(file_path).strip())
        if is_windows_abs_path(file_path) or is_hpc_style_path(file_path):
            logger.info(
                "✅ [Path Resolver] Local/HPC 路径：跳过 UPLOAD_DIR/Results 拼接，原样返回: %s",
                file_path,
            )
            return file_path
        
        original_path = file_path
        attempted_paths = []
        
        # Check 1: Absolute & Exists
        # 🔥 CRITICAL: If absolute and exists, RETURN IMMEDIATELY (Do not touch it!)
        if Path(file_path).is_absolute():
            path_obj = Path(file_path)
            if path_obj.exists():
                resolved = str(path_obj.resolve())
                logger.info(f"✅ [Path Resolver] 绝对路径已存在，直接返回: {original_path} -> {resolved}")
                return resolved
            else:
                attempted_paths.append(str(path_obj.resolve()))
                logger.debug(f"🔍 [Path Resolver] 绝对路径不存在: {file_path}")
        
        # Check 2: Results Directory
        # Construct path = RESULTS_DIR / file_path
        results_dir_path = Path(self.results_dir)
        # Remove leading slash if present
        file_path_clean = file_path.lstrip('/')
        potential_results_path = results_dir_path / file_path_clean
        
        if potential_results_path.exists():
            resolved = str(potential_results_path.resolve())
            logger.info(f"✅ [Path Resolver] 在 Results 目录找到: {original_path} -> {resolved}")
            return resolved
        attempted_paths.append(str(potential_results_path.resolve()))
        
        # Also try with original path (in case it's already relative to results)
        if not Path(file_path).is_absolute():
            potential_results_path2 = results_dir_path / file_path
            if potential_results_path2.exists():
                resolved = str(potential_results_path2.resolve())
                logger.info(f"✅ [Path Resolver] 在 Results 目录找到（原始路径）: {original_path} -> {resolved}")
                return resolved
            attempted_paths.append(str(potential_results_path2.resolve()))
        
        # Check 3: Uploads Directory
        # Construct path = UPLOAD_DIR / file_path
        upload_dir_path = Path(self.upload_dir)
        potential_upload_path = upload_dir_path / file_path_clean
        
        if potential_upload_path.exists():
            resolved = str(potential_upload_path.resolve())
            logger.info(f"✅ [Path Resolver] 在 Uploads 目录找到: {original_path} -> {resolved}")
            return resolved
        attempted_paths.append(str(potential_upload_path.resolve()))
        
        # Try with original relative path
        if not Path(file_path).is_absolute():
            potential_upload_path2 = upload_dir_path / file_path
            if potential_upload_path2.exists():
                resolved = str(potential_upload_path2.resolve())
                logger.info(f"✅ [Path Resolver] 在 Uploads 目录找到（原始路径）: {original_path} -> {resolved}")
                return resolved
            attempted_paths.append(str(potential_upload_path2.resolve()))
        
        # Check 4: Relative to Current Work Dir
        # Try as relative path from current working directory
        if not Path(file_path).is_absolute():
            try:
                cwd_path = Path.cwd() / file_path
                if cwd_path.exists():
                    resolved = str(cwd_path.resolve())
                    logger.info(f"✅ [Path Resolver] 在当前工作目录找到: {original_path} -> {resolved}")
                    return resolved
                attempted_paths.append(str(cwd_path.resolve()))
            except Exception as e:
                logger.debug(f"⚠️ [Path Resolver] 检查当前工作目录失败: {e}")
        
        # Check 5: Filename search (only if it's just a filename)
        filename = Path(file_path).name
        if filename == file_path or '/' not in file_path.replace('\\', '/'):
            # Search in results directory recursively
            logger.debug(f"🔍 [Path Resolver] 搜索文件名: {filename} 在 {self.results_dir}")
            for found_path in results_dir_path.rglob(filename):
                if found_path.is_file():
                    resolved = str(found_path.resolve())
                    logger.info(f"✅ [Path Resolver] 在 Results 目录递归找到: {original_path} -> {resolved}")
                    return resolved
            
            # Search in uploads directory recursively
            logger.debug(f"🔍 [Path Resolver] 搜索文件名: {filename} 在 {self.upload_dir}")
            for found_path in upload_dir_path.rglob(filename):
                if found_path.is_file():
                    resolved = str(found_path.resolve())
                    logger.info(f"✅ [Path Resolver] 在 Uploads 目录递归找到: {original_path} -> {resolved}")
                    return resolved
        
        # Failure: All checks failed
        error_msg = (
            f"Could not resolve path: {original_path}. "
            f"Checked: {attempted_paths[:5]}"  # Limit to first 5 for readability
        )
        logger.error(f"❌ [Path Resolver] {error_msg}")
        raise FileNotFoundError(error_msg)

    def _to_results_url(self, raw_path: Any) -> Any:
        """将工具返回的物理路径规范化为 /results/... URL（用于 SSE 与 state_snapshot 持久化一致性）。"""
        if not isinstance(raw_path, str):
            return raw_path
        p = raw_path.strip()
        if not p:
            return p
        if p.startswith("/results/") or p.startswith("/uploads/") or p.startswith("http://") or p.startswith("https://"):
            return p
        results_prefix = str(self.results_dir).rstrip("/")
        if p.startswith(results_prefix + "/") or p == results_prefix:
            rel = p[len(results_prefix):].lstrip("/")
            return f"/results/{rel}" if rel else "/results"
        if p.startswith("results/"):
            return "/" + p
        if p.startswith("/app/results/"):
            return "/results/" + p[len("/app/results/"):]
        # 其他相对路径按 run 目录兜底（与 server.py 逻辑对齐）
        run_name = os.path.basename(self.output_dir or "").strip()
        if run_name and "/" not in p:
            return f"/results/{run_name}/{p}"
        if p.startswith("/"):
            return p
        return f"/results/{p}"

    def _normalize_result_paths(self, payload: Any) -> Any:
        """递归规范化结果中的图片/下载链接/常见路径字段。"""
        if isinstance(payload, dict):
            out: Dict[str, Any] = {}
            for k, v in payload.items():
                if k in (
                    "path", "plot_path", "output_plot_path", "h5ad_path", "trajectory_data_path",
                    "tmaps_dir", "output_path", "output_file", "file_path", "zip_path",
                    "markdown_file_path",
                ):
                    out[k] = self._to_results_url(v)
                elif k == "download_links" and isinstance(v, list):
                    fixed = []
                    for item in v:
                        if isinstance(item, dict):
                            item = {**item, "path": self._to_results_url(item.get("path"))}
                        fixed.append(item)
                    out[k] = fixed
                elif k == "images" and isinstance(v, list):
                    fixed_imgs = []
                    for item in v:
                        if isinstance(item, dict):
                            fixed_imgs.append({**item, "path": self._to_results_url(item.get("path"))})
                        else:
                            fixed_imgs.append(self._to_results_url(item))
                    out[k] = fixed_imgs
                else:
                    out[k] = self._normalize_result_paths(v)
            return out
        if isinstance(payload, list):
            return [self._normalize_result_paths(x) for x in payload]
        return payload

    @staticmethod
    def _unwrap_optional_annotation(annotation: Any) -> Any:
        """剥离 Optional / Union[..., None]，便于识别 List[float] 等内层类型。"""
        if annotation is inspect.Parameter.empty:
            return annotation
        ann = annotation
        for _ in range(12):
            origin = get_origin(ann)
            args = get_args(ann)
            if not origin or not args or type(None) not in args:
                break
            non_none = [a for a in args if a is not type(None)]
            if len(non_none) != 1:
                break
            ann = non_none[0]
        return ann

    @staticmethod
    def _resolve_callable_for_signature(tool_func: Callable) -> Callable:
        """
        ToolRegistry 的 wrapper 签名为 (*args, **kwargs)，会导致具名参数无法匹配、类型转换被跳过。
        优先沿 __wrapped__ 取原始函数以恢复真实注解。
        """
        f: Callable = tool_func
        seen: set[int] = set()
        for _ in range(8):
            wid = id(f)
            if wid in seen:
                break
            seen.add(wid)
            inner = getattr(f, "__wrapped__", None)
            if inner is None or not callable(inner):
                break
            f = inner
        return f

    @staticmethod
    def _coerce_value(annotation: Any, value: Any) -> Any:
        """根据参数注解将 LLM 输出的值转换为目标类型。转换失败返回 _COERCE_FAILED，由调用方丢弃该参数，避免毒药透传导致底层 TypeError。"""
        if value is None:
            return value
        if annotation is inspect.Parameter.empty:
            return value
        try:
            ann = WorkflowExecutor._unwrap_optional_annotation(annotation)
            origin = get_origin(ann)
            is_list_ann = ann is list or origin is list

            # 前端表单常见：逗号分隔字符串 -> List[float] / List[int] / List[str] / list
            if isinstance(value, str) and is_list_ann:
                items = [x.strip() for x in value.split(",") if x.strip()]
                if not items:
                    return []
                elem_args = get_args(ann) if origin is list else ()
                elem_type = elem_args[0] if elem_args else None
                if elem_type is str:
                    return list(items)
                if elem_type is int:
                    try:
                        return [int(float(x)) for x in items]
                    except (TypeError, ValueError):
                        return _COERCE_FAILED
                if elem_type is float:
                    try:
                        return [float(x) for x in items]
                    except (TypeError, ValueError):
                        return _COERCE_FAILED
                try:
                    return [float(x) for x in items]
                except (TypeError, ValueError):
                    return list(items)

            if ann is float:
                return float(value) if not isinstance(value, (int, float)) else float(value)
            if ann is int:
                return int(value) if not isinstance(value, int) else int(value)
            if ann is bool:
                if isinstance(value, bool):
                    return value
                if isinstance(value, str):
                    return value.strip().lower() in ("1", "true", "yes", "on")
                return bool(value)
            return value
        except (TypeError, ValueError):
            return _COERCE_FAILED

    def _filter_and_coerce_params_to_signature(
        self, params: Dict[str, Any], tool_func: Callable
    ) -> Dict[str, Any]:
        """
        防幻觉参数过滤 + 类型转换：只保留工具函数签名中存在的参数，并对值做类型转换。
        支持 **kwargs：若函数有 VAR_KEYWORD，则未在具名参数中的键也会保留。
        """
        try:
            sig = inspect.signature(self._resolve_callable_for_signature(tool_func))
        except Exception as e:
            logger.debug("无法获取工具签名，跳过过滤: %s", e)
            return params
        allowed = set(sig.parameters.keys())
        has_kwargs = any(
            p.kind == inspect.Parameter.VAR_KEYWORD for p in sig.parameters.values()
        )
        valid_params = {}
        for key, value in params.items():
            if key in sig.parameters:
                param = sig.parameters[key]
                coerced = self._coerce_value(param.annotation, value)
                if coerced is _COERCE_FAILED:
                    logger.warning("丢弃非法推荐参数（类型转换失败）: %s = %s，将使用工具默认值", key, value)
                    continue
                valid_params[key] = coerced
            elif has_kwargs:
                valid_params[key] = value
            else:
                logger.warning("丢弃无效的推荐参数: %s (工具签名中不存在)", key)
        return valid_params

    def execute_step(
        self,
        step_data: Dict[str, Any],
        step_context: Optional[Dict[str, Any]] = None
    ) -> Dict[str, Any]:
        """
        执行单个步骤
        
        Args:
            step_data: 步骤数据，包含 tool_id 和 params
            step_context: 步骤上下文（包含前序步骤的输出等）
        
        Returns:
            步骤执行结果
        """
        step_id = step_data.get("step_id", "unknown")
        tool_id = step_data.get("tool_id")
        params = dict(step_data.get("params", {}))
        # 动态参数推荐与注入：优先使用步骤级 recommended_params，否则使用工作流级（step_context）
        recommended_params = step_data.get("recommended_params") or (step_context or {}).get("recommended_params") or {}
        if not isinstance(recommended_params, dict):
            recommended_params = {}
        if recommended_params:
            for k, v in recommended_params.items():
                if v is not None:
                    params[k] = v
            logger.debug(f"🔄 [Executor] 已合并 recommended_params: {list(recommended_params.keys())}")
        step_name = step_data.get("name", tool_id)
        
        logger.info(f"🔧 执行步骤: {step_id} ({tool_id})")
        
        if not tool_id:
            error_msg = f"步骤 {step_id} 缺少 tool_id"
            logger.error(f"❌ {error_msg}")
            return {
                "status": "error",
                "step_id": step_id,
                "step_name": step_name,
                "error": error_msg,
                "message": error_msg
            }
        
        # 查找工具
        tool_func = registry.get_tool(tool_id)
        if not tool_func:
            error_msg = f"工具 '{tool_id}' 未在注册表中找到"
            logger.error(f"❌ {error_msg}")
            return {
                "status": "error",
                "step_id": step_id,
                "step_name": step_name,
                "error": error_msg,
                "message": error_msg
            }
        
        tool_metadata = registry.get_metadata(tool_id)
        tool_category = tool_metadata.category if tool_metadata else None
        _valid_early = self._tool_accepted_param_names(tool_id, tool_func)
        if (
            "adata_path" in _valid_early
            and self._param_slot_empty(params, "adata_path")
            and "file_path" in params
            and not self._param_slot_empty(params, "file_path")
        ):
            params["adata_path"] = params.pop("file_path")
            logger.info(
                "🔄 [Executor] 参数映射: file_path -> adata_path（反射, 工具: %s）",
                tool_id,
            )
        
        # 验证参数（可选但推荐）
        # 🔥 CRITICAL FIX: 参数验证失败时不应该清空 params，应该保留原始参数
        # 🔥 CRITICAL FIX: 保存关键参数（如 group_column）以防验证后丢失
        critical_params_backup = {}
        tools_requiring_group_column = ["differential_analysis", "metabolomics_plsda", "metabolomics_pathway_enrichment"]
        if tool_id in tools_requiring_group_column and "group_column" in params:
            critical_params_backup["group_column"] = params["group_column"]
        
        try:
            if tool_metadata:
                # 使用 Pydantic schema 验证参数（允许额外字段）
                # 🔥 CRITICAL FIX: 使用 model_validate 并设置 extra='ignore' 来忽略额外字段
                try:
                    validated_params = tool_metadata.args_schema.model_validate(params, strict=False)
                    params = validated_params.model_dump(exclude_unset=False)
                    logger.debug(f"✅ 参数验证通过: {step_id}")
                except Exception as e:
                    # 如果 model_validate 失败，尝试使用 __init__ 但捕获额外字段
                    try:
                        # 只提取模型定义的字段
                        schema_fields = tool_metadata.args_schema.model_fields.keys()
                        filtered_params = {k: v for k, v in params.items() if k in schema_fields}
                        validated_params = tool_metadata.args_schema(**filtered_params)
                        params = validated_params.model_dump(exclude_unset=False)
                        # 保留不在 schema 中的额外参数（如 kwargs）
                        extra_params = {k: v for k, v in params.items() if k not in schema_fields}
                        params.update(extra_params)
                        logger.debug(f"✅ 参数验证通过（保留额外参数）: {step_id}")
                    except Exception as e2:
                        logger.warning(f"⚠️ 参数验证失败（继续执行，保留原始参数）: {e2}")
                        # 保留原始 params
                
                # 🔥 CRITICAL FIX: 恢复关键参数（如果验证后丢失）
                for key, value in critical_params_backup.items():
                    if key not in params:
                        params[key] = value
                        logger.warning(f"⚠️ [Executor] 验证后恢复关键参数 {key}: {value}")
        except Exception as validation_err:
            logger.warning(f"⚠️ 参数验证失败（继续执行，保留原始参数）: {validation_err}")
            # 🔥 CRITICAL: 验证失败时保留原始 params，不因验证失败而中断或清空参数
            # params 保持不变，继续执行
        
        # 🔥 CRITICAL DEBUG: 记录原始参数
        if tool_id in ["metabolomics_plsda", "differential_analysis", "metabolomics_pathway_enrichment"]:
            logger.info(f"🔍 [Executor] {tool_id} 原始参数: {list(params.keys())}")
            if "group_column" in params:
                logger.info(f"✅ [Executor] 原始参数中包含 group_column: {params['group_column']}")
            else:
                logger.warning(f"⚠️ [Executor] 原始参数中缺少 group_column")
        
        # 处理数据流：替换占位符（传递工具类别与当前步骤，供 STED-EC 等注入上一步输出）
        processed_params = self._process_data_flow(
            params,
            step_context,
            tool_category=tool_category,
            current_tool_id=tool_id,
            current_step_id=step_id,
        )
        
        # 🔥 CRITICAL DEBUG: 记录处理后的参数
        if tool_id in ["metabolomics_plsda", "differential_analysis", "metabolomics_pathway_enrichment"]:
            logger.info(f"🔍 [Executor] {tool_id} 处理后参数: {list(processed_params.keys())}")
            if "group_column" in processed_params:
                logger.info(f"✅ [Executor] 处理后参数中包含 group_column: {processed_params['group_column']}")
            else:
                logger.error(f"❌ [Executor] 处理后参数中缺少 group_column！")
        
        # 🔥 CRITICAL FIX: Radiomics/Spatial 等只传入函数签名内参数，避免 file_path/output_dir 等误注入
        if tool_category in ("Radiomics", "Spatial"):
            if "file_path" in processed_params:
                del processed_params["file_path"]
            if "output_dir" in processed_params:
                del processed_params["output_dir"]
            try:
                sig = inspect.signature(tool_func)
                allowed = set(sig.parameters.keys())
                extra = set(processed_params.keys()) - allowed
                if extra:
                    for k in extra:
                        del processed_params[k]
                    logger.debug(f"🔄 [Executor] {tool_category} 工具 {tool_id} 移除非签名参数: {extra}")
            except Exception as e:
                logger.debug("Executor: could not get tool signature for param filter: %s", e)
        # 🔥 CRITICAL FIX: 对于 scRNA-seq 工具，确保移除 file_path 参数（如果存在）
        # 🔥 EXCEPTION: rna_cellranger_count 和 rna_convert_cellranger_to_h5ad 不使用 adata_path
        tools_not_using_adata_path = ["rna_cellranger_count", "rna_convert_cellranger_to_h5ad"]
        if tool_category == "scRNA-seq" and tool_id not in tools_not_using_adata_path:
            if "file_path" in processed_params:
                # 如果已经有 adata_path，移除 file_path
                if "adata_path" in processed_params:
                    del processed_params["file_path"]
                    logger.info(f"🔄 [Executor] 移除多余的 file_path 参数（工具已有 adata_path）")
                else:
                    # 如果没有 adata_path，将 file_path 映射为 adata_path
                    processed_params["adata_path"] = processed_params.pop("file_path")
                    logger.info(f"🔄 [Executor] 参数映射: file_path -> adata_path (工具: {tool_id})")
        
        # 🔥 CRITICAL FIX: 对于 rna_cellranger_count 和 rna_convert_cellranger_to_h5ad，移除 adata_path 参数（如果存在）
        if tool_id in tools_not_using_adata_path and "adata_path" in processed_params:
            logger.warning(f"⚠️ [Executor] 工具 {tool_id} 不接受 adata_path 参数，已移除")
            del processed_params["adata_path"]
        
        # 🔥 CRITICAL FIX: 对于 visualize_volcano，自动注入 diff_results（如果缺失）
        if tool_id == "visualize_volcano":
            # 如果 diff_results 缺失或仍然是占位符，尝试从 differential_analysis 步骤获取
            if "diff_results" not in processed_params or (
                isinstance(processed_params.get("diff_results"), str) and 
                processed_params.get("diff_results", "").startswith("<")
            ):
                # 查找 differential_analysis 步骤的结果
                diff_result = None
                for step_id, step_result in self.step_results.items():
                    if step_id == "differential_analysis" or "differential" in step_id.lower():
                        # step_results 存储的是工具返回的原始结果
                        if isinstance(step_result, dict):
                            diff_result = step_result
                            break
                
                if isinstance(diff_result, dict):
                    processed_params["diff_results"] = diff_result
                    logger.info(f"✅ [Executor] 自动注入 diff_results 到 visualize_volcano")
                else:
                    logger.error(f"❌ [Executor] 无法找到 differential_analysis 结果，visualize_volcano 将失败")
            
            # 过滤掉工具不接受的参数（保留 diff_results）
            allowed_params = {"diff_results", "output_path", "fdr_threshold", "log2fc_threshold"}
            filtered_params = {k: v for k, v in processed_params.items() if k in allowed_params}
            if len(filtered_params) < len(processed_params):
                removed = set(processed_params.keys()) - allowed_params
                logger.warning(f"⚠️ [Executor] 移除 {tool_id} 不接受的参数: {removed}")
            processed_params = filtered_params
        
        # 🔥 CRITICAL REGRESSION FIX: Normalize all path-like parameters before tool execution
        # 🔥 TASK 2: Only resolve paths that are NOT already absolute and existing
        path_params = ["file_path", "adata_path", "h5ad_path", "trajectory_data_path", "image_path", "mask_path", "output_path", "output_file", "fastq_path", "reference_path", "output_h5ad", "features_csv_path", "rad_score_csv_path"]
        for param_name in path_params:
            if param_name in processed_params:
                original_path = processed_params[param_name]
                if original_path and isinstance(original_path, str):
                    # Skip placeholder values and unresolved step refs (e.g. <sted_ec_moscot_trajectory>)
                    if original_path not in ["<待上传数据>", "<PENDING_UPLOAD>", ""] and not (
                        original_path.strip().startswith("<") and original_path.strip().endswith(">")
                    ):
                        # 🔥 CRITICAL: If already absolute and exists, do NOT modify it
                        if Path(original_path).is_absolute() and Path(original_path).exists():
                            logger.debug(f"✅ [Executor] 路径参数 {param_name} 已是绝对路径且存在，不修改: {original_path}")
                            # Keep it as is - do not call _resolve_file_path
                        else:
                            # Only resolve if not absolute or doesn't exist
                            try:
                                resolved_path = self._resolve_file_path(original_path)
                                if resolved_path != original_path:
                                    logger.info(f"🔄 [Executor] 解析路径参数 {param_name}: {original_path} -> {resolved_path}")
                                processed_params[param_name] = resolved_path
                            except FileNotFoundError as e:
                                # If resolution fails, keep original and let tool handle the error
                                logger.warning(f"⚠️ [Executor] 路径解析失败，保留原始路径: {original_path} (错误: {e})")
                                # Keep original_path - tool will handle the error
        
        # 🔥 Security: Verify signature only for user-uploaded files (under upload_dir), not pipeline outputs
        try:
            from ..utils.security import verify_file_signature
            from .security_config import get_signing_public_key
            public_key_b64 = get_signing_public_key()
            upload_dir_path = Path(self.upload_dir or os.getenv("UPLOAD_DIR", "/app/uploads")).resolve()
            if public_key_b64 and upload_dir_path.exists():
                for input_param in ["file_path", "adata_path"]:
                    path_val = processed_params.get(input_param)
                    if path_val and isinstance(path_val, str) and path_val not in ["<待上传数据>", "<PENDING_UPLOAD>", ""]:
                        path_obj = Path(path_val).resolve()
                        if not path_obj.is_file():
                            continue
                        # Only verify files under upload dir (user-uploaded); skip pipeline outputs
                        try:
                            path_obj.relative_to(upload_dir_path)
                        except ValueError:
                            logger.debug("✅ [Executor] 跳过签名校验（非上传文件）: %s", path_val)
                            continue
                        # 仅当存在 .sig 侧车文件时才校验；解压/未签名文件无 .sig 则跳过
                        sig_path = Path(str(path_obj) + ".sig")
                        if not sig_path.exists():
                            logger.debug("✅ [Executor] 跳过签名校验（无 .sig）: %s", path_obj.name)
                            continue
                        if not verify_file_signature(path_obj, public_key_b64):
                            raise SecurityException(
                                "数据完整性校验未通过，无法执行分析。请重新上传已签名的数据。",
                                path=str(path_obj),
                            )
                        logger.debug(f"✅ [Executor] 文件签名校验通过: %s", path_val)
        except SecurityException:
            raise
        except Exception as e:
            logger.debug("Signature verification skipped or failed (non-blocking): %s", e)
        
        # 🔥 ARCHITECTURAL UPGRADE: Phase 3 - Pre-Flight Check & Auto-Correction
        # 对于需要 group_column 的工具，验证列是否存在，如果不存在则使用 semantic_map 自动修正
        tools_requiring_group_column = ["differential_analysis", "metabolomics_plsda", "metabolomics_pathway_enrichment"]
        if tool_id in tools_requiring_group_column and "group_column" in processed_params:
            group_column = processed_params.get("group_column")
            file_path = processed_params.get("file_path")
            
            if file_path and os.path.exists(file_path) and not path_looks_like_medical_imaging(Path(file_path)):
                try:
                    import pandas as pd
                    # Pre-Flight Check: 读取文件检查列是否存在
                    df = pd.read_csv(file_path, index_col=0, nrows=1)  # 只读第一行检查列名
                    
                    # 🔥 修复：先尝试精确匹配，如果失败则尝试模糊匹配（忽略大小写、空格、下划线）
                    group_column_found = group_column in df.columns
                    matched_column = None
                    
                    if not group_column_found:
                        # 尝试模糊匹配：忽略大小写、空格、下划线
                        group_column_normalized = group_column.lower().replace(' ', '').replace('_', '').replace('-', '')
                        for col in df.columns:
                            col_normalized = col.lower().replace(' ', '').replace('_', '').replace('-', '')
                            if col_normalized == group_column_normalized:
                                matched_column = col
                                logger.info(f"🔄 [Executor] 模糊匹配分组列: '{group_column}' -> '{col}'")
                                break
                        
                        if matched_column:
                            group_column = matched_column
                            group_column_found = True
                            processed_params["group_column"] = matched_column
                    
                    if not group_column_found:
                        logger.warning(f"⚠️ [Executor] Pre-Flight Check: 分组列 '{group_column}' 不存在于数据中")
                        
                        # Auto-Correction: 尝试获取 semantic_map
                        semantic_map = None
                        
                        # 方法1: 从 step_context 获取（如果传递了）
                        if step_context:
                            file_metadata = step_context.get("file_metadata")
                            if file_metadata:
                                semantic_map = file_metadata.get("semantic_map", {})
                        
                        # 方法2: 如果未找到，重新读取文件元数据（使用 FileInspector）
                        if not semantic_map:
                            try:
                                from .file_inspector import FileInspector
                                # os 已在文件顶部导入，无需重复导入
                                upload_dir = os.getenv("UPLOAD_DIR", "/app/uploads")
                                inspector = FileInspector(upload_dir)
                                file_metadata = inspector.inspect_file(file_path)
                                if file_metadata.get("status") == "success":
                                    semantic_map = file_metadata.get("semantic_map", {})
                                    logger.info(f"✅ [Executor] 重新读取 file_metadata 获取 semantic_map")
                            except Exception as e:
                                logger.debug(f"⚠️ [Executor] 无法重新读取 file_metadata: {e}")
                        
                        # 如果找到 semantic_map，使用它进行自动修正
                        if semantic_map:
                            group_cols = semantic_map.get("group_cols", [])
                            if len(group_cols) == 1:
                                # 如果 FileInspector 找到了恰好一个分组列，自动使用它
                                auto_group_col = group_cols[0]
                                logger.info(f"✅ [Executor] Auto-Correction: 使用 semantic_map 中的唯一分组列: {auto_group_col}")
                                processed_params["group_column"] = auto_group_col
                            elif len(group_cols) > 1:
                                # 如果有多个，使用第一个
                                auto_group_col = group_cols[0]
                                logger.info(f"✅ [Executor] Auto-Correction: 使用 semantic_map 中的第一个分组列: {auto_group_col}")
                                processed_params["group_column"] = auto_group_col
                            else:
                                logger.error(f"❌ [Executor] semantic_map['group_cols'] 为空，无法自动修正")
                        else:
                            # 回退到启发式检测
                            logger.warning(f"⚠️ [Executor] 未找到 semantic_map，使用启发式检测...")
                            detected_group_column = self._detect_group_column_from_file(file_path)
                            
                            if detected_group_column:
                                logger.info(f"✅ [Executor] 启发式检测到分组列: {detected_group_column}，替换 '{group_column}'")
                                processed_params["group_column"] = detected_group_column
                            else:
                                # 列出所有可用的列名
                                available_cols = list(df.columns)
                                logger.error(f"❌ [Executor] 无法自动检测分组列。可用列: {available_cols[:10]}")
                                # 不立即失败，让工具自己处理错误（工具会返回友好的错误信息）
                except Exception as e:
                    logger.warning(f"⚠️ [Executor] Pre-Flight Check 失败: {e}，继续执行")
            elif file_path and os.path.exists(file_path) and path_looks_like_medical_imaging(Path(file_path)):
                logger.info("⏭️ [Executor] 医学影像文件跳过 group_column 的 pandas Pre-Flight Check")
        
        # 🔥 CRITICAL FIX: 对于需要 group_column 的工具，强制确保参数存在
        # 这是执行前的最后一道防线，确保参数绝对不会丢失
        tools_requiring_group_column = ["differential_analysis", "metabolomics_plsda", "metabolomics_pathway_enrichment"]
        if tool_id in tools_requiring_group_column:
            # 多重检查：从多个来源尝试获取 group_column
            group_column_value = None
            
            # 检查1: processed_params 中是否有
            if "group_column" in processed_params:
                group_column_value = processed_params["group_column"]
                logger.info(f"✅ [Executor] {tool_id} group_column 已存在: {group_column_value}")
            
            # 检查2: 原始 step_data 中是否有
            elif "group_column" in step_data.get("params", {}):
                group_column_value = step_data["params"]["group_column"]
                processed_params["group_column"] = group_column_value
                logger.warning(f"⚠️ [Executor] {tool_id} 从原始步骤数据恢复 group_column: {group_column_value}")
            
            # 检查3: step_context 中的 file_metadata
            elif step_context and step_context.get("file_metadata"):
                semantic_map = step_context["file_metadata"].get("semantic_map", {})
                group_cols = semantic_map.get("group_cols", [])
                if group_cols:
                    group_column_value = group_cols[0]
                    processed_params["group_column"] = group_column_value
                    logger.warning(f"⚠️ [Executor] {tool_id} 从 file_metadata 自动注入 group_column: {group_column_value}")
            
            # 如果所有检查都失败，返回明确错误
            if not group_column_value:
                logger.error(f"❌ [Executor] CRITICAL: {tool_id} 需要 group_column 参数，但所有来源都未找到！")
                logger.error(f"   processed_params: {list(processed_params.keys())}")
                logger.error(f"   原始步骤参数: {list(step_data.get('params', {}).keys())}")
                logger.error(f"   step_context keys: {list(step_context.keys()) if step_context else 'None'}")
                return {
                    "status": "error",
                    "step_id": step_id,
                    "step_name": step_name,
                    "error": f"工具 {tool_id} 需要 group_column 参数，但参数缺失且无法自动获取。请检查 Planner 是否正确生成了 group_column 参数。",
                    "message": f"步骤 {step_name} 执行失败：缺少必需参数 group_column"
                }
        
        # 专家报告（sted_ec_expert_report）等：将 MCP 开关注入 **kwargs，供工具内动态挂载 OpenAI tools
        if step_context and step_context.get("enabled_mcps") is not None:
            try:
                sig = inspect.signature(tool_func)
                has_kw = any(p.kind == inspect.Parameter.VAR_KEYWORD for p in sig.parameters.values())
                if has_kw or "enabled_mcps" in sig.parameters:
                    processed_params["enabled_mcps"] = step_context["enabled_mcps"]
            except Exception:
                pass

        # 防幻觉 + 类型转换：只保留签名内参数并做类型转换（如 resolution="0.8" -> 0.8）
        processed_params = self._filter_and_coerce_params_to_signature(processed_params, tool_func)

        # 执行工具（可选 ContextVar 流式日志，供工具内 emit_tool_log → SSE）
        _tl_sink = None
        if step_context and step_context.get("tool_log_sink") is not None:
            _tl_sink = step_context.get("tool_log_sink")
        elif getattr(self, "_tool_log_sink", None) is not None:
            _tl_sink = self._tool_log_sink
        try:
            # 🔥 CRITICAL DEBUG: 记录所有参数（包括 group_column）
            logger.info(f"🚀 调用工具: {tool_id} with params: {list(processed_params.keys())}")
            if "group_column" in processed_params:
                logger.info(f"✅ [Executor] group_column 参数存在: {processed_params['group_column']}")
            else:
                if tool_id in tools_requiring_group_column:
                    logger.warning(f"⚠️ [Executor] group_column 参数缺失（但已尝试修复）！可用参数: {list(processed_params.keys())}")
            with tool_log_emitter_scope(_tl_sink if callable(_tl_sink) else None):
                if inspect.iscoroutinefunction(tool_func):
                    result = asyncio.run(tool_func(**processed_params))
                else:
                    result = tool_func(**processed_params)
            
            # 确保结果是字典格式
            if not isinstance(result, dict):
                result = {
                    "status": "success",
                    "data": result,
                    "message": f"步骤 {step_name} 执行完成"
                }
            
            # 确保包含 status 字段
            if "status" not in result:
                result["status"] = "success"
            
            # 🔥 根据工具返回的状态记录正确的日志和消息
            tool_status = result.get("status", "success")
            if tool_status == "error":
                error_msg = result.get("error") or result.get("message") or f"步骤 {step_name} 执行失败"
                logger.error(f"❌ 步骤 {step_id} 执行失败: {error_msg}")
                
                # 🔥 全局报错一致性：工具已返回 user_message 时优先使用，否则再走 ErrorFormatter
                from .error_formatter import ErrorFormatter
                tool_user_message = result.get("user_message") or result.get("message")
                if tool_user_message and len(str(tool_user_message).strip()) > 5:
                    user_message = str(tool_user_message).strip()
                    formatted_error = ErrorFormatter.format_error(error_msg, tool_id, step_name)
                    formatted_error["user_message"] = user_message  # 优先保留工具侧用户向文案
                else:
                    formatted_error = ErrorFormatter.format_error(error_msg, tool_id, step_name)
                
                # 存储结果供后续步骤使用（即使失败也存储，用于调试）
                self.step_results[step_id] = result
                return {
                    "status": "error",
                    "step_id": step_id,
                    "step_name": step_name,
                    "tool_id": tool_id,
                    "result": result,
                    "error": error_msg,
                    "message": formatted_error["user_message"],
                    "user_message": formatted_error["user_message"],
                    "error_category": formatted_error["error_category"],
                    "suggestion": formatted_error["suggestion"],
                    "can_skip": formatted_error["can_skip"],
                    "technical_details": formatted_error.get("technical_details", error_msg),  # 保留技术细节
                    "traceback": result.get("traceback", ""),  # 工具返回的 traceback（如有）
                    "debug_info": result.get("traceback", result.get("debug_info", "")),
                }
            else:
                # success 或 skipped
                is_skipped = result.get("status") == "skipped"
                if is_skipped:
                    logger.info(f"⏭️ 步骤 {step_id} 已跳过（依赖未满足）")
                else:
                    logger.info(f"✅ 步骤 {step_id} 执行成功")
                self.step_results[step_id] = result
                msg = result.get("message") or result.get("error") or (f"步骤 {step_name} 已跳过" if is_skipped else f"步骤 {step_name} 执行完成")
                return {
                    "status": result.get("status", "success"),
                    "step_id": step_id,
                    "step_name": step_name,
                    "tool_id": tool_id,
                    "result": result,
                    "message": msg,
                }
        
        except Exception as e:
            error_msg = str(e)
            tb_str = traceback.format_exc()
            logger.error("🔥 [步骤执行崩溃] %s:\n%s", step_name, tb_str)
            
            # 🔥 全局报错一致性：ValueError 等业务异常常含中文说明，优先保留为 user_message
            from .error_formatter import ErrorFormatter
            formatted_error = ErrorFormatter.format_error(error_msg, tool_id, step_name, error_type="exception")
            if isinstance(e, ValueError) and len(error_msg.strip()) >= 15 and any(c in error_msg for c in ["请", "需要", "检查", "上传", "指定"]):
                formatted_error["user_message"] = error_msg
                formatted_error["error_category"] = "data_issue"
            
            return {
                "status": "error",
                "step_id": step_id,
                "step_name": step_name,
                "tool_id": tool_id,
                "error": error_msg,
                "message": formatted_error["user_message"],
                "user_message": formatted_error["user_message"],
                "error_category": formatted_error["error_category"],
                "suggestion": formatted_error["suggestion"],
                "can_skip": formatted_error["can_skip"],
                "technical_details": formatted_error.get("technical_details", str(e)),
                "traceback": tb_str,
                "debug_info": tb_str,
            }
    
    def _build_sted_ec_expert_readonly_context_json(self) -> str:
        """
        从 step_results 组装 STED-EC 只读 JSON（遗留：供独立调用 sted_ec_expert_report 工具时使用）。
        标准流程下专家报告由 BaseAgent._generate_analysis_summary + Orchestrator diagnosis 生成。
        """
        payload: Dict[str, Any] = {
            "steps": [],
            "output_dir": getattr(self, "output_dir", None),
        }
        for sid in sorted(self.step_results.keys()):
            if sid == "sted_ec_expert_report" or not str(sid).startswith("sted_ec_"):
                continue
            raw = self.step_results[sid]
            if not isinstance(raw, dict):
                continue
            entry: Dict[str, Any] = {
                "step_id": sid,
                "status": raw.get("status"),
                "message_excerpt": str(raw.get("message", ""))[:2000],
            }
            summ = raw.get("summary")
            if isinstance(summ, dict):
                try:
                    entry["summary_excerpt"] = json.dumps(summ, ensure_ascii=False, default=str)[:2800]
                except TypeError:
                    entry["summary_excerpt"] = str(summ)[:2000]
            elif isinstance(summ, str):
                entry["summary_excerpt"] = summ[:2000]

            rd = raw.get("report_data") if isinstance(raw.get("report_data"), dict) else {}
            imgs = rd.get("images") or raw.get("images")
            paths: List[Dict[str, str]] = []
            if isinstance(imgs, list):
                for it in imgs[:10]:
                    if isinstance(it, dict):
                        p = it.get("path") or it.get("url") or ""
                        t = str(it.get("title", ""))
                        if p:
                            paths.append({"title": t, "path": str(p)})
                    elif isinstance(it, str):
                        paths.append({"path": it})
            entry["figure_paths"] = paths[:10]
            dls = rd.get("download_links") or raw.get("download_links")
            if isinstance(dls, list) and dls:
                entry["download_links_excerpt"] = [
                    (
                        str(x.get("path", x))
                        if isinstance(x, dict)
                        else str(x)
                    )
                    for x in dls[:5]
                ]
            for k in ("time_key", "cell_type_key", "n_obs", "h5ad_path", "n_cells"):
                if k in raw and raw[k] is not None:
                    entry[k] = raw[k]
            payload["steps"].append(entry)
        return json.dumps(payload, ensure_ascii=False, default=str)

    @staticmethod
    def _extract_ann_data_path_from_tool_result(res: Any) -> Optional[str]:
        """从单步工具返回字典中提取 AnnData 磁盘路径（不依赖具体工具名）。"""
        if not isinstance(res, dict):
            return None
        st = res.get("status")
        if st is not None and st != "success" and st != "skipped":
            return None
        for key in ("adata_path", "h5ad_path", "output_h5ad"):
            v = res.get(key)
            if isinstance(v, str) and v.strip():
                return v.strip()
        op = res.get("output_path")
        if isinstance(op, str) and op.strip().lower().endswith(".h5ad"):
            return op.strip()
        of = res.get("output_file")
        if isinstance(of, str) and of.strip().lower().endswith(".h5ad"):
            return of.strip()
        return None

    def _find_prior_step_ann_data_path(
        self,
        step_context: Optional[Dict[str, Any]],
        current_step_id: Optional[str],
    ) -> Optional[str]:
        """按 steps_order 向前查找最近一步产出的 .h5ad / adata 路径。"""
        if not step_context or not current_step_id:
            return None
        steps_order = step_context.get("steps_order")
        if not isinstance(steps_order, list) or current_step_id not in steps_order:
            return None
        idx = steps_order.index(current_step_id)
        for j in range(idx - 1, -1, -1):
            sid = steps_order[j]
            res = self.step_results.get(sid)
            p = self._extract_ann_data_path_from_tool_result(res)
            if p:
                try:
                    return self._resolve_file_path(p)
                except Exception:
                    return p
        return None

    def _backfill_ann_data_path_params(
        self,
        processed: Dict[str, Any],
        step_context: Optional[Dict[str, Any]],
        tool_category: Optional[str],
        current_tool_id: Optional[str],
        current_step_id: Optional[str],
    ) -> None:
        """Spatial / scRNA-seq / STED_EC：h5ad_path、adata_path 仍为空时，用前序步骤输出或 resolved_input_paths 兜底。"""
        if tool_category not in ("Spatial", "scRNA-seq", "STED_EC"):
            return
        tf = registry.get_tool(current_tool_id) if current_tool_id else None
        valid = self._tool_accepted_param_names(current_tool_id or "", tf)
        prior = self._find_prior_step_ann_data_path(step_context, current_step_id)
        if not prior and step_context:
            for rp in reversed(step_context.get("resolved_input_paths") or []):
                if isinstance(rp, str) and rp.strip().lower().endswith(".h5ad"):
                    try:
                        prior = self._resolve_file_path(rp.strip())
                    except Exception:
                        prior = rp.strip()
                    break
        if not prior:
            return
        for key in ("h5ad_path", "adata_path"):
            if key not in valid:
                continue
            if not self._param_slot_empty(processed, key):
                continue
            processed[key] = prior
            logger.info("🔄 数据流: 兜底写入 %s = %s（前序 AnnData 输出/输入列表）", key, prior)

    def _process_data_flow(
        self,
        params: Dict[str, Any],
        step_context: Optional[Dict[str, Any]] = None,
        tool_category: Optional[str] = None,
        current_tool_id: Optional[str] = None,
        current_step_id: Optional[str] = None,
    ) -> Dict[str, Any]:
        """
        处理数据流：替换占位符（如 <step1_output>）和自动注入文件路径

        Args:
            params: 原始参数
            step_context: 步骤上下文（包含 current_file_path）
            tool_category: 工具类别（用于确定文件路径参数名）
            current_tool_id: 当前步骤的 tool_id（用于 STED-EC 等从预处理步骤注入 time_key/cell_type_key）
            current_step_id: 当前步骤 step_id（供前序 h5ad 链兜底）

        Returns:
            处理后的参数
        """
        processed = {}
        
        # 🔥 CRITICAL FIX: 先备份关键参数，防止在处理过程中丢失
        critical_params_backup = {}
        if "group_column" in params:
            critical_params_backup["group_column"] = params["group_column"]
            logger.debug(f"🔍 [数据流处理] 备份 group_column: {params['group_column']}")
        
        # 🔥 CRITICAL FIX: 先复制所有非占位符参数（包括 group_column），确保不会丢失
        # 这样可以保证即使占位符处理失败，关键参数也不会丢失
        for key, value in params.items():
            # 跳过占位符，稍后处理
            if isinstance(value, str) and value.startswith("<") and value.endswith(">"):
                continue
            # 立即复制非占位符参数（包括 group_column）
            processed[key] = value
            if key == "group_column":
                logger.debug(f"✅ [数据流处理] 已复制 group_column: {value}")
        
        # 🔥 反射：仅当签名中的路径形参在 processed 中仍全部为空时，按优先级注入单一主路径
        current_file_path: Optional[str] = None
        if step_context:
            current_file_path = step_context.get("current_file_path")
            if not current_file_path:
                assets_list = step_context.get("omics_assets") or []
                ctx_asset = self._get_asset_for_tool(tool_category, assets_list)
                if ctx_asset:
                    current_file_path = self._resolve_file_path(
                        ctx_asset.primary_path
                    )
            if not current_file_path:
                om = step_context.get("omics_resolved") or {}
                rf = step_context.get("resolved_input_paths") or []
                current_file_path = self._infer_default_path_for_tool(
                    tool_category, om, rf
                )
        tf_ctx = (
            registry.get_tool(current_tool_id)
            if current_tool_id
            else None
        )
        # 任一步骤参数里仍有 <...> 占位符时，禁止提前注入路径：否则 adata_path 会被填成
        # current_file_path（甚至上传目录），而 spatial_clustering_comparison 等工具使用
        # 「adata_path or h5ad_path」，会优先走错路径，占位符解析出的 h5ad_path 永远不生效。
        has_placeholder_in_params = any(
            isinstance(v, str) and v.startswith("<") and v.endswith(">")
            for v in params.values()
        )
        if (
            current_file_path
            and current_tool_id
            and tf_ctx
            and not has_placeholder_in_params
        ):
            valid_ctx = self._tool_accepted_param_names(current_tool_id, tf_ctx)
            path_keys_in_sig = valid_ctx & self._INJECTABLE_PATH_PARAM_KEYS
            has_any_path = any(
                not self._param_slot_empty(processed, k) for k in path_keys_in_sig
            )
            if not has_any_path and path_keys_in_sig:
                if self._inject_fallback_path_by_signature(
                    processed, current_tool_id, tf_ctx, current_file_path
                ):
                    logger.info(
                        "🔄 数据流: 反射兜底注入路径 = %s", current_file_path
                    )

        # Radiomics：radiomics_pair 资产补全 mask_path（仅当签名接受 mask_path 且槽位仍空）
        if tool_category == "Radiomics" and step_context and current_tool_id and tf_ctx:
            _vr = self._tool_accepted_param_names(current_tool_id, tf_ctx)
            if "mask_path" in _vr and self._param_slot_empty(processed, "mask_path"):
                assets_list = step_context.get("omics_assets") or []
                rad_asset = self._get_asset_for_tool("Radiomics", assets_list)
                if (
                    rad_asset
                    and rad_asset.asset_type == "radiomics_pair"
                    and rad_asset.metadata.get("mask_path")
                ):
                    mp = rad_asset.metadata["mask_path"]
                    processed["mask_path"] = self._resolve_file_path(str(mp))
                    logger.info(
                        "🔄 数据流: 自动注入 mask_path = %s", processed["mask_path"]
                    )

        if current_tool_id and tf_ctx:
            _vm = self._tool_accepted_param_names(current_tool_id, tf_ctx)
            if (
                "adata_path" in _vm
                and self._param_slot_empty(processed, "adata_path")
                and "file_path" in processed
                and not self._param_slot_empty(processed, "file_path")
            ):
                processed["adata_path"] = processed.pop("file_path")
                logger.info("🔄 数据流: 参数映射 file_path -> adata_path（反射）")
        
        # 现在处理占位符（可能会覆盖已复制的值）
        for key, value in params.items():
            if isinstance(value, str) and value.startswith("<") and value.endswith(">"):
                # 占位符，尝试从上下文或步骤结果中获取
                placeholder = value[1:-1]  # 移除 < >
                
                # STED-EC / 模板第一步：h5ad_path 常用 <user_input>，非步骤 ID，须从本次上传路径注入
                if placeholder in ("user_input", "USER_INPUT", "PENDING_UPLOAD", "待上传数据"):
                    if key in (
                        "h5ad_path",
                        "adata_path",
                        "file_path",
                        "trajectory_data_path",
                        "data_path",
                    ):
                        prior = None
                        if step_context:
                            prior = step_context.get("current_file_path")
                            if not prior:
                                for rp in reversed(
                                    step_context.get("resolved_input_paths") or []
                                ):
                                    if not isinstance(rp, str) or not rp.strip():
                                        continue
                                    s = rp.strip()
                                    if s.lower().endswith(".h5ad"):
                                        try:
                                            prior = self._resolve_file_path(s)
                                        except Exception:
                                            prior = s
                                        break
                            if not prior and step_context.get("resolved_input_paths"):
                                rp0 = step_context["resolved_input_paths"][0]
                                if isinstance(rp0, str) and rp0.strip():
                                    try:
                                        prior = self._resolve_file_path(rp0.strip())
                                    except Exception:
                                        prior = rp0.strip()
                        if prior:
                            processed[key] = prior
                            logger.info(
                                "🔄 数据流: %s = <%s> -> %s（用户上传/执行上下文）",
                                key,
                                placeholder,
                                prior,
                            )
                            continue

                # 🔥 CRITICAL FIX: 特殊处理 preprocess_data_output 占位符
                # 用于 PCA、PLS-DA、差异分析等步骤，需要从 preprocess_data 步骤获取输出文件路径
                if placeholder == "preprocess_data_output" or "preprocess" in placeholder.lower():
                    # 查找 preprocess_data 步骤的结果
                    preprocess_result = None
                    for step_id, step_result in self.step_results.items():
                        if step_id == "preprocess_data" or "preprocess" in step_id.lower():
                            # step_results 存储的是工具返回的原始结果
                            if isinstance(step_result, dict):
                                preprocess_result = step_result
                                break
                    
                    if isinstance(preprocess_result, dict):
                        # 提取输出文件路径（优先顺序：output_file, output_path, file_path）
                        output_path = (
                            preprocess_result.get("output_file") or
                            preprocess_result.get("output_path") or
                            preprocess_result.get("file_path")
                        )
                        if output_path:
                            # 🔥 TASK 2: The previous step MUST return an Absolute Path
                            # If it's already absolute and exists, use it directly (Check 1 in _resolve_file_path will handle it)
                            # Only resolve if it's not absolute or doesn't exist
                            if Path(output_path).is_absolute() and Path(output_path).exists():
                                # Already absolute and exists, use directly (no modification)
                                processed[key] = output_path
                                logger.info(f"🔄 数据流: {key} = <{placeholder}> -> {output_path} (绝对路径，直接使用)")
                            else:
                                # Try to resolve (may be relative or non-existent absolute)
                                try:
                                    resolved_path = self._resolve_file_path(output_path)
                                    processed[key] = resolved_path
                                    logger.info(f"🔄 数据流: {key} = <{placeholder}> -> {resolved_path} (已解析)")
                                except FileNotFoundError:
                                    # If resolution fails, use original (let tool handle the error)
                                    processed[key] = output_path
                                    logger.warning(f"⚠️ 数据流: {key} = <{placeholder}> -> {output_path} (解析失败，使用原始路径)")
                            continue
                        else:
                            logger.warning(f"⚠️ preprocess_data 结果中没有找到输出文件路径。可用字段: {list(preprocess_result.keys())}")
                    else:
                        logger.warning(f"⚠️ 未找到 preprocess_data 步骤结果，无法解析 <{placeholder}>。可用步骤: {list(self.step_results.keys())}")
                
                # 🔥 CRITICAL FIX: 特殊处理 differential_analysis_output 占位符
                # 用于 visualize_volcano，需要传递完整的 diff_results 字典
                if placeholder == "differential_analysis_output" or "differential_analysis" in placeholder.lower():
                    # 查找 differential_analysis 步骤的结果
                    diff_result = None
                    for step_id, step_result in self.step_results.items():
                        if step_id == "differential_analysis" or "differential" in step_id.lower():
                            # step_results 存储的是工具返回的原始结果
                            if isinstance(step_result, dict):
                                diff_result = step_result
                                break
                    
                    if isinstance(diff_result, dict):
                        # 对于 visualize_volcano，需要传递完整的 diff_results 字典
                        if key == "diff_results":
                            processed[key] = diff_result
                            logger.info(f"🔄 数据流: {key} = <{placeholder}> -> 完整 diff_results 字典")
                            continue
                        else:
                            # 其他情况，尝试提取特定字段
                            output_path = (
                                diff_result.get("output_file") or
                                diff_result.get("output_path") or
                                diff_result.get("file_path")
                            )
                            if output_path:
                                processed[key] = output_path
                                logger.info(f"🔄 数据流: {key} = <{placeholder}> -> {output_path}")
                                continue
                    else:
                        logger.warning(f"⚠️ 未找到 differential_analysis 步骤结果，无法解析 <{placeholder}>")
                
                # 🔥 特殊处理：从 differential_analysis 结果中提取分组信息
                if placeholder in ["differential_analysis_case_group", "differential_analysis_control_group"]:
                    # 查找 differential_analysis 步骤的结果
                    diff_result = None
                    for step_id, step_result in self.step_results.items():
                        if step_id == "differential_analysis" or "differential" in step_id.lower():
                            # 检查 step_result 的结构（可能是包装在 result 字段中）
                            if isinstance(step_result, dict):
                                actual_result = step_result.get("result", step_result)
                                if isinstance(actual_result, dict):
                                    diff_result = actual_result
                                    break
                    
                    if isinstance(diff_result, dict):
                        # 尝试从多个位置提取分组信息
                        summary = diff_result.get("summary", {})
                        if placeholder == "differential_analysis_case_group":
                            extracted_value = (
                                diff_result.get("case_group") or
                                diff_result.get("group1") or
                                summary.get("case_group") or
                                diff_result.get("experimental_group")
                            )
                        else:  # differential_analysis_control_group
                            extracted_value = (
                                diff_result.get("control_group") or
                                diff_result.get("group2") or
                                summary.get("control_group") or
                                diff_result.get("control_group")
                            )
                        
                        if extracted_value:
                            processed[key] = extracted_value
                            logger.info(f"🔄 数据流: {key} = <{placeholder}> -> {extracted_value}")
                            continue
                        else:
                            logger.warning(f"⚠️ 无法从 differential_analysis 结果中提取 {placeholder}")
                    else:
                        logger.warning(f"⚠️ 未找到 differential_analysis 步骤结果，无法提取 {placeholder}")
                
                # 尝试从 step_results 中获取；若占位符对应步骤被跳过，用前序步骤输出（纯无负担跳过）
                step_result = self.step_results.get(placeholder)
                if not step_result and step_context and isinstance(step_context.get("steps_order"), list):
                    steps_order = step_context["steps_order"]
                    if placeholder in steps_order:
                        idx = steps_order.index(placeholder)
                        for j in range(idx - 1, -1, -1):
                            prev_id = steps_order[j]
                            if prev_id in self.step_results:
                                step_result = self.step_results[prev_id]
                                logger.info(f"🔄 占位符 <{placeholder}> 对应步骤未执行，回退使用前序步骤 {prev_id} 的输出")
                                break
                if step_result is not None:
                    if isinstance(step_result, dict):
                        # 对于 scRNA-seq 工具，优先查找 output_h5ad
                        if tool_category == "scRNA-seq":
                            output_path = (
                                step_result.get("output_h5ad") or  # 🔥 优先使用 output_h5ad
                                step_result.get("output_file") or
                                step_result.get("output_path") or
                                step_result.get("file_path")
                            )
                        else:
                            # 其他工具使用标准字段（含 Radiomics、Spatial：adata_path、h5ad_path、output_path 等）
                            output_path = (
                                step_result.get("adata_path") or
                                step_result.get("h5ad_path") or
                                step_result.get("output_file") or
                                step_result.get("output_path") or
                                step_result.get("file_path") or
                                step_result.get("csv_path") or
                                step_result.get("output_csv_path") or
                                step_result.get("plot_path") or
                                step_result.get("result_path") or
                                step_result.get("preprocessed_file")
                            )
                        # STED-EC：通路富集等步骤的 CSV 必须由「驱动基因」步骤显式给出，避免误用前序 h5ad 路径
                        if (
                            key == "driver_genes_csv"
                            or placeholder == "sted_ec_driver_gene_extraction"
                        ):
                            output_path = (
                                step_result.get("driver_genes_csv")
                                or step_result.get("csv_path")
                                or step_result.get("output_file")
                            )
                            if output_path and not isinstance(
                                output_path, (str, bytes, os.PathLike)
                            ):
                                output_path = None
                        if output_path and isinstance(output_path, (str, bytes, os.PathLike)):
                            processed[key] = output_path
                            logger.info(f"🔄 数据流: {key} = <{placeholder}> -> {output_path}")
                        else:
                            # 路径参数不能赋值为 dict，保留占位符或原值，避免 TypeError
                            if key in (
                                "h5ad_path",
                                "adata_path",
                                "file_path",
                                "trajectory_data_path",
                                "output_path",
                                "image_path",
                                "mask_path",
                                "driver_genes_csv",
                            ):
                                logger.warning(f"⚠️ 占位符 <{placeholder}> 未解析出路径，保留原值")
                                processed[key] = value
                            else:
                                processed[key] = step_result
                    else:
                        if key in (
                            "h5ad_path",
                            "adata_path",
                            "file_path",
                            "trajectory_data_path",
                            "output_path",
                            "image_path",
                            "mask_path",
                            "driver_genes_csv",
                        ):
                            processed[key] = value
                        else:
                            processed[key] = step_result
                elif step_context and placeholder in step_context:
                    processed[key] = step_context[placeholder]
                else:
                    # 占位符未解析，保持原样（可能后续步骤会处理）
                    logger.warning(f"⚠️ 无法解析占位符: {value}")
                    processed[key] = value
            else:
                # 非占位符参数已经在循环开始前复制，这里不需要再次复制
                # 但如果这个 key 不在 processed 中（不应该发生），还是复制一下
                if key not in processed:
                    processed[key] = value
        
        # 🔥 CRITICAL FIX: 强制确保 group_column 等关键参数没有被意外移除
        # 这是最后的保护措施，确保即使前面的逻辑有问题，group_column 也不会丢失
        for key, value in critical_params_backup.items():
            if key not in processed:
                processed[key] = value
                logger.error(f"❌ [数据流处理] CRITICAL: {key} 丢失，已强制恢复: {value}")
            else:
                logger.debug(f"✅ [数据流处理] {key} 已保留: {processed[key]}")
        
        # 双重检查：如果原始 params 中有 group_column，但 processed 中没有，再次恢复
        if "group_column" in params and "group_column" not in processed:
            processed["group_column"] = params["group_column"]
            logger.error(f"❌ [数据流处理] CRITICAL: group_column 在双重检查中丢失，已强制恢复: {params['group_column']}")
        
        # 替换参数字符串中的 <output_dir> 为实际输出目录（如 output_plot_path: "<output_dir>/xxx.png"）
        if step_context and step_context.get("output_dir"):
            out_dir = step_context["output_dir"]
            for k, v in list(processed.items()):
                if isinstance(v, str) and "<output_dir>" in v:
                    processed[k] = v.replace("<output_dir>", out_dir)

        # STED-EC：从数据校验步骤向时间序列标准化、moscot、轨迹可视化 三步注入 time_key / cell_type_key，避免画图步骤只出 1 张图
        if current_tool_id in (
            "sted_ec_time_series_formatting",
            "sted_ec_moscot_trajectory",
            "sted_ec_plot_trajectory",
            "sted_ec_driver_gene_extraction",
        ):
            pre = self.step_results.get("sted_ec_data_validation")
            if isinstance(pre, dict) and pre.get("status") == "success":
                if pre.get("time_key"):
                    processed["time_key"] = pre["time_key"]
                    logger.info(f"🔄 数据流: 从 sted_ec_data_validation 注入 time_key = {pre['time_key']}")
                if pre.get("cell_type_key") is not None:
                    processed["cell_type_key"] = pre["cell_type_key"]
                    logger.info(f"🔄 数据流: 从 sted_ec_data_validation 注入 cell_type_key = {pre['cell_type_key']}")

        self._backfill_ann_data_path_params(
            processed,
            step_context,
            tool_category,
            current_tool_id,
            current_step_id,
        )

        return processed

    def _get_asset_for_tool(
        self,
        tool_category: Optional[str],
        assets: List[DataAsset],
    ) -> Optional[DataAsset]:
        """
        按工具类别匹配结构化资产（与 OmicsAssetManager 分类一致）。
        未命中时返回 None，由 _infer_default_path_for_tool 等回退。
        """
        if not assets:
            return None
        if tool_category in ("scRNA-seq", "Spatial"):
            for a in assets:
                if a.asset_type == "spatial_bundle":
                    return a
            for a in assets:
                if a.asset_type == "10x_bundle":
                    return a
            for a in assets:
                if a.asset_type == "10x_h5_matrix":
                    return a
            for a in assets:
                if a.asset_type == "h5ad_single":
                    return a
            return None
        if tool_category == "Metabolomics":
            for a in assets:
                if a.asset_type == "metabolomics_csv":
                    return a
            return None
        if tool_category == "Radiomics":
            for a in assets:
                if a.asset_type == "radiomics_pair":
                    return a
            for a in assets:
                if a.asset_type == "radiomics_image":
                    return a
            return None
        return None

    # 反射适配器：可注入的路径形参名（与具体 tool_id 无关）
    _INJECTABLE_PATH_PARAM_KEYS = frozenset({
        "file_path",
        "adata_path",
        "data_path",
        "image_path",
        "mask_path",
        "matrix_dir",
        "sequence_or_path",
    })
    _FALLBACK_PATH_PARAM_PRIORITY = (
        "adata_path",
        "matrix_dir",
        "file_path",
        "data_path",
        "image_path",
        "mask_path",
        "sequence_or_path",
    )

    @staticmethod
    def _param_slot_empty(params: Dict[str, Any], key: str) -> bool:
        if key not in params:
            return True
        v = params[key]
        if v is None:
            return True
        if isinstance(v, str) and not v.strip():
            return True
        return False

    def _tool_accepted_param_names(
        self, tool_id: str, tool_func: Optional[Callable]
    ) -> frozenset[str]:
        """
        底层算子真实形参名（跳过 self/cls、*args、**kwargs）。
        registry 包装器签名为 *args/**kwargs，须 unwrap 到原始函数。
        """
        if tool_func is not None:
            try:
                target = inspect.unwrap(tool_func)
                sig = inspect.signature(target)
                names: List[str] = []
                for name, p in sig.parameters.items():
                    if name in ("self", "cls"):
                        continue
                    if p.kind in (
                        inspect.Parameter.VAR_POSITIONAL,
                        inspect.Parameter.VAR_KEYWORD,
                    ):
                        continue
                    names.append(name)
                if names:
                    return frozenset(names)
            except (ValueError, TypeError, AttributeError) as e:
                logger.debug(
                    "签名反射失败 tool_id=%s，回退 schema: %s", tool_id, e
                )
        meta = registry.get_metadata(tool_id)
        if meta is not None and getattr(meta, "args_schema", None) is not None:
            try:
                return frozenset(meta.args_schema.model_fields.keys())
            except Exception as e:
                logger.debug("从 args_schema 取参名失败 tool_id=%s: %s", tool_id, e)
        return frozenset()

    def _maybe_upgrade_radiomics_pair_data_path(
        self,
        params: Dict[str, Any],
        valid: frozenset[str],
        img_resolved: str,
        mask_resolved: Optional[str],
    ) -> bool:
        """
        LLM 仅写入单一路径到 data_path 时，若与配对资产中的原图一致且存在 mask，则升级为逗号拼接。
        不依赖 tool_id，仅依赖 radiomics_pair 语义与签名是否包含 data_path。
        """
        if "data_path" not in valid or not mask_resolved:
            return False
        if self._param_slot_empty(params, "data_path"):
            return False
        raw = params.get("data_path")
        if not isinstance(raw, str):
            return False
        parts = [p.strip() for p in raw.replace(";", ",").split(",") if p.strip()]
        if len(parts) != 1:
            return False
        try:
            one = self._resolve_file_path(parts[0])
        except Exception:
            one = parts[0]
        if one == img_resolved:
            params["data_path"] = f"{img_resolved},{mask_resolved}"
            return True
        return False

    def _inject_paths_from_asset_reflective(
        self,
        params: Dict[str, Any],
        asset: DataAsset,
        tool_id: str,
        tool_func: Optional[Callable],
    ) -> bool:
        """
        基于 DataAsset 类型 + 工具签名的路径注入（开闭原则：无 tool_id 分支表）。
        仅在对应形参槽为空时写入；显式非空值不覆盖。
        """
        valid = self._tool_accepted_param_names(tool_id, tool_func)
        if not valid:
            return False
        did = False
        at = asset.asset_type

        if at == "radiomics_pair":
            img = str(asset.metadata.get("image_path") or asset.primary_path)
            mask = asset.metadata.get("mask_path")
            mask_s = str(mask) if mask else None
            try:
                img_r = self._resolve_file_path(img)
            except Exception:
                img_r = img
            mask_r: Optional[str] = None
            if mask_s:
                try:
                    mask_r = self._resolve_file_path(mask_s)
                except Exception:
                    mask_r = mask_s
            if "data_path" in valid:
                if self._param_slot_empty(params, "data_path"):
                    params["data_path"] = (
                        f"{img_r},{mask_r}" if mask_r else img_r
                    )
                    did = True
                elif self._maybe_upgrade_radiomics_pair_data_path(
                    params, valid, img_r, mask_r
                ):
                    did = True
            if "image_path" in valid and self._param_slot_empty(
                params, "image_path"
            ):
                params["image_path"] = img_r
                did = True
            if (
                "mask_path" in valid
                and self._param_slot_empty(params, "mask_path")
                and mask_r
            ):
                params["mask_path"] = mask_r
                did = True
            return did

        if at == "radiomics_image":
            try:
                p = self._resolve_file_path(asset.primary_path)
            except Exception:
                p = asset.primary_path
            for key in ("image_path", "data_path", "file_path"):
                if key in valid and self._param_slot_empty(params, key):
                    params[key] = p
                    return True
            return did

        if at == "10x_bundle":
            try:
                p = self._resolve_file_path(asset.primary_path)
            except Exception:
                p = asset.primary_path
            for key in ("adata_path", "matrix_dir", "data_path", "file_path"):
                if key in valid and self._param_slot_empty(params, key):
                    params[key] = p
                    return True
            return did

        if at == "spatial_bundle":
            root = asset.metadata.get("visium_root") or asset.primary_path
            try:
                p = self._resolve_file_path(str(root))
            except Exception:
                p = str(root)
            for key in ("adata_path", "matrix_dir", "data_path", "file_path"):
                if key in valid and self._param_slot_empty(params, key):
                    params[key] = p
                    return True
            return did

        if at == "h5ad_single":
            try:
                p = self._resolve_file_path(asset.primary_path)
            except Exception:
                p = asset.primary_path
            for key in ("adata_path", "data_path", "file_path"):
                if key in valid and self._param_slot_empty(params, key):
                    params[key] = p
                    return True
            return did

        if at == "10x_h5_matrix":
            try:
                p = self._resolve_file_path(asset.primary_path)
            except Exception:
                p = asset.primary_path
            for key in ("adata_path", "data_path", "file_path"):
                if key in valid and self._param_slot_empty(params, key):
                    params[key] = p
                    return True
            return did

        if at in ("metabolomics_csv", "protein_fasta"):
            try:
                p = self._resolve_file_path(asset.primary_path)
            except Exception:
                p = asset.primary_path
            for key in ("file_path", "data_path", "sequence_or_path"):
                if key in valid and self._param_slot_empty(params, key):
                    params[key] = p
                    return True
            return did

        if at in (
            "generic_unknown",
            "unknown_file",
            "protein_structure",
            "document",
            "plain_text",
        ):
            try:
                p = self._resolve_file_path(asset.primary_path)
            except Exception:
                p = asset.primary_path
            for key in (
                "file_path",
                "data_path",
                "adata_path",
                "image_path",
                "sequence_or_path",
            ):
                if key in valid and self._param_slot_empty(params, key):
                    params[key] = p
                    return True
            return did

        return did

    def _inject_fallback_path_by_signature(
        self,
        params: Dict[str, Any],
        tool_id: str,
        tool_func: Optional[Callable],
        inject_path: str,
    ) -> bool:
        """无匹配 DataAsset 或资产未写入任何形参时，按签名优先级注入单一主路径。"""
        if not inject_path:
            return False
        valid = self._tool_accepted_param_names(tool_id, tool_func)
        try:
            resolved = self._resolve_file_path(str(inject_path))
        except Exception:
            resolved = str(inject_path)
        for key in self._FALLBACK_PATH_PARAM_PRIORITY:
            if key in valid and self._param_slot_empty(params, key):
                params[key] = resolved
                return True
        return False

    def _seed_workflow_current_path(
        self,
        assets: List[DataAsset],
        omics_resolved: Dict[str, List[str]],
        resolved_file_paths: List[str],
    ) -> Optional[str]:
        """工作流链种子路径：优先 scRNA 资产（10x 目录 / h5ad），否则沿用桶推断。"""
        rna_asset = self._get_asset_for_tool("scRNA-seq", assets)
        if rna_asset:
            return self._resolve_file_path(rna_asset.primary_path)
        return self._infer_default_path_for_tool(
            None, omics_resolved, resolved_file_paths
        )
    
    def _infer_default_path_for_tool(
        self,
        tool_category: Optional[str],
        omics_resolved: Dict[str, List[str]],
        resolved_file_paths: List[str],
    ) -> Optional[str]:
        """
        根据 resolve_omics_paths 结果与工具类别推断单一路径（避免 file_paths[0] 或逗号串误指子文件）。
        """
        if not omics_resolved:
            omics_resolved = {
                "tables": [],
                "h5ad": [],
                "10x_mtx": [],
                "images": [],
                "masks": [],
                "unknown": [],
            }

        def _pick(bucket: List[str]) -> Optional[str]:
            if not bucket:
                return None
            return self._resolve_file_path(bucket[0])

        tx = omics_resolved.get("10x_mtx") or []
        ha = omics_resolved.get("h5ad") or []
        tb = omics_resolved.get("tables") or []
        im = omics_resolved.get("images") or []

        if tool_category == "scRNA-seq":
            for b in (tx, ha, tb):
                p = _pick(b)
                if p:
                    return p
        elif tool_category == "Radiomics":
            for b in (im, tx, ha, tb):
                p = _pick(b)
                if p:
                    return p
        elif tool_category is None:
            # 工作流启动时的通用默认：优先 10x 目录，再 h5ad / 表格 / 影像
            for b in (tx, ha, tb, im):
                p = _pick(b)
                if p:
                    return p
        else:
            for b in (tb, ha, tx, im):
                p = _pick(b)
                if p:
                    return p
        if resolved_file_paths:
            return self._resolve_file_path(resolved_file_paths[0])
        return None

    def _maybe_unwrap_composite_path_params(
        self,
        params: Dict[str, Any],
        tool_category: Optional[str],
        omics_resolved: Dict[str, List[str]],
        resolved_file_paths: List[str],
    ) -> None:
        """将 adata_path/file_path/image_path 中的逗号/分号拼接串替换为推断出的单路径。"""
        keys_to_check: List[str] = []
        if tool_category == "scRNA-seq":
            keys_to_check = ["adata_path", "file_path"]
        elif tool_category == "Radiomics":
            keys_to_check = ["image_path", "mask_path", "file_path"]
        else:
            keys_to_check = ["file_path", "adata_path"]
        inferred = self._infer_default_path_for_tool(
            tool_category, omics_resolved, resolved_file_paths
        )
        if not inferred and tool_category != "Radiomics":
            return
        for k in keys_to_check:
            v = params.get(k)
            if not isinstance(v, str):
                continue
            s = v.strip()
            if "," not in s and ";" not in s:
                continue
            if k == "mask_path":
                masks = omics_resolved.get("masks") or []
                rep = self._resolve_file_path(masks[0]) if masks else inferred
            else:
                rep = inferred
            if rep:
                params[k] = rep
                logger.info(
                    "🔄 [Executor] %s 为逗号/分号拼接，已规范为单一路径: %s",
                    k,
                    rep,
                )

    def _is_image_file(self, file_path: str) -> bool:
        """
        检查文件路径是否为图片文件
        
        Args:
            file_path: 文件路径
        
        Returns:
            如果是图片文件返回 True，否则返回 False
        """
        if not file_path:
            return False
        
        # 支持的图片文件扩展名
        image_extensions = {'.png', '.jpg', '.jpeg', '.gif', '.bmp', '.svg', '.pdf', '.tiff', '.tif'}
        
        # 获取文件扩展名（不区分大小写）
        file_ext = os.path.splitext(file_path.lower())[1]
        
        return file_ext in image_extensions
    
    def _is_10x_format(self, file_path: str) -> bool:
        """
        检测文件路径是否为10x Genomics格式（包含matrix.mtx的目录）
        
        Args:
            file_path: 文件路径（可能是目录或文件）
        
        Returns:
            如果是10x格式返回 True，否则返回 False
        """
        if not file_path:
            return False
        
        path_obj = Path(file_path)
        
        # 如果是目录，检查是否包含10x格式文件
        if path_obj.is_dir():
            try:
                dir_contents = os.listdir(path_obj)
                # 检查是否包含matrix.mtx（压缩或未压缩）
                has_matrix = any(f in dir_contents for f in ['matrix.mtx', 'matrix.mtx.gz'])
                # 检查是否包含barcodes或features文件
                has_barcodes = any('barcodes' in f for f in dir_contents)
                has_features = any(f in dir_contents for f in ['features.tsv', 'features.tsv.gz', 'genes.tsv', 'genes.tsv.gz'])
                
                if has_matrix and (has_barcodes or has_features):
                    logger.info(f"✅ [10x检测] 检测到10x格式目录: {file_path}")
                    return True
                
                # 递归搜索子目录（最多搜索2层）
                for root, dirs, files in os.walk(path_obj):
                    depth = root.replace(str(path_obj), '').count(os.sep)
                    if depth > 2:  # 限制搜索深度
                        break
                    if any(f in files for f in ['matrix.mtx', 'matrix.mtx.gz']):
                        logger.info(f"✅ [10x检测] 在子目录检测到10x格式: {root}")
                        return True
            except Exception as e:
                logger.debug(f"⚠️ [10x检测] 检查目录失败: {e}")
        
        # 如果是.h5ad文件，不是10x格式（已经是处理后的格式）
        if file_path.endswith('.h5ad'):
            return False
        
        return False
    
    def _detect_group_column_from_file(self, file_path: str) -> Optional[str]:
        """
        从文件中自动检测分组列
        
        使用启发式方法：
        1. 优先检查包含分组关键词的非数值列
        2. 检查唯一值 <= 5 的数值列
        3. 检查唯一值 <= 5 的非数值列
        
        Args:
            file_path: 文件路径
        
        Returns:
            检测到的分组列名，如果未找到返回 None
        """
        try:
            if path_looks_like_medical_imaging(Path(file_path)):
                logger.debug("⏭️ [Executor] 医学影像文件跳过 _detect_group_column_from_file")
                return None
            import pandas as pd
            
            # 读取文件（采样读取，避免大文件问题）
            df = pd.read_csv(file_path, index_col=0, nrows=1000)
            
            # 优先级关键词列表（支持部分匹配，如"Muscle loss"会匹配包含"Muscle"或"Loss"的列）
            priority_keywords = ['Diet', 'diet', 'Group', 'group', 'Condition', 'condition', 
                                'Treatment', 'treatment', 'Class', 'class', 'Category', 'category',
                                'Type', 'type', 'Label', 'label', 'Status', 'status',
                                'Muscle', 'muscle', 'Loss', 'loss', 'MuscleLoss', 'muscleloss',
                                'Muscle_loss', 'muscle_loss']  # 添加用户数据相关的关键词
            
            # 识别列类型
            metadata_cols = []
            numeric_cols = []
            
            for col in df.columns:
                if pd.api.types.is_numeric_dtype(df[col]):
                    numeric_cols.append(col)
                else:
                    metadata_cols.append(col)
            
            # 方法1: 检查非数值列（metadata_cols），优先关键词匹配
            # 🔥 改进：支持部分匹配，如"Muscle loss"会匹配"MuscleLoss"或"Muscle_loss"
            for col in metadata_cols:
                col_lower = col.lower().replace(' ', '').replace('_', '').replace('-', '')
                # 检查是否包含任何关键词（忽略大小写、空格、下划线）
                if any(keyword.lower().replace(' ', '').replace('_', '').replace('-', '') in col_lower 
                       for keyword in priority_keywords):
                    unique_count = df[col].nunique()
                    if 2 <= unique_count <= 20:  # 合理的分组数量
                        logger.info(f"✅ [Executor] 检测到分组列（关键词匹配）: {col} ({unique_count} 个唯一值)")
                        return col
            
            # 方法2: 检查数值列，如果唯一值 <= 5，当作分类变量
            for col in numeric_cols:
                unique_count = df[col].nunique()
                if 2 <= unique_count <= 5:
                    # 检查是否包含关键词
                    if any(keyword in col for keyword in priority_keywords):
                        logger.info(f"✅ [Executor] 检测到分组列（数值型关键词匹配）: {col} ({unique_count} 个唯一值)")
                        return col
                    # 或者检查是否是二元变量（0/1）
                    unique_values = sorted(df[col].dropna().unique().tolist())
                    if unique_values == [0, 1] or unique_values == [1, 0]:
                        logger.info(f"✅ [Executor] 检测到分组列（二元数值型）: {col}")
                        return col
            
            # 方法3: 检查非数值列，如果唯一值 <= 5
            for col in metadata_cols:
                unique_count = df[col].nunique()
                if 2 <= unique_count <= 5:
                    logger.info(f"✅ [Executor] 检测到分组列（非数值型，唯一值 <= 5）: {col} ({unique_count} 个唯一值)")
                    return col
            
            # 方法4: 如果都没有，返回第一个 metadata_cols（如果有）
            if metadata_cols:
                first_meta_col = metadata_cols[0]
                unique_count = df[first_meta_col].nunique()
                if 2 <= unique_count <= 20:
                    logger.info(f"✅ [Executor] 检测到分组列（默认 metadata_cols）: {first_meta_col} ({unique_count} 个唯一值)")
                    return first_meta_col
            
            logger.warning("⚠️ [Executor] 未检测到合适的分组列")
            return None
            
        except Exception as e:
            logger.error(f"❌ [Executor] 检测分组列失败: {e}", exc_info=True)
            return None
    
    def execute_workflow(
        self,
        workflow_data: Dict[str, Any],
        file_paths: List[str] = None,
        output_dir: Optional[str] = None,
        agent: Optional[Any] = None,  # 可选的 Agent 实例，用于生成诊断
        enabled_mcps: Optional[List[Any]] = None,
        tool_log_sink: Optional[Callable[[str, str], None]] = None,
    ) -> Dict[str, Any]:
        """
        执行整个工作流
        
        Args:
            workflow_data: 工作流配置（包含 workflow_name 和 steps）
            file_paths: 输入文件路径列表
            output_dir: 输出目录（如果为 None，将自动创建）
            tool_log_sink: 可选 (content, state) -> None，供工具内 emit_tool_log 透传过程日志；未传则忽略
        
        Returns:
            执行报告（符合前端 analysis_report 格式）
        """
        workflow_name = workflow_data.get("workflow_name", "Unknown Workflow")
        steps = workflow_data.get("steps", [])
        
        logger.info("=" * 80)
        logger.info(f"🚀 开始执行工作流: {workflow_name}")
        logger.info(f"📋 步骤数: {len(steps)}")
        logger.info("=" * 80)
        
        # 设置输出目录：与前端/API 一致，统一使用 RESULTS_DIR 下的可写绝对路径，避免相对路径导致权限或 cwd 依赖问题
        if output_dir is None:
            timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
            base = os.path.abspath(self.results_dir)
            output_dir = os.path.join(base, f"run_{timestamp}")
        
        output_path = Path(output_dir)
        try:
            output_path.mkdir(parents=True, exist_ok=True)
        except OSError as e:
            logger.warning("⚠️ [Executor] 默认结果目录不可写 %s，尝试使用当前目录下 results: %s", output_dir, e)
            output_dir = os.path.join(os.getcwd(), "results", f"run_{timestamp}")
            output_path = Path(output_dir)
            output_path.mkdir(parents=True, exist_ok=True)
        self.output_dir = str(output_path.resolve())
        
        logger.info(f"📂 输出目录: {self.output_dir}")

        self._enabled_mcps: List[Any] = list(enabled_mcps or [])
        self._tool_log_sink: Optional[Callable[[str, str], None]] = tool_log_sink
        
        # 初始化步骤结果列表
        steps_details = []
        steps_results = []
        
        # 🔥 上下文链：跟踪当前文件路径，用于自动传递给下一个步骤
        # 🔥 CRITICAL REGRESSION FIX: Resolve file paths to absolute paths
        resolved_file_paths: List[str] = []
        if file_paths:
            for fp in file_paths:
                if fp and isinstance(fp, str):
                    resolved = self._resolve_file_path(fp)
                    resolved_file_paths.append(resolved)
                    if resolved != fp:
                        logger.info(f"🔄 [Executor] 解析输入文件路径: {fp} -> {resolved}")
        omics_resolved: Dict[str, List[str]] = {
            "tables": [],
            "h5ad": [],
            "10x_mtx": [],
            "images": [],
            "masks": [],
            "unknown": [],
        }
        omics_asset_manager = OmicsAssetManager()
        assets: List[DataAsset] = []
        if resolved_file_paths:
            assets = omics_asset_manager.classify_assets(resolved_file_paths)
            omics_resolved = to_legacy_resolved_dict(assets)
            logger.info(
                "🔄 [Executor] OmicsAssetManager: %s 个资产 | 10x_mtx=%s h5ad=%s tables=%s images=%s",
                len(assets),
                bool(omics_resolved.get("10x_mtx")),
                bool(omics_resolved.get("h5ad")),
                bool(omics_resolved.get("tables")),
                bool(omics_resolved.get("images")),
            )
        current_file_path = self._seed_workflow_current_path(
            assets, omics_resolved, resolved_file_paths
        )
        
        # 🔥 TASK 1 FIX: 检测输入文件类型，如果是10x格式，标记需要跳过cellranger步骤
        is_10x_input = bool(
            next((a for a in assets if a.asset_type == "10x_bundle"), None)
        )
        if not is_10x_input and current_file_path:
            is_10x_input = self._is_10x_format(current_file_path)
        if is_10x_input:
            logger.info(
                "✅ [Executor] 检测到输入为 10x 捆绑或 10x 目录，将自动跳过 cellranger / 按需处理 convert 步骤"
            )
        
        # 执行每个步骤
        for i, step in enumerate(steps, 1):
            step_id = step.get("step_id", f"step{i}")
            step_name = step.get("name", step.get("step_name", step_id))
            tool_id = step.get("tool_id", step_id)
            params = step.get("params", {})
            
            # 🔥 前端勾选状态：用户取消勾选的步骤不执行，仅持久化到快照
            if step.get("enabled") is False or step.get("selected") is False:
                logger.info(f"⏭️ [Executor] 跳过步骤 {step_name} ({tool_id})：用户未勾选")
                continue
            
            # 🔥 CRITICAL FIX: 完全移除 visualize_pca 步骤（不在流程中显示）
            if tool_id == "visualize_pca" or step_id == "visualize_pca":
                logger.warning(f"⚠️ [Executor] 完全移除 visualize_pca 步骤（pca_analysis 已包含可视化）")
                continue  # 直接跳过，不添加到步骤详情中
            
            # 🔥 TASK 1 FIX: 如果输入是10x格式，自动跳过cellranger和convert步骤
            if is_10x_input:
                if tool_id in ["rna_cellranger_count", "rna_convert_cellranger_to_h5ad"]:
                    logger.info(f"⏭️ [Executor] 跳过步骤 {step_name} ({tool_id})：输入文件已经是10x格式，无需执行此步骤")
                    # 对于convert步骤，我们需要直接转换10x格式为h5ad
                    if tool_id == "rna_convert_cellranger_to_h5ad":
                        # 直接使用10x目录转换为h5ad
                        from ..tools.rna.upstream import convert_cellranger_to_h5ad
                        output_h5ad = os.path.join(self.output_dir, "converted_from_10x.h5ad")
                        convert_result = convert_cellranger_to_h5ad(
                            cellranger_matrix_dir=current_file_path,
                            output_h5ad_path=output_h5ad
                        )
                        if convert_result.get("status") == "success":
                            current_file_path = convert_result.get("output_path")
                            logger.info(f"✅ [Executor] 10x格式已转换为h5ad: {current_file_path}")
                            # 记录步骤结果
                            self.step_results[step_id] = {
                                "status": "success",
                                "result": convert_result,
                                "message": "10x格式已自动转换为h5ad格式"
                            }
                            # 构建步骤详情并添加到列表
                            step_detail = {
                                "step_id": step_id,
                                "tool_id": tool_id,
                                "name": step_name,
                                "status": "success",
                                "summary": "10x格式已自动转换为h5ad格式",
                                "duration": 0,
                                "step_result": {
                                    "step_name": step_name,
                                    "status": "success",
                                    "logs": "10x格式已自动转换为h5ad格式",
                                    "data": convert_result
                                }
                            }
                            step_detail = sanitize_for_json(step_detail)
                            steps_details.append(step_detail)
                            steps_results.append(step_detail["step_result"])
                            continue
                        else:
                            logger.error(f"❌ [Executor] 10x格式转换失败: {convert_result.get('error')}")
                            # 转换失败，继续执行原步骤（可能会失败，但至少尝试了）
                    else:
                        # cellranger步骤直接跳过
                        self.step_results[step_id] = {
                            "status": "skipped",
                            "result": {"message": "输入文件已经是10x格式，无需执行Cell Ranger计数"},
                            "message": "步骤已跳过：输入文件已经是10x格式"
                        }
                        # 构建步骤详情并添加到列表
                        step_detail = {
                            "step_id": step_id,
                            "tool_id": tool_id,
                            "name": step_name,
                            "status": "skipped",
                            "summary": "步骤已跳过：输入文件已经是10x格式",
                            "duration": 0,
                            "step_result": {
                                "step_name": step_name,
                                "status": "skipped",
                                "logs": "输入文件已经是10x格式，无需执行Cell Ranger计数",
                                "data": {"message": "输入文件已经是10x格式，无需执行Cell Ranger计数"}
                            }
                        }
                        step_detail = sanitize_for_json(step_detail)
                        steps_details.append(step_detail)
                        steps_results.append(step_detail["step_result"])
                        continue
            
            logger.info(f"\n{'=' * 80}")
            logger.info(f"📌 步骤 {i}/{len(steps)}: {step_name} ({step_id})")
            logger.info(f"{'=' * 80}")
            
            # 🔥 智能参数映射：工具类别用于资产选取与逗号串展开，路径注入走反射适配器
            tool_metadata = registry.get_metadata(tool_id)
            tool_category = tool_metadata.category if tool_metadata else None
            tool_func_exec = registry.get_tool(tool_id)

            self._maybe_unwrap_composite_path_params(
                params, tool_category, omics_resolved, resolved_file_paths
            )
            
            # 🔥 CRITICAL FIX: 检查是否有占位符需要处理
            # 占位符（如 <preprocess_data_output>）必须优先于自动注入的 current_file_path
            has_placeholder = any(
                isinstance(v, str) and v.startswith("<") and v.endswith(">")
                for v in params.values()
            )
            
            # 🔥 CRITICAL FIX: 如果存在占位符，DO NOT 自动注入文件路径
            # 占位符会在 execute_step 内部的 _process_data_flow 中解析
            # 只有在没有占位符且参数缺失时，才自动注入
            if not has_placeholder:
                step_asset = self._get_asset_for_tool(tool_category, assets)
                injected_from_asset = False
                if step_asset and tool_func_exec:
                    injected_from_asset = self._inject_paths_from_asset_reflective(
                        params, step_asset, tool_id, tool_func_exec
                    )
                    if injected_from_asset:
                        _pk = self._INJECTABLE_PATH_PARAM_KEYS
                        logger.info(
                            "🔄 反射注入（OmicsAsset %s）: %s",
                            step_asset.asset_type,
                            {k: params[k] for k in params if k in _pk},
                        )
                if not injected_from_asset and tool_func_exec:
                    inject_path = current_file_path or self._infer_default_path_for_tool(
                        tool_category, omics_resolved, resolved_file_paths
                    )
                    if inject_path and self._inject_fallback_path_by_signature(
                        params, tool_id, tool_func_exec, inject_path
                    ):
                        logger.info("🔄 反射兜底注入路径: %s", inject_path)
            else:
                # 有占位符，记录日志但不自动注入
                placeholder_keys = [k for k, v in params.items() if isinstance(v, str) and v.startswith("<") and v.endswith(">")]
                logger.info(f"🔄 检测到占位符 {placeholder_keys}，跳过自动注入，等待占位符解析")
            
            # 构建步骤上下文（含 steps_order 供占位符回退：跳过某步时用前序步骤输出）
            steps_order = [s.get("step_id") or s.get("tool_id") or s.get("id") for s in steps if (s.get("step_id") or s.get("tool_id") or s.get("id"))]
            step_context = {
                "file_paths": file_paths or [],
                "output_dir": self.output_dir,
                "workflow_name": workflow_name,
                "current_file_path": current_file_path,
                "omics_resolved": omics_resolved,
                "resolved_input_paths": resolved_file_paths,
                "omics_assets": assets,
                "recommended_params": workflow_data.get("recommended_params"),
                "steps_order": steps_order,
                "enabled_mcps": getattr(self, "_enabled_mcps", None) or [],
                "tool_log_sink": getattr(self, "_tool_log_sink", None),
            }
            
            # 🔥 参数映射：签名需要 adata_path 且槽位空时，将 file_path 迁入（常见于单细胞工具）
            if tool_func_exec:
                _valid = self._tool_accepted_param_names(tool_id, tool_func_exec)
                if (
                    "adata_path" in _valid
                    and self._param_slot_empty(params, "adata_path")
                    and "file_path" in params
                    and not self._param_slot_empty(params, "file_path")
                ):
                    params["adata_path"] = params.pop("file_path")
                    logger.info("🔄 参数映射: file_path -> adata_path（反射）")
            
            # 🔥 修复：自动检测并替换 group_column 参数
            # 如果工具需要 group_column，但指定的列不存在，尝试自动检测
            if "group_column" in params and current_file_path:
                specified_group_col = params["group_column"]
                detected_group_col = self._detect_group_column_from_file(current_file_path)
                
                if detected_group_col and not path_looks_like_medical_imaging(Path(current_file_path)):
                    # 检查指定的列是否存在
                    try:
                        import pandas as pd
                        df = pd.read_csv(current_file_path, nrows=10)
                        if specified_group_col not in df.columns:
                            # 指定的列不存在，使用检测到的列
                            logger.warning(f"⚠️ [Executor] 指定的分组列 '{specified_group_col}' 不存在，自动使用检测到的分组列: '{detected_group_col}'")
                            params["group_column"] = detected_group_col
                        else:
                            # 指定的列存在，使用指定的列
                            logger.info(f"✅ [Executor] 使用指定的分组列: '{specified_group_col}'")
                    except Exception as e:
                        logger.warning(f"⚠️ [Executor] 无法验证分组列，使用检测到的分组列: '{detected_group_col}'")
                        params["group_column"] = detected_group_col
                else:
                    logger.warning(f"⚠️ [Executor] 无法自动检测分组列，使用指定的分组列: '{specified_group_col}'")
            
            # 如果工具需要 output_dir，也自动注入
            if "output_dir" not in params and self.output_dir:
                params["output_dir"] = self.output_dir
            
            # 更新步骤的 params
            step["params"] = params
            
            # 执行步骤（内部会调用 _process_data_flow 处理占位符）
            tool_start = time.time()
            step_result = self.execute_step(step, step_context)
            logger.info("[Profiler] 工具[%s] 执行耗时: %.2fs", step_name, time.time() - tool_start)
            
            # 🔥 TASK 3 FIX: 更新 current_file_path 供下一个步骤使用
            # 对于 scRNA-seq 工具，所有产生 output_h5ad 的步骤都应该更新 current_file_path
            result_data = step_result.get("result", {})
            if isinstance(result_data, dict):
                tool_id = step.get("tool_id", "")
                tool_metadata = registry.get_metadata(tool_id)
                tool_category = tool_metadata.category if tool_metadata else None
                
                # scRNA-seq / Spatial：任一步产出 AnnData 磁盘路径则更新链式 current_file_path
                if tool_category in ("scRNA-seq", "Spatial"):
                    next_file_path = self._extract_ann_data_path_from_tool_result(result_data)
                    if not next_file_path:
                        next_file_path = (
                            result_data.get("output_file") or
                            result_data.get("output_path") or
                            result_data.get("file_path")
                        )
                    if next_file_path:
                        try:
                            resolved_path = self._resolve_file_path(str(next_file_path))
                        except Exception:
                            resolved_path = str(next_file_path)
                        current_file_path = resolved_path
                        if os.path.exists(resolved_path):
                            logger.info(
                                "✅ [Executor] 更新当前文件路径（%s）: %s",
                                tool_category,
                                current_file_path,
                            )
                        else:
                            logger.warning(
                                "⚠️ [Executor] 输出路径不存在，但会使用: %s",
                                resolved_path,
                            )
                    else:
                        logger.debug(
                            "🔍 [Executor] 步骤 %s 未返回 AnnData 路径，保持当前文件路径",
                            tool_id,
                        )
                else:
                    # 其他工具（如代谢组学）只在 preprocess_data 步骤更新
                    if tool_id == "preprocess_data" or "preprocess" in tool_id.lower():
                        next_file_path = (
                            result_data.get("output_file") or
                            result_data.get("output_path") or
                            result_data.get("file_path") or
                            result_data.get("preprocessed_file")
                        )
                        
                        if next_file_path:
                            resolved_path = self._resolve_file_path(next_file_path)
                            current_file_path = resolved_path
                            if os.path.exists(resolved_path):
                                logger.info(f"✅ [Executor] 更新当前文件路径（来自 {tool_id}）: {current_file_path}")
                            else:
                                logger.warning(f"⚠️ [Executor] 输出路径不存在，但会使用: {resolved_path}")
                        else:
                            logger.debug(f"🔍 [Executor] 步骤 {tool_id} 的输出不更新 current_file_path（后续步骤应使用占位符）")
            
            # 构建步骤详情（符合前端格式）；data 合并顶层 plot_path/output_plot_path 以便前端报告渲染
            result_data = step_result.get("result", {}) or {}
            if isinstance(result_data, dict):
                # 🔥 强制前置路径规范化：SSE 实时流与 state_snapshot 入库前统一为 /results/... URL
                result_data = self._normalize_result_paths(result_data)
            if isinstance(result_data, dict):
                for k in ("plot_path", "output_plot_path", "lollipop_path", "clustermap_path"):
                    if step_result.get(k) and k not in result_data:
                        result_data = {**result_data, k: step_result.get(k)}
                # STED-EC 等工具返回 report_data.images: [{title, path}]，展平为 data.images 供前端与 server 路径标准化使用
                rp_images = (result_data.get("report_data") or {}).get("images") or []
                if rp_images and "images" not in result_data:
                    result_data["images"] = [
                        p.get("path") if isinstance(p, dict) else p for p in rp_images
                    ]
            _tb = step_result.get("traceback", "") or step_result.get("debug_info", "")
            _display_summary = ""
            if isinstance(result_data, dict):
                _rs = result_data.get("summary")
                if isinstance(_rs, str) and _rs.strip():
                    _display_summary = _rs
            if not _display_summary:
                _display_summary = step_result.get("message", "") or ""
            _sr_inner: Dict[str, Any] = {
                "step_name": step_name,
                "status": step_result.get("status", "error"),
                "logs": step_result.get("message", ""),
                "data": result_data,
                "error": step_result.get("error", ""),
                "message": step_result.get("message", ""),
                "user_message": step_result.get("user_message", ""),
                "traceback": _tb,
            }
            # 专家报告等：工具返回 result.report_data（含 references），同步到 step_result 顶层供 SSE/前端「参考文献」Tab
            if isinstance(result_data, dict) and isinstance(result_data.get("report_data"), dict):
                _sr_inner["report_data"] = result_data["report_data"]
            step_detail = {
                "step_id": step_id,
                "tool_id": step.get("tool_id"),
                "name": step_name,
                "status": step_result.get("status", "error"),
                "summary": _display_summary,
                "step_result": _sr_inner,
            }
            if step_result.get("status") == "error":
                step_detail["error"] = step_result.get("error", "")
                step_detail["message"] = step_result.get("message", "")
                step_detail["user_message"] = step_result.get("user_message", "")
                step_detail["error_category"] = step_result.get("error_category", "")
                step_detail["suggestion"] = step_result.get("suggestion", "")
                step_detail["can_skip"] = step_result.get("can_skip", False)
                step_detail["technical_details"] = step_result.get("technical_details", "")
                step_detail["traceback"] = step_result.get("traceback", "") or step_result.get("debug_info", "")
                step_detail["debug_info"] = step_result.get("debug_info", "") or step_result.get("traceback", "")
            
            # 🔥 提取图片路径（如果有）- 优先 result，再顶层；支持 plot_path / output_plot_path
            result_data = step_result.get("result", {}) or {}
            if not isinstance(result_data, dict):
                result_data = {}
            cand = (result_data.get("plot_path") or result_data.get("output_plot_path") or
                    step_result.get("plot_path") or step_result.get("output_plot_path"))
            if not cand or not self._is_image_file(cand):
                out_cand = result_data.get("output_path") or result_data.get("output_file") or result_data.get("file_path")
                if out_cand and self._is_image_file(out_cand):
                    cand = out_cand
            if not cand or not self._is_image_file(cand):
                img_cand = result_data.get("image_path") or result_data.get("lollipop_path") or result_data.get("clustermap_path")
                if img_cand and self._is_image_file(img_cand):
                    cand = img_cand
            plot_path = self._to_results_url(cand)
            path_lower = (plot_path or "").lower()
            if plot_path and self._is_image_file(plot_path) and not path_lower.endswith(".nii.gz") and not path_lower.endswith(".nii") and not path_lower.endswith(".dcm"):
                step_detail["plot"] = plot_path
                logger.info(f"🖼️ 检测到图片文件: {plot_path}")
            elif plot_path:
                logger.debug(f"📄 跳过非可展示文件作为 plot: {plot_path}")
            # 🔥 时光机：单步耗时入快照，供历史还原「耗时 Xs」
            step_detail["duration"] = round(time.time() - tool_start, 1)
            step_detail = sanitize_for_json(step_detail)
            steps_details.append(step_detail)
            steps_results.append(step_detail["step_result"])
            
            # 🔥 CRITICAL REGRESSION FIX: If async job started, stop execution
            if step_result.get("status") == "async_job_started":
                logger.info(f"🚀 [Executor] 检测到异步作业已启动: {step_id}, job_id: {step_result.get('job_id', 'N/A')}")
                logger.info(f"🚀 [Executor] 停止执行后续步骤，等待异步作业完成")
                # Add job_id to step_detail if present
                if "job_id" in step_result:
                    step_detail["job_id"] = step_result["job_id"]
                break  # STOP HERE - Do not execute next steps
            
            # 🔥 CRITICAL FIX: 即使步骤失败，也继续执行后续步骤
            # 这样可以确保前面的步骤结果正常显示，不会因为后续步骤失败而影响前面的显示
            if step_result.get("status") == "error":
                error_msg = step_result.get("error") or step_result.get("message") or "未知错误"
                logger.error(f"❌ 步骤 {step_id} 失败: {error_msg}")
                logger.warning(f"⚠️ 继续执行后续步骤，确保前面的步骤结果正常返回")
                # 不 break，继续执行后续步骤
                # 注意：如果后续步骤依赖于当前失败的步骤，它们可能会失败，但至少会尝试执行
        
        # 确定最终状态：success 与 skipped 均视为非失败（skipped 为依赖未满足时的优雅跳过）
        all_ok = all(
            detail.get("status") in ("success", "skipped")
            for detail in steps_details
        )
        workflow_status = "success" if all_ok else "error"
        
        # 提取最终图片（最后一个成功步骤的图片）
        final_plot = None
        for detail in reversed(steps_details):
            if detail.get("plot"):
                final_plot = detail["plot"]
                break
        
        # 构建执行报告（符合前端格式）
        report_data = {
            "status": workflow_status,
            "workflow_name": workflow_name,
            "steps_details": steps_details,
            "steps_results": steps_results,
            "output_dir": self.output_dir
        }
        
        if final_plot:
            report_data["final_plot"] = final_plot
        
        logger.info("=" * 80)
        logger.info(f"✅ 工作流执行完成: {workflow_name} (状态: {workflow_status})")
        logger.info(f"📊 成功步骤: {sum(1 for d in steps_details if d.get('status') == 'success')}/{len(steps_details)}")
        logger.info("=" * 80)
        
        # 🔥 生成 AI Expert Diagnosis（如果提供了 Agent 实例）
        # 注意：由于 execute_workflow 是同步方法，诊断生成将在后台进行或跳过
        # 实际诊断生成应该在 Agent 的 execute_workflow 方法中调用
        # 这里只标记需要诊断，不实际生成（避免异步/同步混用问题）
        if agent and hasattr(agent, '_generate_analysis_summary'):
            logger.info("📝 [Executor] Agent 实例已提供，诊断将在 Agent 层生成")
        
        # 🔥 清理数据以确保 JSON 序列化安全（处理 Numpy 类型、NaN/Infinity 等）
        logger.info("🧹 清理数据以确保 JSON 序列化安全...")
        sanitized_report = sanitize_for_json(report_data)
        logger.info("✅ 数据清理完成")
        
        return sanitized_report

