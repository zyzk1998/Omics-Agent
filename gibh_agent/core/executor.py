"""
通用工作流执行器 - The Hands

动态执行工作流，不依赖硬编码的工具逻辑。
使用 ToolRegistry 查找和执行工具。
"""
import os
import logging
from typing import Dict, Any, List, Optional
from pathlib import Path
from datetime import datetime

from .tool_registry import registry
from .utils import sanitize_for_json

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
        self.step_results: Dict[str, Any] = {}  # 存储步骤结果，用于数据流传递
    
    def _resolve_file_path(self, file_path: str) -> str:
        """
        🔥 CRITICAL REGRESSION FIX: 解析文件路径为绝对路径
        
        检查顺序：
        1. 如果已经是绝对路径且存在 -> 直接使用
        2. 如果是相对路径 -> 前置 UPLOAD_DIR
        3. 如果只是文件名 -> 在 UPLOAD_DIR 中搜索
        
        Args:
            file_path: 原始文件路径（可能是绝对路径、相对路径或文件名）
            
        Returns:
            解析后的绝对路径
        """
        if not file_path or file_path in ["<待上传数据>", "<PENDING_UPLOAD>", ""]:
            return file_path
        
        original_path = file_path
        upload_dir_path = Path(self.upload_dir)
        
        # Check 1: Is it already absolute and existing?
        if os.path.isabs(file_path):
            path_obj = Path(file_path)
            if path_obj.exists():
                logger.info(f"✅ [Executor] 路径已为绝对路径且存在: {file_path}")
                return str(path_obj.resolve())
            else:
                logger.warning(f"⚠️ [Executor] 绝对路径不存在: {file_path}，尝试在 UPLOAD_DIR 中查找")
        
        # Check 2: Is it relative? (e.g., "guest/.../cow_diet.csv")
        # Try prepending UPLOAD_DIR
        if not os.path.isabs(file_path):
            # Remove leading slash if present (e.g., "/guest/..." -> "guest/...")
            file_path_clean = file_path.lstrip('/')
            potential_path = upload_dir_path / file_path_clean
            
            if potential_path.exists():
                resolved = str(potential_path.resolve())
                logger.info(f"✅ [Executor] 解析相对路径: {original_path} -> {resolved}")
                return resolved
            
            # Try with original relative path
            potential_path2 = upload_dir_path / file_path
            if potential_path2.exists():
                resolved = str(potential_path2.resolve())
                logger.info(f"✅ [Executor] 解析相对路径（原始）: {original_path} -> {resolved}")
                return resolved
        
        # Check 3: Is it just a filename? (e.g., "cow_diet.csv")
        # Search in UPLOAD_DIR recursively
        filename = Path(file_path).name
        if filename == file_path or '/' not in file_path.replace('\\', '/'):
            logger.info(f"🔍 [Executor] 搜索文件名: {filename} 在 {self.upload_dir}")
            for found_path in upload_dir_path.rglob(filename):
                if found_path.is_file():
                    resolved = str(found_path.resolve())
                    logger.info(f"✅ [Executor] 找到文件: {original_path} -> {resolved}")
                    return resolved
        
        # If all checks fail, try to construct absolute path anyway
        if not os.path.isabs(file_path):
            final_path = upload_dir_path / file_path.lstrip('/')
            resolved = str(final_path.resolve())
            logger.warning(f"⚠️ [Executor] 无法验证路径存在，但构造绝对路径: {original_path} -> {resolved}")
            return resolved
        
        # Return original if already absolute (even if doesn't exist)
        logger.warning(f"⚠️ [Executor] 无法解析路径，返回原始值: {original_path}")
        return original_path
    
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
        params = step_data.get("params", {})
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
        
        # 🔥 参数映射：根据工具类别映射文件路径参数
        tool_metadata = registry.get_metadata(tool_id)
        tool_category = tool_metadata.category if tool_metadata else None
        
        # 确定工具期望的文件路径参数名
        if tool_category == "scRNA-seq":
            # RNA 工具使用 adata_path
            file_param_name = "adata_path"
            # 如果提供了 file_path 但没有 adata_path，进行映射
            if "file_path" in params and file_param_name not in params:
                params[file_param_name] = params.pop("file_path")
                logger.info(f"🔄 [Executor] 参数映射: file_path -> {file_param_name} (工具: {tool_id})")
        else:
            # 其他工具（如代谢组学）使用 file_path
            file_param_name = "file_path"
        
        # 验证参数（可选但推荐）
        try:
            if tool_metadata:
                # 使用 Pydantic schema 验证参数
                validated_params = tool_metadata.args_schema(**params)
                params = validated_params.model_dump()
                logger.debug(f"✅ 参数验证通过: {step_id}")
        except Exception as validation_err:
            logger.warning(f"⚠️ 参数验证失败（继续执行）: {validation_err}")
            # 继续执行，不因验证失败而中断
        
        # 处理数据流：替换占位符（传递工具类别信息）
        processed_params = self._process_data_flow(params, step_context, tool_category=tool_category)
        
        # 🔥 CRITICAL FIX: 对于 scRNA-seq 工具，确保移除 file_path 参数（如果存在）
        if tool_category == "scRNA-seq" and "file_path" in processed_params:
            # 如果已经有 adata_path，移除 file_path
            if "adata_path" in processed_params:
                del processed_params["file_path"]
                logger.info(f"🔄 [Executor] 移除多余的 file_path 参数（工具已有 adata_path）")
            else:
                # 如果没有 adata_path，将 file_path 映射为 adata_path
                processed_params["adata_path"] = processed_params.pop("file_path")
                logger.info(f"🔄 [Executor] 参数映射: file_path -> adata_path (工具: {tool_id})")
        
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
        path_params = ["file_path", "adata_path", "output_path", "output_file", "fastq_path", "reference_path", "output_h5ad"]
        for param_name in path_params:
            if param_name in processed_params:
                original_path = processed_params[param_name]
                if original_path and isinstance(original_path, str):
                    # Skip placeholder values
                    if original_path not in ["<待上传数据>", "<PENDING_UPLOAD>", ""]:
                        resolved_path = self._resolve_file_path(original_path)
                        if resolved_path != original_path:
                            logger.info(f"🔄 [Executor] 解析路径参数 {param_name}: {original_path} -> {resolved_path}")
                        processed_params[param_name] = resolved_path
        
        # 🔥 ARCHITECTURAL UPGRADE: Phase 3 - Pre-Flight Check & Auto-Correction
        # 对于需要 group_column 的工具，验证列是否存在，如果不存在则使用 semantic_map 自动修正
        tools_requiring_group_column = ["differential_analysis", "metabolomics_plsda", "metabolomics_pathway_enrichment"]
        if tool_id in tools_requiring_group_column and "group_column" in processed_params:
            group_column = processed_params.get("group_column")
            file_path = processed_params.get("file_path")
            
            if file_path and os.path.exists(file_path):
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
        
        # 执行工具
        try:
            logger.info(f"🚀 调用工具: {tool_id} with params: {list(processed_params.keys())}")
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
                # 存储结果供后续步骤使用（即使失败也存储，用于调试）
                self.step_results[step_id] = result
                return {
                    "status": "error",
                    "step_id": step_id,
                    "step_name": step_name,
                    "tool_id": tool_id,
                    "result": result,
                    "error": error_msg,
                    "message": error_msg  # 🔥 使用错误消息，而不是"执行完成"
                }
            else:
                logger.info(f"✅ 步骤 {step_id} 执行成功")
                # 存储结果供后续步骤使用
                self.step_results[step_id] = result
                return {
                    "status": result.get("status", "success"),
                    "step_id": step_id,
                    "step_name": step_name,
                    "tool_id": tool_id,
                    "result": result,
                    "message": result.get("message", f"步骤 {step_name} 执行完成")
                }
        
        except Exception as e:
            error_msg = f"步骤 {step_id} 执行失败: {str(e)}"
            logger.error(f"❌ {error_msg}", exc_info=True)
            
            return {
                "status": "error",
                "step_id": step_id,
                "step_name": step_name,
                "tool_id": tool_id,
                "error": str(e),
                "message": error_msg
            }
    
    def _process_data_flow(
        self,
        params: Dict[str, Any],
        step_context: Optional[Dict[str, Any]] = None,
        tool_category: Optional[str] = None
    ) -> Dict[str, Any]:
        """
        处理数据流：替换占位符（如 <step1_output>）和自动注入文件路径
        
        Args:
            params: 原始参数
            step_context: 步骤上下文（包含 current_file_path）
            tool_category: 工具类别（用于确定文件路径参数名）
        
        Returns:
            处理后的参数
        """
        processed = {}
        
        # 🔥 根据工具类别确定文件路径参数名
        if tool_category == "scRNA-seq":
            file_param_name = "adata_path"
        else:
            file_param_name = "file_path"
        
        # 🔥 如果文件路径参数缺失，尝试从上下文注入
        if file_param_name not in params and step_context:
            current_file_path = step_context.get("current_file_path")
            if current_file_path:
                processed[file_param_name] = current_file_path
                logger.info(f"🔄 数据流: 自动注入 {file_param_name} = {current_file_path}")
        
        # 🔥 如果提供了 file_path 但工具需要 adata_path，进行映射
        if tool_category == "scRNA-seq" and "file_path" in params and "adata_path" not in params:
            processed["adata_path"] = params.pop("file_path")
            logger.info(f"🔄 数据流: 参数映射 file_path -> adata_path")
        
        for key, value in params.items():
            if isinstance(value, str) and value.startswith("<") and value.endswith(">"):
                # 占位符，尝试从上下文或步骤结果中获取
                placeholder = value[1:-1]  # 移除 < >
                
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
                
                # 尝试从 step_results 中获取
                if placeholder in self.step_results:
                    step_result = self.step_results[placeholder]
                    # 🔥 CRITICAL FIX: 对于 scRNA-seq 工具，优先提取 output_h5ad
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
                            # 其他工具使用标准字段
                            output_path = (
                                step_result.get("output_file") or
                                step_result.get("output_path") or
                                step_result.get("file_path") or
                                step_result.get("plot_path") or
                                step_result.get("result_path") or
                                step_result.get("preprocessed_file")
                            )
                        if output_path:
                            processed[key] = output_path
                            logger.info(f"🔄 数据流: {key} = <{placeholder}> -> {output_path}")
                        else:
                            # 如果没有找到路径，使用整个结果
                            processed[key] = step_result
                    else:
                        processed[key] = step_result
                elif step_context and placeholder in step_context:
                    processed[key] = step_context[placeholder]
                else:
                    # 占位符未解析，保持原样（可能后续步骤会处理）
                    logger.warning(f"⚠️ 无法解析占位符: {value}")
                    processed[key] = value
            else:
                processed[key] = value
        
        return processed
    
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
        agent: Optional[Any] = None  # 可选的 Agent 实例，用于生成诊断
    ) -> Dict[str, Any]:
        """
        执行整个工作流
        
        Args:
            workflow_data: 工作流配置（包含 workflow_name 和 steps）
            file_paths: 输入文件路径列表
            output_dir: 输出目录（如果为 None，将自动创建）
        
        Returns:
            执行报告（符合前端 analysis_report 格式）
        """
        workflow_name = workflow_data.get("workflow_name", "Unknown Workflow")
        steps = workflow_data.get("steps", [])
        
        logger.info("=" * 80)
        logger.info(f"🚀 开始执行工作流: {workflow_name}")
        logger.info(f"📋 步骤数: {len(steps)}")
        logger.info("=" * 80)
        
        # 设置输出目录
        if output_dir is None:
            timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
            output_dir = f"./results/run_{timestamp}"
        
        output_path = Path(output_dir)
        output_path.mkdir(parents=True, exist_ok=True)
        self.output_dir = str(output_path)
        
        logger.info(f"📂 输出目录: {self.output_dir}")
        
        # 初始化步骤结果列表
        steps_details = []
        steps_results = []
        
        # 🔥 上下文链：跟踪当前文件路径，用于自动传递给下一个步骤
        # 🔥 CRITICAL REGRESSION FIX: Resolve file paths to absolute paths
        resolved_file_paths = []
        if file_paths:
            for fp in file_paths:
                if fp and isinstance(fp, str):
                    resolved = self._resolve_file_path(fp)
                    resolved_file_paths.append(resolved)
                    if resolved != fp:
                        logger.info(f"🔄 [Executor] 解析输入文件路径: {fp} -> {resolved}")
        current_file_path = resolved_file_paths[0] if resolved_file_paths else None
        
        # 执行每个步骤
        for i, step in enumerate(steps, 1):
            step_id = step.get("step_id", f"step{i}")
            step_name = step.get("name", step.get("step_name", step_id))
            tool_id = step.get("tool_id", step_id)
            params = step.get("params", {})
            
            # 🔥 CRITICAL FIX: 完全移除 visualize_pca 步骤（不在流程中显示）
            if tool_id == "visualize_pca" or step_id == "visualize_pca":
                logger.warning(f"⚠️ [Executor] 完全移除 visualize_pca 步骤（pca_analysis 已包含可视化）")
                continue  # 直接跳过，不添加到步骤详情中
            
            logger.info(f"\n{'=' * 80}")
            logger.info(f"📌 步骤 {i}/{len(steps)}: {step_name} ({step_id})")
            logger.info(f"{'=' * 80}")
            
            # 🔥 智能参数映射：根据工具类型自动映射文件路径参数
            tool_metadata = registry.get_metadata(tool_id)
            tool_category = tool_metadata.category if tool_metadata else None
            
            # 确定工具期望的文件路径参数名
            if tool_category == "scRNA-seq":
                # RNA 工具使用 adata_path
                file_param_name = "adata_path"
            else:
                # 其他工具（如代谢组学）使用 file_path
                file_param_name = "file_path"
            
            # 自动注入文件路径（如果缺失且我们有当前文件路径）
            if file_param_name not in params and current_file_path:
                params[file_param_name] = current_file_path
                logger.info(f"🔄 自动注入 {file_param_name}: {current_file_path}")
            
            # 🔥 参数映射：如果工具期望 adata_path 但提供了 file_path，进行映射
            if file_param_name == "adata_path" and "file_path" in params and file_param_name not in params:
                params[file_param_name] = params.pop("file_path")
                logger.info(f"🔄 参数映射: file_path -> {file_param_name}")
            
            # 如果工具需要 output_dir，也自动注入
            if "output_dir" not in params and self.output_dir:
                params["output_dir"] = self.output_dir
            
            # 更新步骤的 params
            step["params"] = params
            
            # 构建步骤上下文（包含文件路径等）
            step_context = {
                "file_paths": file_paths or [],
                "output_dir": self.output_dir,
                "workflow_name": workflow_name,
                "current_file_path": current_file_path  # 传递当前文件路径
            }
            
            # 执行步骤
            step_result = self.execute_step(step, step_context)
            
            # 🔥 CRITICAL FIX: 更新 current_file_path 供下一个步骤使用
            # 对于 scRNA-seq 工具，优先使用 output_h5ad
            result_data = step_result.get("result", {})
            if isinstance(result_data, dict):
                # 🔥 对于 scRNA-seq 工具，优先查找 output_h5ad
                tool_metadata = registry.get_metadata(step.get("tool_id", ""))
                tool_category = tool_metadata.category if tool_metadata else None
                
                if tool_category == "scRNA-seq":
                    # scRNA-seq 工具优先使用 output_h5ad
                    next_file_path = (
                        result_data.get("output_h5ad") or  # 🔥 优先使用 output_h5ad
                        result_data.get("output_file") or
                        result_data.get("output_path") or
                        result_data.get("file_path")
                    )
                else:
                    # 其他工具使用标准字段
                    next_file_path = (
                        result_data.get("output_file") or
                        result_data.get("output_path") or
                        result_data.get("file_path") or
                        result_data.get("preprocessed_file")
                    )
                
                if next_file_path:
                    # 🔥 修复：即使文件不存在也更新路径（文件可能稍后创建）
                    current_file_path = next_file_path
                    if os.path.exists(next_file_path):
                        logger.info(f"✅ 更新当前文件路径: {current_file_path}")
                    else:
                        logger.warning(f"⚠️ 输出路径不存在，但会使用: {next_file_path} (文件可能稍后创建)")
            
            # 构建步骤详情（符合前端格式）
            step_detail = {
                "step_id": step_id,
                "tool_id": step.get("tool_id"),
                "name": step_name,
                "status": step_result.get("status", "error"),
                "summary": step_result.get("message", ""),
                "step_result": {
                    "step_name": step_name,
                    "status": step_result.get("status", "error"),
                    "logs": step_result.get("message", ""),
                    "data": step_result.get("result", {})
                }
            }
            
            # 🔥 提取图片路径（如果有）- 严格检查文件类型
            result_data = step_result.get("result", {})
            if isinstance(result_data, dict):
                # 优先检查明确的图片路径字段
                plot_path = result_data.get("plot_path") or result_data.get("image_path")
                
                # 如果没有明确的图片路径字段，检查 output_path 的文件扩展名
                if not plot_path:
                    output_path = result_data.get("output_path") or result_data.get("output_file") or result_data.get("file_path")
                    if output_path and self._is_image_file(output_path):
                        plot_path = output_path
                
                # 只有确认是图片文件才添加到 plot 字段
                if plot_path and self._is_image_file(plot_path):
                    step_detail["plot"] = plot_path
                    logger.info(f"🖼️ 检测到图片文件: {plot_path}")
                elif plot_path:
                    # 如果不是图片文件（如 CSV），记录到 data 字段而不是 plot
                    logger.debug(f"📄 检测到非图片文件: {plot_path}，不添加到 plot 字段")
            
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
            
            # 如果步骤失败，停止执行
            if step_result.get("status") == "error":
                logger.error(f"❌ 步骤 {step_id} 失败，停止工作流执行")
                break
        
        # 确定最终状态
        all_success = all(
            detail.get("status") == "success"
            for detail in steps_details
        )
        workflow_status = "success" if all_success else "error"
        
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

