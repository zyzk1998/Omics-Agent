"""代谢组学智能体（Metabolomics Agent）"""
from typing import Dict, Any, List, AsyncIterator
from ..base_agent import BaseAgent
from ...core.llm_client import LLMClient
from ...core.prompt_manager import PromptManager
from ...tools.metabolomics_tool import MetabolomicsTool
import logging

logger = logging.getLogger(__name__)


class MetabolomicsAgent(BaseAgent):
    """代谢组学智能体"""
    
    def __init__(
        self,
        llm_client: LLMClient,
        prompt_manager: PromptManager,
        metabolomics_config: Dict[str, Any] = None
    ):
        super().__init__(llm_client, prompt_manager, "metabolomics_expert")
        self.metabolomics_config = metabolomics_config or {}
        self.metabolomics_tool = MetabolomicsTool(self.metabolomics_config)
    
    async def process_query(
        self,
        query: str,
        history: List[Dict[str, str]] = None,
        uploaded_files: List[Dict[str, str]] = None,
        **kwargs
    ) -> Dict[str, Any]:
        """
        处理用户查询
        
        Returns:
            处理结果字典，可能包含：
            - chat: 聊天响应（流式）
        """
        query_lower = query.lower().strip()
        file_paths = self.get_file_paths(uploaded_files or [])
        
        # 判断是否是工作流请求
        is_workflow_request = self._is_workflow_request(query_lower, file_paths)
        
        if is_workflow_request:
            # 工作流请求：先检查数据，然后生成工作流配置
            return await self._generate_workflow_config(query, file_paths)
        else:
            # 普通聊天：流式响应
            return {
                "type": "chat",
                "response": self._stream_chat_response(query, file_paths)
            }
    
    def _is_workflow_request(self, query: str, file_paths: List[str]) -> bool:
        """判断是否是工作流请求"""
        workflow_keywords = [
            "规划", "流程", "workflow", "pipeline", "分析", "run",
            "执行", "plan", "做一下", "跑一下", "分析一下", "全流程"
        ]
        
        if any(kw in query for kw in workflow_keywords):
            return True
        
        if file_paths and (not query or len(query) < 5):
            return True
        
        return False
    
    async def _generate_workflow_config(
        self,
        query: str,
        file_paths: List[str]
    ) -> Dict[str, Any]:
        """
        生成工作流配置
        
        流程：
        1. 先检查数据（inspect_data）
        2. 基于检查结果提取参数
        3. 生成工作流配置
        """
        # 强制检查：如果有文件，先检查
        inspection_result = None
        if file_paths:
            input_path = file_paths[0]
            try:
                inspection_result = self.metabolomics_tool.inspect_data(input_path)
                if "error" in inspection_result:
                    logger.warning(f"File inspection failed: {inspection_result.get('error')}")
            except Exception as e:
                logger.error(f"Error inspecting file: {e}", exc_info=True)
        
        # 使用 LLM 提取参数（传入检查结果）
        extracted_params = await self._extract_workflow_params(query, file_paths, inspection_result)
        
        # 构建工作流配置
        workflow_config = {
            "workflow_name": "Metabolomics Analysis Pipeline",
            "steps": [
                {
                    "step_id": "inspect_data",
                    "tool_id": "inspect_data",
                    "desc": "检查数据文件",
                    "params": {"file_path": file_paths[0] if file_paths else ""}
                },
                {
                    "step_id": "preprocess_data",
                    "tool_id": "preprocess_data",
                    "desc": "数据预处理",
                    "params": {
                        "file_path": file_paths[0] if file_paths else "",
                        "missing_threshold": extracted_params.get("missing_threshold", "0.5"),
                        "normalization": extracted_params.get("normalization", "log2"),
                        "scale": extracted_params.get("scale", "true")
                    }
                },
                {
                    "step_id": "pca_analysis",
                    "tool_id": "pca_analysis",
                    "desc": "主成分分析",
                    "params": {
                        "n_components": extracted_params.get("n_components", "10")
                    }
                },
                {
                    "step_id": "differential_analysis",
                    "tool_id": "differential_analysis",
                    "desc": "差异代谢物分析",
                    "params": {
                        "group_column": extracted_params.get("group_column", "Muscle loss"),
                        "method": extracted_params.get("method", "t-test"),
                        "p_value_threshold": extracted_params.get("p_value_threshold", "0.05"),
                        "fold_change_threshold": extracted_params.get("fold_change_threshold", "1.5"),
                        "group1": extracted_params.get("group1"),  # 如果 >2 个分组，需要指定
                        "group2": extracted_params.get("group2")   # 如果 >2 个分组，需要指定
                    }
                },
                {
                    "step_id": "visualize_pca",
                    "tool_id": "visualize_pca",
                    "desc": "PCA 可视化",
                    "params": {
                        "group_column": extracted_params.get("group_column", "Muscle loss"),
                        "pc1": "1",
                        "pc2": "2"
                    }
                },
                {
                    "step_id": "visualize_volcano",
                    "tool_id": "visualize_volcano",
                    "desc": "火山图可视化",
                    "params": {
                        "p_value_threshold": extracted_params.get("p_value_threshold", "0.05"),
                        "fold_change_threshold": extracted_params.get("fold_change_threshold", "1.5")
                    }
                }
            ]
        }
        
        return {
            "type": "workflow_config",
            "workflow_data": workflow_config,
            "file_paths": file_paths
        }
    
    async def _extract_workflow_params(
        self,
        query: str,
        file_paths: List[str],
        inspection_result: Dict[str, Any] = None
    ) -> Dict[str, Any]:
        """
        使用 LLM 提取工作流参数
        
        基于检查结果智能推荐参数
        """
        # 构建包含检查结果的提示
        inspection_info = ""
        if inspection_result and "error" not in inspection_result:
            inspection_info = f"""
【Data Inspection Results】
- Number of samples: {inspection_result.get('n_samples', 'N/A')}
- Number of metabolites: {inspection_result.get('n_metabolites', 'N/A')}
- Missing values: {inspection_result.get('missing_values', {}).get('percentage', 'N/A')}%
- Group column: {inspection_result.get('group_info', {}).get('column', 'N/A')}
- Groups: {inspection_result.get('group_info', {}).get('groups', {})}
"""
        
        prompt = f"""
Extract workflow parameters for metabolomics analysis from the user query.

User Query: {query}
File Paths: {file_paths}

{inspection_info}

Based on the inspection results and user query, extract the following parameters:
- missing_threshold: Threshold for removing metabolites with high missing values (default: 0.5)
- normalization: Normalization method - "log2", "zscore", or "none" (default: "log2")
- scale: Whether to apply StandardScaler (default: true)
- n_components: Number of PCA components (default: 10)
- group_column: Column name for group comparison (default: "Muscle loss" or first metadata column)
- method: Statistical method for differential analysis - "t-test" or "mann-whitney" (default: "t-test")
- p_value_threshold: P-value threshold for significance (default: 0.05)
- fold_change_threshold: Fold change threshold (default: 1.5)

Return JSON only:
{{"missing_threshold": "0.5", "normalization": "log2", "scale": "true", "n_components": "10", "group_column": "Muscle loss", "method": "t-test", "p_value_threshold": "0.05", "fold_change_threshold": "1.5"}}
"""
        
        messages = [
            {"role": "system", "content": "You are a parameter extraction assistant. Return JSON only."},
            {"role": "user", "content": prompt}
        ]
        
        try:
            completion = await self.llm_client.achat(messages, temperature=0.1, max_tokens=256)
            think_content, response = self.llm_client.extract_think_and_content(completion)
            
            # 解析 JSON
            import json
            json_str = response.strip()
            if "```json" in json_str:
                json_str = json_str.split("```json")[1].split("```")[0].strip()
            
            return json.loads(json_str)
        except Exception as e:
            logger.error(f"Error extracting parameters: {e}")
            return {}
    
    async def _stream_chat_response(
        self,
        query: str,
        file_paths: List[str]
    ) -> AsyncIterator[str]:
        """
        流式聊天响应（支持 ReAct 循环和工具调用）
        """
        context = {
            "context": f"Uploaded files: {', '.join(file_paths) if file_paths else 'None'}",
            "available_tools": list(self.metabolomics_tool.tool_map.keys()),
            "tool_descriptions": {
                "inspect_data": "检查代谢组学数据文件，返回数据摘要（样本数、代谢物数、缺失值、分组信息等）",
                "preprocess_data": "预处理数据：处理缺失值、标准化、缩放",
                "pca_analysis": "执行主成分分析 (PCA)",
                "differential_analysis": "执行差异代谢物分析（两组比较）",
                "visualize_pca": "生成 PCA 可视化图",
                "visualize_volcano": "生成火山图（Volcano Plot）"
            }
        }
        
        # 如果有文件，强制先检查（符合 SOP）
        inspection_result = None
        if file_paths:
            input_path = file_paths[0]
            try:
                inspection_result = self.metabolomics_tool.inspect_data(input_path)
                if "error" not in inspection_result:
                    # 将检查结果添加到上下文中
                    inspection_summary = f"""
【Data Inspection Completed】
- Samples: {inspection_result.get('n_samples', 'N/A')}
- Metabolites: {inspection_result.get('n_metabolites', 'N/A')}
- Missing values: {inspection_result.get('missing_values', {}).get('percentage', 'N/A')}%
- Group column: {inspection_result.get('group_info', {}).get('column', 'N/A')}
- Groups: {inspection_result.get('group_info', {}).get('groups', {})}
"""
                    # 先输出检查结果
                    yield f"🔍 **Data Inspection Results:**\n{inspection_summary}\n\n"
                    # 将检查结果添加到查询中，让 LLM 基于此分析
                    query = f"""{query}

{inspection_summary}

Based on the inspection results above, please:
1. Analyze the data characteristics
2. Propose appropriate analysis parameters
3. Ask for confirmation before proceeding with analysis
"""
            except Exception as e:
                logger.error(f"Error inspecting file: {e}", exc_info=True)
                yield f"⚠️ Warning: Could not inspect file: {str(e)}\n\n"
        
        # 构建增强的用户查询，包含工具说明
        enhanced_query = f"""{query}

【Available Tools】
You have access to:
- inspect_data(file_path): Check data file structure (already executed above if files were provided)
- preprocess_data(file_path, missing_threshold, normalization, scale): Preprocess metabolomics data
- pca_analysis(n_components, file_path): Perform PCA analysis
- differential_analysis(group_column, file_path, method, p_value_threshold, fold_change_threshold): Find differential metabolites
- visualize_pca(group_column, pca_file, pc1, pc2): Generate PCA plot
- visualize_volcano(diff_file, p_value_threshold, fold_change_threshold): Generate volcano plot

【Workflow Rule】
- If user provides CSV files: Inspect first (already done above), then analyze and propose parameters
- Before running any analysis, you MUST have inspected the data first
- Now analyze the inspection results and propose parameters.
"""
        
        # 流式输出 LLM 响应
        async for chunk in self.chat(enhanced_query, context, stream=True):
            yield chunk
    
    async def execute_workflow(
        self,
        workflow_config: Dict[str, Any],
        file_paths: List[str],
        output_dir: str
    ) -> Dict[str, Any]:
        """
        执行代谢组学工作流
        
        Args:
            workflow_config: 工作流配置
            file_paths: 文件路径列表
            output_dir: 输出目录
        
        Returns:
            分析报告
        """
        import os
        
        input_path = file_paths[0] if file_paths else None
        if not input_path:
            raise ValueError("No input files provided")
        
        # 设置输出目录
        if not os.path.exists(output_dir):
            os.makedirs(output_dir, exist_ok=True)
        
        # 更新代谢组工具的输出目录
        from pathlib import Path
        self.metabolomics_config["output_dir"] = output_dir
        self.metabolomics_tool.output_dir = Path(output_dir)
        self.metabolomics_tool.output_dir.mkdir(parents=True, exist_ok=True)
        
        # 执行工作流步骤
        steps = workflow_config.get("steps", [])
        steps_details = []
        final_plot = None
        
        try:
            for step in steps:
                step_id = step.get("step_id")
                tool_id = step.get("tool_id")
                params = step.get("params", {})
                
                logger.info(f"🔧 执行步骤: {step_id} ({tool_id})")
                
                if tool_id == "inspect_data":
                    result = self.metabolomics_tool.inspect_data(params.get("file_path", input_path))
                    steps_details.append({
                        "step_id": step_id,
                        "tool_id": tool_id,
                        "name": step.get("desc", step_id),
                        "summary": f"检查完成: {result.get('n_samples', 'N/A')} 个样本, {result.get('n_metabolites', 'N/A')} 个代谢物",
                        "status": result.get("status", "success")
                    })
                
                elif tool_id == "preprocess_data":
                    result = self.metabolomics_tool.preprocess_data(
                        file_path=params.get("file_path", input_path),
                        missing_threshold=float(params.get("missing_threshold", "0.5")),
                        normalization=params.get("normalization", "log2"),
                        scale=params.get("scale", "true").lower() == "true"
                    )
                    steps_details.append({
                        "step_id": step_id,
                        "tool_id": tool_id,
                        "name": step.get("desc", step_id),
                        "summary": result.get("message", "预处理完成"),
                        "status": result.get("status", "success")
                    })
                
                elif tool_id == "pca_analysis":
                    # 如果上一步是预处理，使用预处理后的文件
                    preprocessed_file = None
                    for prev_step in steps_details:
                        if prev_step.get("tool_id") == "preprocess_data":
                            # 从预处理结果中获取文件路径
                            preprocessed_file = os.path.join(output_dir, "preprocessed_data.csv")
                            break
                    
                    result = self.metabolomics_tool.pca_analysis(
                        n_components=int(params.get("n_components", "10")),
                        file_path=preprocessed_file or params.get("file_path", input_path)
                    )
                    steps_details.append({
                        "step_id": step_id,
                        "tool_id": tool_id,
                        "name": step.get("desc", step_id),
                        "summary": result.get("message", "PCA 分析完成"),
                        "status": result.get("status", "success")
                    })
                
                elif tool_id == "differential_analysis":
                    # 使用预处理后的文件
                    preprocessed_file = os.path.join(output_dir, "preprocessed_data.csv")
                    if not os.path.exists(preprocessed_file):
                        preprocessed_file = input_path
                    
                    result = self.metabolomics_tool.differential_analysis(
                        group_column=params.get("group_column", "Muscle loss"),
                        file_path=preprocessed_file,
                        method=params.get("method", "t-test"),
                        p_value_threshold=float(params.get("p_value_threshold", "0.05")),
                        fold_change_threshold=float(params.get("fold_change_threshold", "1.5")),
                        group1=params.get("group1"),
                        group2=params.get("group2")
                    )
                    
                    # 检查是否需要用户选择分组
                    if result.get("status") == "need_selection":
                        return {
                            "status": "need_selection",
                            "message": result.get("message", "需要选择比较的分组"),
                            "groups": result.get("groups", []),
                            "steps_details": steps_details
                        }
                    
                    steps_details.append({
                        "step_id": step_id,
                        "tool_id": tool_id,
                        "name": step.get("desc", step_id),
                        "summary": result.get("message", "差异分析完成"),
                        "status": result.get("status", "success"),
                        "details": result.get("summary", "")
                    })
                
                elif tool_id == "visualize_pca":
                    # 使用 PCA 分析结果文件
                    pca_file = params.get("pca_file") or os.path.join(output_dir, "pca_results.csv")
                    if not os.path.exists(pca_file):
                        # 如果 PCA 文件不存在，尝试从步骤详情中获取
                        for prev_step in steps_details:
                            if prev_step.get("tool_id") == "pca_analysis":
                                pca_file = os.path.join(output_dir, "pca_results.csv")
                                break
                    
                    result = self.metabolomics_tool.visualize_pca(
                        group_column=params.get("group_column", "Muscle loss"),
                        pca_file=pca_file,
                        pc1=int(params.get("pc1", "1")),
                        pc2=int(params.get("pc2", "2"))
                    )
                    plot_path = result.get("plot_path") or result.get("plot_file")
                    if plot_path:
                        # 转换为相对路径（相对于 output_dir）
                        if os.path.isabs(plot_path):
                            plot_path = os.path.relpath(plot_path, output_dir)
                        # 确保路径使用正斜杠（Web 兼容）
                        plot_path = plot_path.replace("\\", "/")
                        final_plot = plot_path
                    steps_details.append({
                        "step_id": step_id,
                        "tool_id": tool_id,
                        "name": step.get("desc", step_id),
                        "summary": "PCA 可视化完成",
                        "status": result.get("status", "success"),
                        "plot": plot_path.replace("\\", "/") if plot_path else None
                    })
                
                elif tool_id == "visualize_volcano":
                    # 使用差异分析结果文件
                    diff_file = params.get("diff_file") or os.path.join(output_dir, "differential_results.csv")
                    if not os.path.exists(diff_file):
                        # 如果差异分析文件不存在，尝试从步骤详情中获取
                        for prev_step in steps_details:
                            if prev_step.get("tool_id") == "differential_analysis":
                                diff_file = os.path.join(output_dir, "differential_results.csv")
                                break
                    
                    result = self.metabolomics_tool.visualize_volcano(
                        diff_file=diff_file,
                        p_value_threshold=float(params.get("p_value_threshold", "0.05")),
                        fold_change_threshold=float(params.get("fold_change_threshold", "1.5"))
                    )
                    plot_path = result.get("plot_path") or result.get("plot_file")
                    if plot_path:
                        # 转换为相对路径（相对于 output_dir）
                        if os.path.isabs(plot_path):
                            plot_path = os.path.relpath(plot_path, output_dir)
                        # 确保路径使用正斜杠（Web 兼容）
                        plot_path = plot_path.replace("\\", "/")
                        if not final_plot:
                            final_plot = plot_path
                    steps_details.append({
                        "step_id": step_id,
                        "tool_id": tool_id,
                        "name": step.get("desc", step_id),
                        "summary": "火山图可视化完成",
                        "status": result.get("status", "success"),
                        "plot": plot_path.replace("\\", "/") if plot_path else None
                    })
            
            return {
                "status": "success",
                "workflow_name": workflow_config.get("workflow_name", "Metabolomics Analysis"),
                "steps_details": steps_details,
                "final_plot": final_plot,
                "output_dir": output_dir
            }
            
        except Exception as e:
            logger.error(f"❌ 工作流执行失败: {e}", exc_info=True)
            return {
                "status": "error",
                "error": str(e),
                "steps_details": steps_details,
                "output_dir": output_dir
            }

