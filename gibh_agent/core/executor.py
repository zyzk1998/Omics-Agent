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
    
    def __init__(self, output_dir: Optional[str] = None):
        """
        初始化工作流执行器
        
        Args:
            output_dir: 输出目录（如果为 None，将在执行时创建）
        """
        self.output_dir = output_dir
        self.step_results: Dict[str, Any] = {}  # 存储步骤结果，用于数据流传递
    
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
        
        # 验证参数（可选但推荐）
        try:
            tool_metadata = registry.get_metadata(tool_id)
            if tool_metadata:
                # 使用 Pydantic schema 验证参数
                validated_params = tool_metadata.args_schema(**params)
                params = validated_params.model_dump()
                logger.debug(f"✅ 参数验证通过: {step_id}")
        except Exception as validation_err:
            logger.warning(f"⚠️ 参数验证失败（继续执行）: {validation_err}")
            # 继续执行，不因验证失败而中断
        
        # 处理数据流：替换占位符
        processed_params = self._process_data_flow(params, step_context)
        
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
        step_context: Optional[Dict[str, Any]] = None
    ) -> Dict[str, Any]:
        """
        处理数据流：替换占位符（如 <step1_output>）
        
        Args:
            params: 原始参数
            step_context: 步骤上下文
        
        Returns:
            处理后的参数
        """
        processed = {}
        
        for key, value in params.items():
            if isinstance(value, str) and value.startswith("<") and value.endswith(">"):
                # 占位符，尝试从上下文或步骤结果中获取
                placeholder = value[1:-1]  # 移除 < >
                
                # 尝试从 step_results 中获取
                if placeholder in self.step_results:
                    step_result = self.step_results[placeholder]
                    # 提取输出路径（如果存在）
                    if isinstance(step_result, dict):
                        # 尝试多种可能的输出路径字段
                        output_path = (
                            step_result.get("output_path") or
                            step_result.get("file_path") or
                            step_result.get("plot_path") or
                            step_result.get("result_path")
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
    
    def execute_workflow(
        self,
        workflow_data: Dict[str, Any],
        file_paths: List[str] = None,
        output_dir: Optional[str] = None
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
        
        # 执行每个步骤
        for i, step in enumerate(steps, 1):
            step_id = step.get("step_id", f"step{i}")
            step_name = step.get("name", step.get("step_name", step_id))
            
            logger.info(f"\n{'=' * 80}")
            logger.info(f"📌 步骤 {i}/{len(steps)}: {step_name} ({step_id})")
            logger.info(f"{'=' * 80}")
            
            # 构建步骤上下文（包含文件路径等）
            step_context = {
                "file_paths": file_paths or [],
                "output_dir": self.output_dir,
                "workflow_name": workflow_name
            }
            
            # 执行步骤
            step_result = self.execute_step(step, step_context)
            
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
            
            # 提取图片路径（如果有）
            result_data = step_result.get("result", {})
            if isinstance(result_data, dict):
                plot_path = (
                    result_data.get("plot_path") or
                    result_data.get("image_path") or
                    result_data.get("output_path")
                )
                if plot_path:
                    step_detail["plot"] = plot_path
            
            steps_details.append(step_detail)
            steps_results.append(step_detail["step_result"])
            
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
        
        return report_data

