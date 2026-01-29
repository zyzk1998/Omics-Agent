"""
错误信息格式化器

将面向程序员的错误信息转换为用户友好的提示。
"""
import logging
from typing import Dict, Any, Optional

logger = logging.getLogger(__name__)


class ErrorFormatter:
    """错误信息格式化器"""
    
    @staticmethod
    def format_error(
        error: str,
        tool_id: str,
        step_name: str,
        error_type: Optional[str] = None
    ) -> Dict[str, Any]:
        """
        格式化错误信息
        
        Args:
            error: 原始错误信息
            tool_id: 工具ID
            step_name: 步骤名称
            error_type: 错误类型（可选）
        
        Returns:
            格式化后的错误信息字典，包含：
            - user_message: 用户友好的错误消息
            - error_category: 错误类别（data_issue, config_issue, skipnable, etc.）
            - suggestion: 优化建议
            - can_skip: 是否可以跳过此步骤
        """
        error_lower = error.lower()
        
        # 1. 数据问题类错误
        if any(keyword in error_lower for keyword in [
            "file not found", "文件不存在", "path does not exist", "no such file",
            "empty", "空", "no data", "没有数据", "missing", "缺失"
        ]):
            return {
                "user_message": f"数据文件问题：{step_name} 步骤无法找到或读取数据文件。",
                "error_category": "data_issue",
                "suggestion": "请检查数据文件路径是否正确，文件是否存在且可读。",
                "can_skip": False
            }
        
        # 2. 配置问题类错误
        if any(keyword in error_lower for keyword in [
            "not installed", "未安装", "import error", "导入错误",
            "module not found", "模块未找到", "missing dependency", "缺少依赖"
        ]):
            return {
                "user_message": f"依赖缺失：{step_name} 步骤需要额外的软件包或工具。",
                "error_category": "config_issue",
                "suggestion": "请联系管理员安装必要的依赖包，或跳过此步骤。",
                "can_skip": True
            }
        
        # 3. 参数问题类错误
        if any(keyword in error_lower for keyword in [
            "invalid parameter", "参数无效", "parameter error", "参数错误",
            "missing parameter", "缺少参数", "required parameter", "必需参数"
        ]):
            return {
                "user_message": f"参数配置问题：{step_name} 步骤的参数配置不正确。",
                "error_category": "config_issue",
                "suggestion": "请检查步骤参数设置，或使用默认参数重试。",
                "can_skip": False
            }
        
        # 4. 内存/资源问题
        if any(keyword in error_lower for keyword in [
            "memory", "内存", "out of memory", "内存不足", "resource", "资源"
        ]):
            return {
                "user_message": f"资源不足：{step_name} 步骤需要更多内存或计算资源。",
                "error_category": "resource_issue",
                "suggestion": "请尝试减小数据规模，或联系管理员增加系统资源。",
                "can_skip": True
            }
        
        # 5. 网络/下载问题（如 CellTypist 模型下载）
        if any(keyword in error_lower for keyword in [
            "download", "下载", "network", "网络", "connection", "连接",
            "timeout", "超时", "无法下载"
        ]):
            return {
                "user_message": f"网络问题：{step_name} 步骤需要下载资源，但网络连接失败。",
                "error_category": "network_issue",
                "suggestion": "请检查网络连接，或稍后重试。此步骤可以跳过，不影响其他分析。",
                "can_skip": True
            }
        
        # 6. 特定工具错误
        
        # CellTypist 相关错误
        if "celltypist" in error_lower or "cell_annotation" in tool_id:
            if "not installed" in error_lower or "import" in error_lower:
                return {
                    "user_message": f"细胞类型注释工具未安装：{step_name} 步骤需要 CellTypist 工具包。",
                    "error_category": "config_issue",
                    "suggestion": "此步骤可以跳过，不影响其他分析结果。如需使用，请联系管理员安装 CellTypist。",
                    "can_skip": True
                }
            elif "model" in error_lower or "下载" in error_lower:
                return {
                    "user_message": f"模型下载失败：{step_name} 步骤需要下载 CellTypist 模型，但下载失败。",
                    "error_category": "network_issue",
                    "suggestion": "此步骤可以跳过，不影响其他分析结果。您可以稍后手动下载模型或使用其他注释方法。",
                    "can_skip": True
                }
        
        # Path 变量错误（annotation.py 中的问题）
        if "path" in error_lower and ("referenced before assignment" in error_lower or "未定义" in error_lower):
            return {
                "user_message": f"内部错误：{step_name} 步骤遇到内部配置问题。",
                "error_category": "internal_error",
                "suggestion": "此错误已修复，请重试。如果问题仍然存在，可以跳过此步骤。",
                "can_skip": True
            }
        
        # 7. 数据质量问题（可忽略）
        if any(keyword in error_lower for keyword in [
            "quality", "质量", "low quality", "低质量", "filter", "过滤"
        ]):
            return {
                "user_message": f"数据质量提示：{step_name} 步骤检测到数据质量问题。",
                "error_category": "data_quality",
                "suggestion": "这可能是正常的数据过滤过程，不影响整体分析。",
                "can_skip": False
            }
        
        # 8. 默认错误（未知错误）
        return {
            "user_message": f"执行失败：{step_name} 步骤执行时遇到问题。",
            "error_category": "unknown_error",
            "suggestion": "请查看详细错误信息，或尝试跳过此步骤继续其他分析。",
            "can_skip": True,
            "technical_details": error  # 保留原始错误信息供调试
        }
    
    @staticmethod
    def format_step_error(
        step_result: Dict[str, Any],
        tool_id: str,
        step_name: str
    ) -> Dict[str, Any]:
        """
        格式化步骤错误信息
        
        Args:
            step_result: 步骤执行结果
            tool_id: 工具ID
            step_name: 步骤名称
        
        Returns:
            格式化后的错误信息
        """
        error = step_result.get("error", "") or step_result.get("message", "")
        if not error:
            return step_result
        
        formatted = ErrorFormatter.format_error(error, tool_id, step_name)
        
        # 更新 step_result
        step_result["user_message"] = formatted["user_message"]
        step_result["error_category"] = formatted["error_category"]
        step_result["suggestion"] = formatted["suggestion"]
        step_result["can_skip"] = formatted["can_skip"]
        if "technical_details" in formatted:
            step_result["technical_details"] = formatted["technical_details"]
        
        return step_result
