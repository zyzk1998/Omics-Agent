"""
错误信息格式化器

将面向程序员的错误信息转换为用户友好的提示。
"""
import logging
import re
from typing import Dict, Any, Optional, Tuple

logger = logging.getLogger(__name__)

# 精确模式 → (用户向文案, 错误类别, 建议, 可跳过)
_ERROR_PATTERN_MAP: Tuple[Tuple[re.Pattern[str], Dict[str, Any]], ...] = (
    (
        re.compile(
            r"(file not found|no such file|errno\s*2|"
            r"image not found|path does not exist|文件不存在|找不到文件)",
            re.I,
        ),
        {
            "user_message": (
                "很抱歉，系统未能定位到所需的文件。请检查您是否上传了正确格式的文件"
                "（如 .nii、.nii.gz、.dcm），并确认文件未损坏。"
            ),
            "error_category": "data_issue",
            "suggestion": "请在工作区重新上传数据，或确认上一步已成功生成输出文件。",
            "can_skip": False,
        },
    ),
    (
        re.compile(
            r"(unsupported format|invalid format|无法识别.*格式|格式不支持|"
            r"not a valid|cannot read.*format)",
            re.I,
        ),
        {
            "user_message": "该文件格式暂不被支持，建议检查数据类型是否与当前分析任务匹配。",
            "error_category": "data_issue",
            "suggestion": "请确认文件扩展名与内容一致（如 NIfTI、DICOM、CSV），必要时转换格式后重试。",
            "can_skip": False,
        },
    ),
    (
        re.compile(r"(memoryerror|out of memory|内存溢出|内存不足)", re.I),
        {
            "user_message": "数据量过大导致内存溢出，建议先对数据进行降采样或分批处理。",
            "error_category": "resource_issue",
            "suggestion": "可尝试缩小 ROI、降低分辨率，或拆分队列后分批分析。",
            "can_skip": True,
        },
    ),
)


class ErrorFormatter:
    """错误信息格式化器"""

    @staticmethod
    def _match_pattern_map(error: str) -> Optional[Dict[str, Any]]:
        for pattern, spec in _ERROR_PATTERN_MAP:
            if pattern.search(error):
                out = dict(spec)
                out["technical_details"] = error
                return out
        return None

    @staticmethod
    def format_error(
        error: str,
        tool_id: str,
        step_name: str,
        error_type: Optional[str] = None,
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
        mapped = ErrorFormatter._match_pattern_map(error)
        if mapped:
            return mapped

        # 0. 宿主环境 / 组学 CLI 缺失（不可静默「跳过」）
        if any(
            tok in error_lower
            for tok in (
                "前置条件",
                "未满足",
                "host_env_prereq",
                "fatal_env_prereq",
                "missing from path",
                "不在注册表",
            )
        ) or any(tok in error for tok in ("GIBH_REF_", "PATH", "宿主环境")):
            return {
                "user_message": error if len(error.strip()) > 12 else f"宿主运行环境不满足：{step_name} 无法执行。",
                "error_category": "host_env_prereq",
                "suggestion": "请为后端 Python 进程配置与运维终端一致的 Conda PATH，并设置 GIBH_REF_HG38 等变量；参见 scripts/install_omics_real_env.sh。",
                "can_skip": False,
                "technical_details": error,
            }
        if any(
            keyword in error_lower
            for keyword in [
                "file not found",
                "文件不存在",
                "path does not exist",
                "no such file",
                "empty",
                "空",
                "no data",
                "没有数据",
                "missing",
                "缺失",
            ]
        ):
            return {
                "user_message": f"数据文件问题：{step_name} 步骤无法找到或读取数据文件。",
                "error_category": "data_issue",
                "suggestion": "请检查数据文件路径是否正确，文件是否存在且可读。",
                "can_skip": False,
                "technical_details": error,
            }

        # 2. 配置问题类错误
        if any(
            keyword in error_lower
            for keyword in [
                "not installed",
                "未安装",
                "import error",
                "导入错误",
                "module not found",
                "模块未找到",
                "missing dependency",
                "缺少依赖",
            ]
        ):
            return {
                "user_message": f"依赖缺失：{step_name} 步骤需要额外的软件包或工具。",
                "error_category": "config_issue",
                "suggestion": "请联系管理员安装必要的依赖包，或跳过此步骤。",
                "can_skip": True,
                "technical_details": error,
            }

        # 3. 参数问题类错误
        if any(
            keyword in error_lower
            for keyword in [
                "invalid parameter",
                "参数无效",
                "parameter error",
                "参数错误",
                "missing parameter",
                "缺少参数",
                "required parameter",
                "必需参数",
            ]
        ):
            return {
                "user_message": f"参数配置问题：{step_name} 步骤的参数配置不正确。",
                "error_category": "config_issue",
                "suggestion": "请检查步骤参数设置，或使用默认参数重试。",
                "can_skip": False,
                "technical_details": error,
            }

        # 4. 内存/资源问题
        if any(
            keyword in error_lower
            for keyword in [
                "memory",
                "内存",
                "out of memory",
                "内存不足",
                "resource",
                "资源",
            ]
        ):
            return {
                "user_message": f"资源不足：{step_name} 步骤需要更多内存或计算资源。",
                "error_category": "resource_issue",
                "suggestion": "请尝试减小数据规模，或联系管理员增加系统资源。",
                "can_skip": True,
                "technical_details": error,
            }

        # 5. 网络/下载问题（如 CellTypist 模型下载）
        if any(
            keyword in error_lower
            for keyword in [
                "download",
                "下载",
                "network",
                "网络",
                "connection",
                "连接",
                "timeout",
                "超时",
                "无法下载",
            ]
        ):
            return {
                "user_message": f"网络问题：{step_name} 步骤需要下载资源，但网络连接失败。",
                "error_category": "network_issue",
                "suggestion": "请检查网络连接，或稍后重试。此步骤可以跳过，不影响其他分析。",
                "can_skip": True,
                "technical_details": error,
            }

        # 6. 特定工具错误

        # CellTypist 相关错误
        if "celltypist" in error_lower or "cell_annotation" in tool_id:
            if "not installed" in error_lower or "import" in error_lower:
                return {
                    "user_message": f"细胞类型注释工具未安装：{step_name} 步骤需要 CellTypist 工具包。",
                    "error_category": "config_issue",
                    "suggestion": "此步骤可以跳过，不影响其他分析结果。如需使用，请联系管理员安装 CellTypist。",
                    "can_skip": True,
                    "technical_details": error,
                }
            if "model" in error_lower or "下载" in error_lower:
                return {
                    "user_message": f"模型下载失败：{step_name} 步骤需要下载 CellTypist 模型，但下载失败。",
                    "error_category": "network_issue",
                    "suggestion": "此步骤可以跳过，不影响其他分析结果。您可以稍后手动下载模型或使用其他注释方法。",
                    "can_skip": True,
                    "technical_details": error,
                }

        # Path 变量错误（annotation.py 中的问题）
        if "path" in error_lower and (
            "referenced before assignment" in error_lower or "未定义" in error_lower
        ):
            return {
                "user_message": f"内部错误：{step_name} 步骤遇到内部配置问题。",
                "error_category": "internal_error",
                "suggestion": "此错误已修复，请重试。如果问题仍然存在，可以跳过此步骤。",
                "can_skip": True,
                "technical_details": error,
            }

        # 7. 数据质量问题（可忽略）
        if any(
            keyword in error_lower
            for keyword in [
                "quality",
                "质量",
                "low quality",
                "低质量",
                "filter",
                "过滤",
            ]
        ):
            return {
                "user_message": f"数据质量提示：{step_name} 步骤检测到数据质量问题。",
                "error_category": "data_quality",
                "suggestion": "这可能是正常的数据过滤过程，不影响整体分析。",
                "can_skip": False,
                "technical_details": error,
            }

        # 8. 分组/标签缺失（PLS-DA、ROC、差异分析等需分组）
        if any(
            keyword in error
            for keyword in [
                "分组",
                "分组标签",
                "至少两个",
                "group",
                "label_col",
                "metadata",
                "样本分组",
                "分类标签",
                "pls-da",
                "pls_da",
                "两组的样本",
                "label column",
            ]
        ) or "需要至少两" in error or "至少 2" in error:
            return {
                "user_message": error
                if len(error.strip()) > 10
                else f"数据分组问题：{step_name} 需要至少两个有效分组才能执行（如 PLS-DA、ROC）。请上传或指定样本分组信息。",
                "error_category": "data_issue",
                "suggestion": "请检查是否上传了样本分组文件 (Metadata)，或在参数中正确指定分组列/标签列。",
                "can_skip": True,
                "technical_details": error,
            }

        # 9. 工具已返回用户向文案（中文句子）：优先保留，避免被 generic 覆盖
        if len(error.strip()) >= 20 and any(c in error for c in ["请", "需要", "检查", "上传", "指定", "正确"]):
            return {
                "user_message": error,
                "error_category": "data_issue",
                "suggestion": "请根据上述提示检查数据或参数后重试，或跳过此步骤继续其他分析。",
                "can_skip": True,
                "technical_details": error,
            }

        # 10. 默认错误（未知错误）
        return {
            "user_message": f"执行失败：{step_name} 步骤执行时遇到问题。",
            "error_category": "unknown_error",
            "suggestion": "请查看详细错误信息，或尝试跳过此步骤继续其他分析。",
            "can_skip": True,
            "technical_details": error,
        }

    @staticmethod
    def format_step_error(
        step_result: Dict[str, Any],
        tool_id: str,
        step_name: str,
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

    @staticmethod
    def format_tool_exception_message(
        exc: BaseException,
        tool_id: str = "",
        step_name: str = "工具执行",
    ) -> Dict[str, Any]:
        """ReAct / 直连工具调用异常 → 用户向文案 + 技术日志。"""
        raw = str(exc).strip() or exc.__class__.__name__
        formatted = ErrorFormatter.format_error(raw, tool_id, step_name, error_type="exception")
        return formatted
