# -*- coding: utf-8 -*-
"""
MCP 插件目录（与 gibh_agent/tools 物理隔离）。

启动时调用 load_mcp_plugins() 以注册 MCP 工具到全局 ToolRegistry。
"""
import logging

logger = logging.getLogger(__name__)


def load_mcp_plugins() -> None:
    """导入 MCP 子模块以触发 @registry.register。"""
    try:
        from . import web_search_mcp  # noqa: F401
        logger.info("✅ [MCP] 已加载模块: web_search_mcp")
    except Exception as e:
        logger.warning("⚠️ [MCP] web_search_mcp 加载失败: %s", e, exc_info=True)
    try:
        from . import ncbi_mcp  # noqa: F401
        logger.info("✅ [MCP] 已加载模块: ncbi_mcp")
    except Exception as e:
        logger.warning("⚠️ [MCP] ncbi_mcp 加载失败: %s", e, exc_info=True)
