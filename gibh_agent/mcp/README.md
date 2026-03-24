# MCP 插件（与 `gibh_agent/tools` 隔离）

- 入口：`load_mcp_plugins()`（在 `server.py` 启动与 `executor` 中调用）
- 工具仍注册到全局 `ToolRegistry`，供编排器在聊天模式下以 OpenAI `tools` 形式挂载（如 `mcp_web_search`）

## 全网搜索

- 默认：`ddgs` 包（零配置，已实弹测试）
- 备选：环境变量见仓库根目录 `.env.example`（Tavily / Serper / Bing）
