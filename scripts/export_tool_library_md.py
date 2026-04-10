#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
将 ToolRegistry 全量导出为仓库根目录 `工具库.md`，并附加超算 MCP 22 工具固定清单
（docs/hpc_mcp_tools_catalog.json，与网关 GET /tools 对齐）。

用法（在仓库根目录）:
    PYTHONPATH=. python3 scripts/export_tool_library_md.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

import gibh_agent.tools  # noqa: F401, E402 — 触发 @registry.register
from gibh_agent.core.tool_registry import registry  # noqa: E402

OUT_MD = REPO_ROOT / "工具库.md"
CATALOG_JSON = REPO_ROOT / "docs" / "hpc_mcp_tools_catalog.json"

ASK_HUMAN_APPENDIX = r"""
## 附录：DeepReAct 虚拟工具（不在 Registry）

`ask_human_for_clarification` — 由 `orchestrator._stream_deep_react_chat` 追加到 `tools` 列表，无独立 Registry 条目。

```json
{
  "type": "function",
  "function": {
    "name": "ask_human_for_clarification",
    "description": "当需要用户确认、选择或补充信息时调用；系统将暂停并等待用户下一轮输入。若有明确选项，请通过 options 传入字符串数组。",
    "parameters": {
      "type": "object",
      "properties": {
        "question": { "type": "string", "description": "要向用户说明或询问的内容" },
        "options": {
          "type": "array",
          "items": { "type": "string" },
          "description": "可选的快捷回复选项，供前端渲染为可点击按钮"
        }
      }
    }
  }
}
```
"""


def main() -> None:
    tools_json = registry.get_all_tools_json()
    n_reg = len(tools_json)

    if not CATALOG_JSON.is_file():
        raise SystemExit(f"缺少 {CATALOG_JSON}，请先创建 MCP 工具清单。")
    with CATALOG_JSON.open("r", encoding="utf-8") as f:
        mcp_catalog = json.load(f)
    if not isinstance(mcp_catalog, list):
        raise SystemExit("hpc_mcp_tools_catalog.json 须为 JSON 数组")
    n_mcp = len(mcp_catalog)

    body_json = json.dumps(tools_json, ensure_ascii=False, indent=2)
    mcp_block = json.dumps(mcp_catalog, ensure_ascii=False, indent=2)

    lines = [
        "# 工具库（ToolRegistry JSON 导出）",
        "",
        "本文档由 `scripts/export_tool_library_md.py` 生成，反映**当前 Python 进程成功 import 并注册**的 Registry 工具。",
        "",
        "**说明：**",
        "",
        "1. **OpenAI / DeepReAct 形态**：各工具的 `args_schema` 为 Pydantic `model_json_schema()` 输出，与 `tool_names_to_openai_tools` 使用的 `parameters` 同源。",
        "2. **MCP 超算工具 `hpc_mcp_*`**：若在运行时由 `HPCMCPManager.connect()` 从网关同步，会出现在 Registry 列表中；未连接网关时 Registry 条数可能少于线上。",
        "   **固定附录** `docs/hpc_mcp_tools_catalog.json` 列出与网关 `tools/list` 对齐的 22 个远端工具及对应 OpenAI 函数名，供检索与 prompt 对齐。",
        "3. **未纳入 Registry 的虚拟工具**（仅 DeepReAct 编排器拼接）：`ask_human_for_clarification`，见文末附录。",
        "",
        f"**Registry 工具总数：** {n_reg}",
        f"**附录 MCP 清单条数（HPC+工作站）：** {n_mcp}",
        "",
        "---",
        "",
        "## 完整 JSON 数组（Registry）",
        "",
        "以下为一组 JSON 对象，每项含 `name`, `description`, `category`, `output_type`, `args_schema`。",
        "",
        "```json",
        body_json,
        "```",
        "",
        "---",
        "",
        "## 附录：超算 / 工作站 MCP 工具固定清单（OpenAI 名 `hpc_mcp_*`）",
        "",
        "与 `docs/hpc_mcp_tools_catalog.json` 一致；动态 schema 以便 `GET {MCP_GATEWAY_URL}/tools` 为准。",
        "",
        "```json",
        mcp_block,
        "```",
        "",
        ASK_HUMAN_APPENDIX.strip(),
        "",
        "*生成命令：`PYTHONPATH=. python3 scripts/export_tool_library_md.py`*",
        "",
    ]

    OUT_MD.write_text("\n".join(lines), encoding="utf-8")
    print(f"Wrote {OUT_MD} (registry={n_reg}, mcp_catalog={n_mcp})")


if __name__ == "__main__":
    main()
