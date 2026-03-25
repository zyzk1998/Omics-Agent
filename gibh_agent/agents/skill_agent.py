"""
技能快车道（Skill Fast-Lane）：轻量级执行器，绕开 SOPPlanner。

通过 Orchestrator 识别 `[Skill_Route: tool_name]` 后，由本类完成：参数抽取 → 工具调用 → Markdown 汇总。
"""
from __future__ import annotations

import inspect
import json
import logging
import re
from typing import Any, AsyncIterator, Callable, Dict, List, Optional

from gibh_agent.core.llm_client import LLMClient
from gibh_agent.core.openai_tools import tool_names_to_openai_tools
from gibh_agent.core.tool_registry import registry

logger = logging.getLogger(__name__)

_FASTA_SUFFIXES = (".fa", ".fasta", ".faa", ".ffn")


def strip_skill_route_tag(query: str) -> str:
    return re.sub(
        r"\[Skill_Route:\s*[a-zA-Z0-9_]+\]\s*",
        "",
        (query or "").strip(),
        count=1,
    ).strip()


class SkillAgent:
    """
    轻量级技能执行器：依赖 ToolRegistry + LLM Function Calling 填参，再调用注册工具。
    """

    def __init__(
        self,
        llm_client: Optional[LLMClient],
        format_sse: Callable[[str, Dict[str, Any]], str],
    ):
        self.llm_client = llm_client
        self._format_sse = format_sse

    def _emit(self, event_type: str, data: Dict[str, Any]) -> str:
        return self._format_sse(event_type, data)

    def _fallback_args(self, tool_name: str, user_query: str, file_paths: List[str]) -> Dict[str, Any]:
        """无 LLM 时的最小启发式（如 FASTA 路径 → sequence_or_path）。"""
        meta = registry.get_metadata(tool_name)
        if not meta:
            return {}
        fields = list(meta.args_schema.model_fields.keys())
        out: Dict[str, Any] = {}
        clean = strip_skill_route_tag(user_query)
        if "sequence_or_path" in fields:
            chosen: Optional[str] = None
            for p in file_paths:
                pl = p.lower()
                if pl.endswith(_FASTA_SUFFIXES):
                    chosen = p
                    break
            out["sequence_or_path"] = chosen or (file_paths[0] if file_paths else clean)
        if "file_path" in fields and file_paths:
            out.setdefault("file_path", file_paths[0])
        return out

    async def _llm_extract_args(
        self,
        tool_name: str,
        user_query: str,
        file_paths: List[str],
        model_name: str,
    ) -> Dict[str, Any]:
        if not self.llm_client:
            return {}
        tools = tool_names_to_openai_tools([tool_name])
        if not tools:
            return {}
        clean = strip_skill_route_tag(user_query)
        fp_block = "\n".join(f"- {p}" for p in file_paths) if file_paths else "(无上传文件)"
        system = (
            "你是生物信息学助手的参数抽取模块。用户已通过 [Skill_Route] 指定工具，"
            "你必须调用给定的 function，并从用户描述与文件路径中推断合法参数。"
            "禁止编造文件路径；若序列表在正文中，应填入 sequence_or_path。"
        )
        user_msg = (
            f"【用户原文（可含 FASTA）】\n{clean}\n\n"
            f"【可用本地文件绝对路径】\n{fp_block}\n"
            f"请调用工具 `{tool_name}` 并给出完整参数。"
        )
        messages = [
            {"role": "system", "content": system},
            {"role": "user", "content": user_msg},
        ]
        try:
            completion = await self.llm_client.achat(
                messages,
                tools=tools,
                tool_choice={"type": "function", "function": {"name": tool_name}},
                temperature=0.1,
                max_tokens=4096,
                model=model_name,
            )
        except Exception as e:
            logger.warning("SkillAgent LLM 填参失败: %s", e, exc_info=True)
            return {}
        msg = completion.choices[0].message
        tcs = getattr(msg, "tool_calls", None) or []
        if not tcs:
            raw = (getattr(msg, "content", None) or "").strip()
            if raw and "```" in raw:
                try:
                    blob = raw.split("```", 2)[1]
                    if blob.lstrip().startswith("json"):
                        blob = blob.lstrip()[4:]
                    return json.loads(blob.strip())
                except Exception:
                    pass
            return {}
        fn = getattr(tcs[0], "function", None)
        raw_args = (getattr(fn, "arguments", None) if fn is not None else None) or "{}"
        try:
            return json.loads(raw_args)
        except json.JSONDecodeError:
            return {}

    def _pick_result_url(self, result: Dict[str, Any], key: str) -> str:
        """顶层或嵌套 data 中的下载链接（兼容不同工具返回形状）。"""
        v = result.get(key)
        if isinstance(v, str) and v.strip():
            return v.strip()
        data = result.get("data")
        if isinstance(data, dict):
            inner = data.get(key)
            if isinstance(inner, str) and inner.strip():
                return inner.strip()
        return ""

    def _format_result_markdown(self, tool_name: str, result: Dict[str, Any]) -> str:
        if not isinstance(result, dict):
            return f"### `{tool_name}`\n\n```\n{result!r}\n```"
        if result.get("status") == "error":
            err = (result.get("message") or result.get("error") or str(result)).strip() or "(无具体说明)"
            quoted = "\n".join(f"> {line}" if line else ">" for line in err.splitlines()) or f"> {err}"
            return (
                "### ⚠️ 技能执行未完成\n\n"
                "**系统提示：**\n\n"
                f"{quoted}\n\n"
                f"_工具：`{tool_name}`_"
            )

        if result.get("status") == "success":
            markdown_lines = [f"### 🧬 技能 `{tool_name}` 执行成功\n"]
            if "message" in result and result.get("message"):
                markdown_lines.append(f"{result['message']}\n")
            markdown_lines.append("#### 📊 结果下载与可视化\n")
            html_u = self._pick_result_url(result, "html_url")
            csv_u = self._pick_result_url(result, "csv_url")
            zip_u = self._pick_result_url(result, "zip_url")
            if html_u:
                markdown_lines.append(f"- [🌐 查看交互式图表 (HTML)]({html_u})")
            if csv_u:
                markdown_lines.append(f"- [📄 下载原始数据 (CSV)]({csv_u})")
            if zip_u:
                markdown_lines.append(f"- [📦 下载完整结果包 (ZIP)]({zip_u})")
            if not (html_u or csv_u or zip_u):
                markdown_lines.append(
                    "- （本次响应未返回 html_url / csv_url / zip_url 直链；若任务仍在排队，请稍后重试或查看下方扩展字段。）"
                )
            extra = {
                k: v
                for k, v in result.items()
                if k not in ("status", "message", "html_url", "csv_url", "zip_url", "data")
            }
            data_block = result.get("data")
            if isinstance(data_block, dict):
                rest = {
                    kk: vv
                    for kk, vv in data_block.items()
                    if kk not in ("html_url", "csv_url", "zip_url", "status", "message")
                }
                if rest:
                    extra = {**extra, "data": rest}
            if extra:
                markdown_lines.append("\n<details><summary>其它字段</summary>\n\n```json\n")
                try:
                    markdown_lines.append(json.dumps(extra, ensure_ascii=False, indent=2))
                except Exception:
                    markdown_lines.append(str(extra))
                markdown_lines.append("\n```\n</details>")
            return "\n".join(markdown_lines)

        lines = [f"### 技能 `{tool_name}` 已完成\n"]
        msg = result.get("message")
        if msg:
            lines.append(f"{msg}\n")
        for label, key in (
            ("HTML 报告", "html_url"),
            ("CSV 数据", "csv_url"),
            ("ZIP 打包", "zip_url"),
        ):
            url = self._pick_result_url(result, key)
            if url:
                lines.append(f"- [🔗 {label}]({url})")
        extra = {k: v for k, v in result.items() if k not in ("status", "message", "html_url", "csv_url", "zip_url")}
        if extra:
            lines.append("\n<details><summary>其它字段</summary>\n\n```json\n")
            try:
                lines.append(json.dumps(extra, ensure_ascii=False, indent=2))
            except Exception:
                lines.append(str(extra))
            lines.append("\n```\n</details>")
        return "\n".join(lines)

    async def execute_skill(
        self,
        tool_name: str,
        user_query: str,
        file_paths: List[str],
        *,
        model_name: str,
    ) -> AsyncIterator[str]:
        """
        异步生成器：依次产出 SSE 字符串（status / message），由 Orchestrator 直接 yield。
        """
        yield self._emit(
            "status",
            {"content": "[ProcessLog] 正在解析技能参数...", "state": "running"},
        )

        tool_fn = registry.get_tool(tool_name)
        if not tool_fn:
            md = f"### 技能失败\n\n未在 ToolRegistry 中找到工具 `{tool_name}`。"
            yield self._emit("message", {"content": md})
            return

        tool_kwargs = await self._llm_extract_args(tool_name, user_query, file_paths, model_name)
        logger.info("🧠 [SkillAgent] LLM 提取的参数: %s", tool_kwargs)
        if not tool_kwargs:
            args = self._fallback_args(tool_name, user_query, file_paths)
            logger.info("SkillAgent 使用回退参数: %s", args)
        else:
            args = tool_kwargs
        logger.info("🧠 [SkillAgent] 最终执行参数: %s", args)

        yield self._emit(
            "status",
            {"content": f"[ProcessLog] 正在执行技能: {tool_name} ...", "state": "running"},
        )

        try:
            if inspect.iscoroutinefunction(tool_fn):
                result = await tool_fn(**args)
            else:
                result = tool_fn(**args)
        except Exception as e:
            logger.exception("SkillAgent 工具调用异常: %s", e)
            result = {"status": "error", "message": str(e)}

        if not isinstance(result, dict):
            result = {"status": "success", "data": result, "message": str(result)}

        md = self._format_result_markdown(tool_name, result)
        yield self._emit("message", {"content": md})
