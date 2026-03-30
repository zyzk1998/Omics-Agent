"""
技能快车道（Skill Fast-Lane）：轻量级执行器，绕开 SOPPlanner。

通过 Orchestrator 识别 `[Skill_Route: tool_name]` 后，由本类完成：参数抽取 → 工具调用 → Markdown 汇总。
"""
from __future__ import annotations

import inspect
import json
import logging
import re
import time
import uuid
from typing import Any, AsyncIterator, Callable, Dict, List, Optional

from gibh_agent.core.llm_client import LLMClient
from gibh_agent.core.openai_tools import tool_names_to_openai_tools
from gibh_agent.core.stream_utils import THINK_CLOSE, THINK_OPEN, stream_from_llm_chunks
from gibh_agent.core.tool_registry import registry
from gibh_agent.core.utils import sanitize_for_json

logger = logging.getLogger(__name__)

_FASTA_SUFFIXES = (".fa", ".fasta", ".faa", ".ffn")


def _extract_inline_fasta_for_fallback(text: str) -> str:
    """无上传文件时从技能模板正文中提取 FASTA（兜底）。"""
    t = (text or "").strip()
    if not t:
        return ""
    if "```" in t:
        chunks = re.split(r"```\w*", t)
        for blob in chunks:
            b = blob.strip()
            if b.startswith(">"):
                return b if b.endswith("\n") else b + "\n"
            seq = re.sub(r"[^AUCGTNaucgtn]", "", b.replace("\n", ""))
            if len(seq) >= 4:
                return f">sequence\n{seq.upper()}\n"
    m = re.search(r">[^\n]+\r?\n[AUCGTNaucgtn\s]+", t, re.I)
    if m:
        return m.group(0).strip() + "\n"
    return ""


def _looks_like_smiles_token(s: str) -> bool:
    t = (s or "").strip()
    if len(t) < 4:
        return False
    if not re.match(r"^[A-Za-z0-9@+\-\[\]()=#\\.\\/,;%:]+$", t):
        return False
    return "(" in t or "=" in t or "[" in t or any(c.isdigit() for c in t) or (len(t) <= 30 and t[0].isupper())


def _extract_inline_smiles_for_fallback(text: str) -> str:
    """从技能模板正文提取 SMILES：反引号块、列表项「：」后、或独立化学式行。"""
    t = text or ""
    for m in re.finditer(r"`([A-Za-z0-9@+\-\[\]()=#\\.\\/,;%:]+)`", t):
        cand = m.group(1).strip()
        if _looks_like_smiles_token(cand):
            return cand.split()[0]
    for m in re.finditer(r"[:：]\s*`?([A-Za-z0-9@+\-\[\]()=#\\.\\/,;%]{5,})", t):
        cand = m.group(1).strip().rstrip("`").strip()
        if _looks_like_smiles_token(cand):
            return cand.split()[0]
    for ln in t.splitlines():
        s = ln.strip()
        if not s or s.startswith("#") or s.startswith("【") or s.startswith("*"):
            continue
        if s.startswith("-"):
            if "`" in s:
                inner = re.search(r"`([^`]+)`", s)
                if inner:
                    cand = inner.group(1).strip()
                    if _looks_like_smiles_token(cand):
                        return cand.split()[0]
            continue
        if len(s) >= 2 and re.match(r"^[A-Za-z0-9@+\-\[\]()=#\.\\/,;:]+$", s):
            if any(c in s for c in "()[]=#") or (len(s) <= 20 and s[0].isalpha()):
                return s.split()[0]
    return ""


def _extract_inline_pdb_for_fallback(text: str) -> str:
    t = text or ""
    if "```" in t:
        for blob in re.split(r"```\w*", t):
            b = blob.strip()
            if "ATOM" in b or "HEADER" in b:
                return b
    return ""


def _extract_inline_csv_table_for_fallback(text: str) -> str:
    lines: List[str] = []
    grab = False
    for ln in (text or "").splitlines():
        low = ln.lower()
        if "gene" in low and "score" in low and "," in ln:
            grab = True
        if grab:
            lines.append(ln)
            if len(lines) >= 12:
                break
    return "\n".join(lines) if lines else ""


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
        self._pending_done_tool: Optional[Dict[str, Any]] = None

    def consume_pending_done_tool(self) -> Optional[Dict[str, Any]]:
        """供 Orchestrator 在 `done` 事件中合并 `tool_result`（取后即清空）。"""
        bundle = self._pending_done_tool
        self._pending_done_tool = None
        return bundle

    def _emit(self, event_type: str, data: Dict[str, Any]) -> str:
        return self._format_sse(event_type, data)

    async def _iter_throttled_skill_thought(
        self,
        step_id: str,
        throttle_parts: List[str],
        st: Dict[str, Any],
        *,
        force: bool,
        max_chars: int = 48000,
    ):
        """按字数/时间节流输出 substep_delta，避免逐 token 数千条 SSE，又能在执行记录里流式看到思考。"""
        full = "".join(throttle_parts)[:max_chars]
        sent = int(st.get("sent", 0))
        if len(full) <= sent:
            return
        tail = full[sent:]
        now = time.monotonic()
        char_th, time_th = 320, 0.22
        if not force and len(tail) < char_th and (now - float(st.get("last", 0))) < time_th:
            return
        st["sent"] = len(full)
        st["last"] = now
        yield self._agentic_emit(
            {
                "action": "substep_delta",
                "step_id": step_id,
                "substep_id": "think_1",
                "delta": tail,
                "state": "running",
            }
        )

    def _agentic_emit(self, payload: Dict[str, Any]) -> str:
        """[AgenticLog] 行通过 status.content 透出，保证整条 SSE 合法。"""
        line = f"[AgenticLog] {json.dumps(payload, ensure_ascii=False, allow_nan=False, default=str)}"
        sse_state = str(payload.get("state") or "running")
        return self._emit("status", {"content": line, "state": sse_state})

    def _strip_json_noise(self, s: str) -> str:
        """去掉 think 残留、Markdown 围栏与常见非 JSON 前缀。"""
        t = (s or "").strip()
        if not t:
            return t
        if THINK_CLOSE in t:
            t = t.split(THINK_CLOSE, 1)[-1].strip()
        t = t.removeprefix("\ufeff")
        if t.startswith("```"):
            lines = t.split("\n")
            if len(lines) >= 2 and lines[0].strip().startswith("```"):
                lines = lines[1:]
            if lines and lines[-1].strip() == "```":
                lines = lines[:-1]
            t = "\n".join(lines).strip()
        for prefix in ("JSON:", "json:", "参数:", "输出:"):
            if t.lower().startswith(prefix.lower()):
                t = t[len(prefix) :].strip()
        return t

    def _json_loose_repair(self, blob: str) -> str:
        """仅做安全、高频修复：尾随逗号。"""
        b = blob.strip()
        b = re.sub(r",(\s*[\]\}])", r"\1", b)
        return b

    def _parse_tool_args_from_text(self, raw: str) -> Dict[str, Any]:
        """将 think 闭合标签之后的模型输出解析为工具参数字典（严格 JSON）。"""
        s = self._strip_json_noise(raw)
        if not s:
            return {}

        def _loads(blob: str) -> Optional[Dict[str, Any]]:
            for candidate in (blob, self._json_loose_repair(blob)):
                try:
                    out = json.loads(candidate)
                    return out if isinstance(out, dict) else None
                except json.JSONDecodeError:
                    continue
            return None

        if "```" in s:
            parts = s.split("```")
            for blob in parts:
                b = blob.strip()
                if b.lower().startswith("json"):
                    b = b[4:].strip()
                if b.startswith("{") and b.endswith("}"):
                    got = _loads(b)
                    if got is not None:
                        return got

        if "{" in s:
            start, end = s.find("{"), s.rfind("}")
            if start != -1 and end > start:
                got = _loads(s[start : end + 1])
                if got is not None:
                    return got

        got = _loads(s)
        return got if got is not None else {}

    def _messages_for_streaming_arg_extract(
        self, tool_name: str, user_query: str, file_paths: List[str]
    ) -> List[Dict[str, str]]:
        """流式填参：纯文本协议（think 标签 + JSON），不再使用 chat tools。"""
        meta = registry.get_metadata(tool_name)
        schema_hint = "{}"
        if meta:
            try:
                schema_hint = json.dumps(meta.args_schema.model_json_schema(), ensure_ascii=False, indent=2)
            except Exception:
                schema_hint = str(list(meta.args_schema.model_fields.keys()))
        clean = strip_skill_route_tag(user_query)
        fp_block = "\n".join(f"- {p}" for p in file_paths) if file_paths else "(无上传文件)"
        system = (
            "你是生物信息学助手的参数抽取模块。用户已通过 [Skill_Route] 指定工具。"
            f"先用 {THINK_OPEN}…{THINK_CLOSE} 写简短推理；闭合标签之后，输出必须且仅为一个符合 RFC 8259 的 JSON 对象，"
            "供程序 json.loads 直接解析（解析失败会导致技能失败）。"
            "JSON 规则：双引号包裹所有键名与字符串值；布尔必须小写 true/false；禁止尾逗号；禁止在 JSON 外再写任何字符。"
            "FASTA 作为字符串值时，换行必须转义为 \\n，写成一行合法 JSON（例如 \"fasta_content\": \">id\\nAUGC...\"）。"
            "禁止编造 file_path：路径必须逐字来自【可用本地文件绝对路径】列表；若列表为空或用户未上传，"
            "对 PySkills 类工具请改用内联字段（勿写 /app/uploads/demo_*.fa 等虚构路径）："
            "rnafold_analysis→fasta_content；sascore_analysis→smiles_text；pymol_analysis→pdb_content；gseapy_analysis→table_content。"
            "有上传文件时优先 file_path 并用列表中的绝对路径。"
            "若用户未单独写一行 SMILES，但正文中用反引号标出参考分子 SMILES（或列表项中带 SMILES），"
            "须将**第一个**可解析的 SMILES 填入 smiles_text，不得留空导致工具失败。"
        )
        user_msg = (
            f"工具名称: `{tool_name}`\n\n"
            f"参数 JSON Schema（供参考）:\n```json\n{schema_hint}\n```\n\n"
            f"【用户原文（可含 FASTA）】\n{clean}\n\n"
            f"【可用本地文件绝对路径】\n{fp_block}\n\n"
            f"输出要求：\n"
            f"1) {THINK_OPEN} 与 {THINK_CLOSE} 之间只写推理，不要在此放 JSON。\n"
            f"2) {THINK_CLOSE} 之后：第一非空字符必须是 `{{`，最后一行必须以 `}}` 结束；"
            "中间为单个 JSON 对象。不要用 ```json 代码块包裹；不要加「最终 JSON：」等前缀。\n"
            "3) 若含 sequence_or_path 且为 FASTA 文本，必须整段（含 > 头行）作为单个 JSON 字符串，内部换行用 \\n。"
        )
        return [
            {"role": "system", "content": system},
            {"role": "user", "content": user_msg},
        ]

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
            if tool_name == "rnafold_analysis":
                chosen: Optional[str] = None
                for p in file_paths:
                    pl = p.lower()
                    if pl.endswith((".fa", ".fasta", ".fna", ".faa", ".ffn")):
                        chosen = p
                        break
                out["file_path"] = chosen or file_paths[0]
            elif tool_name == "sascore_analysis":
                chosen_sm: Optional[str] = None
                for p in file_paths:
                    pl = p.lower()
                    if pl.endswith((".smi", ".txt", ".csv", ".tsv")):
                        chosen_sm = p
                        break
                out["file_path"] = chosen_sm or file_paths[0]
            elif tool_name == "pymol_analysis":
                chosen_pdb: Optional[str] = None
                for p in file_paths:
                    pl = p.lower()
                    if pl.endswith((".pdb", ".cif", ".mmcif", ".ent")):
                        chosen_pdb = p
                        break
                out["file_path"] = chosen_pdb or file_paths[0]
            elif tool_name == "gseapy_analysis":
                chosen_g: Optional[str] = None
                for p in file_paths:
                    pl = p.lower()
                    if pl.endswith((".csv", ".txt", ".tsv", ".gmt")):
                        chosen_g = p
                        break
                out["file_path"] = chosen_g or file_paths[0]
            else:
                out.setdefault("file_path", file_paths[0])

        if tool_name == "rnafold_analysis" and "fasta_content" in fields and not file_paths:
            ex = _extract_inline_fasta_for_fallback(clean)
            if ex:
                out["fasta_content"] = ex
        if tool_name == "sascore_analysis" and "smiles_text" in fields and not file_paths:
            exs = _extract_inline_smiles_for_fallback(clean)
            if exs:
                out["smiles_text"] = exs
        if tool_name == "pymol_analysis" and "pdb_content" in fields and not file_paths:
            exp = _extract_inline_pdb_for_fallback(clean)
            if exp:
                out["pdb_content"] = exp
        if tool_name == "gseapy_analysis" and "table_content" in fields and not file_paths:
            exc = _extract_inline_csv_table_for_fallback(clean)
            if exc:
                out["table_content"] = exc
        return out

    def _augment_args_from_user_text(
        self,
        tool_name: str,
        args: Dict[str, Any],
        user_query: str,
        file_paths: List[str],
    ) -> Dict[str, Any]:
        """LLM 留空内联字段时，用用户原文启发式补全，避免工具因空参失败。"""
        clean = strip_skill_route_tag(user_query or "")
        meta = registry.get_metadata(tool_name)
        fields = set(meta.args_schema.model_fields.keys()) if meta else set()

        def _blank(v: Any) -> bool:
            return v is None or (isinstance(v, str) and not v.strip())

        fp = (args.get("file_path") or "").strip() if isinstance(args.get("file_path"), str) else ""
        has_upload = bool(file_paths) or bool(fp)

        if tool_name == "sascore_analysis" and "smiles_text" in fields:
            if _blank(args.get("smiles_text")) and not has_upload:
                sm = _extract_inline_smiles_for_fallback(clean)
                if sm:
                    args["smiles_text"] = sm
        if tool_name == "pymol_analysis" and "pdb_content" in fields:
            if _blank(args.get("pdb_content")) and not has_upload:
                pdb = _extract_inline_pdb_for_fallback(clean)
                if pdb:
                    args["pdb_content"] = pdb
        if tool_name == "gseapy_analysis" and "table_content" in fields:
            if _blank(args.get("table_content")) and not has_upload:
                tb = _extract_inline_csv_table_for_fallback(clean)
                if tb:
                    args["table_content"] = tb
        return args

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
            "禁止编造 file_path：路径必须来自用户上传列表；若无上传，使用 fasta_content / smiles_text / pdb_content / table_content 等内联字段。"
            "function 参数字段须可被服务端严格校验：字符串用双引号、布尔为小写 true/false；"
            "FASTA/PDB/CSV 整段作为字符串时，换行在参数里须为 \\n。"
            "若用户仅在模板中列出带反引号的参考 SMILES，须任选第一个合法 SMILES 写入 smiles_text。"
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
            json_u = self._pick_result_url(result, "json_url")
            if html_u:
                markdown_lines.append(f"- [🌐 查看交互式图表 (HTML)]({html_u})")
            if csv_u:
                markdown_lines.append(f"- [📄 下载原始数据 (CSV)]({csv_u})")
            if zip_u:
                markdown_lines.append(f"- [📦 下载完整结果包 (ZIP)]({zip_u})")
            if json_u:
                markdown_lines.append(f"- [📋 下载完整结果 (JSON)]({json_u})")
            if not (html_u or csv_u or zip_u or json_u):
                markdown_lines.append(
                    "- （本次响应未返回 html_url / csv_url / zip_url / json_url 直链；若任务仍在排队，请稍后重试或查看下方扩展字段。）"
                )
            gen_seq = result.get("generated_sequence")
            if isinstance(gen_seq, str) and gen_seq.strip():
                markdown_lines.append("\n#### 🧪 生成序列\n\n```\n" + gen_seq.strip() + "\n```\n")
            extra = {
                k: v
                for k, v in result.items()
                if k
                not in (
                    "status",
                    "message",
                    "html_url",
                    "csv_url",
                    "zip_url",
                    "json_url",
                    "data",
                    "generated_sequence",
                )
            }
            data_block = result.get("data")
            if isinstance(data_block, dict):
                rest = {
                    kk: vv
                    for kk, vv in data_block.items()
                    if kk not in ("html_url", "csv_url", "zip_url", "json_url", "status", "message")
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
        异步生成器：产出 SSE（含 [AgenticLog] 结构化 status 与最终 message），由 Orchestrator 直接 yield。
        填参阶段使用 LLMClient.astream + stream_from_llm_chunks：think/reasoning 与 JSON 前的正文
        经节流后以 substep_delta 流式写入执行记录（非逐 token、不入 process_log 爆炸性快照）。
        """
        self._pending_done_tool = None
        step_id = f"skill_{uuid.uuid4().hex[:8]}"
        yield self._agentic_emit(
            {
                "action": "step_start",
                "step_id": step_id,
                "content": f"执行技能: {tool_name}",
                "state": "running",
            }
        )

        tool_fn = registry.get_tool(tool_name)
        if not tool_fn:
            md = f"### 技能失败\n\n未在 ToolRegistry 中找到工具 `{tool_name}`。"
            yield self._emit("message", {"content": md})
            yield self._agentic_emit(
                {
                    "action": "step_complete",
                    "step_id": step_id,
                    "state": "error",
                    "collapsed": True,
                }
            )
            return

        tool_kwargs: Dict[str, Any] = {}
        if self.llm_client:
            yield self._agentic_emit(
                {
                    "action": "substep_start",
                    "step_id": step_id,
                    "substep_id": "think_1",
                    "kind": "thought",
                    "state": "running",
                }
            )
            json_fragments: List[str] = []
            stream_ok = False
            try:
                messages = self._messages_for_streaming_arg_extract(tool_name, user_query, file_paths)
                stream = self.llm_client.astream(
                    messages,
                    temperature=0.1,
                    max_tokens=4096,
                    model=model_name,
                )
                # 思考过程：节流流式 substep_delta（~320 字或 ~0.22s 一批），避免逐 token 数千条 SSE；
                # message 轨仅在首个「{」之前 mirror 到思考区，不把参数 JSON 整段塞进执行记录。
                _throttle_buf: List[str] = []
                _flush_st: Dict[str, Any] = {"sent": 0, "last": time.monotonic()}
                _json_started = False
                async for event_type, data in stream_from_llm_chunks(stream, model_name=model_name):
                    if event_type == "thought":
                        ct = data.get("content") or ""
                        if ct:
                            _throttle_buf.append(ct)
                    elif event_type == "message":
                        mc = data.get("content") or ""
                        json_fragments.append(mc)
                        if not _json_started:
                            full_j = "".join(json_fragments)
                            if "{" not in full_j:
                                _throttle_buf.append(mc)
                            else:
                                bi = full_j.find("{")
                                prem = full_j[:bi]
                                already = len("".join(_throttle_buf))
                                if len(prem) > already:
                                    _throttle_buf.append(prem[already:])
                                _json_started = True
                    async for _sse in self._iter_throttled_skill_thought(
                        step_id, _throttle_buf, _flush_st, force=False
                    ):
                        yield _sse
                stream_ok = True
                async for _sse in self._iter_throttled_skill_thought(
                    step_id, _throttle_buf, _flush_st, force=True
                ):
                    yield _sse
            except Exception as e:
                logger.warning("SkillAgent 流式填参失败，回退 achat: %s", e, exc_info=True)
                yield self._agentic_emit(
                    {
                        "action": "substep_complete",
                        "step_id": step_id,
                        "substep_id": "think_1",
                        "state": "completed",
                    }
                )
                tool_kwargs = await self._llm_extract_args(tool_name, user_query, file_paths, model_name)
            if stream_ok:
                if not "".join(_throttle_buf).strip():
                    yield self._agentic_emit(
                        {
                            "action": "substep_delta",
                            "step_id": step_id,
                            "substep_id": "think_1",
                            "delta": "（模型未返回可展示的推理文本；已直接解析参数 JSON。）\n",
                            "state": "running",
                        }
                    )
                yield self._agentic_emit(
                    {
                        "action": "substep_complete",
                        "step_id": step_id,
                        "substep_id": "think_1",
                        "state": "completed",
                    }
                )
                raw_json = "".join(json_fragments)
                tool_kwargs = self._parse_tool_args_from_text(raw_json)
                if not tool_kwargs:
                    tool_kwargs = await self._llm_extract_args(
                        tool_name, user_query, file_paths, model_name
                    )

        if not tool_kwargs:
            args = self._fallback_args(tool_name, user_query, file_paths)
            logger.info("SkillAgent 使用回退参数: %s", args)
        else:
            args = dict(tool_kwargs)
        args = self._augment_args_from_user_text(tool_name, args, user_query, file_paths)
        logger.info("🧠 [SkillAgent] 最终执行参数: %s", args)

        yield self._agentic_emit(
            {
                "action": "substep_start",
                "step_id": step_id,
                "substep_id": "tool_1",
                "kind": "code",
                "content": "提取参数",
                "tool_args": args,
                "state": "completed",
            }
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

        step_final_state = "error" if result.get("status") == "error" else "completed"
        yield self._agentic_emit(
            {
                "action": "step_complete",
                "step_id": step_id,
                "state": step_final_state,
                "collapsed": True,
            }
        )
        if result.get("status") == "error":
            md_err = self._format_result_markdown(tool_name, result)
            yield self._emit("message", {"content": md_err})
        else:
            tr = sanitize_for_json(result)
            self._pending_done_tool = {"tool_name": tool_name, "tool_result": tr}
            # 独立事件：避免仅靠 done 合并时丢字段；前端可提前渲染右栏
            yield self._emit("skill_tool_result", {"tool_name": tool_name, "tool_result": tr})
            yield self._emit(
                "message",
                {
                    "content": "✅ **技能执行成功，请在右侧工作台查看详细报告。**",
                },
            )
