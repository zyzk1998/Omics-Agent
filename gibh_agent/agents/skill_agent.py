"""
技能快车道（Skill Fast-Lane）：轻量级执行器，绕开 SOPPlanner。

通过 Orchestrator 识别 `[Skill_Route: tool_name]` 后，由本类完成：参数抽取 → 工具调用 → Markdown 汇总。
"""
from __future__ import annotations

import asyncio
import inspect
import json
import logging
import queue
import re
import time
import uuid
from typing import Any, AsyncIterator, Callable, Dict, List, Optional

from gibh_agent.core.llm_client import LLMClient
from gibh_agent.core.openai_tools import tool_names_to_openai_tools
from gibh_agent.core.stream_utils import THINK_CLOSE, THINK_OPEN, stream_from_llm_chunks
from gibh_agent.core.tool_registry import registry
from gibh_agent.core.tool_output_memory import build_tool_output_memory_text
from gibh_agent.core.tool_stream_log import tool_log_emitter_scope
from gibh_agent.core.utils import sanitize_for_json

logger = logging.getLogger(__name__)

# 双保险：快车道解析工具名前，确保 GI 等「侧车模块」已完成 @registry.register
# （与 gibh_agent/tools/chem_rdkit_tools.py 末尾侧车一致；避免部分入口未 import gibh_agent.tools 全包）
try:
    import importlib

    importlib.import_module("gibh_agent.tools.chem_gi_absorption_tools")
    importlib.import_module("gibh_agent.tools.chem_misc_tools")
    importlib.import_module("gibh_agent.tools.chem_rdkit_batch2_tools")
except Exception:
    logger.warning(
        "chem_gi_absorption_tools / chem_misc_tools 预加载失败（若仅使用其它化学工具可忽略）",
        exc_info=True,
    )
try:
    import importlib as _importlib_bioml

    _importlib_bioml.import_module("gibh_agent.tools.crispr_cas9_tool")
    _importlib_bioml.import_module("gibh_agent.tools.bioml_batch3_tools")
except Exception:
    logger.warning(
        "crispr_cas9_tool / bioml_batch3_tools 预加载失败（快车道 CRISPR/第三批生医仿真将不可用）",
        exc_info=True,
    )
try:
    from gibh_agent.skills import bootstrap_skills

    bootstrap_skills()
except Exception:
    logger.warning("gibh_agent/skills BaseSkill 预加载失败（首发核心技能快车道可能不可用）", exc_info=True)

# 首发核心技能 tool_id（与 gibh_agent/skills/skill_*.py · skill_id 一致）
_LAUNCH_SMILES_TOOLS = frozenset(
    {
        "pubchem_smiles_to_cid",
        "rdkit_3d_mol_render",
        "chembl_similar_molecules",
        "rdkit_substructure_search",
    }
)
_LAUNCH_QUERY_TOOLS = frozenset({"chembl_drug_search", "chipatlas_experiment_search"})
_LAUNCH_SEQUENCE_TOOLS = frozenset({"nucleotide_sequence_blast", "protein_sequence_blast"})

try:
    from gibh_agent.skills.launch_skill_demos import (
        LAUNCH_SKILL_DEMO_ARGS,
        apply_launch_demo_defaults as _apply_launch_demo_defaults,
    )
except Exception:
    LAUNCH_SKILL_DEMO_ARGS = {}
    _apply_launch_demo_defaults = None  # type: ignore


def _extract_json_block_from_prompt(text: str) -> Dict[str, Any]:
    """从 prompt_template 内 ```json ... ``` 广场一键体验块解析工具参数。"""
    if not (text or "").strip():
        return {}
    m = re.search(r"```(?:json)?\s*\n(\{.*?\})\s*```", text, re.DOTALL | re.IGNORECASE)
    if not m:
        return {}
    try:
        obj = json.loads(m.group(1))
        return obj if isinstance(obj, dict) else {}
    except json.JSONDecodeError:
        return {}


def _merge_launch_tool_args(tool_name: str, args: Dict[str, Any], user_query: str) -> Dict[str, Any]:
    """合并模板 JSON 演示块 + launch_skill_demos 默认值（不覆盖非空用户参数）。"""
    if tool_name not in LAUNCH_SKILL_DEMO_ARGS or _apply_launch_demo_defaults is None:
        return args
    merged = dict(args)
    for key, val in _extract_json_block_from_prompt(user_query).items():
        if _arg_blank(merged.get(key)):
            merged[key] = val
    return _apply_launch_demo_defaults(tool_name, merged)

# 药物相似性技能（逻辑 ID drug_similarity）：双工具 persona，供填参 LLM 对齐科学叙事
DRUG_SIM_PERSONA_LIPINSKI = (
    "你是计算化学与成药性（drug-likeness）方向的助理：当前任务为 **Lipinski 五规则** 快筛，"
    "侧重口服小分子 **ADME 相关理化描述符**（分子量、logP、氢键供体/受体、可旋转键等）与 **Rule-of-Five 合规性** 解读，"
    "避免展开 PubChem 相似性检索话术。"
)
DRUG_SIM_PERSONA_SIMILARITY = (
    "你是化学信息学与药物发现方向的助理：当前任务为 **高通量结构相似性检索**，"
    "侧重 **分子指纹 / 相似度度量**、命中列表与 **类似物/先导系列** 的发现与聚类叙事，"
    "避免把重点写成单纯 Lipinski 合规判断。"
    "填参时勿随意改小检索面：未写明时保持 **similarity_threshold=0.5**、**max_results_per_db=100**，"
    "与原始 Skills 包 `find_similar.py` → `get_similar_from_all` 默认一致，以免 Top 命中缩水。"
)

_FASTA_SUFFIXES = (".fa", ".fasta", ".faa", ".ffn")


def _arg_blank(v: Any) -> bool:
    return v is None or (isinstance(v, str) and not v.strip())


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


def _extract_multi_smiles_for_matrix_fallback(text: str) -> str:
    """从用户正文、反引号或表格中抽取至少两条疑似 SMILES，供 ``chem_tanimoto_matrix``。"""
    t = text or ""
    hits: List[str] = []
    for m in re.finditer(r"`([A-Za-z0-9@+\-\[\]()=#\\.\\/,;%:]+)`", t):
        cand = m.group(1).strip()
        if _looks_like_smiles_token(cand):
            hits.append(cand.split()[0])
    dedup: List[str] = []
    seen: set[str] = set()
    for h in hits:
        if h not in seen:
            seen.add(h)
            dedup.append(h)
    if len(dedup) >= 2:
        return "\n".join(dedup)

    line_hits: List[str] = []
    for ln in t.splitlines():
        s = ln.strip()
        if not s or s.startswith("#"):
            continue
        if all(c in "-:| " for c in s):
            continue
        if s.startswith("|"):
            cells = [c.strip() for c in s.split("|")]
            for c in cells:
                if c and _looks_like_smiles_token(c):
                    line_hits.append(c.split()[0])
                    break
            continue
        if _looks_like_smiles_token(s):
            line_hits.append(s.split()[0])
    out_ln: List[str] = []
    seen2: set[str] = set()
    for h in line_hits:
        if h not in seen2:
            seen2.add(h)
            out_ln.append(h)
    if len(out_ln) >= 2:
        return "\n".join(out_ln)

    flat = t.replace("\n", " ").strip()
    if "|" in flat:
        parts = [p.strip() for p in flat.split("|") if p.strip()]
        ok = [p.split()[0] for p in parts if _looks_like_smiles_token(p)]
        if len(ok) >= 2:
            return "\n".join(ok)
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


# 人源性评估（OASis）技能广场默认演示序列（与 seed_skills / index.html 一致）
_OASIS_DEMO_HEAVY = (
    "EVQLVESGGGLVQPGGSLRLSCAASGFTFSSYAMSWVRQAPGKGLEWVSAISGSGGSTYYADSVKGRFTISRDNSKNTLYLQMNSLRAEDTAVYYCARDYGDYWGQGTLVTVSS"
)
# Open Babel 广场一键演示（须 len≥4 以满足 SMILES 启发式；勿用 CCO 等过短串）
_OPENBABEL_DEMO_SMILES = "CC(=O)Oc1ccccc1C(=O)O"

_OASIS_DEMO_LIGHT = (
    "DIQMTQSPSSLSASVGDRVTITCRASQDVNTAVAWYQQKPGKAPKLLIYSASFLYSGVPSRFSGSRSGTDFTLTISSLQPEDFATYYCQQHYTTPPTFGQGTKVEIK"
)
_AA_ONLY = re.compile(r"^[ACDEFGHIKLMNPQRSTVWY]+$", re.I)


def _extract_antibody_chains_for_fallback(text: str) -> Dict[str, str]:
    """从用户正文 / 技能模板中解析抗体重链、轻链氨基酸序列。"""
    raw = text or ""
    out: Dict[str, str] = {}

    # 1) 单行 JSON：antibody_sequences_json 或内联 {"heavy_chain":...,"light_chain":...}
    for blob_m in re.finditer(r"\{[^{}]*heavy[^{}]*light[^{}]*\}", raw, re.I | re.S):
        try:
            obj = json.loads(blob_m.group(0))
            hc = (
                obj.get("heavy_chain")
                or obj.get("heavy")
                or obj.get("H")
                or obj.get("heavyChain")
                or ""
            )
            lc = (
                obj.get("light_chain")
                or obj.get("light")
                or obj.get("L")
                or obj.get("lightChain")
                or ""
            )
            if isinstance(hc, str) and isinstance(lc, str):
                hc, lc = hc.strip().upper(), lc.strip().upper()
                if len(hc) >= 20 and len(lc) >= 20 and _AA_ONLY.match(hc) and _AA_ONLY.match(lc):
                    return {"heavy_chain": hc, "light_chain": lc}
        except (json.JSONDecodeError, TypeError):
            continue

    # 2) 显式键值：heavy_chain / light_chain（含反引号、中英文标签）
    kv_patterns = [
        (r"(?:heavy[_\s-]*chain|重链|H\s*链)\s*[:：=]\s*`?([ACDEFGHIKLMNPQRSTVWY]{20,})`?", "heavy_chain"),
        (r"(?:light[_\s-]*chain|轻链|L\s*链)\s*[:：=]\s*`?([ACDEFGHIKLMNPQRSTVWY]{20,})`?", "light_chain"),
    ]
    for pat, key in kv_patterns:
        m = re.search(pat, raw, re.I)
        if m:
            out[key] = m.group(1).strip().upper()

    if out.get("heavy_chain") and out.get("light_chain"):
        return out

    # 3) FASTA：>heavy / >light 或链类型在 header
    current: Optional[str] = None
    buf_h: List[str] = []
    buf_l: List[str] = []
    for ln in raw.splitlines():
        s = ln.strip()
        if not s:
            continue
        if s.startswith(">"):
            header = s[1:].lower()
            if any(x in header for x in ("heavy", "h_chain", "vh", "重链")):
                current = "h"
            elif any(x in header for x in ("light", "l_chain", "vl", "轻链", "kappa", "lambda")):
                current = "l"
            else:
                current = None
            continue
        seq = re.sub(r"[^A-Za-z]", "", s).upper()
        if len(seq) < 20 or not _AA_ONLY.match(seq):
            continue
        if current == "h":
            buf_h.append(seq)
        elif current == "l":
            buf_l.append(seq)
    if buf_h and buf_l:
        return {"heavy_chain": "".join(buf_h), "light_chain": "".join(buf_l)}

    # 4) 技能模板内嵌的两条长肽段（按长度降序取前两条）
    candidates: List[str] = []
    for m in re.finditer(r"`([ACDEFGHIKLMNPQRSTVWY]{80,})`", raw, re.I):
        candidates.append(m.group(1).upper())
    for m in re.finditer(r"\b([ACDEFGHIKLMNPQRSTVWY]{80,})\b", raw):
        candidates.append(m.group(1).upper())
    uniq: List[str] = []
    for c in candidates:
        if c not in uniq:
            uniq.append(c)
    if len(uniq) >= 2:
        uniq.sort(key=len, reverse=True)
        return {"heavy_chain": uniq[0], "light_chain": uniq[1]}

    norm = raw.replace(" ", "").upper()
    if _OASIS_DEMO_HEAVY in norm and _OASIS_DEMO_LIGHT in norm:
        return {"heavy_chain": _OASIS_DEMO_HEAVY, "light_chain": _OASIS_DEMO_LIGHT}

    return out


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

    def _drug_similarity_persona_prefix(self, tool_name: str) -> str:
        if tool_name == "lipinski_druglikeness":
            return DRUG_SIM_PERSONA_LIPINSKI + "\n\n"
        if tool_name == "drug_similarity_search":
            return DRUG_SIM_PERSONA_SIMILARITY + "\n\n"
        return ""

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
        persona = self._drug_similarity_persona_prefix(tool_name)
        system = (
            persona
            + "你是生物信息学助手的参数抽取模块。用户已通过 [Skill_Route] 指定工具。"
            f"先用 {THINK_OPEN}…{THINK_CLOSE} 写简短推理；闭合标签之后，输出必须且仅为一个符合 RFC 8259 的 JSON 对象，"
            "供程序 json.loads 直接解析（解析失败会导致技能失败）。"
            "JSON 规则：双引号包裹所有键名与字符串值；布尔必须小写 true/false；禁止尾逗号；禁止在 JSON 外再写任何字符。"
            "FASTA 作为字符串值时，换行必须转义为 \\n，写成一行合法 JSON（例如 \"fasta_content\": \">id\\nAUGC...\"）。"
            "禁止编造 file_path：路径必须逐字来自【可用本地文件绝对路径】列表；若列表为空或用户未上传，"
            "对 PySkills 类工具请改用内联字段（勿写 /app/uploads/demo_*.fa 等虚构路径）："
            "rnafold_analysis→fasta_content；sascore_analysis→smiles_text；lipinski_druglikeness→smiles_text；"
            "drug_similarity_search→smiles_text；chem_openbabel→smiles_text 或 file_path（二选一，禁止同时填）；"
            "广场默认演示：smiles_text=模板中阿司匹林 SMILES、compute_properties=true；格式转换则 compute_properties=false 且必填 out_format；"
            "pymol_analysis→pdb_content；gseapy_analysis→table_content。"
            "crispr_cas9_simulation→guides_text、target_sequence、cell_line（可选）、result_format、random_seed（可选）；"
            "无上传文件时从正文提取 gRNA 与靶序列填入 guides_text/target_sequence，禁止编造 file_path。"
            "immune_cell_isolation_simulation→simulation_json、output_format、random_seed（可选）。"
            "rna_secondary_structure_analysis→structure_json 或 file_path（二选一）、result_format。"
            "circular_dichroism_analysis→spectrum_json、output_format。"
            "itc_binding_thermodynamic_analysis→protein_conc、ligand_conc、itc_injections_json 或 file_path、cell_volume（可选）、temperature（可选）、volume_unit（可选）；"
            "cell_cycle_phase_duration_estimation→flow_cytometry_json、initial_estimates_json（可选 JSON，可留空使用默认初值）；"
            "protease_kinetics_analysis→kinetics_input_json；enzyme_kinetics_analysis→kinetics_input_json；"
            "protein_sequence_conservation_analysis→alignment_json（含 protein_sequences 数组）。"
            "protein_homology_structure_assessment→uniprot_id、blast_evalue（可选）、struct_threshold、ptm_threshold、top_n、test_mode（可选布尔）；"
            "antibody_humanness_oasis_evaluation→**必填** heavy_chain 与 light_chain（各为一条氨基酸串，无空格）；"
            "或 antibody_sequences_json 单行 JSON {\"heavy_chain\":\"...\",\"light_chain\":\"...\"}；"
            "禁止留空；用户未提供序列时从正文/FASTA/模板示例中提取，仍无则使用技能模板内演示 H/L 序列。"
            "antibody_name、scheme、cdr_definition、threshold 可选。"
            "有上传文件时优先 file_path 并用列表中的绝对路径。"
            "若用户未单独写一行 SMILES，但正文中用反引号标出参考分子 SMILES（或列表项中带 SMILES），"
            "须将**第一个**可解析的 SMILES 填入 smiles_text，不得留空导致工具失败。"
            "首发核心技能：chembl_drug_search→query；pubchem_smiles_to_cid/rdkit_3d_mol_render/chembl_similar_molecules→smiles；"
            "rdkit_substructure_search→substructure_smiles+target_smiles；rdkit_mol_format_convert→input_text；"
            "nucleotide_sequence_blast/protein_sequence_blast→sequence_text 或 sequence_or_path；"
            "mhc_epitope_search→mhc_allele 和/或 peptide_sequence；chipatlas_experiment_search→query 或 expid。"
            "首发技能广场模板内含 ```json``` 一键体验块：须优先解析该 JSON 填入工具参数；"
            "用户未改需求且字段仍空时，使用模板 JSON 或 launch_skill_demos 默认值，禁止留空必填项。"
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
            elif tool_name in (
                "sascore_analysis",
                "lipinski_druglikeness",
                "drug_similarity_search",
                "chem_aromaticity_perception",
                "chem_lipinski_five",
                "chem_functional_groups",
                "chem_kekulization",
                "chem_pattern_fingerprint",
            ):
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
            elif tool_name == "chem_openbabel":
                chosen_ob: Optional[str] = None
                for p in file_paths:
                    pl = p.lower()
                    if pl.endswith(
                        (".mol", ".sdf", ".smi", ".smiles", ".pdb", ".xyz", ".cml", ".mol2", ".cdx")
                    ):
                        chosen_ob = p
                        break
                out["file_path"] = chosen_ob or file_paths[0]
            elif tool_name == "chem_tanimoto_matrix":
                chosen_tm: Optional[str] = None
                for p in file_paths:
                    pl = p.lower()
                    if pl.endswith((".smi", ".txt", ".csv", ".tsv")):
                        chosen_tm = p
                        break
                out["file_path"] = chosen_tm or file_paths[0]
            else:
                out.setdefault("file_path", file_paths[0])

        if tool_name == "rnafold_analysis" and "fasta_content" in fields and not file_paths:
            ex = _extract_inline_fasta_for_fallback(clean)
            if ex:
                out["fasta_content"] = ex
        if tool_name in (
            "sascore_analysis",
            "lipinski_druglikeness",
            "drug_similarity_search",
            "chem_aromaticity_perception",
            "chem_lipinski_five",
            "chem_functional_groups",
            "chem_kekulization",
            "chem_pattern_fingerprint",
        ) and "smiles_text" in fields and not file_paths:
            exs = _extract_inline_smiles_for_fallback(clean)
            if exs:
                out["smiles_text"] = exs
        if tool_name == "chem_openbabel":
            if "compute_properties" in fields:
                out.setdefault("compute_properties", True)
            if "in_format" in fields:
                out.setdefault("in_format", "smi")
            if "smiles_text" in fields and not file_paths:
                exo = _extract_inline_smiles_for_fallback(clean)
                out["smiles_text"] = exo or _OPENBABEL_DEMO_SMILES
        if tool_name == "chem_tanimoto_matrix" and "smiles_text" in fields and not file_paths:
            exm = _extract_multi_smiles_for_matrix_fallback(clean)
            if exm:
                out["smiles_text"] = exm
        if tool_name == "pymol_analysis" and "pdb_content" in fields and not file_paths:
            exp = _extract_inline_pdb_for_fallback(clean)
            if exp:
                out["pdb_content"] = exp
        if tool_name == "gseapy_analysis" and "table_content" in fields and not file_paths:
            exc = _extract_inline_csv_table_for_fallback(clean)
            if exc:
                out["table_content"] = exc
        if tool_name == "rna_secondary_structure_analysis" and "file_path" in fields and file_paths:
            chosen_rna: Optional[str] = None
            for p in file_paths:
                pl = p.lower()
                if pl.endswith((".fa", ".fasta", ".fna", ".txt", ".db", ".dot")):
                    chosen_rna = p
                    break
            if chosen_rna:
                out["file_path"] = chosen_rna
            elif not out.get("file_path"):
                out["file_path"] = file_paths[0]
        if tool_name == "antibody_humanness_oasis_evaluation":
            pairs = _extract_antibody_chains_for_fallback(clean)
            if pairs:
                out["heavy_chain"] = pairs["heavy_chain"]
                out["light_chain"] = pairs["light_chain"]
            else:
                out["heavy_chain"] = _OASIS_DEMO_HEAVY
                out["light_chain"] = _OASIS_DEMO_LIGHT
            out.setdefault("scheme", "kabat")
            out.setdefault("cdr_definition", "kabat")
            out.setdefault("threshold", "relaxed")
            out.setdefault("antibody_name", "Antibody")
        if tool_name in _LAUNCH_SMILES_TOOLS and "smiles" in fields and not file_paths:
            ex_sm = _extract_inline_smiles_for_fallback(clean)
            if ex_sm:
                out["smiles"] = ex_sm
        if tool_name == "rdkit_substructure_search":
            if "substructure_smiles" in fields and _arg_blank(out.get("substructure_smiles")):
                out.setdefault("substructure_smiles", "c1ccccc1")
            if "target_smiles" in fields and _arg_blank(out.get("target_smiles")):
                ex_t = _extract_inline_smiles_for_fallback(clean)
                if ex_t:
                    out["target_smiles"] = ex_t
        if tool_name in _LAUNCH_QUERY_TOOLS and "query" in fields and _arg_blank(out.get("query")):
            qline = re.sub(r"\[Skill_Route:[^\]]+\]", "", clean).strip()
            if qline:
                out["query"] = qline[:200]
        if tool_name == "chipatlas_experiment_search" and "expid" in fields:
            m_acc = re.search(r"\b(SRX|GSM|SRP)\d+\b", clean, re.I)
            if m_acc and _arg_blank(out.get("expid")):
                out["expid"] = m_acc.group(0)
        if tool_name in _LAUNCH_SEQUENCE_TOOLS:
            if "sequence_text" in fields and _arg_blank(out.get("sequence_text")):
                seq = re.sub(r"[^A-Za-z]", "", clean)
                if len(seq) >= 8:
                    out["sequence_text"] = seq[:5000]
        if tool_name == "rdkit_mol_format_convert" and "input_text" in fields:
            if _arg_blank(out.get("input_text")):
                ex_in = _extract_inline_smiles_for_fallback(clean)
                if ex_in:
                    out["input_text"] = ex_in
        if tool_name == "mhc_epitope_search" and "mhc_allele" in fields:
            m_hla = re.search(r"HLA-[A-Z0-9\*:\-]+", clean, re.I)
            if m_hla and _arg_blank(out.get("mhc_allele")):
                out["mhc_allele"] = m_hla.group(0)
        return _merge_launch_tool_args(tool_name, out, clean)

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

        if tool_name in (
            "sascore_analysis",
            "lipinski_druglikeness",
            "drug_similarity_search",
            "chem_aromaticity_perception",
            "chem_lipinski_five",
            "chem_functional_groups",
            "chem_kekulization",
            "chem_pattern_fingerprint",
        ) and "smiles_text" in fields:
            if _blank(args.get("smiles_text")) and not has_upload:
                sm = _extract_inline_smiles_for_fallback(clean)
                if sm:
                    args["smiles_text"] = sm
        if tool_name == "chem_openbabel":
            # LLM 偶发使用 smiles 而非 smiles_text
            sm_alias = args.get("smiles")
            if isinstance(sm_alias, str) and sm_alias.strip() and _blank(args.get("smiles_text")):
                args["smiles_text"] = sm_alias.strip()
            args.pop("smiles", None)
            fp_ob = (args.get("file_path") or "").strip() if isinstance(args.get("file_path"), str) else ""
            st_ob = (args.get("smiles_text") or "").strip() if isinstance(args.get("smiles_text"), str) else ""
            if _blank(fp_ob) and _blank(st_ob):
                if file_paths:
                    chosen_ob: Optional[str] = None
                    for p in file_paths:
                        pl = p.lower()
                        if pl.endswith(
                            (".mol", ".sdf", ".smi", ".smiles", ".pdb", ".xyz", ".cml", ".mol2")
                        ):
                            chosen_ob = p
                            break
                    args["file_path"] = chosen_ob or file_paths[0]
                else:
                    sm2 = _extract_inline_smiles_for_fallback(clean)
                    args["smiles_text"] = sm2 or _OPENBABEL_DEMO_SMILES
            if "compute_properties" in fields:
                cp = args.get("compute_properties")
                if cp is None or (isinstance(cp, str) and not str(cp).strip()):
                    args["compute_properties"] = True
                elif isinstance(cp, str):
                    args["compute_properties"] = cp.strip().lower() in ("1", "true", "yes", "on")
            fp_ob2 = (args.get("file_path") or "").strip() if isinstance(args.get("file_path"), str) else ""
            st_ob2 = (args.get("smiles_text") or "").strip() if isinstance(args.get("smiles_text"), str) else ""
            if st_ob2 and fp_ob2:
                args.pop("file_path", None)
            if "in_format" in fields and _blank(args.get("in_format")) and st_ob2 and not fp_ob2:
                args["in_format"] = "smi"
        if tool_name == "chem_tanimoto_matrix":
            fp_tm = (args.get("file_path") or "").strip() if isinstance(args.get("file_path"), str) else ""
            st_tm = (args.get("smiles_text") or "").strip() if isinstance(args.get("smiles_text"), str) else ""
            if _blank(fp_tm) and _blank(st_tm) and file_paths:
                chosen_tm: Optional[str] = None
                for p in file_paths:
                    pl = p.lower()
                    if pl.endswith((".smi", ".txt", ".csv", ".tsv")):
                        chosen_tm = p
                        break
                args["file_path"] = chosen_tm or file_paths[0]
            elif _blank(fp_tm) and _blank(st_tm) and not file_paths:
                mm = _extract_multi_smiles_for_matrix_fallback(clean)
                if mm:
                    args["smiles_text"] = mm
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
        if tool_name == "crispr_cas9_simulation":
            if "cell_line" in fields and _blank(args.get("cell_line")):
                args["cell_line"] = "HEK293"
            if "result_format" in fields and _blank(args.get("result_format")):
                args["result_format"] = "markdown"
        if tool_name == "immune_cell_isolation_simulation":
            if "output_format" in fields and _blank(args.get("output_format")):
                args["output_format"] = "text"
        if tool_name == "rna_secondary_structure_analysis":
            if "result_format" in fields and _blank(args.get("result_format")):
                args["result_format"] = "json"
        if tool_name == "circular_dichroism_analysis":
            if "output_format" in fields and _blank(args.get("output_format")):
                args["output_format"] = "text"
        if tool_name == "antibody_humanness_oasis_evaluation":
            js = (args.get("antibody_sequences_json") or "").strip()
            if js and isinstance(js, str):
                try:
                    obj = json.loads(js)
                    if _blank(args.get("heavy_chain")):
                        args["heavy_chain"] = (
                            obj.get("heavy_chain") or obj.get("heavy") or obj.get("H") or ""
                        )
                    if _blank(args.get("light_chain")):
                        args["light_chain"] = (
                            obj.get("light_chain") or obj.get("light") or obj.get("L") or ""
                        )
                except json.JSONDecodeError:
                    pass
            for alias_h, alias_l in (
                ("heavy", "light"),
                ("Heavy_chain", "Light_chain"),
                ("H_chain", "L_chain"),
            ):
                if _blank(args.get("heavy_chain")) and isinstance(args.get(alias_h), str):
                    args["heavy_chain"] = args[alias_h]
                if _blank(args.get("light_chain")) and isinstance(args.get(alias_l), str):
                    args["light_chain"] = args[alias_l]
            pairs = _extract_antibody_chains_for_fallback(clean)
            if pairs:
                if _blank(args.get("heavy_chain")):
                    args["heavy_chain"] = pairs["heavy_chain"]
                if _blank(args.get("light_chain")):
                    args["light_chain"] = pairs["light_chain"]
            if _blank(args.get("heavy_chain")) or _blank(args.get("light_chain")):
                args["heavy_chain"] = _OASIS_DEMO_HEAVY
                args["light_chain"] = _OASIS_DEMO_LIGHT
            if "scheme" in fields and _blank(args.get("scheme")):
                args["scheme"] = "kabat"
            if "cdr_definition" in fields and _blank(args.get("cdr_definition")):
                args["cdr_definition"] = "kabat"
            if "threshold" in fields and _blank(args.get("threshold")):
                args["threshold"] = "relaxed"
            if "antibody_name" in fields and _blank(args.get("antibody_name")):
                args["antibody_name"] = "Antibody"
        if tool_name in _LAUNCH_SMILES_TOOLS and "smiles" in fields and _blank(args.get("smiles")) and not has_upload:
            sm_l = _extract_inline_smiles_for_fallback(clean)
            if sm_l:
                args["smiles"] = sm_l
        if tool_name == "rdkit_mol_format_convert" and "input_text" in fields and _blank(args.get("input_text")):
            ex_l = _extract_inline_smiles_for_fallback(clean)
            if ex_l:
                args["input_text"] = ex_l
        if tool_name in _LAUNCH_SEQUENCE_TOOLS and "sequence_text" in fields and _blank(args.get("sequence_text")):
            fa = _extract_inline_fasta_for_fallback(clean)
            if fa:
                args["sequence_text"] = re.sub(r"[^A-Za-z]", "", fa)[:5000]
        if tool_name == "mhc_epitope_search" and "mhc_allele" in fields and _blank(args.get("mhc_allele")):
            m_hla2 = re.search(r"HLA-[A-Z0-9\*:\-]+", clean, re.I)
            if m_hla2:
                args["mhc_allele"] = m_hla2.group(0)
        return _merge_launch_tool_args(tool_name, args, clean)

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
        persona = self._drug_similarity_persona_prefix(tool_name)
        system = (
            persona
            + "你是生物信息学助手的参数抽取模块。用户已通过 [Skill_Route] 指定工具，"
            "你必须调用给定的 function，并从用户描述与文件路径中推断合法参数。"
            "禁止编造 file_path：路径必须来自用户上传列表；若无上传，使用 fasta_content / smiles_text / pdb_content / table_content 等内联字段。"
            "lipinski_druglikeness 与 drug_similarity_search 与 sascore 相同：优先 file_path，否则 smiles_text。"
            "function 参数字段须可被服务端严格校验：字符串用双引号、布尔为小写 true/false；"
            "FASTA/PDB/CSV 整段作为字符串时，换行在参数里须为 \\n。"
            "若用户仅在模板中列出带反引号的参考 SMILES，须任选第一个合法 SMILES 写入 smiles_text。"
            "chem_openbabel：无上传时 smiles_text 必填（可取模板演示阿司匹林 SMILES），广场一键体验默认 compute_properties=true。"
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
        step_banner = f"执行技能: {tool_name}"
        if tool_name == "lipinski_druglikeness":
            step_banner = "执行技能: Lipinski 类药性快筛（逻辑技能 drug_similarity / PK & Rule-of-Five）"
        elif tool_name == "drug_similarity_search":
            step_banner = "执行技能: 结构相似性检索（逻辑技能 drug_similarity / 聚类与类似物发现）"
        yield self._agentic_emit(
            {
                "action": "step_start",
                "step_id": step_id,
                "content": step_banner,
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

        _remote_blast = tool_name in _LAUNCH_SEQUENCE_TOOLS and bool(
            args.get("use_remote", True)
        )
        if _remote_blast:
            yield self._emit(
                "status",
                {
                    "content": (
                        "📡 正在向 NCBI 全球公共集群提交比对任务，因官方队列波动，"
                        "通常需要 1~5 分钟，请耐心等待..."
                    ),
                    "state": "running",
                },
            )

        try:
            if inspect.iscoroutinefunction(tool_fn):
                result = await tool_fn(**args)
            else:
                _log_q: queue.Queue = queue.Queue(maxsize=256)

                def _tool_log_sink(content: str, state: str = "running") -> None:
                    try:
                        _log_q.put_nowait({"content": content, "state": state})
                    except queue.Full:
                        pass

                def _run_sync_tool() -> Dict[str, Any]:
                    with tool_log_emitter_scope(_tool_log_sink):
                        out = tool_fn(**args)
                        return out if isinstance(out, dict) else {"status": "success", "data": out}

                _tool_task = asyncio.create_task(asyncio.to_thread(_run_sync_tool))
                _last_hb = time.monotonic()
                _hb_interval = 30.0
                while not _tool_task.done():
                    await asyncio.sleep(0.05)
                    while True:
                        try:
                            _line = _log_q.get_nowait()
                        except queue.Empty:
                            break
                        yield self._emit(
                            "status",
                            {
                                "content": str(_line.get("content") or ""),
                                "state": str(_line.get("state") or "running"),
                            },
                        )
                    if _remote_blast and (time.monotonic() - _last_hb) >= _hb_interval:
                        yield self._emit(
                            "status",
                            {
                                "content": (
                                    "⏳ NCBI 远程排队中...（官方集群负载波动属正常现象，"
                                    "请勿关闭页面）"
                                ),
                                "state": "running",
                            },
                        )
                        _last_hb = time.monotonic()
                result = await _tool_task
                while True:
                    try:
                        _line = _log_q.get_nowait()
                    except queue.Empty:
                        break
                    yield self._emit(
                        "status",
                        {
                            "content": str(_line.get("content") or ""),
                            "state": str(_line.get("state") or "running"),
                        },
                    )
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
            mem = ""
            try:
                mem = build_tool_output_memory_text(tool_name, tr)
            except Exception as mem_exc:  # noqa: BLE001
                logger.warning("工具输出记忆注入构建失败（忽略）: %s", mem_exc, exc_info=True)
            ux = "✅ **技能执行成功，请在右侧工作台查看详细报告。**"
            full_text = f"{ux}\n\n{mem}" if mem else ux
            yield self._emit(
                "message",
                {
                    "content": full_text,
                    **({"tool_output_memory": mem} if mem else {}),
                },
            )
