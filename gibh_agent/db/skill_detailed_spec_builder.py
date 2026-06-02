# -*- coding: utf-8 -*-
"""
从 seed_skills + launch_skill_demos + prompt_template 批量生成 detailed_spec。

排除 main_category=多模态组学；已存在于 SKILL_DETAILED_SPECS_BY_TOOL_ID 的 tool_id 不覆盖。
"""
from __future__ import annotations

import json
import re
from typing import Any, Dict, List, Optional, Tuple

from gibh_agent.db.seed_skills import get_all_system_skills_list
from gibh_agent.skills.launch_skill_demos import LAUNCH_SKILL_DEMO_ARGS

_ROUTE_RE = re.compile(r"\[Skill_Route:\s*([a-zA-Z0-9_]+)\s*\]", re.I)


def extract_tool_id_from_prompt(prompt_template: Optional[str]) -> str:
    if not prompt_template:
        return ""
    m = _ROUTE_RE.search(prompt_template)
    return m.group(1) if m else ""
_JSON_BLOCK_RE = re.compile(r"```json\s*([\s\S]*?)```", re.I)

_PROMPT_SOFT_TOOL_IDS = frozenset(
    {
        "ppt_outline",
        "mindmap_gen",
        "weekly_report_writer",
        "tailored_resume",
        "academic_poster_generator",
        "academic_abstract_refiner",
        "email_manager",
        "pdf_extractor",
        "deep_research",
        "blueprint_drafter",
        "diff_expr_interpreter",
        "go_kegg_narrative",
        "single_cell_checklist",
        "journal_cover_letter",
        "review_rebuttal_outline",
        "acmg_variant_interpretation",
        "pipeline_selection_memo",
        "omics_metadata_review",
        "qpcr_primer_design_guide",
        "oral_presentation_outline",
        "clinical_trial_protocol_skeleton",
        "drug_repositioning_memo",
        "protein_function_hypothesis",
        "experiment_failure_postmortem",
        "literature_matrix_notes",
        "multi_omics_storyline",
        "ethics_consent_checklist",
        "stats_method_advisor",
    }
)

_PARAM_TYPE_HINTS: Dict[str, str] = {
    "query": "string",
    "limit": "integer",
    "user_request": "string",
    "context": "string",
    "smiles": "string",
    "substructure_smiles": "string",
    "target_smiles": "string",
    "sequence_text": "string",
    "sequence_or_path": "string",
    "file_path": "string",
    "input_text": "string",
    "input_format": "string",
    "output_format": "string",
    "record_id": "string",
    "database": "string",
    "use_remote": "boolean",
    "max_target_seqs": "integer",
    "mhc_allele": "string",
    "expid": "string",
    "max_cids": "integer",
    "match_threshold": "number",
    "similarity_threshold": "integer",
    "search_type": "string",
    "num_tokens": "integer",
    "temperature": "number",
    "data_path": "string",
    "input_dir": "string",
    "matrix_dir": "string",
}


def _param_label(key: str) -> str:
    labels = {
        "user_request": "用户请求 / 主题",
        "context": "补充材料",
        "query": "检索词",
        "sequence_text": "序列文本",
        "sequence_or_path": "序列或文件路径",
        "smiles": "SMILES",
        "file_path": "文件路径",
        "input_text": "输入结构文本",
    }
    return labels.get(key, key)


def _param_desc(key: str, tool_id: str) -> str:
    common = {
        "user_request": "任务主题、对比设计或输出结构要求。",
        "context": "待分析表格、摘要、metadata 或前序材料。",
        "query": "数据库检索关键词、基因名、rsID 或 CHEMBL ID。",
        "limit": "返回条数上限。",
        "sequence_text": "FASTA 或裸序列文本。",
        "sequence_or_path": "基因符号或 FASTA 绝对路径。",
        "smiles": "小分子 SMILES 字符串。",
        "file_path": "上传文件的绝对路径（若适用）。",
    }
    if key in common:
        return common[key]
    return f"见技能 `{tool_id}` 实现与广场一键体验 JSON。"


def _infer_kind(tool_id: str, main_cat: str, sub_cat: str) -> str:
    tid = tool_id.lower()
    sc = (sub_cat or "").strip()
    mc = (main_cat or "").strip()
    if tool_id in _PROMPT_SOFT_TOOL_IDS or sc == "文本处理":
        return "prompt_soft"
    if sc == "信息检索" or tid.endswith("_query") or "fetch" in tid or "crossref" in tid:
        return "db_query"
    if "blast" in tid:
        return "blast"
    if mc == "化学" or tid.startswith("chem_") or "rdkit" in tid or "lipinski" in tid:
        return "chem"
    if sc == "数据可视化" or "render" in tid or "image" in tid or "pymol" in tid:
        return "viz"
    if "simulation" in tid or sc == "临床应用":
        return "simulation"
    if sc == "预测与建模":
        return "predict"
    return "analysis"


def _viz_html(kind: str, name: str, tool_id: str) -> str:
    safe_name = (name or tool_id).replace("<", "").replace(">", "")[:40]
    if kind == "db_query":
        return f"""
<table style="width:100%;border-collapse:collapse;font-family:system-ui;font-size:11px;">
  <thead><tr style="background:#f1f5f9;">
  <th style="padding:6px 8px;text-align:left;border:1px solid #e2e8f0;">字段</th>
  <th style="padding:6px 8px;text-align:left;border:1px solid #e2e8f0;">示例值</th></tr></thead>
  <tbody>
  <tr><td style="padding:6px 8px;border-bottom:1px solid #e2e8f0;">命中记录</td>
  <td style="padding:6px 8px;border-bottom:1px solid #e2e8f0;">{safe_name} 检索结果摘要</td></tr>
  </tbody>
</table>"""
    if kind == "blast":
        return """
<table style="width:100%;border-collapse:collapse;font-family:system-ui;font-size:11px;">
  <thead><tr style="background:#1e293b;color:#f8fafc;">
  <th style="padding:6px 8px;">Accession</th><th style="padding:6px 8px;">E-value</th></tr></thead>
  <tbody><tr><td style="padding:6px 8px;border-bottom:1px solid #e2e8f0;">NM_001234</td>
  <td style="padding:6px 8px;border-bottom:1px solid #e2e8f0;">2.3e-45</td></tr></tbody>
</table>"""
    if kind == "chem":
        return f"""
<div style="font-family:system-ui;font-size:12px;padding:12px;background:#f8fafc;border-radius:8px;">
  <div style="font-weight:600;margin-bottom:6px;color:#334155;">{safe_name}</div>
  <div style="padding:8px;background:#fff;border:1px solid #e2e8f0;border-radius:6px;font-family:ui-monospace;font-size:11px;">
  SMILES → 性质 / 结构 / 检索结果
  </div>
</div>"""
    if kind == "viz":
        return f"""
<div style="font-family:system-ui;font-size:12px;padding:14px;background:#eff6ff;border-radius:8px;text-align:center;color:#1e40af;">
  <div style="font-weight:600;">{safe_name}</div>
  <div style="margin-top:8px;font-size:11px;color:#64748b;">结构图 / 3D 渲染 / 图表预览（示意）</div>
</div>"""
    if kind == "simulation":
        return f"""
<div style="font-family:system-ui;font-size:12px;padding:12px;background:#fff;border:1px solid #e2e8f0;border-radius:8px;">
  <div style="font-weight:600;margin-bottom:8px;">{safe_name} · 流程示意</div>
  <div style="display:flex;gap:6px;flex-wrap:wrap;">
  <span style="padding:4px 10px;background:#dbeafe;border-radius:4px;font-size:11px;">输入校验</span>
  <span style="color:#94a3b8;">→</span>
  <span style="padding:4px 10px;background:#d1fae5;border-radius:4px;font-size:11px;">核心计算</span>
  <span style="color:#94a3b8;">→</span>
  <span style="padding:4px 10px;background:#fef3c7;border-radius:4px;font-size:11px;">结果报告</span>
  </div>
</div>"""
    if kind == "prompt_soft":
        return f"""
<div style="font-family:system-ui;font-size:12px;padding:12px;background:#fff;border:1px solid #e2e8f0;border-radius:8px;">
  <div style="font-weight:600;margin-bottom:6px;">{safe_name}</div>
  <div style="color:#475569;line-height:1.55;">Markdown 结构化交付 · 可复制调用示例一键体验</div>
</div>"""
    # analysis / predict default
    return f"""
<div style="font-family:system-ui;font-size:12px;padding:12px;background:#f8fafc;border:1px solid #e2e8f0;border-radius:8px;">
  <div style="font-weight:600;color:#334155;">{safe_name}</div>
  <div style="margin-top:6px;font-size:11px;color:#64748b;">分析结果表 / 指标卡 / Markdown 报告（示意）</div>
</div>"""


def _usage_hint(kind: str, tool_id: str, main_cat: str) -> str:
    if kind == "prompt_soft":
        return "本技能为 Prompt 驱动，无独立表单。请在对话框粘贴材料，或使用下方「对话调用示例」一键复制后发送。"
    if kind == "db_query":
        return "需 launch-skills 联网访问外部数据库；未改写检索词时使用广场一键体验 JSON。"
    if kind == "blast":
        return "远程 BLAST 通常需 1–3 分钟；请粘贴序列或提供 FASTA 路径。"
    if kind == "chem":
        return "化学技能由 launch-skills / RDKit 执行；SMILES 与格式字段须符合规范。"
    if kind == "viz":
        return "输出含可视化文件或 HTML；大图请在对话工作台查看完整结果。"
    if kind == "simulation":
        return "仿真类技能根据您提供的参数生成流程日志与结果文件说明；不替代真实湿实验。"
    if main_cat == "化学":
        return "请确认已启动 launch-skills 服务；部分任务需上传数据文件。"
    return "在首页对话框发送调用示例，或由智能体根据 [Skill_Route] 自动路由执行。"


def _outputs_for_kind(kind: str) -> List[Dict[str, Any]]:
    if kind == "prompt_soft":
        return [{"name": "markdown", "type": "Markdown", "description": "结构化报告 / 大纲 / 清单。"}]
    if kind == "viz":
        return [
            {"name": "artifact_url", "type": "URL / 文件", "description": "图片、HTML 或 Mol 块下载链接。"},
            {"name": "markdown", "type": "Markdown", "description": "结果说明与参数摘要。"},
        ]
    if kind == "db_query":
        return [{"name": "markdown", "type": "Markdown + 表格", "description": "命中记录与字段摘要。"}]
    if kind == "blast":
        return [{"name": "markdown", "type": "Markdown + 表格", "description": "BLAST 命中与比对统计。"}]
    return [
        {"name": "markdown", "type": "Markdown", "description": "分析报告、表格或日志。"},
        {"name": "json_url", "type": "JSON（可选）", "description": "完整结果下载链接（若条目较多）。"},
    ]


def _inputs_for_kind(kind: str, demo: Dict[str, Any]) -> List[Dict[str, Any]]:
    if kind == "prompt_soft":
        return [
            {"name": "任务描述", "type": "自然语言", "required": True, "description": "写入 user_request。"},
            {"name": "补充材料", "type": "文本 / 表格", "required": False, "description": "写入 context。"},
        ]
    if kind in ("db_query", "blast"):
        return [{"name": "检索条件", "type": "string / 序列", "required": True, "description": "见参数表 query 或 sequence 字段。"}]
    if kind == "chem":
        return [{"name": "分子表示", "type": "SMILES / 文件", "required": True, "description": "SMILES 文本或结构文件路径。"}]
    if demo.get("file_path") is not None or demo.get("data_path") is not None:
        return [{"name": "数据文件", "type": "文件路径", "required": True, "description": "上传后填入 file_path / data_path。"}]
    return [{"name": "分析输入", "type": "文本 / 文件", "required": True, "description": "按参数表提供序列、表格或配置。"}]


def _build_params_table(tool_id: str, demo: Dict[str, Any]) -> List[Dict[str, Any]]:
    if not demo:
        if tool_id in _PROMPT_SOFT_TOOL_IDS:
            return [
                {"name": "user_request", "type": "string", "required": True, "description": _param_desc("user_request", tool_id)},
                {"name": "context", "type": "string", "required": False, "description": _param_desc("context", tool_id)},
            ]
        return [
            {"name": "user_request", "type": "string", "required": True, "description": "任务说明（若技能支持）。"},
        ]
    rows = []
    for key in demo.keys():
        rows.append(
            {
                "name": key,
                "type": _PARAM_TYPE_HINTS.get(key, "string"),
                "required": key in ("query", "user_request", "smiles", "sequence_text", "substructure_smiles", "input_text"),
                "description": _param_desc(key, tool_id),
            }
        )
    return rows


def _extract_query_example(name: str, prompt: str, tool_id: str, demo: Dict[str, Any]) -> str:
    m = _JSON_BLOCK_RE.search(prompt or "")
    if m:
        block = m.group(1).strip()
        intro = f"您好。请使用技能「{name}」完成下列任务。\n\n```json\n{block}\n```"
        return intro
    if demo:
        body = json.dumps(demo, ensure_ascii=False, indent=0)
        return f"您好。请执行「{name}」。\n\n```json\n{body}\n```"
    return f"您好。请调用技能 [Skill_Route: {tool_id}] 完成相关分析。"


def build_detailed_spec_from_seed(skill: Dict[str, Any]) -> Optional[Dict[str, Any]]:
    """由单条 seed 技能字典生成 detailed_spec。"""
    pt = skill.get("prompt_template") or ""
    tool_id = extract_tool_id_from_prompt(pt)
    if not tool_id:
        return None
    name = (skill.get("name") or tool_id).strip()
    main_cat = (skill.get("main_category") or "").strip()
    sub_cat = (skill.get("sub_category") or "").strip()
    desc = (skill.get("description") or "").strip()
    kind = _infer_kind(tool_id, main_cat, sub_cat)
    demo = dict(LAUNCH_SKILL_DEMO_ARGS.get(tool_id) or {})

    description_long = desc
    if description_long and not description_long.endswith("。"):
        description_long += "。"
    if len(description_long) < 20:
        description_long = f"「{name}」：{desc or '科研分析技能'}。通过对话调用 [Skill_Route: {tool_id}] 执行。"

    return {
        "tool_id": tool_id,
        "description_long": description_long,
        "usage_hint": _usage_hint(kind, tool_id, main_cat),
        "inputs": _inputs_for_kind(kind, demo),
        "outputs": _outputs_for_kind(kind),
        "parameters_table": _build_params_table(tool_id, demo),
        "query_examples": [_extract_query_example(name, pt, tool_id, demo)],
        "demo_visualization": _viz_html(kind, name, tool_id),
    }


def build_bulk_detailed_specs(
    *,
    exclude_main_categories: Optional[Tuple[str, ...]] = ("多模态组学",),
    skip_existing_keys: Optional[set] = None,
) -> Dict[str, Dict[str, Any]]:
    """生成尚未登记的技能 detailed_spec 字典（按 tool_id 去重，保留首条）。"""
    existing = set(skip_existing_keys or [])
    out: Dict[str, Dict[str, Any]] = {}
    for skill in get_all_system_skills_list():
        mc = (skill.get("main_category") or "").strip()
        if exclude_main_categories and mc in exclude_main_categories:
            continue
        spec = build_detailed_spec_from_seed(skill)
        if not spec:
            continue
        tid = spec["tool_id"]
        if tid in existing:
            continue
        if tid in out:
            continue
        out[tid] = spec
    return out
