# -*- coding: utf-8 -*-
"""
技能广场「详情说明书」首批数据（Phase 1：10 项核心技能）。

主键为 tool_id（与 `[Skill_Route: <tool_id>]` 一致）。
后续批次在此文件追加，或通过 Skill.detailed_spec 列覆盖。
"""
from __future__ import annotations

import re
from typing import Any, Dict, Optional

from gibh_agent.db.skill_detail_demo_visualizations import (
    DEMO_VIZ_ABSTRACT,
    DEMO_VIZ_BLUEPRINT_RNASEQ,
    DEMO_VIZ_CHEMBL,
    DEMO_VIZ_DEG_TABLE,
    DEMO_VIZ_DEEP_RESEARCH,
    DEMO_VIZ_GO_KEGG,
    DEMO_VIZ_MINDMAP,
    DEMO_VIZ_PPT_OUTLINE,
    DEMO_VIZ_SC_CHECKLIST,
    DEMO_VIZ_WEEKLY_REPORT,
)

_SKILL_ROUTE_RE = re.compile(r"\[Skill_Route:\s*([a-zA-Z0-9_]+)\s*\]", re.I)

# 标准 detailed_spec 结构见 docs/技能详情 Demo 页面规范.md
SKILL_DETAILED_SPECS_BY_TOOL_ID: Dict[str, Dict[str, Any]] = {
    "ppt_outline": {
        "tool_id": "ppt_outline",
        "description_long": (
            "将研究主题、方法要点与核心结论整理为**分页科研汇报 PPT 大纲**。"
            "输出为结构化 JSON（总页数、每页标题与 3–6 条要点），可直接导入 PowerPoint / Keynote，"
            "或作为组会汇报的口述提纲。"
        ),
        "usage_hint": "本技能无独立表单界面，请在首页对话框粘贴研究材料或通过下方示例一键复制后发送。",
        "inputs": [
            {
                "name": "研究主题与材料",
                "type": "自然语言 / 可选前序分析摘要",
                "required": True,
                "description": "汇报主题、科学问题、主要发现；若有单细胞/组学分析摘要可一并提供。",
            }
        ],
        "outputs": [
            {
                "name": "ppt_outline",
                "type": "JSON",
                "description": "含 theme（总标题）、slides[]（page/title/bullets）、total_pages。",
            },
            {
                "name": "markdown",
                "type": "Markdown",
                "description": "工作台可读的汇报大纲摘要（由编排器渲染）。",
            },
        ],
        "parameters_table": [
            {
                "name": "user_request",
                "type": "string",
                "required": True,
                "description": "汇报主题、受众（组会/答辩）与希望覆盖的章节重点。",
            },
            {
                "name": "context",
                "type": "string",
                "required": False,
                "description": "已有分析结论、图表说明或论文摘要片段。",
            },
        ],
        "query_examples": [
            (
                "您好。请根据下列内容生成**科研汇报 PPT 大纲**（结构化 JSON）。\n\n"
                "```json\n"
                '{"user_request": "单细胞时空动力学分析：研究背景、最优传输轨迹推断、'
                '驱动基因与通路富集、核心结论与展望", "context": ""}\n'
                "```"
            ),
        ],
        "demo_visualization": DEMO_VIZ_PPT_OUTLINE,
    },
    "mindmap_gen": {
        "tool_id": "mindmap_gen",
        "description_long": (
            "从分析逻辑、实验流程或决策树中提取层级关系，生成 **Mermaid 思维导图**。"
            "适合梳理组学分析流水线、机制假说或项目里程碑，并在工作台可视化展示。"
        ),
        "usage_hint": "在对话中描述「从 A 到 B 的分析步骤」或粘贴已有流程说明；技能将输出 Mermaid 源码。",
        "inputs": [
            {
                "name": "分析逻辑描述",
                "type": "自然语言",
                "required": True,
                "description": "步骤链、分支条件或模块依赖（如 QC → 聚类 → 注释 → 富集）。",
            }
        ],
        "outputs": [
            {
                "name": "mermaid",
                "type": "Mermaid 文本",
                "description": "可渲染的思维导图源码（graph TD / flowchart）。",
            },
            {
                "name": "markdown",
                "type": "Markdown",
                "description": "附带简要说明的交付块。",
            },
        ],
        "parameters_table": [
            {
                "name": "user_request",
                "type": "string",
                "required": True,
                "description": "要绘制的逻辑主题与层级要点（箭头连接步骤）。",
            },
            {
                "name": "context",
                "type": "string",
                "required": False,
                "description": "补充材料、已有流程图文字版或项目背景。",
            },
        ],
        "query_examples": [
            (
                "您好。请将下列分析逻辑整理为 **Mermaid 思维导图**。\n\n"
                "```json\n"
                '{"user_request": "单细胞时空动力学：数据校验→时序标准化→轨迹推断→'
                '驱动基因→通路富集→专家解读", "context": ""}\n'
                "```"
            ),
        ],
        "demo_visualization": DEMO_VIZ_MINDMAP,
    },
    "diff_expr_interpreter": {
        "tool_id": "diff_expr_interpreter",
        "description_long": (
            "将用户粘贴的**差异表达（DEG）结果表**解读为组会/论文可用的生物学叙事："
            "关键基因功能、通路关联、与实验设计的对应关系，以及后续 qPCR / Western / 功能实验验证建议。"
            "不重新计算统计量，仅基于您提供的数据解读。"
        ),
        "usage_hint": "请粘贴 gene、log2FC、padj（或 pvalue）列；说明对比组与筛选阈值。",
        "inputs": [
            {
                "name": "DEG 结果表",
                "type": "TSV / 文本表",
                "required": True,
                "description": "至少包含基因名与效应量、校正 p 值；可附火山图文字摘要。",
            }
        ],
        "outputs": [
            {
                "name": "markdown",
                "type": "Markdown 报告",
                "description": "分节叙事 + 重点基因表 + 验证建议清单。",
            },
        ],
        "parameters_table": [
            {
                "name": "user_request",
                "type": "string",
                "required": True,
                "description": "对比设计（如处理 vs 对照）、阈值与解读侧重点。",
            },
            {
                "name": "context",
                "type": "string",
                "required": True,
                "description": "完整 DEG 表或 Top 基因列表（制表符/逗号分隔）。",
            },
        ],
        "query_examples": [
            (
                "您好。请将我提供的**差异表达（DEG）结果**解读为组会可用的生物学叙事报告。\n\n"
                "```json\n"
                '{"user_request": "解读下列 Demo DEG 结果并给出验证建议", "context": '
                '"对比：处理组 vs 对照组；阈值 padj<0.05, |log2FC|>1\\n'
                'TP53\\t2.1\\t1.2e-8\\nCDKN1A\\t1.8\\t3.4e-6\\nMYC\\t-1.5\\t0.002"}\n'
                "```"
            ),
        ],
        "demo_visualization": DEMO_VIZ_DEG_TABLE,
    },
    "weekly_report_writer": {
        "tool_id": "weekly_report_writer",
        "description_long": (
            "根据周期范围、工作笔记与上一份周报，合成**双视角周报**："
            "个人状态同步（进展/阻塞/学习）与团队对外同步（里程碑/风险/协作需求），"
            "并继承上一份中的高优先级待办（🟥/🟨）。"
        ),
        "usage_hint": "在 user_request 中写明日期范围；context 中放日志路径或本周要点摘录。",
        "inputs": [
            {
                "name": "周期与笔记",
                "type": "自然语言 / 路径说明",
                "required": True,
                "description": "起止日期、项目名、本周完成与未完成事项。",
            }
        ],
        "outputs": [
            {
                "name": "markdown",
                "type": "Markdown",
                "description": "完整周报正文（个人 + 团队章节）。",
            },
        ],
        "parameters_table": [
            {
                "name": "user_request",
                "type": "string",
                "required": True,
                "description": "周期、关注项目与输出结构要求。",
            },
            {
                "name": "context",
                "type": "string",
                "required": False,
                "description": "上一份周报路径、Obsidian/日志摘录或待办继承说明。",
            },
        ],
        "query_examples": [
            (
                "请为我起草本周**周报**（个人 + 团队双视角）。\n\n"
                "```json\n"
                '{"user_request": "周期：2026-05-19 至 2026-05-25；按 SKILL 工作流合成周报", '
                '"context": "上一份报告路径与本周实验/分析进展摘要"}\n'
                "```"
            ),
        ],
        "demo_visualization": DEMO_VIZ_WEEKLY_REPORT,
    },
    "deep_research": {
        "tool_id": "deep_research",
        "description_long": (
            "对指定科研或临床主题进行**多视角深度调研**：执行摘要、分主题发现、"
            "共识与争议、证据强度提示与引用占位。适合开题前的文献脉络梳理；"
            "若未提供原文，须在结论中标注「需数据库核实」。"
        ),
        "usage_hint": "主题尽量具体；可在 context 粘贴已有摘要或笔记，禁止编造未给出的文献。",
        "inputs": [
            {
                "name": "研究问题",
                "type": "自然语言",
                "required": True,
                "description": "如「间歇性禁食对代谢健康的影响」。",
            }
        ],
        "outputs": [
            {
                "name": "markdown",
                "type": "Markdown 长报告",
                "description": "Executive Summary + 分节论述 + 参考文献占位。",
            },
        ],
        "parameters_table": [
            {
                "name": "user_request",
                "type": "string",
                "required": True,
                "description": "调研主题与希望覆盖的子问题。",
            },
            {
                "name": "context",
                "type": "string",
                "required": False,
                "description": "用户已有论文摘要、课程笔记或约束（如仅近 5 年）。",
            },
        ],
        "query_examples": [
            (
                "您好。请对下列主题做**深度调研**（摘要 + 发现 + 共识/争议）。\n\n"
                "```json\n"
                '{"user_request": "间歇性禁食对代谢健康的影响：益处、风险与适用人群", "context": ""}\n'
                "```"
            ),
        ],
        "demo_visualization": DEMO_VIZ_DEEP_RESEARCH,
    },
    "blueprint_drafter": {
        "tool_id": "blueprint_drafter",
        "description_long": (
            "根据业务流程描述生成**扁平工程蓝图 HTML**（完整 <!DOCTYPE html>，"
            "系统字体、无外部 CDN），用于 RNA-seq / 蛋白组等流水线示意或系统架构草图。"
            "节点标签应为用户业务术语，而非平台内部代号。"
        ),
        "usage_hint": "用箭头描述步骤顺序；生成后在工作台 iframe 预览，可导出或截图。",
        "inputs": [
            {
                "name": "流程描述",
                "type": "自然语言",
                "required": True,
                "description": "从原始数据到交付物的步骤链（含质控、分析、可视化）。",
            }
        ],
        "outputs": [
            {
                "name": "html",
                "type": "HTML 文档",
                "description": "可 iframe 渲染的单文件蓝图。",
            },
        ],
        "parameters_table": [
            {
                "name": "user_request",
                "type": "string",
                "required": True,
                "description": "要绘制的流程或架构（中文节点名）。",
            },
            {
                "name": "context",
                "type": "string",
                "required": False,
                "description": "数据来源、物种、额外分支说明。",
            },
        ],
        "query_examples": [
            (
                "您好。请生成**扁平工程蓝图 HTML**。\n\n"
                "```json\n"
                '{"user_request": "绘制 RNA-seq 差异分析流程蓝图：原始 FASTQ → 质控 → '
                '表达矩阵 → 差异分析 → 可视化", "context": "中文节点标签；仅描述业务数据流"}\n'
                "```"
            ),
        ],
        "demo_visualization": DEMO_VIZ_BLUEPRINT_RNASEQ,
    },
    "go_kegg_narrative": {
        "tool_id": "go_kegg_narrative",
        "description_long": (
            "基于用户提供的 **GO/KEGG 富集结果表**（Term、padj、GeneRatio 等），"
            "撰写「机制—疾病关联—验证方向」连贯叙事，突出主题通路。"
            "不编造表中未出现的 p 值或基因数。"
        ),
        "usage_hint": "粘贴富集结果前几行即可体验；说明希望突出的生物学主题（如免疫、代谢）。",
        "inputs": [
            {
                "name": "富集结果表",
                "type": "TSV / 文本",
                "required": True,
                "description": "GO ID 或 KEGG 通路名、校正 p 值、GeneRatio 等。",
            }
        ],
        "outputs": [
            {
                "name": "markdown",
                "type": "Markdown",
                "description": "分节机制叙事 + 后续实验/文献检索建议。",
            },
        ],
        "parameters_table": [
            {
                "name": "user_request",
                "type": "string",
                "required": True,
                "description": "叙事角度、疾病背景或希望强调的通路主题。",
            },
            {
                "name": "context",
                "type": "string",
                "required": True,
                "description": "富集表内容（制表符分隔）。",
            },
        ],
        "query_examples": [
            (
                "您好。请根据下列富集结果撰写**GO/KEGG 叙事报告**。\n\n"
                "```json\n"
                '{"user_request": "撰写富集结果叙事，突出免疫与代谢主题", "context": '
                '"GO:0006955 immune response\\tpadj=2.1e-5\\tGeneRatio=12/180"}\n'
                "```"
            ),
        ],
        "demo_visualization": DEMO_VIZ_GO_KEGG,
    },
    "single_cell_checklist": {
        "tool_id": "single_cell_checklist",
        "description_long": (
            "根据研究问题与实验设计文字描述，输出 **scRNA-seq 实验设计检查清单**："
            "样本量、批次、双胞、参考面板、质控阈值、对照设置等可勾选条目，"
            "并标注常见混杂风险。"
        ),
        "usage_hint": "说明物种、平台（如 10x 3' v3）、分组与生物学重复数。",
        "inputs": [
            {
                "name": "实验设计描述",
                "type": "自然语言",
                "required": True,
                "description": "疾病模型、处理因素、取样时间点与关注细胞类型。",
            }
        ],
        "outputs": [
            {
                "name": "markdown",
                "type": "Markdown Checklist",
                "description": "分级检查项与风险备注。",
            },
        ],
        "parameters_table": [
            {
                "name": "user_request",
                "type": "string",
                "required": True,
                "description": "课题背景、平台与分组设计。",
            },
            {
                "name": "context",
                "type": "string",
                "required": False,
                "description": "预期细胞数、已有预实验或文献参考。",
            },
        ],
        "query_examples": [
            (
                "您好。请生成 **scRNA-seq 实验设计检查清单**。\n\n"
                "```json\n"
                '{"user_request": "人肺腺癌免疫治疗前后 scRNA-seq；10x 3\' v3", '
                '"context": "每组 6 例重复；关注 T 细胞耗竭"}\n'
                "```"
            ),
        ],
        "demo_visualization": DEMO_VIZ_SC_CHECKLIST,
    },
    "academic_abstract_refiner": {
        "tool_id": "academic_abstract_refiner",
        "description_long": (
            "将用户提供的长摘要或结论段落精炼为**中英双语单段 SCI 风格摘要**"
            "（Summary_Report 格式，无小标题）。保持数据与结论与原文一致，不虚构实验结果。"
        ),
        "usage_hint": "将待精炼英文/中文草稿放入 context；说明目标期刊风格（可选）。",
        "inputs": [
            {
                "name": "待精炼文稿",
                "type": "纯文本",
                "required": True,
                "description": "摘要、Introduction 末段或 Discussion 结论段。",
            }
        ],
        "outputs": [
            {
                "name": "markdown",
                "type": "Markdown",
                "description": "中英双语摘要 + 可选修改说明。",
            },
        ],
        "parameters_table": [
            {
                "name": "user_request",
                "type": "string",
                "required": True,
                "description": "精炼要求（字数、时态、是否保留方法细节）。",
            },
            {
                "name": "context",
                "type": "string",
                "required": True,
                "description": "完整待精炼原文。",
            },
        ],
        "query_examples": [
            (
                "您好。请将下列文稿精炼为**中英双语 SCI 风格摘要**。\n\n"
                "```json\n"
                '{"user_request": "精炼为双语单段摘要，无小标题", "context": "（在此粘贴待精炼长文）"}\n'
                "```"
            ),
        ],
        "demo_visualization": DEMO_VIZ_ABSTRACT,
    },
    "chembl_drug_search": {
        "tool_id": "chembl_drug_search",
        "description_long": (
            "在 **ChEMBL** 数据库中检索药物或分子记录：支持药物名、同义词或 CHEMBL ID。"
            "返回结构化表格与 Markdown 摘要，便于先导化合物调研与靶点关联浏览。"
            "需 launch-skills 微服务联网执行。"
        ),
        "usage_hint": "未改写检索词时使用 Demo JSON；有明确药物名时替换 query 字段。",
        "inputs": [
            {
                "name": "检索词",
                "type": "string",
                "required": True,
                "description": "药物通用名、别名或 CHEMBL 分子 ID。",
            }
        ],
        "outputs": [
            {
                "name": "markdown",
                "type": "Markdown + 表格",
                "description": "命中记录摘要与关键字段。",
            },
            {
                "name": "json_url",
                "type": "JSON（可选）",
                "description": "完整结果下载链接（若条目过多）。",
            },
        ],
        "parameters_table": [
            {
                "name": "query",
                "type": "string",
                "required": True,
                "description": "检索关键词，如 aspirin。",
            },
            {
                "name": "search_type",
                "type": "string",
                "required": False,
                "description": "drug（默认）或 molecule（按 CHEMBL ID）。",
            },
            {
                "name": "limit",
                "type": "integer",
                "required": False,
                "description": "返回条数上限，默认 10。",
            },
        ],
        "query_examples": [
            (
                "您好。我需要在 ChEMBL 中检索药物记录。\n\n"
                "```json\n"
                '{"query": "aspirin", "search_type": "drug", "limit": 10}\n'
                "```\n\n"
                "请按上述 JSON 检索 ChEMBL 中包含 aspirin 的药物，最多 10 条。"
            ),
        ],
        "demo_visualization": DEMO_VIZ_CHEMBL,
    },
}

# Phase 2：合并第二批 10 项（PubMed / BLAST / BepiPred 等）
from gibh_agent.db.skill_detailed_specs_phase2 import SKILL_DETAILED_SPECS_PHASE2  # noqa: E402

SKILL_DETAILED_SPECS_BY_TOOL_ID.update(SKILL_DETAILED_SPECS_PHASE2)

# 批量：除多模态组学外其余已实现 [Skill_Route] 技能（约 70+ 项，不覆盖 Phase1/2 手工精修条目）
from gibh_agent.db.skill_detailed_specs_bulk import SKILL_DETAILED_SPECS_BULK  # noqa: E402

SKILL_DETAILED_SPECS_BY_TOOL_ID.update(SKILL_DETAILED_SPECS_BULK)


def extract_tool_id_from_prompt(prompt_template: Optional[str]) -> str:
    """从 prompt_template 解析 Skill_Route tool_id。"""
    if not prompt_template:
        return ""
    m = _SKILL_ROUTE_RE.search(prompt_template)
    return m.group(1) if m else ""


def resolve_detailed_spec(
    prompt_template: Optional[str],
    db_detailed_spec: Any = None,
    tool_id_hint: Optional[str] = None,
) -> Optional[Dict[str, Any]]:
    """
    解析技能的 detailed_spec：DB 列优先，否则静态注册表。
    tool_id_hint 用于前端 Featured 卡片（id 即 tool_id）。
    """
    if isinstance(db_detailed_spec, dict) and db_detailed_spec:
        out = dict(db_detailed_spec)
        if not out.get("tool_id"):
            tid = tool_id_hint or extract_tool_id_from_prompt(prompt_template)
            if tid:
                out["tool_id"] = tid
        return out
    tid = (tool_id_hint or "").strip() or extract_tool_id_from_prompt(prompt_template)
    if tid and tid in SKILL_DETAILED_SPECS_BY_TOOL_ID:
        return dict(SKILL_DETAILED_SPECS_BY_TOOL_ID[tid])
    return None


def enrich_skill_plaza_item(item: Dict[str, Any], db_detailed_spec: Any = None) -> Dict[str, Any]:
    """为广场列表项附加 tool_name 与 detailed_spec（若有）。"""
    pt = item.get("prompt_template") or ""
    hint = str(item.get("tool_name") or item.get("id") or "").strip()
    if hint and not _SKILL_ROUTE_RE.search(str(hint)):
        tool_id = hint if hint in SKILL_DETAILED_SPECS_BY_TOOL_ID else extract_tool_id_from_prompt(pt)
    else:
        tool_id = extract_tool_id_from_prompt(pt) or (
            hint if hint in SKILL_DETAILED_SPECS_BY_TOOL_ID else ""
        )
    if tool_id:
        item["tool_name"] = tool_id
    spec = resolve_detailed_spec(pt, db_detailed_spec=db_detailed_spec, tool_id_hint=tool_id or hint)
    if spec:
        item["detailed_spec"] = spec
    return item
