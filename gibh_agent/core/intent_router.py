# -*- coding: utf-8 -*-
"""
Master LLM Intent Router — LUI 自然语言入口的全局路由网关。

在 /api/chat 进入 AgentOrchestrator 之前，用轻量 LLM 调用将用户 message
映射到 GLOBAL_INTENT_REGISTRY 中的 workflow / skill，或触发澄清 / 闲聊分支。
"""
from __future__ import annotations

import json
import logging
import os
import re
from enum import Enum
from typing import Any, Dict, List, Optional

from pydantic import BaseModel, Field

from gibh_agent.core.llm_client import LLMClient
from gibh_agent.core.llm_json_extract import extract_json_object_from_llm_text

logger = logging.getLogger(__name__)


class IntentRouteKind(str, Enum):
    workflow = "workflow"
    skill = "skill"


class IntentRegistryEntry(BaseModel):
    """全局 Tool/Skill 注册表条目（注入 LLM Prompt）。"""

    id: str
    name: str
    aliases: List[str] = Field(default_factory=list)
    description: str = ""
    route_kind: IntentRouteKind = IntentRouteKind.workflow
    target_domain: Optional[str] = None
    required_parameters: List[str] = Field(default_factory=list)


class IntentRouterResult(BaseModel):
    """LLM 意图探针输出（exact | missing_param | clarify | chitchat）。"""

    status: str  # exact | missing_param | clarify | chitchat
    target_id: Optional[str] = None
    question: Optional[str] = None
    candidates: List[Dict[str, str]] = Field(default_factory=list)
    rationale_short: str = ""


# 技能必填参数（注入 LLM 注册表，供 MISSING_PARAM 槽位校验）
_SKILL_REQUIRED_PARAMETERS: Dict[str, List[str]] = {
    "lipinski_druglikeness": [
        "目标化合物：SMILES 结构式、可唯一映射的常用药物/小分子名称，或已挂载的 SMILES 文本文件",
    ],
    "drug_similarity_search": [
        "目标化合物：SMILES 结构式、可唯一映射的化合物俗名/药物名，或已挂载的分子/SMILES 文件",
    ],
    "bepipred3_epitope_prediction": [
        "蛋白/抗原氨基酸序列（FASTA 文本或 sequence_or_path 文件）",
    ],
    # 首发核心技能（gibh_agent/skills/ BaseSkill，与 [Skill_Route] 一致）
    "chembl_drug_search": [
        "药物名称关键词或 ChEMBL ID（query）",
    ],
    "rdkit_substructure_search": [
        "子结构 SMILES/SMARTS（substructure_smiles）与待检分子 SMILES 或 SMILES 列表文件",
    ],
    "nucleotide_sequence_blast": [
        "核酸查询序列（sequence_text）或 FASTA 文件（sequence_or_path）",
    ],
    "protein_sequence_blast": [
        "蛋白氨基酸序列（sequence_text）或 FASTA 文件（sequence_or_path）",
    ],
    "pubchem_smiles_to_cid": [
        "化合物 SMILES 字符串（smiles）",
    ],
    "mhc_epitope_search": [
        "MHC 等位基因（mhc_allele）和/或肽段序列（peptide_sequence）",
    ],
    "chipatlas_experiment_search": [
        "实验 accession（expid）或检索关键词（query）",
    ],
    "rdkit_3d_mol_render": [
        "分子 SMILES（smiles）",
    ],
    "chembl_similar_molecules": [
        "查询 SMILES 或 ChEMBL 分子 ID（smiles / chembl_id）",
    ],
    "rdkit_mol_format_convert": [
        "待转换结构文本（input_text）及 input_format / output_format",
    ],
    # 阵列一 Batch1 · 纯 Prompt 生信软技能（文本粘贴即可，无需大文件）
    "diff_expr_interpreter": [
        "DEG 结果表或火山图/Top 基因文字摘要（粘贴至 context 即可）",
    ],
    "go_kegg_narrative": [
        "GO/KEGG 富集结果表或通路条目（粘贴至 context）",
    ],
    "single_cell_checklist": [
        "研究问题、物种、平台与分组设计（user_request / context 文本）",
    ],
    "journal_cover_letter": [
        "论文题目、核心亮点与目标期刊名称",
    ],
    "review_rebuttal_outline": [
        "审稿意见原文或要点列表（粘贴至 context）",
    ],
    "acmg_variant_interpretation": [
        "变异位点、基因、HGVS、人群频率与功能预测（粘贴至 context）",
    ],
    "pipeline_selection_memo": [
        "分析目标、组学类型、样本量与计算资源约束（文本）",
    ],
    "omics_metadata_review": [
        "样本 metadata 表或字段列表（粘贴至 context）",
    ],
    "qpcr_primer_design_guide": [
        "扩增目标、样本类型与可选引物序列（文本）",
    ],
    "oral_presentation_outline": [
        "报告主题、时长、听众；论文结构或摘要（context）",
    ],
    "clinical_trial_protocol_skeleton": [
        "试验设计类型、人群、主要/次要终点描述",
    ],
    "drug_repositioning_memo": [
        "疾病—靶点—药物与已知证据（文字）",
    ],
    "protein_function_hypothesis": [
        "蛋白结构域、修饰、定位等注释（粘贴至 context）",
    ],
    "experiment_failure_postmortem": [
        "实验目标、失败现象与已尝试排查（文本）",
    ],
    "literature_matrix_notes": [
        "多篇文献摘要或笔记（粘贴至 context）",
    ],
    "multi_omics_storyline": [
        "各组学层级关键发现文字摘要",
    ],
    "ethics_consent_checklist": [
        "研究类型、干预风险与人群特征",
    ],
    "stats_method_advisor": [
        "研究设计、分组、结局变量与样本量结构",
    ],
    "pubmed_query": [
        "PubMed 检索关键词或主题词（query）",
    ],
    "uniprot_query": [
        "蛋白名、基因名或 UniProt accession（query）",
    ],
    "smiles_to_cid": [
        "化合物 SMILES 字符串（smiles）",
    ],
    "seq_format_converter": [
        "核酸/蛋白序列文本（sequence_text）或 FASTA 文件（sequence_or_path）",
    ],
    "calc_molecular_weight": [
        "小分子 SMILES 结构式（smiles）",
    ],
    "clinvar_query": [
        "基因名、rsID 或 HGVS 变异描述（query）",
    ],
    "dbsnp_query": [
        "rsID、基因名或染色体位置（query）",
    ],
    "gwas_catalog_query": [
        "rsID 或疾病/性状名称（query）",
    ],
    "reactome_query": [
        "通路名或关键词（query，如 apoptosis）",
    ],
    "interpro_query": [
        "结构域名、IPR 编号或关键词（query）",
    ],
    "geo_query": [
        "关键词、GSE accession 或 MeSH 主题（query）",
    ],
    "gene_protein_info_query": [
        "人类基因符号（query，如 TP53）",
    ],
    "mrna_sequence_fetch": [
        "人类基因符号（sequence_or_path 或 query）",
    ],
    "drug_id_crossref": [
        "药物名称（query）或 SMILES（smiles）",
    ],
    "refmet_metabolite_search": [
        "代谢物名称（query）",
    ],
    "rnacentral_query": [
        "RNAcentral URS 编号或 RNA 关键词（query）",
    ],
    "ensembl_go_descendants": [
        "GO 术语 ID（query，如 GO:0006915）",
    ],
    "opentargets_chembl_hierarchy": [
        "ChEMBL 药物/分子 ID（query，如 CHEMBL25）",
    ],
    "fda_drug_label_search": [
        "药品品牌名或通用名（query）",
    ],
}

# 工作流必填参数（按注册表 id）
_WORKFLOW_REQUIRED_PARAMETERS: Dict[str, List[str]] = {
    "genomics_pipeline": [
        "测序原始数据或比对结果（如 FASTQ、BAM、VCF 等）的文件路径或已上传挂载",
    ],
    "rna_scrna": [
        "表达矩阵或单细胞 count 矩阵目录（matrix_dir）或已上传表达数据文件",
    ],
    "metabolomics": [
        "代谢组原始数据或特征表文件/目录（input_dir 或 file_path）",
    ],
    "spatial_transcriptomics": [
        "空间转录组表达矩阵与坐标信息（matrix_dir 或配套挂载文件）",
    ],
    "radiomics": [
        "医学影像文件对（image_path 与 mask_path，或 DICOM/NIfTI 挂载）",
    ],
    "proteomics": [
        "蛋白组质谱原始数据或鉴定结果表（file_path / input_dir）",
    ],
    "epigenomics": [
        "表观组测序数据（如 BAM、bed、bigWig 等）文件或目录",
    ],
    "sted_ec_trajectory": [
        "单细胞表达矩阵与细胞注释（matrix_dir 或等效挂载）",
    ],
    "spatiotemporal_dynamics": [
        "时空动力学输入矩阵或 h5ad 等挂载文件",
    ],
}

# 精选生物医药单点技能（Registry tool_id → 展示名）
_CURATED_SKILL_ENTRIES: List[IntentRegistryEntry] = [
    IntentRegistryEntry(
        id="lipinski_druglikeness",
        name="Lipinski 类药性筛查",
        aliases=["五规则", "类药性", "口服小分子", "lipinski", "成药潜势"],
        description="基于 Lipinski 五规则对 SMILES 做口服小分子成药潜势快筛（本地计算）。",
        route_kind=IntentRouteKind.skill,
        required_parameters=_SKILL_REQUIRED_PARAMETERS["lipinski_druglikeness"],
    ),
    IntentRegistryEntry(
        id="drug_similarity_search",
        name="药物结构相似性搜索",
        aliases=["分子相似", "pubchem", "chembl", "结构类似物", "tanimoto", "药物相似性评估"],
        description="在 PubChem/ChEMBL 等库检索结构相似化合物并生成 HTML 报告（需联网）。",
        route_kind=IntentRouteKind.skill,
        required_parameters=_SKILL_REQUIRED_PARAMETERS["drug_similarity_search"],
    ),
    IntentRegistryEntry(
        id="bepipred3_epitope_prediction",
        name="BepiPred 表位预测",
        aliases=["表位", "抗体表位", "bepipred", "免疫表位"],
        description="基于 BepiPred-3 对蛋白序列进行线性 B 细胞表位预测。",
        route_kind=IntentRouteKind.skill,
        required_parameters=_SKILL_REQUIRED_PARAMETERS["bepipred3_epitope_prediction"],
    ),
]

# 工作流 id → orchestrator target_domain（与 orchestrator._map 一致）
_WORKFLOW_TARGET_DOMAIN: Dict[str, str] = {
    "genomics_pipeline": "genomics",
    "rna_scrna": "rna",
    "metabolomics": "metabolomics",
    "spatial_transcriptomics": "spatial",
    "radiomics": "radiomics",
    "proteomics": "proteomics",
    "epigenomics": "epigenomics",
    "sted_ec_trajectory": "sted_ec_trajectory",
    "spatiotemporal_dynamics": "spatiotemporal_dynamics",
}


def _workflow_registry_entries() -> List[IntentRegistryEntry]:
    entries: List[IntentRegistryEntry] = []
    try:
        from gibh_agent.core.workflows import WorkflowRegistry

        wr = WorkflowRegistry()
        for domain_key, desc in (wr.list_workflows() or {}).items():
            dk = (domain_key or "").strip()
            if not dk:
                continue
            lid = dk.lower().replace(" ", "_")
            if dk == "RNA":
                reg_id, aliases = "rna_scrna", ["转录组", "单细胞", "scrna", "rna-seq", "bulk rna"]
            elif dk == "Metabolomics":
                reg_id, aliases = "metabolomics", ["代谢组", "lc-ms", "gc-ms", "代谢物"]
            elif dk == "Spatial":
                reg_id, aliases = "spatial_transcriptomics", ["空间转录", "visium", "空间组学"]
            elif dk == "Radiomics":
                reg_id, aliases = "radiomics", ["影像组学", "医学影像", "dicom", "nii"]
            elif dk == "genomics":
                reg_id, aliases = "genomics_pipeline", ["基因组", "wgs", "wes", "胚系", "变异检测"]
            elif dk == "proteomics":
                reg_id, aliases = "proteomics", ["蛋白组", "质谱蛋白"]
            elif dk == "epigenomics":
                reg_id, aliases = "epigenomics", ["表观组", "甲基化", "chip-seq"]
            elif dk == "STED_EC":
                reg_id, aliases = "sted_ec_trajectory", ["轨迹推断", "sted-ec", "细胞轨迹"]
            elif dk == "SPATIOTEMPORAL_DYNAMICS":
                reg_id, aliases = "spatiotemporal_dynamics", ["时空动力学", "单细胞时空"]
            else:
                reg_id, aliases = lid, [dk.lower()]

            td = _WORKFLOW_TARGET_DOMAIN.get(reg_id) or dk.lower()
            req_params = _WORKFLOW_REQUIRED_PARAMETERS.get(reg_id) or [
                "与本流程匹配的组学/影像原始数据文件或目录（须已上传挂载或给出可解析路径）",
            ]
            entries.append(
                IntentRegistryEntry(
                    id=reg_id,
                    name=dk if dk not in ("RNA", "Metabolomics") else (
                        "单细胞/转录组全流程" if dk == "RNA" else "代谢组学全流程"
                    ),
                    aliases=aliases,
                    description=(desc or "").strip()[:400],
                    route_kind=IntentRouteKind.workflow,
                    target_domain=td,
                    required_parameters=req_params,
                )
            )
    except Exception as e:
        logger.warning("intent_router: WorkflowRegistry 加载失败: %s", e)
    return entries


def _dynamic_skill_registry_entries() -> List[IntentRegistryEntry]:
    """从 gibh_agent/skills/ 动态发现的 BaseSkill（Skill Factory）。"""
    try:
        from gibh_agent.core.skill_registry import (
            discover_and_register_skills,
            get_discovered_skill_metas,
        )

        if not get_discovered_skill_metas():
            discover_and_register_skills()
        entries: List[IntentRegistryEntry] = []
        for m in get_discovered_skill_metas():
            aliases = list(m.aliases or [])
            if m.tool_chain_key and m.tool_chain_key not in aliases:
                aliases.append(m.tool_chain_key)
            req = list(m.required_parameters or [])
            if not req:
                req = [
                    "与本技能匹配的数据文件或序列（file_path / sequence_or_path 等，须已上传挂载）",
                ]
            entries.append(
                IntentRegistryEntry(
                    id=m.skill_id,
                    name=m.display_name or m.skill_id,
                    aliases=aliases[:12],
                    description=(m.description or "")[:400],
                    route_kind=IntentRouteKind.skill,
                    required_parameters=req[:8],
                )
            )
        return entries
    except Exception as e:
        logger.warning("intent_router: 动态技能注册表加载失败: %s", e)
        return []


def build_global_intent_registry() -> List[IntentRegistryEntry]:
    """合并工作流 SOP、精选技能与动态技能目录，供 Prompt 注入。"""
    seen: set = set()
    out: List[IntentRegistryEntry] = []
    for e in (
        _workflow_registry_entries()
        + _CURATED_SKILL_ENTRIES
        + _dynamic_skill_registry_entries()
    ):
        if e.id in seen:
            continue
        seen.add(e.id)
        out.append(e)
    return out


def registry_entry_by_id(entry_id: str) -> Optional[IntentRegistryEntry]:
    eid = (entry_id or "").strip()
    if not eid:
        return None
    for e in build_global_intent_registry():
        if e.id == eid:
            return e
    return None


def build_master_intent_router_system_prompt(registry: List[IntentRegistryEntry]) -> str:
    """领域无关的全局 ReAct 意图 Prompt：注册表注入 + 四态 JSON + 动态缺参示例。"""
    compact = [
        {
            "id": e.id,
            "name": e.name,
            "aliases": e.aliases[:8],
            "description": (e.description or "")[:200],
            "route_kind": e.route_kind.value,
            "required_parameters": e.required_parameters[:6],
        }
        for e in registry
    ]
    registry_json = json.dumps(compact, ensure_ascii=False, indent=2)
    return f"""你是 Omics Agent 的「全局意图泛化路由引擎」（Universal ReAct Router）。你的唯一依据是下方【工具/技能注册表】及用户上下文（message、has_uploaded_files、conversation_history）。

ReAct 流程（内部完成，禁止输出思考过程）：(1) 解析用户要做什么；(2) 在注册表中**指哪打哪**匹配至多一个执行目标；(3) 按 required_parameters 判断槽位是否齐备；(4) 输出**一个** JSON（禁止 Markdown 围栏与多余字段）。

【工具/技能注册表 — 每项含 required_parameters】
{registry_json}

【输出契约 — 纯 JSON，四态之一】

1) EXACT — 已唯一命中注册表 id，且 required_parameters 在**当前句、上传标记、或对话历史**中可执行（含经知识映射后的等效参数）：
{{"status":"exact","target_id":"<注册表 id>","rationale_short":"<≤40字>"}}

2) MISSING_PARAM — 已唯一命中 id，但**完全**不足以启动该工具（无文件、无路径、无可映射实体、无格式正确的输入）：
{{"status":"missing_param","target_id":"<注册表 id>","question":"<专业且友好的中文追问>","candidates":[{{"id":"<必须=target_id>","name":"<动态示例1>","type":"<skill|workflow>"}},{{"id":"<target_id>","name":"<动态示例2>","type":"..."}},...],"rationale_short":"<≤40字>"}}
- **必须**输出 question + **2～3 个** candidates；candidates **只能**由你根据**当前 target_id 的 required_parameters 与 description 动态撰写**，禁止照搬无关领域示例。
- 每个 candidate.id **必须等于** target_id；name 为可点击短句（如化合物工具可用「以 <小分子俗名> 为例」；基因工具可用「以 <基因符号> 序列为例」；组学流程可用「使用内置演示数据集」「查看输入格式要求」等——须与**该工具**一致）。
- 若用户已给出可唯一映射的常识实体 → **禁止** missing_param，应 exact。

3) CLARIFY — 合理命中 ≥2 个注册表 id，用户未明确择一：
{{"status":"clarify","question":"<澄清问句>","candidates":[{{"id":"<注册表 id>","name":"<展示名>"}},...],"rationale_short":"<≤40字>"}}
- candidates 2～5 个，id 均来自注册表。

4) CHITCHAT / 通用文件处理（GENERAL_QA）— 下列情形**必须** chitchat，**禁止** exact 到任何组学 workflow 或生信 skill：
- (a) 纯闲聊、百科、写作润色、翻译、找代码 Bug、总结普通文档等**无**注册表工具执行诉求；
- (b) 用户挂载的是**通用**文本/代码/表格（.txt/.md/.py/.csv 等）且需求为「总结」「翻译」「解释」「改错」「写文章」等**非**专业生信分析 —— 一律 chitchat，**不得**因 has_uploaded_files=true 就拉起 genomics/RNA/代谢等工作流；
- (c) 需求**超出**注册表能力（如分子对接、AlphaFold、未列出的专有软件），须委婉说明当前平台不支持，**禁止**编造 target_id。
{{"status":"chitchat","rationale_short":"<说明不支持、通用问答或转为科普>"}}

【全局知识映射与参数转换原则 (Universal Knowledge Mapping)】
注册表覆盖组学流程、小分子/蛋白序列技能等；各工具 required_parameters 可能要求 FASTA、SMILES、表达矩阵路径、影像文件对等严格格式。
1. **智能推断**：用户提供常识性自然语言实体（化合物俗名、基因符号如 TP53、蛋白名、细胞类型、疾病语境下的标准分析对象等）时，用生物医学常识判断能否**唯一**映射到该工具所需输入。能映射 → 视为参数已齐备 → exact；**不得**在已可映射时仍索要底层格式。
2. **拒绝幻觉**：target_id、candidates[].id **必须**来自注册表；禁止为用户未提及的数据编造具体 SMILES/突变/文件路径并 exact；禁止将无关工具强行关联（节外生枝）。
3. **指哪打哪**：命中工具后不得改判无关 id；仅用户明确更换分析类型或 clarify 选型后可切换。

【MISSING_PARAM 动态泛化要求】
question 须说明该工具**缺什么**（对照 required_parameters）。candidates 的 name 须是**该工具领域**的一键示例（由你生成，非固定模板）：
- 小分子/化合物类工具 → 可选不同俗名或「上传 SMILES 文件」类短句；
- 序列/表位类 → 不同基因/蛋白示例或「粘贴 FASTA」；
- 组学/影像 workflow → 演示数据集、路径占位、格式说明类短句。
type 字段与注册表 route_kind 一致（skill / workflow）。

【多轮上下文 — 任务延续】
结合 conversation_history：上一轮 missing_param/clarify 后，用户补充实体、换对象、或点击式短句 → 延续同一分析意图与 target_id（除非用户明确改换）；补充后若可映射 → exact。

【槽位齐备（领域无关）】
- 结构化字符串（SMILES、FASTA、路径、矩阵目录名）已出现 → 齐备。
- 可唯一映射的常识实体 → 齐备。
- has_uploaded_files=true 且与当前工具数据类型相符 → 齐备。
- 仅「分析一下」「帮我跑」无任何对象 → missing_param（动态 candidates）或 clarify（多工具）。

【AMBIGUOUS 典型模式（仍须 clarify，非 exhaustive）】
- 化合物「评估/看看」未区分类药性 vs 相似检索 → 对应注册表多项。
- 「组学/测序/图像/蛋白/轨迹」泛称未指明流程 → 从注册表选 2～5 个最相关 workflow/skill。

【总原则】
- 超出注册表 → chitchat，不 hallucinate 工具。
- 通用文件 + 非生信诉求 → chitchat（FILE_PROCESSING 语义），绝不误触组学流水线。
- 多工具 → clarify；单工具缺参 → missing_param（**LLM 自生成** 2～3 candidates）；单工具够参 → exact。
- rationale_short ≤40 字。
"""


def _extract_history_text(content: Any) -> str:
    if content is None:
        return ""
    if isinstance(content, str):
        return content.strip()
    if isinstance(content, dict):
        for key in ("text", "message", "content"):
            val = content.get(key)
            if isinstance(val, str) and val.strip():
                return val.strip()
        clar = content.get("clarification")
        if isinstance(clar, dict) and clar.get("question"):
            return str(clar["question"]).strip()
        try:
            return json.dumps(content, ensure_ascii=False)[:600]
        except Exception:
            return str(content)[:600]
    return str(content).strip()


def trim_history_for_intent_router(
    history: Optional[List[Dict[str, Any]]],
    *,
    max_messages: int = 10,
) -> List[Dict[str, str]]:
    """保留最近若干条 user/assistant 轮次，供 MasterIntentRouter 多轮推理。"""
    trimmed: List[Dict[str, str]] = []
    for h in history or []:
        if not isinstance(h, dict):
            continue
        role = (h.get("role") or "").strip().lower()
        if role == "agent":
            role = "assistant"
        if role not in ("user", "assistant"):
            continue
        text = _extract_history_text(h.get("content") or h.get("message"))
        if not text:
            continue
        trimmed.append({"role": role, "content": text[:800]})
    return trimmed[-max_messages:]


def _serialize_user_context(
    message: str,
    has_files: bool,
    history: Optional[List[Dict[str, Any]]] = None,
) -> str:
    return json.dumps(
        {
            "message": (message or "").strip(),
            "has_uploaded_files": bool(has_files),
            "conversation_history": trim_history_for_intent_router(history),
        },
        ensure_ascii=False,
        indent=2,
    )


def master_intent_router_enabled() -> bool:
    raw = (os.environ.get("GIBH_MASTER_INTENT_ROUTER") or "1").strip().lower()
    return raw not in ("0", "false", "no", "off")


def should_bypass_master_intent_router(req: Any) -> bool:
    """与快车道、工作流执行、澄清闭环等约定对齐的绕过条件。"""
    if getattr(req, "intent_override_id", None):
        return True
    if getattr(req, "target_domain", None):
        td = (req.target_domain or "").strip()
        if td:
            return True
    if getattr(req, "workflow_data", None):
        return True
    if getattr(req, "test_dataset_id", None):
        ts = (req.test_dataset_id or "").strip()
        if ts:
            return True
    if getattr(req, "local_tool_response", None):
        return True
    msg = (getattr(req, "message", None) or "").strip()
    if msg.startswith("⚡ 启动工作流："):
        return True
    if re.search(r"\[Skill_Route:\s*\w+", msg, re.I):
        return True
    if "[Omics_Route:" in msg:
        return True
    return False


def apply_intent_dispatch(
    *,
    message: str,
    target_id: str,
) -> Dict[str, Any]:
    """
    将 intent_override_id / exact 的 target_id 转为 orchestrator 可消费的 kwargs 补丁。

    Returns:
        dict 可能含 target_domain、message（注入 Skill_Route）
    """
    entry = registry_entry_by_id(target_id)
    if not entry:
        return {}
    out: Dict[str, Any] = {}
    if entry.route_kind == IntentRouteKind.workflow and entry.target_domain:
        out["target_domain"] = entry.target_domain
    elif entry.route_kind == IntentRouteKind.skill:
        base = (message or "").strip()
        prefix = f"[Skill_Route: {entry.id}]"
        if not re.search(r"\[Skill_Route:\s*" + re.escape(entry.id), base, re.I):
            out["message"] = f"{prefix}\n我选择：{entry.name}\n{base}".strip()
    return out


class MasterIntentRouter:
    """领域无关 LLM 全局路由：注册表注入 + 历史感知 + 四态 JSON（缺参示例由 LLM 动态生成）。"""

    CONFIDENCE_FALLBACK_CLARIFY = "无法可靠解析意图，请选择最接近的分析类型。"

    def __init__(self, llm_client: LLMClient) -> None:
        self._llm = llm_client
        self._registry = build_global_intent_registry()

    async def classify(
        self,
        message: str,
        *,
        has_files: bool = False,
        history: Optional[List[Dict[str, Any]]] = None,
    ) -> IntentRouterResult:
        system_prompt = build_master_intent_router_system_prompt(self._registry)
        user_blob = _serialize_user_context(message, has_files, history)
        messages = [
            {"role": "system", "content": system_prompt},
            {
                "role": "user",
                "content": (
                    "请结合 conversation_history 做 ReAct 意图分类，仅输出 JSON：\n\n"
                    + user_blob
                ),
            },
        ]
        try:
            raw = await self._llm.achat(
                messages,
                temperature=0.1,
                max_tokens=768,
            )
            text = ""
            if hasattr(raw, "choices") and raw.choices:
                text = (raw.choices[0].message.content or "").strip()
            else:
                text = str(raw or "").strip()
        except Exception as e:
            logger.warning("MasterIntentRouter LLM 调用失败: %s", e, exc_info=True)
            return IntentRouterResult(
                status="clarify",
                question=self.CONFIDENCE_FALLBACK_CLARIFY,
                candidates=[
                    {"id": c.id, "name": c.name}
                    for c in self._registry[:5]
                ],
                rationale_short="路由模型不可用，降级为澄清",
            )

        parsed = extract_json_object_from_llm_text(text)
        if not parsed or not isinstance(parsed, dict):
            logger.warning("MasterIntentRouter JSON 解析失败: %s", text[:300])
            return IntentRouterResult(
                status="clarify",
                question=self.CONFIDENCE_FALLBACK_CLARIFY,
                candidates=[{"id": c.id, "name": c.name} for c in self._registry[:4]],
                rationale_short="JSON 解析失败",
            )

        status = (parsed.get("status") or "").strip().lower()
        if status not in ("exact", "missing_param", "clarify", "chitchat"):
            status = "clarify"

        parsed_candidates: List[Dict[str, str]] = []
        for c in (parsed.get("candidates") or []):
            if not isinstance(c, dict):
                continue
            cid = str(c.get("id", "")).strip()
            cname = str(c.get("name", "")).strip()
            if not cid or not cname:
                continue
            item: Dict[str, str] = {"id": cid, "name": cname}
            ctype = str(c.get("type", "")).strip()
            if ctype:
                item["type"] = ctype
            parsed_candidates.append(item)

        result = IntentRouterResult(
            status=status,
            target_id=(parsed.get("target_id") or "").strip() or None,
            question=(parsed.get("question") or "").strip() or None,
            candidates=parsed_candidates,
            rationale_short=(parsed.get("rationale_short") or "")[:200],
        )

        if result.status in ("exact", "missing_param") and result.target_id:
            if not registry_entry_by_id(result.target_id):
                logger.warning(
                    "MasterIntentRouter 未知 target_id=%s，改 clarify", result.target_id
                )
                result.status = "clarify"
                result.target_id = None
                result.question = result.question or "未能唯一匹配，请选择您需要的分析类型："
                result.candidates = result.candidates or [
                    {"id": c.id, "name": c.name} for c in self._registry[:4]
                ]

        if result.status == "missing_param" and not result.question:
            entry = registry_entry_by_id(result.target_id or "")
            if entry and entry.required_parameters:
                req_hint = "、".join(entry.required_parameters[:2])
                result.question = (
                    f"请提供执行「{entry.name}」所需的参数：{req_hint}。"
                )
            else:
                result.question = "请补充执行任务所需的参数或数据文件。"

        if result.status == "missing_param" and result.target_id:
            result.candidates = sanitize_missing_param_candidates(
                result.target_id,
                result.candidates,
            )
            if len(result.candidates) < 2:
                logger.warning(
                    "MasterIntentRouter missing_param 候选不足（%s 条），"
                    "须由 LLM 动态生成 2～3 个；target_id=%s",
                    len(result.candidates),
                    result.target_id,
                )

        if result.status == "clarify" and not result.candidates:
            result.candidates = [{"id": c.id, "name": c.name} for c in self._registry[:4]]

        logger.info(
            "🧭 [MasterIntentRouter] status=%s target_id=%s candidates=%s reason=%s",
            result.status,
            result.target_id,
            len(result.candidates),
            result.rationale_short,
        )
        return result


def sanitize_missing_param_candidates(
    target_id: str,
    candidates: List[Dict[str, str]],
) -> List[Dict[str, str]]:
    """
    missing_param 候选的结构化清洗：id 对齐 target_id、去重、补全 type。
    不注入任何领域硬编码示例——示例文案 100% 由 LLM 生成。
    """
    tid = (target_id or "").strip()
    if not tid:
        return []
    entry = registry_entry_by_id(tid)
    route_type = entry.route_kind.value if entry else "skill"
    out: List[Dict[str, str]] = []
    seen_names: set = set()
    for raw in candidates or []:
        if not isinstance(raw, dict):
            continue
        name = (raw.get("name") or "").strip()
        if not name or name in seen_names:
            continue
        seen_names.add(name)
        ctype = (raw.get("type") or route_type).strip() or route_type
        out.append({"id": tid, "name": name, "type": ctype})
    return out[:5]


def enrich_clarification_candidates(
    candidates: List[Dict[str, str]],
) -> List[Dict[str, str]]:
    """为 SSE 澄清候选补充 route_kind（skill/workflow），供前端填充快车道 Prompt。"""
    out: List[Dict[str, str]] = []
    for raw in candidates or []:
        if not isinstance(raw, dict):
            continue
        cid = (raw.get("id") or "").strip()
        if not cid:
            continue
        entry = registry_entry_by_id(cid)
        out.append(
            {
                "id": cid,
                "name": (raw.get("name") or "").strip()
                or (entry.name if entry else cid),
                "type": (
                    (raw.get("type") or "").strip()
                    or (
                        entry.route_kind.value
                        if entry
                        else (
                            "skill"
                            if cid in {e.id for e in _CURATED_SKILL_ENTRIES}
                            else "workflow"
                        )
                    )
                ),
            }
        )
    return out


def format_clarification_sse_payload(result: IntentRouterResult) -> Dict[str, Any]:
    """供 event: clarification 使用的 data 载荷。"""
    return {
        "type": "clarification",
        "question": result.question or "请选择您需要的分析类型：",
        "candidates": enrich_clarification_candidates(result.candidates),
        "rationale_short": result.rationale_short,
    }
