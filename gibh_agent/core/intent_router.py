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


class IntentRouterResult(BaseModel):
    """LLM 意图探针输出（与产品约定的三种状态对齐）。"""

    status: str  # exact | clarify | chitchat
    target_id: Optional[str] = None
    question: Optional[str] = None
    candidates: List[Dict[str, str]] = Field(default_factory=list)
    rationale_short: str = ""


# 精选生物医药单点技能（Registry tool_id → 展示名）
_CURATED_SKILL_ENTRIES: List[IntentRegistryEntry] = [
    IntentRegistryEntry(
        id="lipinski_druglikeness",
        name="Lipinski 类药性筛查",
        aliases=["五规则", "类药性", "口服小分子", "lipinski", "成药潜势"],
        description="基于 Lipinski 五规则对 SMILES 做口服小分子成药潜势快筛（本地计算）。",
        route_kind=IntentRouteKind.skill,
    ),
    IntentRegistryEntry(
        id="drug_similarity_search",
        name="药物结构相似性搜索",
        aliases=["分子相似", "pubchem", "chembl", "结构类似物", "tanimoto"],
        description="在 PubChem/ChEMBL 等库检索结构相似化合物并生成 HTML 报告（需联网）。",
        route_kind=IntentRouteKind.skill,
    ),
    IntentRegistryEntry(
        id="bepipred3_epitope_prediction",
        name="BepiPred 表位预测",
        aliases=["表位", "抗体表位", "bepipred", "免疫表位"],
        description="基于 BepiPred-3 对蛋白序列进行线性 B 细胞表位预测。",
        route_kind=IntentRouteKind.skill,
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
                )
            )
    except Exception as e:
        logger.warning("intent_router: WorkflowRegistry 加载失败: %s", e)
    return entries


def build_global_intent_registry() -> List[IntentRegistryEntry]:
    """合并工作流 SOP 与精选技能，供 Prompt 注入。"""
    seen: set = set()
    out: List[IntentRegistryEntry] = []
    for e in _workflow_registry_entries() + _CURATED_SKILL_ENTRIES:
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
    """向 LLM 注入工具注册表并规定 JSON 输出契约。"""
    compact = [
        {
            "id": e.id,
            "name": e.name,
            "aliases": e.aliases[:8],
            "description": (e.description or "")[:200],
            "route_kind": e.route_kind.value,
        }
        for e in registry
    ]
    registry_json = json.dumps(compact, ensure_ascii=False, indent=2)
    return f"""你是 Omics Agent 的「全局意图路由网关」。用户将用自然语言描述分析需求（含专业名词、模糊表述或白话）。
你只能根据下方【工具/技能注册表】做意图分类，并输出**一个** JSON 对象（禁止 Markdown 围栏、禁止解释文字）。

【工具/技能注册表】
{registry_json}

【输出契约 — 三选一，键名必须严格一致】

1) 极大概率唯一匹配（EXACT_MATCH）：
{{"status":"exact","target_id":"<注册表 id>","rationale_short":"<一句中文理由>"}}

2) 匹配多个或意图模糊（AMBIGUOUS）：
{{"status":"clarify","question":"<面向用户的中文澄清问句>","candidates":[{{"id":"<id>","name":"<展示名>"}},...],"rationale_short":"<理由>"}}
- candidates 必须从注册表中选取 2～5 个最相关项，id/name 与注册表一致。

3) 与平台工具无关的闲聊、百科、天气、写代码示例、打招呼（CHITCHAT）：
{{"status":"chitchat","rationale_short":"<理由>"}}

【判定规则 — 总原则】
- **宁可 clarify，不可误 exact**：只要合理对应注册表中 **≥2 个** id，且用户未明确点名某一流程/工具，必须判 clarify，禁止 exact。
- 仅当用户**明确且唯一**指向某一注册表 id（或无可争议的单一域，如「跑 Lipinski 五规则筛查」）才判 exact。
- 仅概念问答（如「什么是 WGS」「科普一下转录组」）无「帮我跑/分析我的数据」执行诉求：chitchat。
- 禁止编造注册表外的 id；candidates 的 id 必须来自注册表。
- rationale_short 控制在约 40 字以内。

【必须判 clarify 的典型模糊模式（非 exhaustive）】
- **药物/化合物**：「评估/分析一下这个化合物/小分子/药物」但未区分类药性(Lipinski) vs 结构相似检索 → candidates 含 lipinski_druglikeness 与 drug_similarity_search。
- **组学泛称**：「跑个组学/做组学分析/上机组学」未指明基因组/转录/蛋白/代谢/表观/空间 → candidates 至少含 genomics_pipeline、rna_scrna、metabolomics、proteomics 中多项。
- **变异/测序文件**：「BAM/VCF/FASTQ + 找突变/变异」未指明胚系 WGS、RNA、表观或蛋白层面 → candidates 含 genomics_pipeline，并搭配 rna_scrna 或 epigenomics 等合理备选项。
- **转录/表达**：「测序下机了想分析一下」未说明 bulk RNA、单细胞还是空间 → candidates 含 rna_scrna 与 spatial_transcriptomics（及/或 metabolomics）。
- **影像/图像**：「医学图像/CT/MRI/病理切片」未说明影像组学特征 vs 空间组学 → candidates 含 radiomics 与 spatial_transcriptomics。
- **蛋白/免疫**：「蛋白序列分析一下」可能为表位预测或蛋白组流程 → candidates 含 bepipred3_epitope_prediction 与 proteomics。
- **细胞/轨迹**：「细胞运动/轨迹」可能 STED_EC 或时空动力学 → candidates 含 sted_ec_trajectory 与 spatiotemporal_dynamics。
- **代谢**：「质谱数据想看一下」未指明靶向/非靶向代谢组 → 若仅对应 metabolomics 一项可 exact；若同时像蛋白组则 clarify（metabolomics + proteomics）。

【禁止误判】
- 上述模糊模式**禁止**判 chitchat（用户有执行/analysis 诉求）。
- 上述模糊模式**禁止**只返回 1 个 candidate（必须 2～5 个）。
"""


def _serialize_user_context(message: str, has_files: bool) -> str:
    return json.dumps(
        {"message": (message or "").strip(), "has_uploaded_files": bool(has_files)},
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
    """轻量 LLM 全局路由：注册表注入 + 三态 JSON 分类。"""

    CONFIDENCE_FALLBACK_CLARIFY = "无法可靠解析意图，请选择最接近的分析类型。"

    def __init__(self, llm_client: LLMClient) -> None:
        self._llm = llm_client
        self._registry = build_global_intent_registry()

    async def classify(
        self,
        message: str,
        *,
        has_files: bool = False,
    ) -> IntentRouterResult:
        system_prompt = build_master_intent_router_system_prompt(self._registry)
        user_blob = _serialize_user_context(message, has_files)
        messages = [
            {"role": "system", "content": system_prompt},
            {
                "role": "user",
                "content": "请对以下用户输入做意图分类，仅输出 JSON：\n\n" + user_blob,
            },
        ]
        try:
            raw = await self._llm.achat(
                messages,
                temperature=0.1,
                max_tokens=512,
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
        if status not in ("exact", "clarify", "chitchat"):
            status = "clarify"

        result = IntentRouterResult(
            status=status,
            target_id=(parsed.get("target_id") or "").strip() or None,
            question=(parsed.get("question") or "").strip() or None,
            candidates=[
                {"id": str(c.get("id", "")), "name": str(c.get("name", ""))}
                for c in (parsed.get("candidates") or [])
                if isinstance(c, dict) and c.get("id")
            ],
            rationale_short=(parsed.get("rationale_short") or "")[:200],
        )

        if result.status == "exact" and result.target_id:
            if not registry_entry_by_id(result.target_id):
                logger.warning("MasterIntentRouter 未知 target_id=%s，改 clarify", result.target_id)
                result.status = "clarify"
                result.target_id = None
                result.question = result.question or "未能唯一匹配，请选择您需要的分析类型："
                result.candidates = result.candidates or [
                    {"id": c.id, "name": c.name} for c in self._registry[:4]
                ]

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
                    entry.route_kind.value
                    if entry
                    else ("skill" if cid in {e.id for e in _CURATED_SKILL_ENTRIES} else "workflow")
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
