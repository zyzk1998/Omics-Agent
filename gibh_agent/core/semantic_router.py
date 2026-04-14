# -*- coding: utf-8 -*-
"""
语义路由中枢（Semantic Router）— Phase 1 沙盒模块

与编排器、server、执行器解耦：仅根据用户 query 与上下文元数据，调用统一 LLMClient
产出结构化路由决策（JSON），供后续独立接入。

设计要点：
- 强规则写在 System Prompt 中，由模型遵守；解析失败或低置信度走 Fail-safe，避免误路由到重型任务。
"""
from __future__ import annotations

import json
import logging
import re
from enum import Enum
from typing import Any, Dict, Optional, Union

from pydantic import BaseModel, Field, field_validator

from gibh_agent.core.llm_client import LLMClient

logger = logging.getLogger(__name__)


class RouteKind(str, Enum):
    """路由目标：与产品约定的五类出口。"""

    task = "task"  # 需要走工具/工作流执行的分析类任务
    hpc = "hpc"  # 超算 / 调度器相关（如作业状态、队列）
    chat = "chat"  # 闲聊、通用知识、无执行诉求
    clarify = "clarify"  # 缺前提、需用户补充（如无文件却要跑流程）
    skill_fast_lane = "skill_fast_lane"  # 明确技能快捷通道（预留）


class RouterInput(BaseModel):
    """
    路由输入：聚合「用户说了什么」与「当前会话/系统状态」。

    file_status 允许 dict 或 str：
    - dict：如 {"has_files": true, "types": [".fastq"]}，便于程序构造；
    - str：兼容上游直接把描述文本塞进来。
    """

    query: str = Field(..., description="用户自然语言")
    file_status: Union[Dict[str, Any], str] = Field(
        default_factory=dict,
        description="是否有文件、类型列表或自由文本描述",
    )
    mcp_status: Dict[str, Any] = Field(
        default_factory=dict,
        description="MCP 能力开关，例如 compute_scheduler 是否可用",
    )
    session_flags: Dict[str, Any] = Field(
        default_factory=dict,
        description="会话级标记，如 target_domain、实验上下文等",
    )
    recent_history: str = Field(
        default="",
        description="最近 2～3 轮对话的压缩文本（用户/助手轮替），供路由理解追问与跨轮语境；无则空字符串",
    )


class RouterOutput(BaseModel):
    """路由输出：必须可被 JSON 反序列化并校验。"""

    route: RouteKind
    confidence: float = Field(..., ge=0.0, le=1.0, description="0~1，低于阈值触发重试或 clarify")
    rationale_short: str = Field(
        ...,
        max_length=512,
        description="简短判别理由，便于日志与排障",
    )

    @field_validator("rationale_short")
    @classmethod
    def strip_rationale(cls, v: str) -> str:
        return (v or "").strip()


# ---------------------------------------------------------------------------
# System Prompt：业务红线写死在此，模型必须优先遵守（与代码 Fail-safe 双保险）
# ---------------------------------------------------------------------------
SEMANTIC_ROUTER_SYSTEM_PROMPT = """你是 GIBH 智能体的「语义路由」判定器。你只输出一个 JSON 对象，不要 markdown、不要代码围栏、不要多余解释。
JSON 必须严格包含三个键：
- "route": 必须是 task、hpc、chat、clarify、skill_fast_lane 之一
- "confidence": 0.0 到 1.0 的小数
- "rationale_short": 极短说明原因（<120字）

【最高优先级：反向红线与意图剥离 (Negative Guardrails)】
!!! 警告：你极易因看到“基因、测序”或遇到“没见过的句式”而误判为 task。必须强制遵守以下规则 !!!
1. 事实检索与生活常识免疫：如果用户的意图是获取外部知识、天气、新闻、穿衣建议，或查询生物数据库（如 TP53 基因功能）。即使句中包含复杂的专有名词，也绝对属于 Chat 范畴，必须且只能路由到 `chat`！绝对禁止判为 `task`。
2. Task 的物理边界：`task` 的唯一合法场景是「用户要对自己的数据文件进行加工、跑生信流水线或分析（DAG）」。

【常规判定规则（优先级递减）】
1) clarify（缺前提）
   - 若用户明确要「跑分析管线」，但 file_status 表明无可用数据，且问题不是纯概念咨询，选 clarify。
2) hpc（超算）
   - 仅当 mcp_status.compute_scheduler 为 true 时允许。涉及超算运维、作业状态、pestat 等选 hpc。
3) chat（对话/检索）
   - 普通问候、生活闲聊/咨询（如“怎么穿衣”）、实时联网检索（天气/新闻），或查询医学百科级知识，选 chat。
4) task（管线分析）
   - 用户明确要跑组学分析、流程等，且 file_status 有匹配文件就绪时，选 task。
5) skill_fast_lane（技能快车道：无上传文件也可选）
   - 当用户意图是**已落地的 PySkills 化学/成药技能**（见下），且不是「泛泛了解 Lipinski 概念」式的百科问答时，选 `skill_fast_lane`。
   - **药物相似性 / 药物筛选 / 类药性 / Lipinski / 五规则 / 口服成药 / 分子相似度 / 结构相似性 / Tanimoto / 类似物发现 / PubChem 检索 / ChEMBL 相似化合物** 等表述，只要语境是「算一算、筛一筛、跑工具」而非纯概念课，一律判为 `skill_fast_lane`。
   - 此时 `rationale_short` **必须**包含字面片段 **`skill_id=drug_similarity`**（便于编排器与审计识别逻辑技能包）；可另附简短中文，例如：`skill_id=drug_similarity，用户要 Lipinski 快筛` 或 `skill_id=drug_similarity，结构相似性检索`。
   - 若同时像闲聊（例如仅问「Lipinski 五规则是什么」）且无执行诉求，选 `chat`。

【上下文感知原则】（若 recent_history 非空则必须结合判断）
- 若 recent_history 显示刚完成某条生信工作流，且当前 query 是含糊追问（如「进一步做交叉分析」），优先判 task，继续走分析链。
- 若 recent_history 为闲聊/检索，按 query 本身判 chat 或 clarify。

【强制分类防幻觉样本 (Few-Shot Guardrails)】
- "广州出门穿什么合适？" -> route: "chat" (生活常识/天气检索)
- "查询 TP53 基因的功能和染色体位置" -> route: "chat" (调库检索，非数据管线)
- "你能帮我搜一下关于 CRISPR 的文献吗" -> route: "chat" (文献检索)
- "帮我对这两个 fastq 文件做差异表达分析" -> route: "task" (针对文件的分析管线)
- "进一步做跨组学交叉分析" (且 recent_history 刚完成转录组) -> route: "task" (基于有效上下文继续执行)
- "帮我对这个 SMILES 做一下 Lipinski 五规则筛查" -> route: "skill_fast_lane", rationale_short 含 skill_id=drug_similarity
- "在 PubChem 里找和这个分子结构相似的化合物" -> route: "skill_fast_lane", rationale_short 含 skill_id=drug_similarity
- "Lipinski 五规则是哪五条？"（纯概念） -> route: "chat"
"""


def _serialize_router_context(inp: RouterInput) -> str:
    """将结构化上下文压成一段 JSON 字符串，便于模型一次性阅读。"""
    payload = {
        "query": inp.query,
        "file_status": inp.file_status,
        "mcp_status": inp.mcp_status,
        "session_flags": inp.session_flags,
        "recent_history": (inp.recent_history or "").strip(),
    }
    return json.dumps(payload, ensure_ascii=False, indent=2)


def _extract_json_object(text: str) -> Optional[Dict[str, Any]]:
    """
    从模型原文中抠出 JSON 对象。

    部分模型仍可能包一层 ```json ... ```，此处做宽松剥离；失败则返回 None。
    """
    if not text or not text.strip():
        return None
    raw = text.strip()
    fence = re.search(r"```(?:json)?\s*([\s\S]*?)\s*```", raw, re.IGNORECASE)
    if fence:
        raw = fence.group(1).strip()
    raw = raw.strip()
    try:
        obj = json.loads(raw)
        if isinstance(obj, dict):
            return obj
    except json.JSONDecodeError:
        pass
    # 再尝试从文本中截取第一个 { ... } 块（弱兜底）
    start = raw.find("{")
    end = raw.rfind("}")
    if start != -1 and end != -1 and end > start:
        try:
            obj = json.loads(raw[start : end + 1])
            if isinstance(obj, dict):
                return obj
        except json.JSONDecodeError:
            return None
    return None


class SemanticRouter:
    """
    语义路由引擎：封装 LLMClient.achat，强制 JSON 形态输出并做阈值与重试 Fail-safe。

    不在此模块内调用 orchestrator / server；仅依赖 LLMClient。
    """

    # 置信度阈值：低于此值视为「不可靠」，触发一次重试；两次仍不足则 clarify
    CONFIDENCE_THRESHOLD = 0.7
    MAX_LLM_ATTEMPTS = 2

    def __init__(self, llm_client: LLMClient) -> None:
        self._llm = llm_client

    def _clarify_fallback(self, rationale_short: str) -> RouterOutput:
        """统一兜底：避免 silent failure 把用户导向错误执行链。"""
        return RouterOutput(
            route=RouteKind.clarify,
            confidence=0.0,
            rationale_short=rationale_short[:512],
        )

    async def decide_route(self, input_data: RouterInput) -> RouterOutput:
        """
        调用大模型做路由判定；解析/校验失败或置信度不足时最多再试 1 次，仍失败则 clarify。

        重试原因：
        - 模型偶发输出截断、非 JSON、或 confidence 偏低；
        - 一次重试成本远低于误路由到 HPC/任务执行带来的风险。
        """
        user_blob = _serialize_router_context(input_data)
        messages = [
            {"role": "system", "content": SEMANTIC_ROUTER_SYSTEM_PROMPT},
            {
                "role": "user",
                "content": (
                    "请根据以下上下文输出路由 JSON（仅 JSON 对象）：\n\n" + user_blob
                ),
            },
        ]

        last_issue = ""

        for attempt in range(self.MAX_LLM_ATTEMPTS):
            try:
                completion = await self._llm.achat(
                    messages,
                    temperature=0.1,
                    max_tokens=512,
                )
                choice0 = completion.choices[0]
                content = (choice0.message.content or "").strip()
            except Exception as e:
                last_issue = f"LLM 调用异常: {e!s}"
                logger.warning("SemanticRouter achat 失败 attempt=%s: %s", attempt, e)
                continue

            parsed = _extract_json_object(content)
            if parsed is None:
                last_issue = "模型返回无法解析为 JSON 对象"
                logger.warning(
                    "SemanticRouter JSON 解析失败 attempt=%s raw_prefix=%r",
                    attempt,
                    content[:200],
                )
                continue

            try:
                out = RouterOutput.model_validate(parsed)
            except Exception as e:
                last_issue = f"Pydantic 校验失败: {e!s}"
                logger.warning("SemanticRouter 输出结构非法 attempt=%s: %s", attempt, e)
                continue

            if out.confidence >= self.CONFIDENCE_THRESHOLD:
                return out

            last_issue = f"置信度 {out.confidence:.2f} 低于阈值 {self.CONFIDENCE_THRESHOLD}"
            logger.info("SemanticRouter 低置信度 attempt=%s: %s", attempt, last_issue)

        # 两次尝试仍无法得到可信结构化结果 → 强制 clarify，把原因写入 rationale 便于审计
        fb_msg = "路由判定失败或置信不足，已降级为澄清。"
        if last_issue:
            fb_msg = f"{fb_msg} 详情: {last_issue}"
        return self._clarify_fallback(fb_msg)
