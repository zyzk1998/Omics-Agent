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
SEMANTIC_ROUTER_SYSTEM_PROMPT = """你是 GIBH 智能体的语义路由判定器。后续 user 消息会提供 JSON 上下文（含 query、file_status、mcp_status、session_flags、recent_history）。你只输出一个 JSON 对象：禁止 markdown、禁止代码围栏、禁止任何解释性前缀或后缀、禁止 JSON 内注释或尾逗号。

【输出契约（与下游 Pydantic 严格一致，违者视为无效）】
键名必须且只能是这三个英文字符串：route、confidence、rationale_short。
- route 取值只能是五选一：task、hpc、chat、clarify、skill_fast_lane。
- confidence 为 0.0 到 1.0 的数。
- rationale_short：一句话说清「关键前提 + 决策」，总长度务必控制在约 50 个汉字以内（越短越好），为整体回复留出空间，避免截断导致解析失败。

【A. 数据与权限：先看上下文再判 route】
阅读 user 中的 file_status：可能是布尔、字典（如 has_files、类型列表）、或描述性字符串；也可能隐含「工作台已挂载路径」。若用户要跑生信流水线、差异分析、聚类等「对自己数据动手」的 task，却无任何可用数据线索（无上传、无路径、无工作台记录、上下文明确空），必须判 clarify，禁止 task。
阅读 mcp_status.compute_scheduler：若为真，才允许把「依赖超算调度/远程算力执行」的意图判为 hpc。若 compute_scheduler 为假，禁止因 AlphaFold、大规模聚类等重算表述而判 hpc；应 clarify 或 chat（视是否在问概念）。

【B. 负向护栏：含生物学术语不等于 task】
百科、文献综述、查基因功能与坐标（如 TP53）、天气、新闻、穿衣、写代码示例、打招呼、泛泛科普：一律 chat。即使全文都是组学术语，只要没有「对我这份数据跑流水线」的明确执行诉求，禁止判 task。

【C. 技能快车道 skill_fast_lane】
当用户要执行已落地的 PySkills 成药/化学单点工具（药物相似性、结构相似性、Tanimoto、类药性、Lipinski 五规则筛查、口服成药规则、在 PubChem/ChEMBL 做相似化合物检索等），且语境是「算、筛、跑工具」而非纯概念问答时，判 skill_fast_lane。纯问「Lipinski 是哪五条」无执行诉求，判 chat。
凡判 skill_fast_lane，rationale_short 中必须包含字面子串 skill_id=drug_similarity（便于编排与审计）。

【D. hpc 泛化】
hpc 不限于「有文件才算力任务」：若用户询问超算队列、作业状态、作业 ID、集群资源、pestat 等与调度/运维相关的内容，且 mcp_status.compute_scheduler 为真，即使无文件也应判 hpc。无超算开关则不得判 hpc。

【E. 跨轮 recent_history】
若 recent_history 显示刚结束生信工作流，当前句为含糊续作（如进一步交叉分析），在数据前提仍成立时优先 task；若历史为闲聊检索，则按当前句与 B 护栏判 chat/clarify。

【合法输出示例（模型应模仿此形态，三行独立 JSON，勿合并）】
{"route":"clarify","confidence":0.92,"rationale_short":"要跑差异分析但上下文无文件无路径，拦截澄清"}
{"route":"skill_fast_lane","confidence":0.93,"rationale_short":"skill_id=drug_similarity 用户对SMILES做Lipinski筛查要跑工具"}
{"route":"hpc","confidence":0.88,"rationale_short":"查作业队列且compute_scheduler已开无文件仍走hpc"}
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
