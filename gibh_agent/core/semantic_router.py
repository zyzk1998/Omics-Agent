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

JSON 必须严格包含三个键，且仅这三个键：
- "route": 字符串，必须是以下之一：task、hpc、chat、clarify、skill_fast_lane
- "confidence": 0 到 1 之间的小数，表示你对该分类的确信程度
- "rationale_short": 一句中文，极短说明为何选该 route（不超过 120 字）

判定规则（按优先级思考，冲突时以更保守、更安全的为准）：

1) clarify（缺前提 / 需用户补充）
   - 若用户意图明显是要「执行」生信/组学/数据分析类任务，但 file_status 表明当前没有可用数据文件（或无路径、无上传），且用户问题不是「纯概念咨询」（例如不是只问「QC 是什么」），则必须选 clarify，让系统向用户索要文件或明确数据位置。
   - 若信息严重不足、无法区分是执行还是闲聊，也选 clarify。

2) hpc（超算 / 调度）
   - 仅当 mcp_status 中 compute_scheduler 为 true（或明确开启）时，才允许选 hpc。
   - 在此前提下，若用户自然语言涉及超算运维、作业状态、队列、节点资源、常见命令语义（如 pestat、squeue、作业 ID、排队原因等），选 hpc。

3) chat
   - 普通问候、感谢、闲聊，或明显是通用百科/教材级知识且不要求在本系统内执行分析或访问用户数据时，选 chat。

4) task
   - 用户要跑分析、流程、差异分析、比对、组装等，且 file_status 已表明有匹配类型的文件或数据已就绪时，选 task。

5) skill_fast_lane
   - 仅当用户意图极度明确属于某个已命名的「技能」快捷路径、且不需要 clarify 时选用；否则不要用。

输出前自检：若选了 hpc，务必确认 mcp_status.compute_scheduler 为真；若应 clarify 却选了 task，属于严重错误。
"""


def _serialize_router_context(inp: RouterInput) -> str:
    """将结构化上下文压成一段 JSON 字符串，便于模型一次性阅读。"""
    payload = {
        "query": inp.query,
        "file_status": inp.file_status,
        "mcp_status": inp.mcp_status,
        "session_flags": inp.session_flags,
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
