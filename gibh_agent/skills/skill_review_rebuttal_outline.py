# -*- coding: utf-8 -*-
"""审稿意见 Rebuttal 要点提纲 — 纯 Prompt → 编号列表 Markdown。"""
from __future__ import annotations

from typing import Any, Dict

from gibh_agent.skills._prompt_skill_markdown import execute_markdown_prompt_skill
from gibh_agent.skills.base_skill import BaseSkill

_SYSTEM = """你是生物医学论文返修（Revision）与审稿回复（Rebuttal）战略顾问，熟悉 Major/Minor revision 应答结构。

【任务】将审稿人意见逐条拆解为可执行的回复策略与补充实验优先级，输出编号列表体例。

【Markdown 结构 — 必须严格按序输出】
# 审稿意见 Rebuttal 要点提纲

> 稿件题目 + 修订版本建议（1 句）

## 总体策略
1. 用 3–5 条编号概括本次修订的叙事主线（态度谦逊、证据导向）

## 分条回复提纲
（每条审稿意见一个三级块，**必须**使用有序编号）

### 意见 1：[复述审稿人原意 — 用户未给则用 Demo 示意]
1. **态度与理解**：一句话承认关切
2. **回复要点**：2–4 条编号子点（数据/文献/统计）
3. **文稿修改位置**：Methods / Results / Discussion 建议
4. **补充实验优先级**：P0 / P1 / P2 + 预估周期 `[待补充：周期]`
5. **若无法完全满足**：替代论证策略（诚实边界）

### 意见 2：…
（按用户提供的意见条数扩展；Demo 至少 3 条意见）

## 修订后图表清单
1. …
2. …

## 时间线与分工
1. …
2. …"""


class ReviewRebuttalOutlineSkill(BaseSkill):
    __abstractskill__ = False

    skill_id = "review_rebuttal_outline"
    display_name = "审稿意见 Rebuttal 要点提纲"
    description = (
        "将审稿意见映射为分条回复策略、文稿修改位置、补充实验优先级（P0/P1/P2）"
        "与时间线；输出编号列表 Markdown。"
    )
    category = "生物医药"
    sub_category = "文本处理"
    aliases = ["rebuttal", "审稿回复", "大修回复", "revision letter", "审稿意见"]
    required_parameters = [
        "审稿意见原文或要点列表（粘贴至 context）",
    ]
    output_type = "markdown"
    __dependencies__: list[str] = []

    def execute(
        self,
        user_request: str = "",
        context: str = "",
        **kwargs: Any,
    ) -> Dict[str, Any]:
        return execute_markdown_prompt_skill(
            skill_id=self.skill_id,
            system_prompt=_SYSTEM,
            user_request=(user_request or "").strip(),
            context=(context or "").strip(),
            deliver_message="Rebuttal 要点提纲已生成，请在右侧工作台查看",
        )
