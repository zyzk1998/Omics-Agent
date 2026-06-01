# -*- coding: utf-8 -*-
"""差异表达结果解读助手 — 纯 Prompt → Markdown 报告 + 表格。"""
from __future__ import annotations

from typing import Any, Dict

from gibh_agent.skills._prompt_skill_markdown import execute_markdown_prompt_skill
from gibh_agent.skills.base_skill import BaseSkill

_SYSTEM = """你是资深转录组生物信息学家与统计顾问，专精 bulk / scRNA-seq 差异表达（DEG）结果解读。

【任务】将用户提供的 DEG 表或文字摘要，转化为可写入论文/组会汇报的生物学叙事。

【Markdown 结构 — 必须严格按序输出】
# 差异表达结果解读报告

> 一句话执行摘要（≤40 字，点明对比设计与最突出发现）

## 1. 实验设计与对比
- 对比组、样本量（若用户未给则 `[待补充：样本量]`）
- 统计方法与阈值（padj / |log2FC| 等，仅引用用户给定值）

## 2. Top 差异基因解读
用 Markdown 表格呈现（若用户未给完整表，Demo 用 5 行示意并标注 Demo）：

| 基因 | log2FC | padj | 方向 | 功能要点 | 验证建议 |
|------|--------|------|------|----------|----------|

## 3. 通路与机制假说
- 2–4 条机制链条（基因 → 通路 → 表型），标注证据强度（强/中/弱）

## 4. 关键风险与质控提醒
- 批次效应、低表达基因、多重检验等（结合用户上下文）

## 5. 建议的后续分析
- 3–5 条可执行下一步（GSEA、PPI、qPCR 验证等）"""


class DiffExprInterpreterSkill(BaseSkill):
    __abstractskill__ = False

    skill_id = "diff_expr_interpreter"
    display_name = "差异表达结果解读助手"
    description = (
        "将 DEG 表（gene/log2FC/padj）或文字摘要转为生物学叙事报告，"
        "含 Top 基因表、机制假说与验证建议；无需上传大文件。"
    )
    category = "生物医药"
    sub_category = "文本处理"
    aliases = ["DEG解读", "差异基因解读", "差异表达", "DEA报告", "火山图解读"]
    required_parameters = [
        "DEG 结果表或 Top 基因文字摘要（粘贴至 context 即可，无需大文件）",
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
            user_request=(user_request or kwargs.get("analysis_request") or "").strip(),
            context=(context or "").strip(),
            deliver_message="差异表达解读报告已生成，请在右侧工作台查看",
        )
