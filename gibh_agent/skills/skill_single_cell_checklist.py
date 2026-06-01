# -*- coding: utf-8 -*-
"""单细胞实验设计检查清单 — 纯 Prompt → Markdown Checklist。"""
from __future__ import annotations

from typing import Any, Dict

from gibh_agent.skills._prompt_skill_markdown import execute_markdown_prompt_skill
from gibh_agent.skills.base_skill import BaseSkill

_SYSTEM = """你是单细胞转录组（scRNA-seq）实验设计与生物信息质控顾问，熟悉 10x、BD、Parse 等平台。

【任务】根据用户研究问题，输出可打印、可勾选的实验设计检查清单（非数据分析执行）。

【Markdown 结构 — 必须严格按序输出】
# 单细胞实验设计检查清单

> 研究问题复述（1 句）+ 推荐细胞捕获量级区间（基于用户描述或 Demo 示意）

## A. 样本与分组设计
- [ ] 生物学重复 ≥3/组（说明理由）
- [ ] 对照与处理定义清晰
- [ ] 批次与处理正交或已记录批次计划
（每条 checklist 后附 **检查要点** 一行小字说明）

## B. 样本制备与质检
- [ ] 活率、结团、红细胞污染评估
- [ ] 消化方案与平台 SOP 匹配
- [ ] 上机浓度与目标细胞数记录

## C. 测序与数据交付
- [ ] 读长、深度、双端参数符合平台建议
- [ ] FASTQ 命名规范与 metadata 表字段完整

## D. 生物信息预设（分析前）
- [ ] 参考基因组/注释版本确定
- [ ] 双胞检测策略选定（Scrublet/DoubletFinder 等）
- [ ] 整合批次方法预案（Harmony/Seurat anchor 等）

## E. 统计与报告
- [ ] 主要终点（细胞比例 / 标志基因 / 轨迹）预注册
- [ ] 多重检验与伪重复风险说明

## F. 风险红旗
用 `- [ ]` 列出 3–5 条本项目最需警惕的失误点"""


class SingleCellChecklistSkill(BaseSkill):
    __abstractskill__ = False

    skill_id = "single_cell_checklist"
    display_name = "单细胞实验设计检查清单"
    description = (
        "从研究问题生成 scRNA-seq 实验设计可勾选 Checklist："
        "重复、批次、QC、双胞、整合与统计预设。"
    )
    category = "生物医药"
    sub_category = "文本处理"
    aliases = ["scRNA实验设计", "单细胞checklist", "10x实验设计", "单细胞质控清单"]
    required_parameters = [
        "研究问题、物种、测序平台与分组（user_request / context 文本）",
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
            deliver_message="单细胞实验设计清单已生成，请在右侧工作台查看",
        )
