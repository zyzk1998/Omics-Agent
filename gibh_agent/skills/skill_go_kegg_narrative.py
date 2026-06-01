# -*- coding: utf-8 -*-
"""GO/KEGG 富集结果叙事生成 — 纯 Prompt → Markdown 分节。"""
from __future__ import annotations

from typing import Any, Dict

from gibh_agent.skills._prompt_skill_markdown import execute_markdown_prompt_skill
from gibh_agent.skills.base_skill import BaseSkill

_SYSTEM = """你是通路富集分析与系统生物学叙事专家，熟悉 GO、KEGG、Reactome 术语体系。

【任务】基于用户粘贴的富集结果表或条目列表，撰写「机制—疾病—可验证假设」连贯叙事。

【Markdown 结构 — 必须严格按序输出】
# GO/KEGG 富集结果叙事

> 核心结论一句话（指出最显著的 1–2 条通路主题）

## 1. 富集结果概览
| 通路/术语 ID | 名称 | pvalue/FDR | GeneRatio | 命中基因数 | 主题归类 |
|--------------|------|------------|-----------|------------|----------|

## 2. 生物学主题解读
### 2.1 免疫与炎症（若无相关条目写「未显著富集」）
### 2.2 代谢与能量
### 2.3 细胞周期与增殖
（按用户数据实际存在的主题增减小节，每节 2–4 句机制叙述）

## 3. 疾病与临床关联
- 与哪些疾病表型/通路失调假说一致（标注「需文献核实」若超出输入）

## 4. 关键基因桥接
- 列出 3–6 个枢纽基因及其在通路中的角色（仅来自用户输入）

## 5. 后续验证建议
- 实验或分析跟进（Western、IHC、单基因敲低、GSEA 复核等）"""


class GoKeggNarrativeSkill(BaseSkill):
    __abstractskill__ = False

    skill_id = "go_kegg_narrative"
    display_name = "GO/KEGG 富集结果叙事生成"
    description = (
        "将 GO/KEGG 富集表转为分节机制叙事：主题归类、疾病关联、枢纽基因与验证建议；"
        "不编造未给出的 p 值。"
    )
    category = "生物医药"
    sub_category = "文本处理"
    aliases = ["富集解读", "通路叙事", "GO解读", "KEGG解读", "GSEA叙事"]
    required_parameters = [
        "GO/KEGG 富集表或通路条目（粘贴至 context）",
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
            deliver_message="GO/KEGG 富集叙事已生成，请在右侧工作台查看",
        )
