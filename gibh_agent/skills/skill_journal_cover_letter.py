# -*- coding: utf-8 -*-
"""期刊 Cover Letter 草稿 — 纯 Prompt → Markdown。"""
from __future__ import annotations

from typing import Any, Dict

from gibh_agent.skills._prompt_skill_markdown import execute_markdown_prompt_skill
from gibh_agent.skills.base_skill import BaseSkill

_SYSTEM = """你是生物医学 SCI 期刊投稿顾问，熟悉 Nature 子刊、Cell Press、Frontiers、MDPI 等 Cover Letter 体裁差异。

【任务】根据论文题目、亮点、目标期刊与作者信息（用户给定或 Demo），生成可直接润色后提交的投稿信 Markdown。

【Markdown 结构 — 必须严格按序输出】
# Cover Letter

**Date:** [待补充：日期]  
**Journal:** [期刊全称]  
**Manuscript title:** [论文题目]  
**Article type:** [Research Article / Short communication 等]

---

Dear Editor-in-Chief,

（正文 3–5 段，每段 3–5 句，商务正式英文；若用户要求中文则全文中文）

1. **研究重要性与领域空白**（Why now / Why this journal）
2. **核心发现与创新点**（量化结果仅引用用户给定数据）
3. **方法学可信度**（队列、模型、验证层级）
4. **与期刊 scope 的匹配**（具体栏目/读者群）
5. **声明段落**：无利益冲突模板句 + 建议审稿人可选段（若用户未要求则省略）

Sincerely,  
[待补充：通讯作者姓名]  
[待补充：单位与邮箱]

---

## 投稿策略备注（编辑内部备忘，非信件正文）
- 3 条标题优化建议
- 2 条可能被审稿人质疑的点及预备回应角度"""


class JournalCoverLetterSkill(BaseSkill):
    __abstractskill__ = False

    skill_id = "journal_cover_letter"
    display_name = "期刊 Cover Letter 草稿"
    description = (
        "根据题目、亮点与目标期刊生成 Cover Letter Markdown 草稿（中英可选），"
        "含投稿策略备忘；不连接投稿系统。"
    )
    category = "生物医药"
    sub_category = "文本处理"
    aliases = ["投稿信", "cover letter", "期刊投稿", "编辑信", "投稿附信"]
    required_parameters = [
        "论文题目、核心亮点与目标期刊名称",
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
            deliver_message="Cover Letter 草稿已生成，请在右侧工作台查看",
        )
