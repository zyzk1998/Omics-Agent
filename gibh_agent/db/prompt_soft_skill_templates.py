# -*- coding: utf-8 -*-
"""
Prompt-Engineered 软技能 prompt_template（与 seed_skills / 广场 Featured 同源）。

首行 [Skill_Route: <tool_id>] 须与 gibh_agent/skills/skill_*.py 中 skill_id 完全一致。
SKILL 正文规范见 gibh_agent/skills/prompt_specs/<spec_id>/SKILL.md
"""
from __future__ import annotations

PROMPT_SOFT_SKILL_TEMPLATES: dict[str, str] = {
    "科研汇报 PPT 大纲生成": """[Skill_Route: ppt_outline]
您好。我需要将下列研究内容整理为**科研汇报 PPT 大纲**（结构化 JSON，含页数、每页标题与要点）。

**广场一键体验（助手侧：用户未改写主题时，须将下列 JSON 原样写入工具调用）**
```json
{"user_request": "单细胞时空动力学分析：研究背景、最优传输轨迹推断、驱动基因与通路富集、核心结论与展望", "context": ""}
```

（助手侧：必填 `user_request`；有前序分析摘要时写入 `context`。）
""",
    "分析逻辑思维导图生成": """[Skill_Route: mindmap_gen]
您好。请将下列分析逻辑整理为**Mermaid 思维导图**（仅逻辑树，不写 Markdown 列表）。

**广场一键体验（助手侧：用户未改写主题时，须将下列 JSON 原样写入工具调用）**
```json
{"user_request": "单细胞时空动力学：数据校验→时序标准化→轨迹推断→驱动基因→通路富集→专家解读", "context": ""}
```

（助手侧：必填 `user_request`；补充材料写入 `context`。）
""",
    "周报撰写助手": """[Skill_Route: weekly_report_writer]
您好。请为我起草本周**周报**（个人状态同步 + 团队对外同步）。

**广场一键体验（助手侧：须写入下列 JSON）**
```json
{"user_request": "周期：2026-05-19 至 2026-05-25；笔记库已挂载；请按 SKILL 工作流合成周报，待办继承上一份报告中的 🟥/🟨 项", "context": "上一份报告路径：reports/2026-05-12-weekly.md；重点关注：项目 A 里程碑、实验 B 失败复盘"}
```

（助手侧：必填 `user_request`；日期、vault、上一份报告路径写入 `context` 或 `user_request`。）
""",
    "定制化简历生成": """[Skill_Route: tailored_resume]
您好。请根据职位描述为我生成**定制化 Markdown 简历**（含优势/差距建议）。

**广场一键体验（助手侧：须写入下列 JSON）**
```json
{"user_request": "职位：Senior Data Analyst；要求 SQL/Python/可视化/A/B 测试；偏好医疗行业", "context": "候选人：RetailCo 数据分析师 5 年；Tableau/Power BI；HealthPlus 实习 1 年；Business Analytics 学士"}
```

（助手侧：`user_request` 放 JD；`context` 放履历要点。）
""",
    "学术会议海报生成": """[Skill_Route: academic_poster_generator]
您好。请根据论文材料生成**学术海报 Markdown 故事板**（Intro/Methods/Results/Conclusions + ≥3 配图方案）。

**广场一键体验（助手侧：须写入下列 JSON）**
```json
{"user_request": "主题：单细胞轨迹推断方法；需要会议海报故事板与配图说明", "context": "摘要：基于最优传输的细胞命运推断；结果：UMAP 轨迹与驱动基因列表（见附件）"}
```

（助手侧：有 PDF 时 `file_path` 指向论文；否则摘要写入 `context`。）
""",
    "学术摘要精炼": """[Skill_Route: academic_abstract_refiner]
您好。请将下列文稿精炼为**中英双语单段摘要**（Summary_Report 格式）。

**广场一键体验（助手侧：须写入下列 JSON）**
```json
{"user_request": "将下方研究综述精炼为 SCI 风格无小标题摘要", "context": "（在此粘贴待精炼的长文摘要或结论段落，勿编造未出现的数据）"}
```

（助手侧：源文放 `context`；禁止虚构实验数据。）
""",
    "邮件管理助手": """[Skill_Route: email_manager]
您好。请按邮件管理 SKILL 为我**撰写邮件文稿**（含主题行，2–3 个语气选项若适用）。

**广场一键体验（助手侧：须写入下列 JSON）**
```json
{"user_request": "write email to client about project delay — professional tone, include milestones and revised date"}
```

（助手侧：不访问邮箱；仅输出可复制 Markdown。）
""",
    "PDF 内容提取助手": """[Skill_Route: pdf_extractor]
您好。请从 PDF **提取并总结**文本结构（章节大纲、表格摘要、后续建议命令）。

**广场一键体验（助手侧：须写入下列 JSON）**
```json
{"user_request": "提取全文结构并总结关键表格与结论", "file_path": ""}
```

（助手侧：有上传 PDF 时将绝对路径写入 `file_path`；无文件时请用户在 `context` 粘贴摘录。）
""",
    "深度调研": """[Skill_Route: deep_research]
您好。请对下列主题做**深度调研**（Executive Summary + 引用 + 共识/争议）。

**广场一键体验（助手侧：须写入下列 JSON）**
```json
{"user_request": "间歇性禁食对代谢健康的影响：益处、风险与适用人群", "context": "可引用用户已提供的论文摘要或笔记；勿编造未给出的文献"}
```

（助手侧：`user_request` 为研究问题；`context` 为已有材料。）
""",
    "工程蓝图制图": """[Skill_Route: blueprint_drafter]
您好。请生成**扁平工程蓝图 HTML**（完整 <!DOCTYPE html>，系统字体，无 CDN）。

**广场一键体验（助手侧：须写入下列 JSON）**
```json
{"user_request": "绘制 GIBH 智能体技能快车道架构图：用户 → 技能广场 → SkillAgent → ToolRegistry → 工具/LLM → 右栏可视化", "context": "包含 api-server、launch-skills 委托分支；中文标签"}
```

（助手侧：必填 `user_request` 描述要画的系统/流程。）
""",
}

PROMPT_SOFT_SKILL_SEEDS: list[dict[str, str]] = [
    {
        "name": "科研汇报 PPT 大纲生成",
        "main_category": "其他技能",
        "sub_category": "文本处理",
        "description": (
            "将研究主题与结论整理为分页 PPT 大纲（JSON：总页数、每页标题与核心要点），"
            "便于导入演示文稿或组会汇报。"
        ),
    },
    {
        "name": "分析逻辑思维导图生成",
        "main_category": "其他技能",
        "sub_category": "文本处理",
        "description": (
            "从分析逻辑中提取层级关系，生成 Mermaid 思维导图并在工作台可视化展示。"
        ),
    },
    {
        "name": "周报撰写助手",
        "main_category": "其他技能",
        "sub_category": "文本处理",
        "description": (
            "按日期范围与笔记库日志合成周报：个人状态同步 + 团队对外同步，"
            "支持待办继承与 Wiki 链接溯源（Obsidian 友好）。"
        ),
    },
    {
        "name": "定制化简历生成",
        "main_category": "其他技能",
        "sub_category": "文本处理",
        "description": "根据 JD 与背景生成 ATS 友好 Markdown 简历及求职策略建议。",
    },
    {
        "name": "学术会议海报生成",
        "main_category": "其他技能",
        "sub_category": "文本处理",
        "description": "从论文内容生成会议海报分节要点、配图方案与质量检查清单（Markdown 故事板）。",
    },
    {
        "name": "学术摘要精炼",
        "main_category": "其他技能",
        "sub_category": "文本处理",
        "description": "长文→中英双语单段 SCI 风格摘要，输出 Summary_Report Markdown。",
    },
    {
        "name": "邮件管理助手",
        "main_category": "其他技能",
        "sub_category": "文本处理",
        "description": "撰写/润色/跟进/冷邮件等多场景邮件文稿；不连接邮箱，仅生成文本。",
    },
    {
        "name": "PDF 内容提取助手",
        "main_category": "其他技能",
        "sub_category": "文本处理",
        "description": "PDF 文本与结构提取摘要（可选 file_path + pypdf）；输出 Markdown 大纲。",
    },
    {
        "name": "深度调研",
        "main_category": "其他技能",
        "sub_category": "信息检索",
        "description": (
            "多视角深度调研报告：执行摘要、关键发现、分主题分析、共识/争议与编号引用来源。"
        ),
    },
    {
        "name": "工程蓝图制图",
        "main_category": "其他技能",
        "sub_category": "数据可视化",
        "description": (
            "生成扁平高数据墨水比工程/架构示意图 HTML，系统字体、无外链 CDN，工作台 iframe 预览。"
        ),
    },
]
