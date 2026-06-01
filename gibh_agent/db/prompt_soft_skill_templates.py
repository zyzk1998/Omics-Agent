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
{"user_request": "绘制 RNA-seq 差异分析流程蓝图：原始 FASTQ → 质控 → 表达矩阵 → 差异分析 → 可视化", "context": "中文节点标签；仅描述用户业务数据流"}
```

（助手侧：必填 `user_request` 描述要画的业务流程，禁止写入平台内部架构代号。）
""",
    "差异表达结果解读助手": """[Skill_Route: diff_expr_interpreter]
您好。请将我提供的**差异表达（DEG）结果**解读为组会/论文可用的生物学叙事报告（Markdown + 基因表）。

**广场一键体验（助手侧：用户未改写时，须将下列 JSON 原样写入工具调用）**
```json
{"user_request": "解读下列 Demo DEG 结果并给出验证建议", "context": "对比：处理组 vs 对照组；阈值 padj<0.05, |log2FC|>1\\nTP53\\t2.1\\t1.2e-8\\nCDKN1A\\t1.8\\t3.4e-6\\nMYC\\t-1.5\\t0.002\\nIL6\\t1.2\\t0.01\\nSTAT3\\t0.9\\t0.03"}
```

（助手侧：完整 DEG 表或火山图摘要放 `context`；对比设计写入 `user_request`。）
""",
    "GO/KEGG 富集结果叙事生成": """[Skill_Route: go_kegg_narrative]
您好。请将**GO/KEGG 富集结果表**整理为机制—疾病—验证方向的连贯叙事（Markdown 分节）。

**广场一键体验（助手侧：须将下列 JSON 原样写入工具调用）**
```json
{"user_request": "撰写富集结果叙事，突出免疫与代谢主题", "context": "GO:0006955 immune response\\tpadj=2.1e-5\\tGeneRatio=12/180\\nKEGG:hsa04668 TNF signaling\\tpadj=8.3e-4\\tGeneRatio=9/180\\nGO:0006091 energy metabolism\\tpadj=0.012\\tGeneRatio=7/180"}
```

（助手侧：富集表粘贴至 `context`；研究背景一句写入 `user_request`。）
""",
    "单细胞实验设计检查清单": """[Skill_Route: single_cell_checklist]
您好。请为我的课题生成**单细胞实验设计可勾选 Checklist**（样本、QC、双胞、整合、统计）。

**广场一键体验（助手侧：须将下列 JSON 原样写入工具调用）**
```json
{"user_request": "人肺腺癌免疫治疗前后 scRNA-seq；比较响应者与非响应者；计划 10x 3' v3", "context": "预期每组 6 例生物学重复；关注 T 细胞耗竭与髓系抑制细胞"}
```

（助手侧：平台、物种、分组写入 `user_request`；样本量与关注点写入 `context`。）
""",
    "期刊 Cover Letter 草稿": """[Skill_Route: journal_cover_letter]
您好。请根据下列信息起草**期刊 Cover Letter**（Markdown，可直接润色投稿）。

**广场一键体验（助手侧：须将下列 JSON 原样写入工具调用）**
```json
{"user_request": "目标期刊：Nature Communications；Article type: Research Article；语言：英文", "context": "Title: Single-cell atlas reveals T cell exhaustion dynamics under PD-1 blockade\\nHighlights: (1) New exhaustion trajectory (2) Cross-cohort validation (3) Proposed biomarker panel"}
```

（助手侧：期刊与体裁放 `user_request`；题目与亮点放 `context`。）
""",
    "审稿意见 Rebuttal 要点提纲": """[Skill_Route: review_rebuttal_outline]
您好。请将下列**审稿意见**拆解为分条 Rebuttal 回复策略与补充实验优先级（编号列表 Markdown）。

**广场一键体验（助手侧：须将下列 JSON 原样写入工具调用）**
```json
{"user_request": "稿件：scRNA-seq 免疫治疗研究；拟回复 Major revision", "context": "Reviewer 1: (1) 缺少独立队列验证 (2) 双胞阈值未说明 (3) 部分标志基因未做 IHC\\nReviewer 2: 统计方法需澄清伪重复"}
```

（助手侧：审稿原文放 `context`；稿件题目与修订类型放 `user_request`。）
""",
    "变异临床意义解读草案（ACMG 风格）": """[Skill_Route: acmg_variant_interpretation]
您好。请根据下列变异信息生成 **ACMG 风格证据条目草稿**（须标注需数据库核实）。

**广场一键体验（助手侧：须将下列 JSON 原样写入工具调用）**
```json
{"user_request": "BRCA1 胚系变异临床意义草稿；GRCh38", "context": "c.5266dupC (p.Gln1756Profs*74); gnomAD 频率未提供; 家族乳腺癌史; 预测：蛋白截断"}
```

（助手侧：HGVS、人群频率、预测工具结果放 `context`。）
""",
    "生物信息 Pipeline 选型备忘录": """[Skill_Route: pipeline_selection_memo]
您好。请为下列目标撰写 **Pipeline 选型对比备忘录**（不执行计算）。

**广场一键体验（助手侧：须将下列 JSON 原样写入工具调用）**
```json
{"user_request": "人肿瘤组织 bulk RNA-seq；有 20 对肿瘤/癌旁；需差异表达与通路分析", "context": "无 WGS；计算资源：32 核 128G；偏好开源工具"}
```

（助手侧：数据类型、样本量、资源约束写入 `user_request` / `context`。）
""",
    "组学样本 Metadata 规范审查": """[Skill_Route: omics_metadata_review]
您好。请审查下列样本 metadata 表的规范性与混杂风险。

**广场一键体验（助手侧：须将下列 JSON 原样写入工具调用）**
```json
{"user_request": "RNA-seq 肿瘤/癌旁对照研究 metadata 审查", "context": "SampleID\\tGroup\\tBatch\\tRIN\\nS1\\tTumor\\tB1\\t8.2\\nS2\\tNormal\\tB1\\t\\nS3\\tTumor\\tB2\\t7.5"}
```

（助手侧：完整样本表粘贴至 `context`。）
""",
    "引物/qPCR 实验设计说明": """[Skill_Route: qpcr_primer_design_guide]
您好。请生成 **引物 / qPCR 实验设计说明**（原则、对照、循环参数）。

**广场一键体验（助手侧：须将下列 JSON 原样写入工具调用）**
```json
{"user_request": "人源 IL6 mRNA 相对定量；SYBR Green；组织 RNA", "context": "内参候选：GAPDH、ACTB；无现成引物序列"}
```

（助手侧：目标基因、样本类型、已有引物序列放 `context`。）
""",
    "会议口头报告讲稿大纲": """[Skill_Route: oral_presentation_outline]
您好。请将下列内容整理为 **10–15 分钟口头报告分页提纲**（JSON，含建议口述时长）。

**广场一键体验（助手侧：须将下列 JSON 原样写入工具调用）**
```json
{"user_request": "组会口头报告 12 分钟；听众：实验室研究生", "context": "题目：PD-1 阻断下 T 细胞耗竭轨迹；结构：背景→方法→scRNA 结果→验证→展望"}
```

（助手侧：时长与听众放 `user_request`；论文结构放 `context`。）
""",
    "临床注册方案骨架生成（CONSORT 导向）": """[Skill_Route: clinical_trial_protocol_skeleton]
您好。请生成 **临床注册方案骨架**（CONSORT 导向 Markdown）。

**广场一键体验（助手侧：须将下列 JSON 原样写入工具调用）**
```json
{"user_request": "II 期单臂免疫联合化疗；晚期 NSCLC；主要终点 ORR", "context": "人群：EGFR 野生型；次要终点 PFS、安全性；需随机化说明框架"}
```

（助手侧：设计类型与终点放 `user_request`。）
""",
    "药物重定位假说备忘录": """[Skill_Route: drug_repositioning_memo]
您好。请撰写 **药物重定位假说备忘录** 与文献检索关键词。

**广场一键体验（助手侧：须将下列 JSON 原样写入工具调用）**
```json
{"user_request": "二甲双胍用于 IPF 抗纤维化重定位假说", "context": "已知：AMPK 激活；动物模型部分有效；靶点 TGF-β/Smad 通路推测"}
```

（助手侧：疾病—靶点—药物证据放 `context`。）
""",
    "蛋白质功能假说推演": """[Skill_Route: protein_function_hypothesis]
您好。请基于下列蛋白注释推演 **功能假说与验证路线**（含 Mermaid）。

**广场一键体验（助手侧：须将下列 JSON 原样写入工具调用）**
```json
{"user_request": "TP53 功能假说推演；肿瘤抑制语境", "context": "结构域：DNA 结合域；修饰：乙酰化；定位：核；已知：细胞周期检查点"}
```

（助手侧：仅解读给定注释，勿编造未提供的实验结果。）
""",
    "实验失败复盘报告": """[Skill_Route: experiment_failure_postmortem]
您好。请整理下列实验失败为 **组会复盘报告**。

**广场一键体验（助手侧：须将下列 JSON 原样写入工具调用）**
```json
{"user_request": "Western Blot 目标条带弱且无特异性", "context": "一抗 1:1000 4°C 过夜；样本：冻融 3 次组织；预期 55 kDa；实际多条非特异"}
```

（助手侧：现象与已尝试步骤放 `context`。）
""",
    "文献矩阵精读笔记": """[Skill_Route: literature_matrix_notes]
您好。请将下列摘要整理为 **文献矩阵精读笔记**。

**广场一键体验（助手侧：须将下列 JSON 原样写入工具调用）**
```json
{"user_request": "免疫治疗预测生物标志物；纳入 3 篇代表文献", "context": "文献A：scRNA 队列 n=24；文献B：bulk 转录组 meta；文献C：空间转录组试点"}
```

（助手侧：摘要原文放 `context`。）
""",
    "多组学整合分析故事线": """[Skill_Route: multi_omics_storyline]
您好。请撰写 **多组学整合故事线** 与 Figure 编排建议（含 Mermaid，不运算）。

**广场一键体验（助手侧：须将下列 JSON 原样写入工具调用）**
```json
{"user_request": "肺癌免疫治疗耐药机制；基因组+转录组+蛋白组", "context": "WGS：TP53 突变；RNA：耗竭 T 细胞模块；蛋白组：STAT3 磷酸化升高"}
```

（助手侧：各层发现摘要放 `context`。）
""",
    "伦理与知情同意要点清单": """[Skill_Route: ethics_consent_checklist]
您好。请生成 **伦理审查 Q&A 与知情同意书章节 Checklist**。

**广场一键体验（助手侧：须将下列 JSON 原样写入工具调用）**
```json
{"user_request": "前瞻性观察性生物样本库；可识别健康人；仅抽血", "context": "存储：5 年；可能基因测序；无干预"}
```

（助手侧：干预风险与人群特征放 `user_request` / `context`。）
""",
    "统计方法选择建议书": """[Skill_Route: stats_method_advisor]
您好。请根据下列研究设计给出 **统计方法选择建议书**。

**广场一键体验（助手侧：须将下列 JSON 原样写入工具调用）**
```json
{"user_request": "两组独立肿瘤 vs 正常；连续结局（表达量）；非正态", "context": "n=18 vs 20；可能有批次；需多重比较校正"}
```

（助手侧：设计、样本量、结局类型放 `user_request` / `context`。）
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
    {
        "name": "差异表达结果解读助手",
        "main_category": "生物医药",
        "sub_category": "文本处理",
        "description": (
            "将 DEG 表或摘要转为生物学叙事：Top 基因表、机制假说、质控提醒与验证建议（Markdown）。"
        ),
    },
    {
        "name": "GO/KEGG 富集结果叙事生成",
        "main_category": "生物医药",
        "sub_category": "文本处理",
        "description": (
            "富集表→分节机制叙事：主题归类、疾病关联、枢纽基因与后续验证建议。"
        ),
    },
    {
        "name": "单细胞实验设计检查清单",
        "main_category": "生物医药",
        "sub_category": "文本处理",
        "description": (
            "按课题生成 scRNA-seq 可勾选实验设计清单：重复、批次、QC、双胞与整合预案。"
        ),
    },
    {
        "name": "期刊 Cover Letter 草稿",
        "main_category": "生物医药",
        "sub_category": "文本处理",
        "description": (
            "根据题目、亮点与目标期刊生成 Cover Letter Markdown 草稿及投稿策略备忘。"
        ),
    },
    {
        "name": "审稿意见 Rebuttal 要点提纲",
        "main_category": "生物医药",
        "sub_category": "文本处理",
        "description": (
            "审稿意见→分条回复策略、文稿修改位置、补充实验 P0/P1/P2 与时间线提纲。"
        ),
    },
    {
        "name": "变异临床意义解读草案（ACMG 风格）",
        "main_category": "生物医药",
        "sub_category": "文本处理",
        "description": (
            "变异位点与人群频率文字→ACMG 证据条目草稿；须标注需 ClinVar/gnomAD 核实。"
        ),
    },
    {
        "name": "生物信息 Pipeline 选型备忘录",
        "main_category": "生物医药",
        "sub_category": "文本处理",
        "description": (
            "RNA-seq/WGS/蛋白组等目标→2–3 套工具链对比表与推荐（不执行计算）。"
        ),
    },
    {
        "name": "组学样本 Metadata 规范审查",
        "main_category": "生物医药",
        "sub_category": "文本处理",
        "description": (
            "审查样本 metadata 表：缺失字段、分组/批次混杂与修订优先级。"
        ),
    },
    {
        "name": "引物/qPCR 实验设计说明",
        "main_category": "生物医药",
        "sub_category": "文本处理",
        "description": (
            "引物设计原则、对照体系、反应参数与验证实验清单。"
        ),
    },
    {
        "name": "会议口头报告讲稿大纲",
        "main_category": "生物医药",
        "sub_category": "文本处理",
        "description": (
            "10–15 分钟口述分页提纲（JSON：标题、要点、建议时长），复用 PPT 大纲组件。"
        ),
    },
    {
        "name": "临床注册方案骨架生成（CONSORT 导向）",
        "main_category": "生物医药",
        "sub_category": "文本处理",
        "description": (
            "CONSORT 导向的方案骨架：目的、入排、终点与统计假设框架。"
        ),
    },
    {
        "name": "药物重定位假说备忘录",
        "main_category": "生物医药",
        "sub_category": "文本处理",
        "description": (
            "疾病—靶点—药物→可验证假说与文献检索关键词（不调用外部库）。"
        ),
    },
    {
        "name": "蛋白质功能假说推演",
        "main_category": "生物医药",
        "sub_category": "文本处理",
        "description": (
            "基于用户给定结构域/修饰/定位推演功能假说、验证路线与 Mermaid 示意。"
        ),
    },
    {
        "name": "实验失败复盘报告",
        "main_category": "生物医药",
        "sub_category": "文本处理",
        "description": (
            "失败现象→原因树、排查步骤与组会 takeaway（Markdown 复盘稿）。"
        ),
    },
    {
        "name": "文献矩阵精读笔记",
        "main_category": "生物医药",
        "sub_category": "文本处理",
        "description": (
            "3–10 篇摘要→对比矩阵表（设计、样本、结论、局限）与跨研究启示。"
        ),
    },
    {
        "name": "多组学整合分析故事线",
        "main_category": "生物医药",
        "sub_category": "文本处理",
        "description": (
            "多组学层级逻辑、Figure 编排建议与 Discussion 叙事（含 Mermaid）。"
        ),
    },
    {
        "name": "伦理与知情同意要点清单",
        "main_category": "生物医药",
        "sub_category": "文本处理",
        "description": (
            "伦理审查常见 Q&A 与知情同意书章节可勾选 Checklist。"
        ),
    },
    {
        "name": "统计方法选择建议书",
        "main_category": "生物医药",
        "sub_category": "文本处理",
        "description": (
            "按研究设计推荐检验方法、前提假设、误区与报告规范提醒。"
        ),
    },
]
