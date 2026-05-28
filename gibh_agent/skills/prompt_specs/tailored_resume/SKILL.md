# Tailored Resume Generator

## When to Use

- Applying for a specific job position
- Customizing resume for different industries or roles
- Career transitions; ATS optimization
- Multiple resume versions for different applications

## What This Skill Does

1. **Analyzes Job Descriptions** — requirements, skills, keywords
2. **Identifies Priorities** — what employers value most
3. **Tailors Content** — reorganizes experience and achievements
4. **Optimizes Keywords** — ATS-friendly natural keywords
5. **Formats Professionally** — clean Markdown resume layout
6. **Provides Recommendations** — gaps, interview hooks, cover letter hooks

## Workflow (Agent Must Follow)

### 1. Gather Information（70 分边界）

自评完整度：**<40%** 才追问（≤2 问）；**≥70%** 直接成稿，次要信息用 `[待补充：…]`。**无输入** → 输出 Demo 简历（Senior 数据/生信方向范例）。

理想输入（非全部必需）：

- Full job description, company, job title
- Candidate background: work history, education, skills, metrics, certifications

### 2. Analyze Job Requirements

Extract must-haves, key skills, soft skills, industry knowledge, ATS keywords. Map Priority 1/2/3.

### 3. Map Experience to Requirements

Match or transferable skills; note gaps; highlight unique strengths.

### 4. Structure Resume (Markdown)

- **Professional Summary** (3–4 lines): years, top skills, industry, value proposition
- **Technical/Core Skills**: grouped by JD categories; exact JD terminology
- **Professional Experience**: action verbs, quantified bullets; most relevant first
- **Education** + optional Certifications/Projects

### 5. ATS Optimization

Standard headings; exact keywords; no tables/graphics in body; acronyms + full terms.

### 6. After Resume — Strategic Recommendations

Strengths analysis, gap analysis, interview prep tips, cover letter hooks.

## Rules

**Do:** truthful, quantify, concise, tailor each time.  
**Don't:** fabricate, personal info (photo/age), first-person pronouns, generic template.

## Output Format

Deliver complete **Markdown resume** first, then optional sections:

- `<center><h1>姓名</h1></center>` + 居中联系方式
- 技能：**Markdown 表格**横向对比（技能 | 熟练度 | 代表项）
- 经历：**加粗** `**公司 · 职位 · 年月**` + 量化 bullets
- 可选：`## 优势分析` / `## 差距与建议` / `## 面试与求职信提示`
