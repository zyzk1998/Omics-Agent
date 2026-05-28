# Prompt-Engineered 软技能规范（SKILL.md）

本目录为 `gibh_agent/skills/prompt_specs/<spec_id>/SKILL.md` 的文档镜像，供阅读与评审。

| spec_id | 广场名称 | tool_id |
|---------|----------|---------|
| `weekly_report_writer` | 周报撰写助手 | `weekly_report_writer` |
| `tailored_resume` | 定制化简历生成 | `tailored_resume` |
| `academic_poster_generator` | 学术会议海报生成 | `academic_poster_generator` |
| `academic_abstract_refiner` | 学术摘要精炼 | `academic_abstract_refiner` |
| `email_manager` | 邮件管理助手 | `email_manager` |
| `pdf_extractor` | PDF 内容提取助手 | `pdf_extractor` |
| `deep_research` | 深度调研 | `deep_research` |
| `blueprint_drafter` | 工程蓝图制图 | `blueprint_drafter` |
| （另见 `ppt_outline` / `mindmap_gen`） | PPT 大纲 / 思维导图 | 专用 `skill_*.py` 实现 |

同步到运行时目录：

```bash
cd /home/ubuntu/GIBH-AGENT-V2
./scripts/sync_prompt_skill_specs.sh
```

存量库写入：

```bash
PYTHONPATH=. python3 scripts/patch_prompt_soft_skills.py
```

**执行路径**：上述 10 项均在 **api-server** 进程内调用 LLM（`LAUNCH_ISOLATED_TOOL_IDS` 不含它们，不会误委托 `launch-skills`）。须配置 `.env` 中 `SILICONFLOW_API_KEY`（或项目默认 LLM 提供商密钥）。
