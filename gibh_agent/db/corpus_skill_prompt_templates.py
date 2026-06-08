# -*- coding: utf-8 -*-
"""
科学语料数据加工技能 prompt_template（与 seed_skills / patch_corpus_skill 同源）。
"""
from __future__ import annotations

CORPUS_SKILL_NAME = "科学语料数据加工"

CORPUS_SKILL_PROMPT_TEMPLATE = """[Skill_Route: skill_corpus_data_processing]
您好。我将引导您完成**科研语料专家标注**：上传图像或文本后，系统会打开内嵌 Label Studio 工作台，由您完成可视化打标；标注结束后自动清洗为 **SFT 微调 JSON**（`instruction` / `input` / `output`）。

**支持上传**
- **图像**：`.png` / `.jpg` / `.jpeg` / `.webp` / `.tif`（框选区域 + 语料类型）
- **文本**：`.txt` / `.md`（填写 instruction / input / output 或 QA 对）
- **任务清单**：`.json`（Label Studio 任务数组，每项含 `data` 字段）

**广场一键体验（助手侧：用户未上传文件时，须在工具调用中保留空路径并提示上传）**
```json
{"image_path": "", "file_path": "", "project_title": "科学语料标注 · 演示项目"}
```

**助手侧纪律**
1. 收到有效 `image_path` 或 `file_path` 后**立即**进入 Label Studio 挂起，勿再向用户追问无关参数。
2. 路径参数仅使用 `image_path` / `file_path`（架构宪法词汇表）。
3. 用户点击「已完成标注」唤醒后，输出 `/results/.../sft_corpus_*.json` 下载链接与预览。
"""

CORPUS_SKILL_SEED: dict[str, str] = {
    "name": CORPUS_SKILL_NAME,
    "main_category": "特色科研流程",
    "sub_category": "语料加工",
    "description": (
        "支持上传科研图像或文本数据，通过内嵌 Label Studio 工作台进行可视化专家标注，"
        "最终一键导出符合大模型微调标准（SFT/RLHF）的高质量 JSON 语料。"
    ),
}
