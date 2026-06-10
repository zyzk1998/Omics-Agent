# -*- coding: utf-8 -*-
"""科学语料数据加工 · 用户可见说明与画板操作指南（前后端同源文案）。"""
from __future__ import annotations

CORPUS_SKILL_DESCRIPTION_SHORT = (
    "【🔬 技能用途】\n"
    "可视化人工标注工具！上传未注释的聚类图、病理图等数据，系统将自动拉起专业标注画板。\n"
    "您框选的结果将一键导出为 SFT 高质量微调语料 JSON。"
)

CORPUS_SKILL_USER_GUIDE = """【🎨 画板操作指南】
发送图片后，您将看到一个深色的专业标注工作台：

👁️ 定位图片：界面正中央即为您上传的待标注图像。

🏷️ 选择标签：在界面右侧的 Labels 面板，点击选中您想标注的类型（如 T-Cell）。

✏️ 框选目标：鼠标移至中央图片上，按住左键并拖动，画出矩形框圈出细胞群。

💾 提交结果：全部画完后，点击界面上方或右侧醒目的蓝色 [Submit (提交)] 按钮。

🔙 返回系统：最后，点击我们系统页面下方的 [✅ 我已标注完成] 按钮，系统将自动回收您的标注数据！"""

CORPUS_SKILL_PROMPT_TEMPLATE = f"""[Skill_Route: skill_corpus_data_processing]
{CORPUS_SKILL_DESCRIPTION_SHORT}

{CORPUS_SKILL_USER_GUIDE}

**支持上传**：图像 .png/.jpg/.webp/.tif · 文本 .txt/.md · 任务清单 .json

```json
{{"image_path": "/assets/images/demos/corpus/test_corpus_umap.png", "project_title": "科学语料标注 · UMAP 演示"}}
```
"""
