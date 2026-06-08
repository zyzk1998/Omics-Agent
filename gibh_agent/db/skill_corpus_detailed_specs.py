# -*- coding: utf-8 -*-
"""科学语料数据加工 · detailed_spec 与 demo_visualization（原图 → SFT JSON）。"""
from __future__ import annotations

from typing import Any, Dict, List

_CORPUS_IMG = "/assets/images/demos/corpus/test_corpus_umap.png"

DEMO_VIZ_CORPUS_PROCESSING = f"""
<div class="skill-viz-corpus-flow" style="font-family:system-ui,-apple-system,sans-serif;">
  <div style="font-size:13px;font-weight:600;color:#334155;margin-bottom:12px;">原图 → 专家标注 → SFT 语料库</div>
  <div style="display:grid;grid-template-columns:1fr auto 1fr;gap:12px;align-items:start;">
    <div style="background:#f8fafc;border:1px solid #e2e8f0;border-radius:10px;padding:10px;">
      <div style="font-size:11px;font-weight:600;color:#64748b;margin-bottom:8px;">📥 输入 · 未注释 UMAP</div>
      <img src="{_CORPUS_IMG}" alt="未注释单细胞 UMAP 聚类图" style="width:100%;border-radius:6px;display:block;" loading="lazy" />
      <div style="font-size:10px;color:#94a3b8;margin-top:6px;">test_corpus_umap.png · 4 个未命名 Cluster</div>
    </div>
    <div style="display:flex;flex-direction:column;align-items:center;justify-content:center;padding-top:48px;color:#6366f1;font-size:22px;">→</div>
    <div style="background:#f0fdf4;border:1px solid #bbf7d0;border-radius:10px;padding:10px;">
      <div style="font-size:11px;font-weight:600;color:#166534;margin-bottom:8px;">📤 输出 · SFT JSON 语料</div>
      <pre style="margin:0;font-size:9.5px;line-height:1.45;background:#fff;border:1px solid #d1fae5;border-radius:6px;padding:8px;overflow:auto;max-height:220px;color:#14532d;"><code>[
  {{
    "instruction": "请识别并标注该单细胞 UMAP 聚类图中的核心细胞类群。",
    "input": "/uploads/test_corpus_umap.png",
    "output": "检测到 4 个主要类群：左上区域 (x≈120, y≈80) 为 T-Cell；右下区域为 B-Cell；中部偏左为 Macrophage；右上为 Unknown 过渡态。"
  }}
]</code></pre>
      <div style="font-size:10px;color:#16a34a;margin-top:6px;">唤醒后自动导出 /results/corpus_hitl/sft_corpus_*.json</div>
    </div>
  </div>
  <div style="margin-top:12px;padding:10px 12px;background:#eff6ff;border-left:3px solid #3b82f6;border-radius:0 8px 8px 0;font-size:11px;color:#1e40af;">
    <strong>🛠️ 处理链路</strong>：上传图像 → 程序化拉起 Label Studio → 专家框选 T-Cell / B-Cell 等区域 → 一键清洗为 instruction / input / output 三元组
  </div>
</div>
"""

_SPEC_CORPUS: Dict[str, Any] = {
    "tool_id": "skill_corpus_data_processing",
    "description_long": (
        "【🔬 核心能力】\n"
        "将杂乱的原始科研图像/文本，通过专家可视化复核，一键转化为可直接用于训练大模型的标准微调语料。\n\n"
        "【✨ 效果预览】\n"
        "📥 **输入**：一张未注释的单细胞聚类图 (.png) 或病理切片 / 文本语料。\n"
        "🛠️ **处理**：自动拉起 Label Studio 工作台，专家在图上框选 T-Cell、B-Cell、Macrophage 等区域并填写 SFT 字段。\n"
        "📤 **输出**：标准 SFT JSON 数组（instruction / input / output），可直接接入 LoRA / 全参微调流水线。\n\n"
        "本技能采用 **硬 HITL** 挂起：收到有效文件后无需复杂 Prompt，系统自动创建 LS 标注项目；"
        "专家完成标注并点击「继续生成报告」后，由 LLM 清洗导出语料文件至工作区。"
    ),
    "usage_hint": (
        "请上传需要进行打标加工的图像或数据文件（点击下方附件图标）。"
        "无需输入复杂指令，发送后系统将自动为您拉起 Label Studio 标注工作台。"
    ),
    "demo_visualization": DEMO_VIZ_CORPUS_PROCESSING,
    "inputs": [
        {
            "name": "image_path",
            "type": "image (.png/.jpg/.webp/.tif)",
            "required": False,
            "description": "待标注科研图像；与 file_path 二选一。",
        },
        {
            "name": "file_path",
            "type": "file (.txt/.md/.json)",
            "required": False,
            "description": "文本语料或 Label Studio 任务 JSON 清单。",
        },
    ],
    "outputs": [
        {
            "name": "sft_corpus_json",
            "type": "JSON",
            "description": "SFT 微调语料数组，每条含 instruction / input / output。",
        },
        {
            "name": "hitl_workspace",
            "type": "Label Studio 项目",
            "description": "内嵌 LS 工作台链接，供专家可视化打标。",
        },
    ],
    "parameters_table": [
        {
            "name": "image_path",
            "type": "path",
            "required": False,
            "description": "单张待标注图像（/uploads/ 或 /results/ 相对路径）。",
        },
        {
            "name": "file_path",
            "type": "path",
            "required": False,
            "description": "文本或任务 JSON 文件路径。",
        },
        {
            "name": "project_title",
            "type": "string",
            "required": False,
            "description": "Label Studio 项目名称（默认「科学语料专家标注」）。",
        },
    ],
    "query_examples": [
        {
            "tier": "zero_shot",
            "label": "一键打标（推荐）",
            "hint": "上传 test_corpus_umap.png 或任意科研图像后直接发送，无需额外指令。",
            "prompt": (
                "[Skill_Route: skill_corpus_data_processing]\n"
                "请上传需要进行打标加工的图像或数据文件（点击下方附件图标）。"
                "无需输入复杂指令，发送后系统将自动为您拉起 Label Studio 标注工作台。"
            ),
        },
        {
            "tier": "dynamic_routing",
            "label": "指定项目名称",
            "hint": "可选自定义 LS 项目标题，便于多批次语料区分。",
            "prompt": (
                '[Skill_Route: skill_corpus_data_processing]\n'
                '```json\n'
                '{"image_path": "/uploads/test_corpus_umap.png", "project_title": "单细胞 UMAP 语料批次-01"}\n'
                '```'
            ),
        },
    ],
    "workflow_highlights": [
        "Label Studio",
        "硬 HITL 挂起",
        "SFT 语料导出",
        "instruction/input/output",
    ],
}

CORPUS_SPECS_BY_TOOL_ID: Dict[str, Dict[str, Any]] = {
    "skill_corpus_data_processing": _SPEC_CORPUS,
}

SKILL_NAME_TO_CORPUS_TOOL_ID: Dict[str, str] = {
    "科学语料数据加工": "skill_corpus_data_processing",
}
