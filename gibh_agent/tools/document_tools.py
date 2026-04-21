"""
文档二次加工工具：由大模型生成 Markdown 幻灯片大纲与思维导图结构，写入结果目录供下载。
"""
from __future__ import annotations

import logging
import os
import uuid
from typing import Any, Dict, Optional

from gibh_agent.core.llm_client import LLMClientFactory
from gibh_agent.core.llm_cloud_providers import get_default_chat_model
from gibh_agent.core.tool_registry import registry

logger = logging.getLogger(__name__)

_RESULTS_SUBDIR = "documents"


def _strip_code_fence(content: str) -> str:
    s = (content or "").strip()
    if not s.startswith("```"):
        return s
    lines = s.split("\n")
    if len(lines) >= 2 and lines[0].strip().startswith("```"):
        lines = lines[1:]
    if lines and lines[-1].strip() == "```":
        lines = lines[:-1]
    return "\n".join(lines).strip()


def _llm_text(system: str, user: str, max_tokens: int = 4096) -> str:
    client = LLMClientFactory.create_for_model(get_default_chat_model())
    completion = client.chat(
        messages=[
            {"role": "system", "content": system},
            {"role": "user", "content": user},
        ],
        temperature=0.3,
        max_tokens=max_tokens,
    )
    choices = getattr(completion, "choices", None) or []
    if not choices:
        raise RuntimeError("LLM 返回空 choices")
    msg = choices[0].message
    text = (getattr(msg, "content", None) or "").strip()
    if not text:
        raise RuntimeError("LLM 正文为空")
    return _strip_code_fence(text)


def _save_md(content: str, prefix: str) -> str:
    results_dir = os.getenv("RESULTS_DIR", "/app/results").rstrip("/")
    out_dir = os.path.join(results_dir, _RESULTS_SUBDIR)
    os.makedirs(out_dir, exist_ok=True)
    fname = f"{prefix}_{uuid.uuid4().hex[:12]}.md"
    path = os.path.join(out_dir, fname)
    with open(path, "w", encoding="utf-8") as f:
        f.write(content)
    return os.path.abspath(path)


def _tool_success_payload(abs_path: str, md: str, short_label: str) -> Dict[str, Any]:
    preview = md.strip()
    if len(preview) > 2400:
        preview = preview[:2400] + "\n\n…"
    return {
        "status": "success",
        "message": f"{short_label}已保存，可通过下方链接下载。",
        "summary": preview,
        "markdown_file_path": abs_path,
        "report_data": {
            "download_links": [
                {
                    "path": abs_path,
                    "title": "📄 下载 Markdown 文档",
                    "kind": "markdown",
                }
            ],
        },
    }


@registry.register(
    name="generate_presentation_md",
    description=(
        "根据用户的分析/汇报需求，调用大模型生成科研汇报用 Markdown 幻灯片大纲；"
        "每页幻灯片之间使用单独一行的 --- 分隔。将内容写入 .md 文件并返回下载信息。"
    ),
    category="Document",
    output_type="mixed",
)
def generate_presentation_md(
    analysis_request: str,
    context: Optional[str] = None,
) -> Dict[str, Any]:
    """
    Args:
        analysis_request: 用户希望体现在幻灯片中的分析主题、结论与结构要求。
        context: 可选的补充上下文（如前序步骤摘要、关键指标）。
    """
    analysis_request = (analysis_request or "").strip()
    if not analysis_request:
        return {
            "status": "error",
            "error": "缺少 analysis_request",
            "message": "请说明需要生成汇报大纲的主题与要点。",
        }
    ctx = (context or "").strip()
    user = f"分析/汇报需求：\n{analysis_request}\n"
    if ctx:
        user += f"\n补充上下文：\n{ctx}\n"
    system = (
        "你是生物医学科研写作助手。输出**仅**为 Markdown 正文，不要前言后语。\n"
        "生成「科研汇报 PPT 大纲」：包含研究背景、方法要点（若上下文有则写）、核心发现、结论与展望等。\n"
        "每一页幻灯片用一级或二级标题标明页题；**页与页之间必须单独一行只写 ---**（三个连字符）。\n"
        "每页下列 3–6 条要点（无序列表即可），语言严谨、可放在学术汇报中使用。\n"
        "不要使用 HTML，不要 ``` 代码围栏包裹全文。"
    )
    try:
        md = _llm_text(system, user, max_tokens=6144)
    except Exception as e:
        logger.exception("generate_presentation_md LLM 失败")
        return {
            "status": "error",
            "error": str(e),
            "message": f"生成 PPT 大纲失败：{e}",
        }
    path = _save_md(md, "presentation")
    return _tool_success_payload(path, md, "PPT 大纲 Markdown ")


@registry.register(
    name="generate_mindmap_md",
    description=(
        "根据用户需求，调用大模型生成用于思维导图的层级 Markdown 列表（使用 - 与缩进表示层级）。"
        "写入 .md 文件并返回下载路径与摘要预览。"
    ),
    category="Document",
    output_type="mixed",
)
def generate_mindmap_md(
    analysis_request: str,
    context: Optional[str] = None,
) -> Dict[str, Any]:
    """
    Args:
        analysis_request: 希望导图涵盖的主题、逻辑链或关键词。
        context: 可选的补充材料（如分析结论、基因列表）。
    """
    analysis_request = (analysis_request or "").strip()
    if not analysis_request:
        return {
            "status": "error",
            "error": "缺少 analysis_request",
            "message": "请说明思维导图要整理的主题与范围。",
        }
    ctx = (context or "").strip()
    user = f"整理需求：\n{analysis_request}\n"
    if ctx:
        user += f"\n补充信息：\n{ctx}\n"
    system = (
        "你是生物医学科研助理。输出**仅**为 Markdown，不要解释性开场白。\n"
        "生成**思维导图用层级列表**：根节点用一行粗体标题 `# 主题`；其下全部使用 `- ` 无序列表，"
        "子层级用两个空格缩进一级（Markdown 标准嵌套列表）。\n"
        "条目简洁，体现因果与并列关系；可含「核心逻辑」「关键驱动基因/通路」「验证要点」等分支（按用户任务取舍）。\n"
        "不要使用 HTML，不要用 ``` 包裹全文。"
    )
    try:
        md = _llm_text(system, user, max_tokens=6144)
    except Exception as e:
        logger.exception("generate_mindmap_md LLM 失败")
        return {
            "status": "error",
            "error": str(e),
            "message": f"生成思维导图 Markdown 失败：{e}",
        }
    path = _save_md(md, "mindmap")
    return _tool_success_payload(path, md, "思维导图 Markdown ")
