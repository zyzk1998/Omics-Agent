# -*- coding: utf-8 -*-
"""科学语料数据加工 — 特色科研流程技能（HITL 挂起由 CorpusProcessingAgent 编排）。"""
from __future__ import annotations

from typing import Any, Dict

from gibh_agent.skills.base_skill import BaseSkill


class CorpusDataProcessingSkill(BaseSkill):
    __abstractskill__ = False

    skill_id = "skill_corpus_data_processing"
    display_name = "科学语料数据加工"
    description = (
        "支持上传科研图像或文本数据，通过内嵌 Label Studio 工作台进行可视化专家标注，"
        "最终一键导出符合大模型微调标准（SFT/RLHF）的高质量 JSON 语料。"
    )
    category = "特色科研流程"
    sub_category = "语料加工"
    aliases = ["语料打标", "SFT语料", "Label Studio语料", "corpus", "语料加工"]
    required_parameters = []
    __dependencies__: list[str] = []

    def execute(
        self,
        image_path: str = "",
        file_path: str = "",
        project_title: str = "",
        **kwargs: Any,
    ) -> Dict[str, Any]:
        """由 CorpusProcessingAgent.stream_skill_flow 承担主路径；此处供 ToolRegistry 与冒烟注册。"""
        from gibh_agent.agents.corpus_processing_agent import CorpusProcessingAgent
        from gibh_agent.core.prompt_manager import create_default_prompt_manager

        agent = CorpusProcessingAgent(llm_client=None, prompt_manager=create_default_prompt_manager())
        if not any([(image_path or "").strip(), (file_path or "").strip()]):
            return {
                "status": "error",
                "message": "请上传 image_path 或 file_path；推荐从技能广场发起以走完整 HITL 编排。",
            }
        return agent.trigger_corpus_hitl(
            image_path=image_path,
            file_path=file_path,
            project_title=project_title,
        )
