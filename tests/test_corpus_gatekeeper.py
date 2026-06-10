# -*- coding: utf-8 -*-
"""Corpus Gatekeeper：模态嗅探与入库前语料就绪。"""
from __future__ import annotations

import asyncio
import json
from pathlib import Path

from gibh_agent.core.corpus_gatekeeper import (
    CORPUS_MODALITY_NLP,
    CORPUS_MODALITY_VLM,
    corpus_archive_is_ready,
    sniff_corpus_modality,
)
from gibh_agent.agents.corpus_processing_agent import CorpusProcessingAgent
from gibh_agent.core.prompt_manager import create_default_prompt_manager


def test_sniff_modality_vlm_when_image_paths_present():
    assert sniff_corpus_modality(hitl_image_paths=["data:image/png;base64,abc"]) == CORPUS_MODALITY_VLM


def test_sniff_modality_nlp_when_no_images():
    assert (
        sniff_corpus_modality(
            steps_details=[{"step_name": "翻译", "status": "success"}],
            hitl_annotations=None,
        )
        == CORPUS_MODALITY_NLP
    )


def test_corpus_archive_is_ready_requires_dataset_json(tmp_path):
    root = tmp_path / "corpus_archive" / "s1"
    root.mkdir(parents=True)
    assert not corpus_archive_is_ready(root)
    (root / "dataset.json").write_text("[]", encoding="utf-8")
    assert corpus_archive_is_ready(root)


def test_generate_corpus_archive_nlp_no_base64_in_dataset(tmp_path, monkeypatch):
    monkeypatch.setenv("RESULTS_DIR", str(tmp_path))
    agent = CorpusProcessingAgent(llm_client=None, prompt_manager=create_default_prompt_manager())
    archive_dir = tmp_path / "corpus_archive" / "nlp-sess"

    async def _run():
        return await agent.generate_corpus_archive_bundle(
            summary_context={
                "expert_report_draft": "纯文本检索任务：基因 SYMBOL 映射与文献摘要。",
            },
            steps_results=[],
            steps_details=[{"step_name": "文献检索", "status": "success"}],
            archive_dir=archive_dir,
            expert_report_md="纯文本检索任务结论。",
        )

    bundle = asyncio.run(_run())
    assert bundle["status"] == "success"
    assert bundle["modality"] == CORPUS_MODALITY_NLP
    dataset = json.loads((archive_dir / "dataset.json").read_text(encoding="utf-8"))
    assert dataset[0]["instruction"]
    assert "image" not in dataset[0]
    raw = json.dumps(dataset)
    assert "base64" not in raw.lower()


def test_generate_corpus_archive_vlm_uses_image_relative_path(tmp_path, monkeypatch):
    monkeypatch.setenv("RESULTS_DIR", str(tmp_path))
    agent = CorpusProcessingAgent(llm_client=None, prompt_manager=create_default_prompt_manager())
    archive_dir = tmp_path / "corpus_archive" / "vlm-sess"

    async def _run():
        return await agent.generate_corpus_archive_bundle(
            summary_context={
                "expert_report_draft": "UMAP 显示 3 群细胞。",
                "hitl_image_paths": ["data:image/png;base64,xyz"],
            },
            steps_results=[],
            steps_details=[],
            archive_dir=archive_dir,
            expert_report_md="UMAP 报告",
        )

    bundle = asyncio.run(_run())
    assert bundle["status"] == "success"
    assert bundle["modality"] == CORPUS_MODALITY_VLM
    dataset = json.loads((archive_dir / "dataset.json").read_text(encoding="utf-8"))
    assert dataset[0]["image"].startswith("images/")
    assert not str(dataset[0]["image"]).startswith("data:")
    assert (archive_dir / dataset[0]["image"]).is_file()
