# -*- coding: utf-8 -*-
"""
HPC 文件扫描与清洗（三层过滤网）：
1) 白名单后缀速筛
2) 正则特征保留
3) LLM 兜底鉴定（deepseek-v4-pro + mode=standard）
"""
from __future__ import annotations

import json
import logging
import re
from typing import List

from gibh_agent.core.llm_client import (
    LLMClientFactory,
    reset_llm_api_mode,
    set_llm_api_mode,
)

logger = logging.getLogger(__name__)


class HpcFileScanner:
    _WHITELIST_SUFFIXES = (
        ".gz",
        ".tar",
        ".zip",
        ".fa",
        ".fasta",
        ".fq",
        ".fastq",
        ".bam",
        ".sam",
        ".cram",
        ".csv",
        ".tsv",
        ".txt",
        ".mtx",
        ".h5",
        ".h5ad",
        ".rds",
        ".rdata",
    )
    _FEATURE_PATTERNS = (
        re.compile(r"_R1_", re.IGNORECASE),
        re.compile(r"_R2_", re.IGNORECASE),
        re.compile(r"_counts", re.IGNORECASE),
        re.compile(r"_expression", re.IGNORECASE),
        re.compile(r"_aligned", re.IGNORECASE),
        re.compile(r"_filtered", re.IGNORECASE),
    )
    _DROP_PATTERNS = (
        re.compile(r"\.log$", re.IGNORECASE),
        re.compile(r"\.out$", re.IGNORECASE),
        re.compile(r"\.err$", re.IGNORECASE),
        re.compile(r"\.sh$", re.IGNORECASE),
        re.compile(r"\.py$", re.IGNORECASE),
        re.compile(r"\.r$", re.IGNORECASE),
        re.compile(r"(^|/)\.", re.IGNORECASE),
    )
    _LLM_UNKNOWN_LIMIT = 50

    async def filter_files(self, raw_file_list: List[str]) -> List[str]:
        normalized = self._normalize_file_list(raw_file_list)
        if not normalized:
            return []

        kept: List[str] = []
        unknown: List[str] = []
        for name in normalized:
            lower = name.lower()
            if self._matches_drop_patterns(name):
                continue
            if lower.endswith(self._WHITELIST_SUFFIXES):
                kept.append(name)
                continue
            if any(p.search(name) for p in self._FEATURE_PATTERNS):
                kept.append(name)
                continue
            unknown.append(name)

        if 0 < len(unknown) <= self._LLM_UNKNOWN_LIMIT:
            llm_kept = await self._llm_identify_files(unknown)
            if llm_kept:
                kept.extend(llm_kept)

        return self._dedupe_preserve_order(kept)

    @staticmethod
    def _normalize_file_list(raw_file_list: List[str]) -> List[str]:
        if not isinstance(raw_file_list, list):
            return []
        out: List[str] = []
        for item in raw_file_list:
            s = str(item or "").strip()
            if not s:
                continue
            if s.lower().startswith("total "):
                continue
            if s.endswith("/"):
                continue
            out.append(s)
        return HpcFileScanner._dedupe_preserve_order(out)

    @staticmethod
    def _dedupe_preserve_order(items: List[str]) -> List[str]:
        seen = set()
        out: List[str] = []
        for x in items:
            if x in seen:
                continue
            seen.add(x)
            out.append(x)
        return out

    def _matches_drop_patterns(self, file_name: str) -> bool:
        return any(p.search(file_name) for p in self._DROP_PATTERNS)

    async def _llm_identify_files(self, unknown_files: List[str]) -> List[str]:
        prompt = (
            "作为生信专家，请从以下文件列表中提取出可能是原始测序数据、组学矩阵或关键分析结果的文件名。"
            "严格忽略日志(log/out)、脚本(sh/py/R)、隐藏文件和临时缓存。"
            "请仅返回一个合法的 JSON 数组，包含你认为保留的文件名字符串，不要输出任何 markdown 标记或其他解释文本。\n\n"
            "文件列表：\n"
            + json.dumps(unknown_files, ensure_ascii=False)
        )
        token = set_llm_api_mode("standard")
        try:
            llm_client = LLMClientFactory.create_for_model("deepseek-v4-pro")
            completion = await llm_client.achat(
                messages=[
                    {"role": "system", "content": "你是严谨的生信文件清洗助手。"},
                    {"role": "user", "content": prompt},
                ],
                model="deepseek-v4-pro",
                temperature=0.0,
                max_tokens=512,
                stream=False,
            )
            content = self._extract_content(completion)
            parsed = self._parse_json_array(content)
            if not parsed:
                return []
            unknown_set = set(unknown_files)
            return [x for x in parsed if x in unknown_set and not self._matches_drop_patterns(x)]
        except Exception as e:  # noqa: BLE001
            logger.warning("[HpcFileScanner] LLM 兜底鉴定失败: %s", e, exc_info=True)
            return []
        finally:
            reset_llm_api_mode(token)

    @staticmethod
    def _extract_content(completion) -> str:
        try:
            choices = getattr(completion, "choices", None) or []
            if choices:
                msg = getattr(choices[0], "message", None)
                if msg is not None:
                    c = getattr(msg, "content", "")
                    if isinstance(c, str):
                        return c.strip()
        except Exception:  # noqa: BLE001
            pass
        return ""

    @staticmethod
    def _parse_json_array(text: str) -> List[str]:
        raw = (text or "").strip()
        if not raw:
            return []
        try:
            obj = json.loads(raw)
            if isinstance(obj, list):
                return [str(x).strip() for x in obj if str(x).strip()]
        except Exception:  # noqa: BLE001
            pass
        m = re.search(r"\[[\s\S]*\]", raw)
        if not m:
            return []
        try:
            obj = json.loads(m.group(0))
            if isinstance(obj, list):
                return [str(x).strip() for x in obj if str(x).strip()]
        except Exception:  # noqa: BLE001
            return []
        return []

