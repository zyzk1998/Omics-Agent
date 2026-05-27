# -*- coding: utf-8 -*-
"""IEDB MHC 关联表位检索 — 首发核心技能。"""
from __future__ import annotations

from typing import Any, Dict, List, Optional

import httpx

from gibh_agent.skills._skill_common import DEFAULT_HTTP_TIMEOUT, err, ok
from gibh_agent.skills.base_skill import BaseSkill
from gibh_agent.skills.launch_skill_demos import apply_launch_demo_defaults

IEDB_MHC_SEARCH_URL = "https://query-api.iedb.org/mhc_search"


class MhcEpitopeSearchSkill(BaseSkill):
    __abstractskill__ = False

    """
    MHC 关联表位检索（IEDB Query API）。

    在 IEDB mhc_search 表中检索与 MHC 限制性、肽段序列相关的实验记录，
    用于发现表位 structure_id 及免疫学上下文。

    参数:
        mhc_allele: MHC 等位基因关键词（如 HLA-A*02:01），模糊匹配 mhc_restriction。
        peptide_sequence: 精确肽段序列（可选）。
        limit: 返回行数，默认 20，最大 100。
        select: PostgREST select 字段列表，默认精简字段以加速响应。
    """

    skill_id = "mhc_epitope_search"
    display_name = "MHC关联表位检索"
    description = (
        "基于 IEDB Query API 检索 MHC 相关表位实验记录，"
        "支持按 MHC 等位基因或肽段序列筛选。"
    )
    category = "生物医药"
    sub_category = "数据分析"
    aliases = ["IEDB", "MHC表位", "mhc", "免疫表位"]
    required_parameters = []
    tool_chain_key = "mhc"
    __dependencies__ = ["pip:httpx"]

    def execute(
        self,
        mhc_allele: str = "",
        peptide_sequence: str = "",
        limit: int = 20,
        select: str = "",
        **kwargs: Any,
    ) -> Dict[str, Any]:
        filled = apply_launch_demo_defaults(
            self.skill_id,
            {
                "mhc_allele": (mhc_allele or kwargs.get("allele") or "").strip(),
                "peptide_sequence": (peptide_sequence or kwargs.get("sequence") or "").strip(),
                "limit": limit,
            },
        )
        allele = str(filled.get("mhc_allele") or "").strip()
        peptide = str(filled.get("peptide_sequence") or "").strip()
        if not allele and not peptide:
            return err("请至少提供 mhc_allele 或 peptide_sequence 之一")

        lim = max(1, min(int(filled.get("limit") or 20), 100))
        fields = (select or "structure_id,linear_sequence,mhc_restriction,source_organism_name,qualitative_measure").strip()

        params: Dict[str, Any] = {
            "limit": lim,
            "select": fields,
        }
        if allele:
            params["mhc_restriction"] = f"ilike.*{allele}*"
        if peptide:
            params["linear_sequence"] = f"eq.{peptide}"

        headers = {"Accept": "application/json"}
        try:
            with httpx.Client(timeout=DEFAULT_HTTP_TIMEOUT, follow_redirects=True) as client:
                resp = client.get(IEDB_MHC_SEARCH_URL, params=params, headers=headers)
                resp.raise_for_status()
                rows = resp.json()
        except httpx.HTTPError as exc:
            return err(f"IEDB API 请求失败: {exc}")

        if not isinstance(rows, list):
            return err("IEDB 返回格式异常", raw_type=type(rows).__name__)

        return ok(
            f"IEDB MHC 检索完成，共 {len(rows)} 条",
            filters={"mhc_allele": allele or None, "peptide_sequence": peptide or None},
            hits=rows,
            count=len(rows),
        )
