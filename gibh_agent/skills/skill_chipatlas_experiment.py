# -*- coding: utf-8 -*-
"""ChIP-Atlas 实验检索 — 首发核心技能。"""
from __future__ import annotations

import re
from typing import Any, Dict, List, Optional

import httpx

from gibh_agent.skills._skill_common import DEFAULT_HTTP_TIMEOUT, err, ok
from gibh_agent.skills.base_skill import BaseSkill
from gibh_agent.skills.launch_skill_demos import apply_launch_demo_defaults

CHIPATLAS_BASE = "https://chip-atlas.org"
ACCESSION_RE = re.compile(r"^(SRX|GSM|SRP|ERP|DRX|ERX)\d+", re.I)


class ChipatlasExperimentSearchSkill(BaseSkill):
    __abstractskill__ = False

    """
    ChIP-Atlas 实验元数据检索。

    - 若 query 为 SRX/GSM 等 accession，直接拉取单条实验元数据。
    - 否则在实验列表元数据中进行关键词匹配（受 limit 与扫描行数约束）。

    参数:
        query: 检索关键词或实验 ID（如 SRX018625、CTCF K-562）。
        expid: 实验 ID（与 query 等价，优先使用）。
        limit: 最大返回条数，默认 20。
    """

    skill_id = "chipatlas_experiment_search"
    display_name = "ChIPAtlas实验获取工具"
    description = (
        "检索 ChIP-Atlas 实验元数据：支持 accession 精确查询或关键词匹配，"
        "返回 SRX/GEO、基因组、抗原与细胞类型等信息。"
    )
    category = "生物医药"
    sub_category = "数据分析"
    aliases = ["ChIP-Atlas", "ChIPAtlas", "chipatlas", "表观遗传实验"]
    required_parameters = []
    tool_chain_key = "chipatlas"
    __dependencies__ = ["pip:httpx"]

    def _get_experiment(self, expid: str) -> Dict[str, Any]:
        url = f"{CHIPATLAS_BASE}/data/exp_metadata.json"
        with httpx.Client(timeout=DEFAULT_HTTP_TIMEOUT, follow_redirects=True) as client:
            resp = client.get(url, params={"expid": expid})
            resp.raise_for_status()
            return resp.json()

    def _search_tab_metadata(self, keyword: str, limit: int) -> List[Dict[str, str]]:
        """在 experimentList.tab 中流式关键词匹配（有扫描上限）。"""
        url = f"{CHIPATLAS_BASE}/data/metadata/experimentList.tab"
        kw = keyword.lower()
        hits: List[Dict[str, str]] = []
        max_scan_lines = int(
            __import__("os").environ.get("CHIPATLAS_SEARCH_MAX_LINES", "80000")
        )

        with httpx.Client(timeout=120.0, follow_redirects=True) as client:
            with client.stream("GET", url) as resp:
                resp.raise_for_status()
                header: Optional[List[str]] = None
                for i, raw_line in enumerate(resp.iter_lines()):
                    if i > max_scan_lines:
                        break
                    line = raw_line.decode("utf-8", errors="replace") if isinstance(raw_line, bytes) else raw_line
                    if not line or line.startswith("#"):
                        continue
                    parts = line.split("\t")
                    if header is None:
                        header = [h.strip() for h in parts]
                        continue
                    row = {header[j]: parts[j] for j in range(min(len(header), len(parts)))}
                    blob = "\t".join(parts).lower()
                    if kw in blob:
                        hits.append(row)
                        if len(hits) >= limit:
                            break
        return hits

    def execute(
        self,
        query: str = "",
        expid: str = "",
        limit: int = 20,
        **kwargs: Any,
    ) -> Dict[str, Any]:
        filled = apply_launch_demo_defaults(
            self.skill_id,
            {
                "expid": (expid or "").strip(),
                "query": (query or kwargs.get("experiment_id") or "").strip(),
                "limit": limit,
            },
        )
        expid = str(filled.get("expid") or "").strip()
        query = str(filled.get("query") or "").strip()
        q = expid or query
        if not q:
            return err("请提供 query 或 expid（关键词或 SRX/GSM accession）")

        lim = max(1, min(int(filled.get("limit") or 20), 100))

        if ACCESSION_RE.match(q):
            try:
                meta = self._get_experiment(q)
            except httpx.HTTPError as exc:
                return err(f"ChIP-Atlas 元数据获取失败: {exc}")
            return ok(
                f"已获取实验 {q} 的元数据",
                expid=q,
                experiment=meta,
                count=1,
            )

        try:
            hits = self._search_tab_metadata(q, lim)
        except httpx.HTTPError as exc:
            return err(f"ChIP-Atlas 实验列表检索失败: {exc}")

        if not hits:
            return ok(
                f"未在扫描范围内匹配到「{q}」相关实验，可尝试更具体的 accession（如 SRX…）",
                query=q,
                hits=[],
                count=0,
            )

        return ok(
            f"ChIP-Atlas 关键词检索完成，共 {len(hits)} 条",
            query=q,
            hits=hits,
            count=len(hits),
        )
