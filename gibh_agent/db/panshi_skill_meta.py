"""
磐石 ScienceOne ToolChain 官方文案对齐（第三批生物医药/化学技能）。

对照源：`gibh_agent/data/scienceone_toolchain_tools.json`
（由 `scripts/refresh_scienceone_toolchain_names.py` 从 ToolChain bundle 解析）

广场大类/小类仍以本仓库 `seed_skills` 为准（生物医药 / 化学 + 数据分析、预测与建模等）；
本模块仅覆盖 **display_name → long_description** 的官方描述，避免种子与技能库文案漂移。
"""
from __future__ import annotations

import json
from functools import lru_cache
from pathlib import Path
from typing import Any

_SCIENCEONE_JSON = Path(__file__).resolve().parents[1] / "data" / "scienceone_toolchain_tools.json"

# 第三批及同批置顶技能：按磐石 display_name 对齐描述（不改变 main/sub_category）
PANSHI_ALIGNED_DISPLAY_NAMES: frozenset[str] = frozenset(
    {
        "CRISPR-Cas9 基因编辑工具",
        "蛋白同源结构评估器",
        "人源性评估",
        "圆二色谱分析工具",
        "蛋白质序列保守性分析工具",
        "ITC结合热力学分析工具",
        "蛋白酶动力学分析工具",
        "酶动力学分析工具",
        "RNA二级结构分析工具",
        "基因集富集分析工具",
        "免疫细胞分离与纯化",
        "估计细胞周期各阶段持续时间",
        "合成可行性分析工具",
    }
)

# display_name → ScienceOne tool_chain key（文档/排障用，不入库）
PANSHI_TOOL_CHAIN_KEY_BY_DISPLAY_NAME: dict[str, str] = {
    "CRISPR-Cas9 基因编辑工具": "CRISPRCas9GeneEditing",
    "蛋白同源结构评估器": "ProteinAssess",
    "人源性评估": "Humanness",
    "圆二色谱分析工具": "CircularDichroismAnalysis",
    "蛋白质序列保守性分析工具": "ProteinSequenceConservationAnalysis",
    "ITC结合热力学分析工具": "ITCAnalysis",
    "蛋白酶动力学分析工具": "ProteaseKineticsAnalysis",
    "酶动力学分析工具": "EnzymeKineticsAnalysis",
    "RNA二级结构分析工具": "RNASecondaryStructureAnalysis",
    "基因集富集分析工具": "GeneSetEnrichmentAnalysis",
    "免疫细胞分离与纯化": "ImmuneCellIsolationAndPurification",
    "估计细胞周期各阶段持续时间": "CellCyclePhaseDurationEstimation",
    "合成可行性分析工具": "SyntheticFeasibilityAnalysis",
}


@lru_cache(maxsize=1)
def load_panshi_descriptions_by_display_name() -> dict[str, str]:
    """从 scienceone_toolchain_tools.json 读取官方 long_description。"""
    if not _SCIENCEONE_JSON.is_file():
        return {}
    try:
        payload = json.loads(_SCIENCEONE_JSON.read_text(encoding="utf-8"))
    except (json.JSONDecodeError, OSError):
        return {}
    tools = payload.get("tools") if isinstance(payload, dict) else payload
    if not isinstance(tools, list):
        return {}
    out: dict[str, str] = {}
    for item in tools:
        if not isinstance(item, dict):
            continue
        name = (item.get("display_name") or "").strip()
        desc = (item.get("long_description") or item.get("description") or "").strip()
        if name and desc:
            out[name] = desc
    return out


def apply_panshi_official_descriptions(skills: list[dict[str, Any]]) -> list[dict[str, Any]]:
    """对白名单技能用磐石 long_description 覆盖 description（原地修改并返回同一列表）。"""
    official = load_panshi_descriptions_by_display_name()
    if not official:
        return skills
    for skill in skills:
        name = (skill.get("name") or "").strip()
        if name not in PANSHI_ALIGNED_DISPLAY_NAMES:
            continue
        desc = official.get(name)
        if desc:
            skill["description"] = desc
        key = PANSHI_TOOL_CHAIN_KEY_BY_DISPLAY_NAME.get(name)
        if key:
            skill["scienceone_tool_chain_key"] = key
    return skills
