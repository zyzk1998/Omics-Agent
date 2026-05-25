#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
生成「生物信息分析工具」磐石 ToolChain ↔ 技能广场对照表 Markdown。

数据源（优先级）：
  1. gibh_agent/data/scienceone_skill_plaza_catalog.json（用户粘贴的技能广场全量，见 parse_scienceone_skill_plaza_paste.py）
  2. gibh_agent/data/scienceone_toolchain_tools.json（ToolChain bundle 约 120 条）
  - gibh_agent/db/seed_skills.py（技能广场种子）
  - gibh_agent/db/panshi_skill_meta.py（已对齐技能的 tool_chain_key）

用法（仓库根目录）:
    PYTHONPATH=. python3 scripts/parse_scienceone_skill_plaza_paste.py
    PYTHONPATH=. python3 scripts/generate_bioinformatics_tools_matrix_md.py
    PYTHONPATH=. python3 scripts/generate_bioinformatics_tools_matrix_md.py --csv docs/生物信息分析工具.csv
    PYTHONPATH=. python3 scripts/generate_bioinformatics_tools_matrix_md.py --include-implemented --csv docs/磐石技能广场_生物医药化学_全量.csv
"""
from __future__ import annotations

import argparse
import csv
import json
import re
import sys
from collections import Counter
from pathlib import Path
from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional, Set

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

TOOLS_JSON = ROOT / "gibh_agent/data/scienceone_toolchain_tools.json"
CATALOG_JSON = ROOT / "gibh_agent/data/scienceone_skill_plaza_catalog.json"

_DATA_SOURCE: str = "bundle"

EXCLUDE_TOOL_CHAIN_KEYS = frozenset(
    {
        "DeepResearch",
        "Writing",
        "AppSet",
        "SmartPPT",
        "DrawIO",
        "DataFormulator",
        "BRIARMBG",
        "Mermaid",
    }
)

OFFICE_DISP_RE = re.compile(r"写作|PPT|流程图|办公|翻译|润色|会议", re.I)

# skills类型：数据库调用 / 纯prompt / 脚本 / 模型
SKILLS_TYPE_DB_RE = re.compile(
    r"Database\b|DB\b|Browser\b|query_|ReactomeDB|UniProt|PubMed|GEODatabase|"
    r"GWASCatalog|dbSNPDatabase|InterProDB|ClinVarDatabase|STRINGDatabase|"
    r"GeneSearch|GSEASupported|AlphaFoldDB|structureMaterial|chemicalElementQuery|"
    r"UCSCGenomeBrowser|pubMedDatabase|信息检索|数据库查询|数据库条目",
    re.I,
)
SKILLS_TYPE_MODEL_RE = re.compile(
    r"AlphaFold|ESMFold|ESM3|ESM-Variants|ESMVariants|Evo2|evo2|DiffDock|DiffSBDD|"
    r"ProtGPT|RFdiffusion|BioGPT|LigandMPNN|OligoFormer|GenMol|BepiPred|ADMET|"
    r"MatterGen|GraphCast|ClimateGPT|ShapeE|LeanDojo|DSAgent|Tora|GenMol|"
    r"基础模型|深度学习|扩散模型|语言模型|蛋白语言模型|生成模型|神经网络|"
    r"Transformer|预训练模型",
    re.I,
)
ROUTE_RE = re.compile(r"\[(?:Skill_Route|Omics_Route):\s*(\w+)\s*\]", re.I)

_REGISTRY_BY_NAME: Optional[Dict[str, Dict[str, Any]]] = None


def _norm(s: str) -> str:
    return re.sub(r"\s+", "", (s or "").strip().lower())


def _escape_md_cell(text: str) -> str:
    return (text or "").replace("|", "\\|").replace("\n", " ").strip()


def load_panshi_tools(*, prefer_catalog: bool = True) -> List[Dict[str, Any]]:
    """加载磐石工具列表；默认优先技能广场用户粘贴全量目录。"""
    global _DATA_SOURCE
    if prefer_catalog and CATALOG_JSON.is_file():
        try:
            payload = json.loads(CATALOG_JSON.read_text(encoding="utf-8"))
            entries = payload.get("entries") if isinstance(payload, dict) else None
            if isinstance(entries, list) and entries:
                tools: List[Dict[str, Any]] = []
                for e in entries:
                    if not isinstance(e, dict):
                        continue
                    tools.append(
                        {
                            "tool_chain_key": (e.get("tool_chain_key") or "").strip(),
                            "display_name": (e.get("display_name") or "").strip(),
                            "long_description": (e.get("long_description") or "").strip(),
                            "source_block": e.get("source_block") or "skill_plaza_user_paste",
                            "main_category": (e.get("main_category") or "").strip(),
                            "skills_type": (e.get("skills_type") or "").strip(),
                            "invoke_status": (e.get("invoke_status") or "").strip(),
                        }
                    )
                _DATA_SOURCE = "skill_plaza_catalog"
                return tools
        except (json.JSONDecodeError, OSError):
            pass
    if not TOOLS_JSON.is_file():
        raise SystemExit(
            f"缺少 {CATALOG_JSON} 或 {TOOLS_JSON}；请先运行 parse_scienceone_skill_plaza_paste.py "
            "或 refresh_scienceone_toolchain_names.py"
        )
    payload = json.loads(TOOLS_JSON.read_text(encoding="utf-8"))
    tools = payload.get("tools") if isinstance(payload, dict) else []
    _DATA_SOURCE = "toolchain_bundle"
    return [t for t in tools if isinstance(t, dict)]


def classify_discipline(tool: Dict[str, Any]) -> str:
    mc = (tool.get("main_category") or "").strip()
    if mc == "化学":
        return "化学与分子信息学"
    if mc == "生物医药":
        return "生物医药"
    key = tool.get("tool_chain_key") or ""
    disp = tool.get("display_name") or ""
    desc = tool.get("long_description") or ""
    blob = f"{key} {disp} {desc}"
    if tool.get("source_block") == "chemistryTools" or re.search(
        r"molecular|Mol|RDKit|OpenBabel|drug|chem|Babel|lipinski|tanimoto|SMILES",
        blob,
        re.I,
    ):
        return "化学与分子信息学"
    if re.search(
        r"组学|基因组|转录|蛋白组|代谢|表观|空间|影像|Radiomics|Visium|FASTQ|BAM|VCF",
        blob,
        re.I,
    ):
        return "多模态组学"
    if re.search(
        r"蛋白|抗体|免疫|酶|RNA|DNA|基因|细胞|CRISPR|表位|Fold|Dock|序列",
        blob,
        re.I,
    ):
        return "生物医药"
    return "生物信息学基础"


def is_bioinformatics_tool(tool: Dict[str, Any]) -> bool:
    if tool.get("source_block") == "skill_plaza_user_paste":
        return (tool.get("main_category") or "") in ("生物医药", "化学")
    key = (tool.get("tool_chain_key") or "").strip()
    if key in EXCLUDE_TOOL_CHAIN_KEYS:
        return False
    blob = "".join(
        [
            key,
            tool.get("display_name") or "",
            tool.get("long_description") or "",
            tool.get("science_tool_dispatch_case") or "",
        ]
    )
    if OFFICE_DISP_RE.search(blob):
        return False
    if tool.get("source_block") == "chemistryTools":
        return True
    if (tool.get("science_tool_dispatch_case") or "").strip():
        return True
    disc = classify_discipline(tool)
    if disc != "生物信息学基础":
        return True
    return bool(
        re.search(
            r"蛋白|基因|RNA|DNA|分子|药物|化学|细胞|免疫|代谢|组学|表位|抗体|酶|激酶|"
            r"Protein|Gene|RNA|DNA|Mol|Drug|Cell|Pathway|Genome|Bepi|Epitope|Antibod|"
            r"CRISPR|Cas9|Fold|Dock|GSEA|Humanness|ITC|Kinase|UniProt|PubMed|GEO|GWAS|"
            r"dbSNP|Reactome|InterPro|AlphaFold|ESM|Oligo|ProtGPT|RNASecondary",
            blob,
            re.I,
        )
    )


def _load_registry_by_name() -> Dict[str, Dict[str, Any]]:
    global _REGISTRY_BY_NAME
    if _REGISTRY_BY_NAME is not None:
        return _REGISTRY_BY_NAME
    import gibh_agent.tools  # noqa: F401
    from gibh_agent.core.tool_registry import registry

    _REGISTRY_BY_NAME = {
        t["name"]: t for t in registry.get_all_tools_json() if t.get("name")
    }
    return _REGISTRY_BY_NAME


def infer_skills_type(
    *,
    tool: Dict[str, Any],
    seed: Optional[Dict[str, Any]],
    route_tool: str,
    skill_name: str,
    sub_category: str,
    intro: str,
) -> str:
    """
    推断技能实现形态（与技能广场 / Registry 对齐）：
    - 数据库调用：公共库/API 检索，无本地重计算管线
    - 纯prompt：仅 Prompt 模板或占位，无 [Skill_Route]/[Omics_Route] 可执行链
    - 脚本：Registry 算子 / Worker 子进程 / RDKit 等可执行脚本链
    - 模型：需调用深度学习基础模型或远端推理服务（AlphaFold/ESM/Evo2 等）
    """
    pt = (seed.get("prompt_template") or "") if seed else ""
    key = (tool.get("tool_chain_key") or "").strip()
    blob = " ".join(
        [
            key,
            skill_name,
            sub_category,
            intro,
            tool.get("display_name") or "",
            tool.get("long_description") or "",
            route_tool,
        ]
    )

    has_route = bool(route_tool)
    has_omics = bool(pt and "[Omics_Route:" in pt)

    if (sub_category or "").strip() == "信息检索":
        return "数据库调用"
    if SKILLS_TYPE_DB_RE.search(blob):
        return "数据库调用"

    reg = _load_registry_by_name().get(route_tool) if route_tool else None
    if reg and (reg.get("category") or "") == "Public API":
        return "数据库调用"
    if route_tool.startswith("query_") or route_tool.startswith("mcp_ncbi"):
        return "数据库调用"

    if SKILLS_TYPE_MODEL_RE.search(blob):
        return "模型"
    if reg:
        desc = (reg.get("description") or "") + (reg.get("category") or "")
        if SKILLS_TYPE_MODEL_RE.search(desc):
            return "模型"

    if has_route or has_omics:
        return "脚本"

    if tool.get("source_block") in ("chemistryTools", "materialTools"):
        return "脚本"

    if seed:
        try:
            from gibh_agent.core.skill_plaza_utils import is_prompt_placeholder_only

            if is_prompt_placeholder_only(pt):
                return "纯prompt"
        except Exception:
            pass

    if not has_route and not has_omics:
        if not seed or not pt.strip():
            return "纯prompt"
        if not ROUTE_RE.search(pt):
            return "纯prompt"

    return "纯prompt"


def infer_sub_category(tool: Dict[str, Any], discipline: str, seed: Optional[Dict[str, Any]]) -> str:
    if seed and (seed.get("sub_category") or "").strip():
        return str(seed["sub_category"]).strip()
    if tool.get("source_block") == "chemistryTools":
        return "化学信息学"
    if discipline == "多模态组学":
        return "流程编排"
    disp = (tool.get("display_name") or "") + (tool.get("tool_chain_key") or "")
    if re.search(r"Fold|Dock|GPT|Diff|预测|生成|Variants", disp, re.I):
        return "预测与建模"
    return "数据分析"


@dataclass
class ImplementedIndex:
    """本仓库技能广场已对接（[Skill_Route]/[Omics_Route] + Registry 可解析）的指纹集合。"""

    seed_names_norm: Set[str] = field(default_factory=set)
    panshi_keys_norm: Set[str] = field(default_factory=set)
    route_tools: Set[str] = field(default_factory=set)


def build_implemented_index() -> ImplementedIndex:
    """与 skill_plaza_utils / seed_skills / Registry 对齐，供表格删减已实现条目。"""
    from gibh_agent.core.skill_plaza_utils import infer_skill_implemented_from_prompt
    from gibh_agent.db.panshi_skill_meta import PANSHI_TOOL_CHAIN_KEY_BY_DISPLAY_NAME
    from gibh_agent.db.seed_skills import get_all_system_skills_list

    import gibh_agent.tools  # noqa: F401
    from gibh_agent.core.tool_registry import registry

    reg_names = {t["name"] for t in registry.get_all_tools_json() if t.get("name")}
    idx = ImplementedIndex()

    for seed in get_all_system_skills_list():
        name = (seed.get("name") or "").strip()
        pt = seed.get("prompt_template") or ""
        if not name or not infer_skill_implemented_from_prompt(pt):
            continue

        idx.seed_names_norm.add(_norm(name))

        panshi_key = PANSHI_TOOL_CHAIN_KEY_BY_DISPLAY_NAME.get(name)
        if panshi_key:
            idx.panshi_keys_norm.add(_norm(panshi_key))

        for m in ROUTE_RE.finditer(pt):
            tag = (m.group(0) or "").lower()
            tool_id = (m.group(1) or "").strip()
            if not tool_id:
                continue
            if tag.startswith("[omics_route:"):
                idx.route_tools.add(tool_id)
                idx.panshi_keys_norm.add(_norm(tool_id))
                continue
            if tool_id in reg_names:
                idx.route_tools.add(tool_id)
                idx.panshi_keys_norm.add(_norm(tool_id))

    return idx


def is_already_implemented_in_agent(
    *,
    tool: Dict[str, Any],
    seed: Optional[Dict[str, Any]],
    route_tool: str,
    skill_name: str,
    index: ImplementedIndex,
) -> bool:
    """该行是否已在技能广场 + Registry（+ 可选磐石 key）完成对接。"""
    key = (tool.get("tool_chain_key") or "").strip()
    if _norm(skill_name) in index.seed_names_norm:
        return True
    if key and _norm(key) in index.panshi_keys_norm:
        return True
    if route_tool and route_tool in index.route_tools:
        return True

    if seed:
        from gibh_agent.core.skill_plaza_utils import infer_skill_implemented_from_prompt

        if infer_skill_implemented_from_prompt(seed.get("prompt_template") or ""):
            return True

    return False


def build_matrix_rows(
    *,
    exclude_implemented: bool = True,
    prefer_catalog: bool = True,
) -> List[Dict[str, str]]:
    from gibh_agent.db.panshi_skill_meta import PANSHI_TOOL_CHAIN_KEY_BY_DISPLAY_NAME
    from gibh_agent.db.seed_skills import get_all_system_skills_list

    tools: List[Dict[str, Any]] = load_panshi_tools(prefer_catalog=prefer_catalog)
    seeds = get_all_system_skills_list()
    by_name = {(s.get("name") or "").strip(): s for s in seeds}
    by_norm = {_norm(s.get("name") or ""): s for s in seeds}
    impl_index = build_implemented_index() if exclude_implemented else None

    rows: List[Dict[str, str]] = []
    skipped_impl = 0
    for tool in sorted(
        [t for t in tools if is_bioinformatics_tool(t)],
        key=lambda x: (classify_discipline(x), x.get("display_name") or ""),
    ):
        key = (tool.get("tool_chain_key") or "").strip()
        disp = (tool.get("display_name") or "").strip()
        desc = (tool.get("long_description") or "").strip()
        discipline = classify_discipline(tool)

        seed = by_name.get(disp) or by_norm.get(_norm(disp)) or by_norm.get(_norm(key))
        if not seed:
            rev_panshi = {v: k for k, v in PANSHI_TOOL_CHAIN_KEY_BY_DISPLAY_NAME.items()}
            seed = by_name.get(rev_panshi.get(key, ""))
        if not seed:
            for s in seeds:
                pt = s.get("prompt_template") or ""
                sname = (s.get("name") or "").strip()
                if key and (
                    key.lower() in pt.lower()
                    or _norm(key) == _norm(sname)
                    or key.lower() in sname.lower()
                ):
                    seed = s
                    break

        skill_name = (seed.get("name") if seed else disp) or key
        sub = infer_sub_category(tool, discipline, seed)
        intro = (seed.get("description") if seed else desc) or (
            f"磐石 ToolChain 工具「{disp}」（{key}）。"
        )
        if len(intro) > 240:
            intro = intro[:237] + "…"

        route_tool = ""
        if seed:
            m = ROUTE_RE.search(seed.get("prompt_template") or "")
            if m:
                route_tool = m.group(1)

        if impl_index and is_already_implemented_in_agent(
            tool=tool,
            seed=seed,
            route_tool=route_tool,
            skill_name=skill_name,
            index=impl_index,
        ):
            skipped_impl += 1
            continue

        skills_type = (tool.get("skills_type") or "").strip() or infer_skills_type(
            tool=tool,
            seed=seed,
            route_tool=route_tool,
            skill_name=skill_name,
            sub_category=sub,
            intro=intro,
        )
        if tool.get("invoke_status"):
            intro = f"{intro}（磐石：{tool['invoke_status']}）"

        rows.append(
            {
                "学科": discipline,
                "技能名称": skill_name,
                "子分类": sub,
                "skills类型": skills_type,
                "简介": intro,
                "tool_chain_key": key,
                "_panshi_display_name": disp,
                "_registry_skill_route": route_tool,
                "_implemented": "是" if route_tool else "广场占位",
            }
        )
    if exclude_implemented and skipped_impl:
        print(
            f"已剔除本仓库已实现技能/算子条目: {skipped_impl}（保留待对接 {len(rows)}）",
            file=sys.stderr,
        )
    return rows


def render_registry_category_table() -> str:
    import gibh_agent.tools  # noqa: F401 — 触发注册
    from gibh_agent.core.tool_registry import registry

    tools = registry.get_all_tools_json()
    counts = Counter(t.get("category") or "Unknown" for t in tools)
    lines = [
        "## Registry 分类统计（本仓库底层算子）",
        "",
        f"当前 Registry 共 **{len(tools)}** 个工具（含运行时 MCP 同步条目时数量可能变化）。",
        "",
        "| Registry category | 数量 | 典型域 |",
        "| --- | ---: | --- |",
    ]
    domain_hints = {
        "Genomics": "基因组胚系/变异",
        "scRNA-seq": "单细胞/转录组",
        "Metabolomics": "代谢组",
        "Proteomics": "蛋白组",
        "Epigenomics": "表观组",
        "Spatial": "空间转录",
        "Radiomics": "影像组学",
        "MedicinalChemistry": "药物化学/RDKit",
        "Biomedicine": "生物医药技能",
        "Bioinformatics": "序列/注释",
        "Public API": "基因/通路公共库",
        "MCP": "联网检索",
    }
    for cat, n in sorted(counts.items(), key=lambda x: (-x[1], x[0])):
        hint = domain_hints.get(cat, "")
        lines.append(f"| {cat} | {n} | {hint} |")
    lines.append("")
    return "\n".join(lines)


def _catalog_coverage_lines() -> List[str]:
    """从 JSON 元数据生成「条数漏斗」说明。"""
    lines = ["", "**数据覆盖说明（条数漏斗）**："]
    if _DATA_SOURCE == "skill_plaza_catalog" and CATALOG_JSON.is_file():
        try:
            payload = json.loads(CATALOG_JSON.read_text(encoding="utf-8"))
        except (json.JSONDecodeError, OSError):
            payload = {}
        total = payload.get("entries_total") or 0
        by_cat = payload.get("entries_by_main_category") or {}
        lines.append(
            f"- 当前数据源：**技能广场用户粘贴全量**（`scienceone_skill_plaza_catalog.json`），"
            f"共 **{total}** 条（生物医药 **{by_cat.get('生物医药', 0)}** + 化学 **{by_cat.get('化学', 0)}**）"
        )
        lines.append(
            "- 更新粘贴文本后请运行：`PYTHONPATH=. python3 scripts/parse_scienceone_skill_plaza_paste.py`"
        )
    if TOOLS_JSON.is_file():
        try:
            payload = json.loads(TOOLS_JSON.read_text(encoding="utf-8"))
            tools = payload.get("tools") if isinstance(payload, dict) else []
            bundle_n = len(tools) if isinstance(tools, list) else 0
            bio_n = sum(1 for t in tools if isinstance(t, dict) and is_bioinformatics_tool(t))
            lines.append(
                f"- 对照：ToolChain 前端 bundle 静态解析 **{bundle_n}** 条，生信/药化筛选 **{bio_n}** 条"
            )
            note = (payload.get("catalog_coverage_note") or "").strip()
            if note and _DATA_SOURCE != "skill_plaza_catalog":
                lines.append(f"- {note}")
        except (json.JSONDecodeError, OSError):
            pass
    return lines


def render_bio_matrix_markdown(
    rows: List[Dict[str, str]],
    *,
    exclude_implemented: bool = True,
) -> str:
    source_page = "https://www.scienceone.cn/toolchain/#/science/tool"
    data_note = (
        "`gibh_agent/data/scienceone_skill_plaza_catalog.json`（技能广场粘贴全量）"
        if _DATA_SOURCE == "skill_plaza_catalog"
        else "`gibh_agent/data/scienceone_toolchain_tools.json`（ToolChain bundle）"
    )
    lines = [
        "## 生物信息分析工具（磐石技能广场 / ToolChain 对照表 · 待对接）",
        "",
        f"对照源：[ScienceOne 技能广场](https://www.scienceone.cn/toolchain/#/skillPlaza) / "
        f"[ToolChain 科学工具]({source_page})；元数据：{data_note}。",
        "",
        "下表用于技能广场**待扩展**能力选型；**已剔除**本仓库 `seed_skills` 中已含 "
        "`[Skill_Route]` / `[Omics_Route]` 且 Registry 可解析的条目（与 GET /api/skills "
        "`is_implemented` 逻辑一致）。",
        "",
        f"- 表中条目（待对接）：**{len(rows)}**",
    ]
    lines.extend(_catalog_coverage_lines())
    if exclude_implemented:
        lines.append(
            "- 已实现条目：已从表中移除（见 `技能库.md` / 技能广场置顶项，勿重复规划）"
        )
    else:
        n_impl = sum(1 for r in rows if r.get("_registry_skill_route"))
        lines.append(f"- 含已实现 `[Skill_Route]`：**{n_impl}**（`--include-implemented` 模式）")
    lines.extend(
        [
        "",
            "**skills类型** 说明：`数据库调用`（公共库/API 检索）｜`纯prompt`（仅模板/占位）｜"
            "`脚本`（Registry/Worker 可执行算子）｜`模型`（深度学习基础模型或远端推理）。",
            "",
            "| 学科 | 技能名称 | 子分类 | skills类型 | 简介 | tool_chain_key |",
            "| --- | --- | --- | --- | --- | --- |",
        ]
    )
    type_counts = Counter(r.get("skills类型", "") for r in rows)
    for r in rows:
        lines.append(
            "| "
            + " | ".join(
                _escape_md_cell(r[k])
                for k in ("学科", "技能名称", "子分类", "skills类型", "简介", "tool_chain_key")
            )
            + " |"
        )
    lines.append("")
    lines.append(
        "skills类型分布："
        + "；".join(f"{k} **{v}**" for k, v in sorted(type_counts.items()) if k)
    )
    lines.extend(
        [
            "",
            "### 维护命令",
            "",
            "```bash",
            "cd /home/ubuntu/GIBH-AGENT-V2",
            "# 1) 将技能广场粘贴文本放入 docs/raw/scienceone_skill_plaza_user_paste.txt 后解析",
            "PYTHONPATH=. python3 scripts/parse_scienceone_skill_plaza_paste.py",
            "# 2) 可选：刷新 ToolChain bundle（约 120 条，需网络）",
            "PYTHONPATH=. python3 scripts/refresh_scienceone_toolchain_names.py",
            "# 3) 重生成本节对照表",
            "PYTHONPATH=. python3 scripts/generate_bioinformatics_tools_matrix_md.py",
            "# 4) 写回根目录 工具库.md（含 Registry JSON + 本表）",
            "PYTHONPATH=. python3 scripts/export_tool_library_md.py",
            "# 5) 可选：含已实现条目（全量对照）",
            "PYTHONPATH=. python3 scripts/generate_bioinformatics_tools_matrix_md.py --include-implemented --csv docs/磐石技能广场_生物医药化学_全量.csv",
            "# 6) 可选：仅用 ToolChain bundle（约 120 条）",
            "PYTHONPATH=. python3 scripts/generate_bioinformatics_tools_matrix_md.py --bundle-only",
            "```",
            "",
        ]
    )
    return "\n".join(lines)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--csv",
        type=Path,
        help="可选：另存 UTF-8-SIG CSV（Excel/WPS）",
    )
    ap.add_argument(
        "--out",
        type=Path,
        default=ROOT / "docs" / "generated" / "bioinformatics_tools_matrix.md",
        help="输出 Markdown 片段路径",
    )
    ap.add_argument(
        "--include-implemented",
        action="store_true",
        help="保留本仓库已实现的技能条目（默认剔除）",
    )
    ap.add_argument(
        "--bundle-only",
        action="store_true",
        help="仅使用 ToolChain bundle（忽略技能广场粘贴 catalog）",
    )
    args = ap.parse_args()

    exclude = not args.include_implemented
    rows = build_matrix_rows(
        exclude_implemented=exclude,
        prefer_catalog=not args.bundle_only,
    )
    md = render_bio_matrix_markdown(rows, exclude_implemented=exclude)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(md, encoding="utf-8")
    print(f"Wrote {args.out} ({len(rows)} rows)")

    if args.csv:
        args.csv.parent.mkdir(parents=True, exist_ok=True)
        fields = ["学科", "技能名称", "子分类", "skills类型", "简介", "tool_chain_key"]
        with args.csv.open("w", encoding="utf-8-sig", newline="") as f:
            w = csv.DictWriter(f, fieldnames=fields, extrasaction="ignore")
            w.writeheader()
            for r in rows:
                w.writerow({k: r.get(k, "") for k in fields})
        print(f"Wrote {args.csv}")


if __name__ == "__main__":
    main()
