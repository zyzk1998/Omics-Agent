#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
磐石 675 工具清单 → 4-D 矩阵评分 → 首发 10–15 精选技能。

用法（仓库根目录）:
  PYTHONPATH=. python scripts/select_priority_skills_top15.py
  PYTHONPATH=. python scripts/select_priority_skills_top15.py --csv docs/生物信息分析工具.csv --top 15
"""
from __future__ import annotations

import argparse
import re
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional, Set, Tuple

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

# --- 硬淘汰：命中任一则不参与首发（算力/依赖/长尾）---
HARD_EXCLUDE_NAME_DESC = re.compile(
    r"CP2K|LAMMPS|量子化学|分子动力学模拟|MPI\b|从头组装|全基因组组装|de\s*novo\s*assembly|"
    r"AlphaFold2?|ESMFold|RFdiffusion|DiffDock|DiffSBDD|GraphCast|ClimateGPT|ProtGPT|"
    r"MatterGen|LeanDojo|DSAgent|OligoFormer|GenMol|细胞周期各阶段持续时间|"
    r"写作|SmartPPT|DrawIO|翻译|润色|会议|DeepResearch|DataFormulator|"
    r"组装安装信息|非活性成分|处置信息|安全处理警告|实验室检测干扰|环境警告|"
    r"商品名|剂型与规格|专利查询|交叉检索|词汇过滤|GraphQL|"
    r"4DN|BioStudies|ENCODE|cBioPortal|GTEx|CELLxGENE_census版本|版本获取|"
    r"标识符交叉|协方差模型|种子比对获取|树数据|XML格式|STITCH",
    re.I,
)

# 冗余 ChEMBL/药物微查询（保留聚合级「ChEMBL药物检索」类）
CHEMBL_MICRO_QUERY = re.compile(
    r"根据ChEMBL|依据ChEMBL|通过ChEMBL|ChEMBL ID获取|ChEMBL ID查询|"
    r"获取药物(同义词|详情|通用名|适应症|不良反应)|查询药物(成分|商品)|检索药物实验室",
    re.I,
)

# 首发必选（产品委员会锁定；须在 CSV 中存在）
MUST_INCLUDE_NAMES: List[str] = [
    "转录组学标准全流程",
    "差异标志物发现",
    "基因富集分析",
    "RNA二级结构可视化工具",
    "分子格式转换工具",
    "二维分子结构可视化",
    "检索相似小分子",
    "ChEMBL药物检索",
    "ADMET性质预测工具",
    "核酸序列比对",
    "蛋白质序列比对",
    "子结构搜索化合物",
    "通过SMILES获取CID",
    "MHC关联表位检索",
    "3D分子结构渲染工具",
]

# 低价值 DB 检索（再高分也不进首发）
LOW_VALUE_DB = re.compile(
    r"项目检索|相似实体|查询转换器|版本获取|标识符|元数据|文献检索|组织信息|"
    r"细胞系信息|审批状态|不良反应事件|黑框警告",
    re.I,
)

# --- 维度 1：业务刚需（加分词）---
FREQ_SCORE_RULES: List[Tuple[re.Pattern, float]] = [
    (re.compile(r"质控|QC|过滤|预处理", re.I), 2.5),
    (re.compile(r"差异表达|差异分析|差异标志|标志物", re.I), 3.0),
    (re.compile(r"富集分析|GSEA|通路富集|基因集", re.I), 3.0),
    (re.compile(r"单细胞|scRNA|转录组|降维|聚类|UMAP|Scanpy", re.I), 3.0),
    (re.compile(r"格式转换|分子格式|SMILES|MOL|SDF", re.I), 3.0),
    (re.compile(r"二维|2D|结构可视化|分子图像|molPreview|结构渲染", re.I), 2.5),
    (re.compile(r"相似|Similarity|Tanimoto|子结构", re.I), 2.5),
    (re.compile(r"ADMET|类药性|Lipinski|五规则|成药性", re.I), 2.5),
    (re.compile(r"表位|Epitope|BepiPred|MHC", re.I), 2.0),
    (re.compile(r"RNA二级|RNAfold|ViennaRNA", re.I), 2.5),
    (re.compile(r"序列比对|比对", re.I), 2.0),
    (re.compile(r"ChEMBL药物检索|化合物检索|PubChem", re.I), 2.0),
    (re.compile(r"全流程|标准流程", re.I), 3.0),
]

# --- 维度 2：环境依赖（从 10 扣分）---
DEP_HEAVY_PENALTY = re.compile(
    r"深度学习|神经网络|预训练模型|扩散模型|GPU集群|超算|CUDA|训练模型|"
    r"conda install.*gcc|特殊编译|Java依赖|MATLAB",
    re.I,
)
DEP_LIGHT_BOOST = re.compile(
    r"RDKit|rdkit|API|数据库|检索|HTTP|REST|GraphQL(?!.*评分)|纯Python|pip install|"
    r"脚本|轻量|快速",
    re.I,
)
DEP_CLI_OK = re.compile(r"samtools|openbabel|obabel|muscle|blast|bedtools", re.I)

# --- 维度 3：算力耗时（从 10 扣分）---
SLOW_PENALTY = re.compile(
    r"全基因组|WGS|WES|从头|组装|比对至参考基因组|长时间|数天|TB级|"
    r"大规模训练|亿级|whole genome",
    re.I,
)
FAST_BOOST = re.compile(r"秒|分钟|快速|实时|轻量|预览|检索|查询|API", re.I)

# --- 维度 4：生态协同（加分）---
SYNERGY_RULES: List[Tuple[re.Pattern, float]] = [
    (re.compile(r"单细胞|时空|轨迹|STED|空间|scRNA|h5ad", re.I), 3.0),
    (re.compile(r"Lipinski|类药性|SMILES|Tanimoto|相似|ADMET|RDKit|分子格式", re.I), 3.0),
    (re.compile(r"BepiPred|表位|抗体|免疫", re.I), 2.5),
    (re.compile(r"富集|GSEA|差异|标志物|Scanpy|蛋白组|基因组|代谢", re.I), 2.5),
    (re.compile(r"CSV|TSV|表格|文本|格式转换", re.I), 2.0),
]

# 与 GIBH 已落地能力强对齐（额外 +5，用于打破平局）
GIBH_SYNERGY_NAMES = frozenset(MUST_INCLUDE_NAMES) | frozenset(
    {
        "三维分子结构可视化",
        "分子可视化工具",
        "相似分子检索",
        "相似性搜索化合物",
        "子结构搜索",
    }
)

# 首发编制：每池目标数量（组学 / 化学 / 通用）
POOL_TARGETS = {"omics": 5, "chem": 5, "general": 4}
POOL_MAX = 15


@dataclass
class ScoredTool:
    row_index: int
    name: str
    discipline: str
    sub_category: str
    skills_type: str
    description: str
    tool_chain_key: str
    pool: str
    score_freq: float = 0.0
    score_dep: float = 0.0
    score_speed: float = 0.0
    score_synergy: float = 0.0
    score_total: float = 0.0
    excluded: bool = False
    exclude_reason: str = ""
    selection_reason: str = ""

    def to_dict(self) -> Dict[str, Any]:
        return {
            "学科": self.discipline,
            "技能名称": self.name,
            "子分类": self.sub_category,
            "skills类型": self.skills_type,
            "简介": self.description,
            "tool_chain_key": self.tool_chain_key,
            "score_frequency": round(self.score_freq, 2),
            "score_dependency": round(self.score_dep, 2),
            "score_speed": round(self.score_speed, 2),
            "score_synergy": round(self.score_synergy, 2),
            "score_total": round(self.score_total, 2),
            "pool": self.pool,
            "Selection_Reason": self.selection_reason,
        }


def _text_blob(name: str, desc: str) -> str:
    return f"{name} {desc}"


def classify_pool(discipline: str, name: str, desc: str) -> str:
    t = _text_blob(name, desc)
    if re.search(r"转录组|单细胞|scRNA|RNA二级|富集分析|差异标志|蛋白.*比对|核酸.*比对", t, re.I):
        return "omics"
    if discipline == "化学与分子信息学" or re.search(
        r"SMILES|分子|药物|化合物|ChEMBL|PubChem|RDKit|ADMET|化学|molPreview", t, re.I
    ):
        return "chem"
    if re.search(r"序列比对|SMILES获取|格式|表格|文本|CID|表位", t, re.I):
        return "general"
    return "omics"


def score_frequency(name: str, desc: str) -> float:
    t = _text_blob(name, desc)
    s = 0.0
    for pat, w in FREQ_SCORE_RULES:
        if pat.search(t):
            s += w
    return min(10.0, s)


def score_dependency(name: str, desc: str, skills_type: str) -> float:
    t = _text_blob(name, desc)
    s = 8.0
    if DEP_HEAVY_PENALTY.search(t):
        s -= 4.0
    if DEP_CLI_OK.search(t):
        s -= 1.0
    if DEP_LIGHT_BOOST.search(t) or (skills_type or "") in ("脚本", "数据库调用"):
        s += 1.5
    if re.search(r"预测与建模", desc) and not DEP_LIGHT_BOOST.search(t):
        s -= 2.0
    return max(0.0, min(10.0, s))


def score_speed(name: str, desc: str) -> float:
    t = _text_blob(name, desc)
    s = 7.0
    if SLOW_PENALTY.search(t):
        s -= 5.0
    if FAST_BOOST.search(t):
        s += 2.0
    if re.search(r"全流程", t):
        s -= 1.0  # 流程可分钟–小时级，仍优于 HPC 组装
    return max(0.0, min(10.0, s))


def score_synergy(name: str, desc: str) -> float:
    t = _text_blob(name, desc)
    s = 0.0
    for pat, w in SYNERGY_RULES:
        if pat.search(t):
            s += w
    if name in GIBH_SYNERGY_NAMES:
        s += 5.0
    return min(10.0, s)


def should_exclude(name: str, desc: str) -> Optional[str]:
    if name in MUST_INCLUDE_NAMES:
        return None
    t = _text_blob(name, desc)
    if HARD_EXCLUDE_NAME_DESC.search(t):
        return "重型/长尾/办公类或纯元数据微查询"
    if LOW_VALUE_DB.search(t):
        return "低频数据库/元数据检索，非分析算子"
    if CHEMBL_MICRO_QUERY.search(name):
        return "ChEMBL 细粒度重复查询（首发保留聚合检索即可）"
    if re.search(r"^获取|^查询|^检索", name) and not re.search(
        r"富集|比对|相似|格式|可视化|预测|全流程|ADMET|表位|RNA二级|ChEMBL药物", name
    ):
        if name not in GIBH_SYNERGY_NAMES:
            return "低频数据库单字段查询"
    return None


REASON_BY_NAME: Dict[str, str] = {
    "转录组学标准全流程": "单细胞 QC→聚类→注释→差异/富集主链路，与 Scanpy 快车道及 STED-EC 下游完全同构；依赖以 Scanpy 为主、分钟–小时级。",
    "差异标志物发现": "代谢/转录标志物筛选为日常高频；对接现有代谢组与 scRNA marker 工具，无 HPC 组装。",
    "基因富集分析": "通路富集是解释差异结果的标配；与 GSEA/富集节点及转录组全流程末段天然串联。",
    "RNA二级结构可视化工具": "序列→结构快速预览；可走 ViennaRNA/PySkills，秒–分钟级，补全核酸分析拼图。",
    "分子格式转换工具": "SMILES/MOL/SDF 互转是化学流水线入口；RDKit/OpenBabel 成熟易部署，支撑药性评估前后格式统一。",
    "二维分子结构可视化": "SMILES→2D 图是药物化学日报级需求；纯 RDKit 渲染，秒级出图。",
    "检索相似小分子": "结构类似物发现核心；与已上线 Tanimoto/药物相似性快车道一致，联网+指纹双模式。",
    "ChEMBL药物检索": "药物发现信息入口；REST API 轻依赖，可与相似性检索、ADMET 串联。",
    "ADMET性质预测工具": "成药性评估高频；与 Lipinski 五规则快车道互补，规则/模型均可快速返回。",
    "核酸序列比对": "序列分析基石；CLI 成熟(blast/muscle 类)，分钟级，可接入 FASTA 挂载链路。",
    "蛋白质序列比对": "蛋白注释与同源分析基础；与保守性/表位技能上下游衔接。",
    "子结构搜索化合物": "药效团与骨架筛选常用；RDKit 子结构匹配，秒–分钟级。",
    "通过SMILES获取CID": "PubChem 标识符桥接；API 调用、秒级，连通相似检索与数据库查询。",
    "MHC关联表位检索": "免疫表位发现；与 BepiPred 表位预测形成「预测+库检」闭环。",
    "3D分子结构渲染工具": "分子 3D 展示与沟通；molPreview/RDKit 轻量可视化，支撑化学报告右栏。",
}


def build_selection_reason(st: ScoredTool) -> str:
    if st.name in REASON_BY_NAME:
        return REASON_BY_NAME[st.name]
    parts: List[str] = []
    if st.name in GIBH_SYNERGY_NAMES:
        parts.append("与现有组学快车道/化学 RDKit·BepiPred 栈高度对齐")
    if st.score_freq >= 7:
        parts.append("高频基石能力")
    if st.score_dep >= 7:
        parts.append("依赖轻(RDKit/API/标准CLI)")
    if st.score_speed >= 7:
        parts.append("可分钟级内响应")
    if st.score_synergy >= 7:
        parts.append("可串联单细胞·药性·表格处理主线")
    if not parts:
        parts.append("四维综合得分领先同池候选")
    return "；".join(parts[:3]) + "。"


def score_dataframe(df: pd.DataFrame) -> List[ScoredTool]:
    scored: List[ScoredTool] = []
    seen_names: Set[str] = set()

    for idx, row in df.iterrows():
        name = str(row.get("技能名称", "")).strip()
        if not name or name in seen_names:
            continue
        seen_names.add(name)

        desc = str(row.get("简介", "")).strip()
        discipline = str(row.get("学科", "")).strip()
        sub = str(row.get("子分类", "")).strip()
        stype = str(row.get("skills类型", "")).strip()
        tkey = str(row.get("tool_chain_key", "")).strip()
        if tkey.lower() == "nan":
            tkey = ""

        ex = should_exclude(name, desc)
        pool = classify_pool(discipline, name, desc)

        st = ScoredTool(
            row_index=int(idx),
            name=name,
            discipline=discipline,
            sub_category=sub,
            skills_type=stype,
            description=desc,
            tool_chain_key=tkey,
            pool=pool,
        )
        if ex:
            st.excluded = True
            st.exclude_reason = ex
            scored.append(st)
            continue

        st.score_freq = score_frequency(name, desc)
        st.score_dep = score_dependency(name, desc, stype)
        st.score_speed = score_speed(name, desc)
        st.score_synergy = score_synergy(name, desc)
        # 4-D 加权：刚需与协同权重略高
        st.score_total = (
            st.score_freq * 0.30
            + st.score_dep * 0.25
            + st.score_speed * 0.20
            + st.score_synergy * 0.25
        )
        scored.append(st)
    return scored


def select_top_by_pool(
    scored: List[ScoredTool],
    df: pd.DataFrame,
    top_n: int = 15,
) -> List[ScoredTool]:
    by_name = {s.name: s for s in scored}
    picked: List[ScoredTool] = []

    # 第一轮：产品锁定必选
    for name in MUST_INCLUDE_NAMES:
        if name in by_name and not by_name[name].excluded:
            picked.append(by_name[name])
        else:
            row = df[df["技能名称"].astype(str).str.strip() == name]
            if row.empty:
                continue
            r = row.iloc[0]
            st = ScoredTool(
                row_index=int(row.index[0]),
                name=name,
                discipline=str(r.get("学科", "")),
                sub_category=str(r.get("子分类", "")),
                skills_type=str(r.get("skills类型", "")),
                description=str(r.get("简介", "")),
                tool_chain_key=str(r.get("tool_chain_key", "")).replace("nan", ""),
                pool=classify_pool(
                    str(r.get("学科", "")),
                    name,
                    str(r.get("简介", "")),
                ),
            )
            st.score_freq = score_frequency(name, st.description)
            st.score_dep = score_dependency(name, st.description, st.skills_type)
            st.score_speed = score_speed(name, st.description)
            st.score_synergy = score_synergy(name, st.description)
            st.score_total = (
                st.score_freq * 0.30
                + st.score_dep * 0.25
                + st.score_speed * 0.20
                + st.score_synergy * 0.25
            )
            picked.append(st)

    # 去重保序
    seen: Set[str] = set()
    unique_picked: List[ScoredTool] = []
    for st in picked:
        if st.name in seen:
            continue
        seen.add(st.name)
        unique_picked.append(st)
    picked = unique_picked[:top_n]

    # 若必选不足 top_n，从候选池补足（排除相似重复）
    if len(picked) < top_n:
        candidates = [s for s in scored if not s.excluded and s.name not in seen]
        candidates.sort(key=lambda x: x.score_total, reverse=True)
        for st in candidates:
            if len(picked) >= top_n:
                break
            # 避免三个「相似检索」重复
            if st.name in ("相似分子检索", "相似性搜索化合物") and any(
                p.name in ("检索相似小分子", "相似分子检索", "相似性搜索化合物") for p in picked
            ):
                continue
            picked.append(st)
            seen.add(st.name)

    for st in picked:
        st.selection_reason = build_selection_reason(st)
    return picked[:top_n]


def main() -> int:
    parser = argparse.ArgumentParser(description="4-D 矩阵筛选首发技能 Top N")
    parser.add_argument(
        "--csv",
        type=Path,
        default=ROOT / "docs" / "生物信息分析工具.csv",
    )
    parser.add_argument(
        "--script-only",
        action="store_true",
        help="仅保留 skills类型=脚本 且学科为生物医药/化学与分子信息学",
    )
    parser.add_argument(
        "--disciplines",
        default="生物医药,化学与分子信息学",
        help="逗号分隔学科过滤（与 --script-only 联用）",
    )
    parser.add_argument("--out", type=Path, default=ROOT / "docs" / "priority_skills_top15.csv")
    parser.add_argument("--top", type=int, default=15)
    args = parser.parse_args()

    df = pd.read_csv(args.csv, encoding="utf-8-sig").fillna("")
    if args.script_only:
        allowed = {d.strip() for d in args.disciplines.split(",") if d.strip()}
        df = df[
            df["skills类型"].astype(str).str.strip() == "脚本"
            & df["学科"].astype(str).str.strip().isin(allowed)
        ].copy()
        logger.info(
            "脚本+学科过滤后候选: %d 条（学科=%s）",
            len(df),
            ",".join(sorted(allowed)),
        )
    scored = score_dataframe(df)
    selected = select_top_by_pool(
        scored, df, top_n=min(max(args.top, 10), 15)
    )

    out_df = pd.DataFrame([s.to_dict() for s in selected])
    args.out.parent.mkdir(parents=True, exist_ok=True)
    out_df.to_csv(args.out, index=False, encoding="utf-8-sig")

    print(f"\n✅ 已写入 {args.out}（{len(selected)} 项）\n")
    print("=" * 72)
    print("首发精选技能（请总指挥确认）")
    print("=" * 72)
    for i, st in enumerate(selected, 1):
        print(f"{i:2d}. [{st.pool:7s}] {st.name}")
        print(f"    理由: {st.selection_reason}")
        print(
            f"    4-D: 刚需={st.score_freq:.1f} 依赖={st.score_dep:.1f} "
            f"速度={st.score_speed:.1f} 协同={st.score_synergy:.1f} 综合={st.score_total:.1f}"
        )
    print("=" * 72)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
