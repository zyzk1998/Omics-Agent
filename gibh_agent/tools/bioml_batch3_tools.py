# -*- coding: utf-8 -*-
"""
第三批生物医药轻量仿真工具：子进程调用 `gibh_agent/assets/bioml_batch3/*.py`（标准库为主）。

与 `crispr_cas9_tool` 并列；须由 `skill_agent` 与 `gibh_agent.tools` 显式预加载，避免快车道下 ToolRegistry 未注册。
"""
from __future__ import annotations

import asyncio
import json
import logging
import os
import sys
import uuid
from pathlib import Path
from typing import Any, Dict, Optional

from ..core.executor import WorkflowExecutor
from ..core.tool_registry import registry
from ..core.utils import safe_tool_execution

logger = logging.getLogger(__name__)


def _package_dir() -> Path:
    return Path(__file__).resolve().parent.parent


def _results_dir() -> Path:
    raw = (os.getenv("RESULTS_DIR") or "").strip()
    if raw:
        return Path(raw).expanduser().resolve()
    docker_results = Path("/app/results")
    if docker_results.parent.is_dir():
        return docker_results.resolve()
    return (_package_dir().parent / "results").resolve()


def _public_path_url(rel_path: str) -> str:
    rel = rel_path if rel_path.startswith("/") else "/" + rel_path.lstrip("/")
    base = (os.getenv("PUBLIC_RESULTS_BASE_URL") or "").rstrip("/")
    return f"{base}{rel}" if base else rel


def _asset_script(name: str) -> Path:
    return _package_dir() / "assets" / "bioml_batch3" / name


def _resolve_uploaded_path(file_path: str) -> str:
    if not (file_path or "").strip():
        raise ValueError("file_path 为空")
    ex = WorkflowExecutor()
    return ex._resolve_file_path(file_path.strip())


async def _run_script(
    script: Path,
    argv: list[str],
    *,
    cwd: Optional[Path] = None,
    timeout: float = 120.0,
) -> tuple[int, str, str]:
    if not script.is_file():
        raise FileNotFoundError(str(script))
    py = (sys.executable or "python3").strip()
    cmd = [py, str(script), *argv]
    env = {**os.environ, "PYTHONUNBUFFERED": "1"}
    logger.info("bioml_batch3 子进程: %s", " ".join(cmd[:10]) + (" ..." if len(cmd) > 10 else ""))
    proc = await asyncio.create_subprocess_exec(
        *cmd,
        cwd=str(cwd or script.parent),
        env=env,
        stdout=asyncio.subprocess.PIPE,
        stderr=asyncio.subprocess.PIPE,
    )
    try:
        out_b, err_b = await asyncio.wait_for(proc.communicate(), timeout=timeout)
    except asyncio.TimeoutError:
        proc.kill()
        await proc.communicate()
        raise TimeoutError(f"子进程超时（>{timeout}s）") from None
    out = (out_b or b"").decode("utf-8", errors="replace").strip()
    err = (err_b or b"").decode("utf-8", errors="replace").strip()
    return int(proc.returncode or 0), out, err


def _build_rna_expert_markdown(data: Dict[str, Any], *, json_url: str = "") -> str:
    """
    将 rna_structure.py 的 JSON 摘要渲染为专家级 Markdown（非仅代码块）。
    """
    lines: list[str] = []
    lines.append("### RNA 二级结构解析报告")
    lines.append("")

    if not data.get("valid", False):
        err = (data.get("error_message") or "结构校验失败").strip()
        lines.append("- **结构校验**: ❌ 不合法")
        lines.append(f"- **说明**: {err}")
        lines.append("")
        if json_url:
            lines.append(f"- [下载完整 JSON]({json_url})")
        return "\n".join(lines).strip()

    struct = data.get("structure") or ""
    seq_eff = data.get("sequence") or ""
    seq_in = data.get("sequence_original") or seq_eff
    n = int(data.get("length") or len(struct))
    align_note = (data.get("alignment_note") or "").strip()

    lines.append(f"- **点括号结构** (`{n}` nt): `{struct}`")
    lines.append(f"- **输入序列**（原始）: `{seq_in}`" + (f"（`{len(seq_in)}` nt）" if seq_in else "（未提供）"))
    if seq_eff and seq_eff != seq_in:
        lines.append(f"- **用于配对的序列**（与结构等长）: `{seq_eff}`")
    if align_note:
        lines.append(f"- **长度对齐说明**: {align_note}")
    lines.append("- **结构校验**: ✅ 合法")
    lines.append("")

    up = data.get("unpaired_percentage")
    if up is not None:
        lines.append(f"- **未配对比例**: {float(up):.1f}%")
    if seq_eff and data.get("gc_content") is not None:
        lines.append(f"- **GC 含量**: {data.get('gc_content')}%")
    if seq_eff and data.get("estimated_free_energy") is not None:
        lines.append(f"- **估计自由能（简化 Turner）**: {data.get('estimated_free_energy')} kcal/mol（仅供参考）")
    lines.append("")

    pairs = data.get("base_pairs") or []
    lines.append("#### 配对详情（Base pairs）")
    lines.append("")
    if not pairs:
        lines.append("*（未解析到碱基对）*")
    else:
        lines.append("| 索引 i | 索引 j | 位点碱基 (i–j) | 配对类型 |")
        lines.append("| ---: | ---: | --- | --- |")
        for p in pairs[:200]:
            i = int(p.get("i", -1))
            j = int(p.get("j", -1))
            pt = (p.get("pair_type") or "").strip() or "—"
            if seq_eff and 0 <= i < len(seq_eff) and 0 <= j < len(seq_eff):
                bi, bj = seq_eff[i], seq_eff[j]
                bases = f"{bi}-{bj}"
            else:
                bases = "—"
            lines.append(f"| {i} | {j} | `{bases}` | **{pt}** |")
        if len(pairs) > 200:
            lines.append("")
            lines.append(f"*（仅展示前 200 对，共 {len(pairs)} 对；完整列表见 JSON。）*")
    lines.append("")

    stems = data.get("stems") or []
    lines.append(f"#### 茎区（Stems）概览 — 共 **{len(stems)}** 条")
    lines.append("")
    if stems:
        for s in stems[:30]:
            sid = s.get("id")
            le = s.get("length")
            ee = s.get("estimated_energy")
            sc = s.get("stability_score")
            extra = ""
            if ee is not None and sc is not None:
                extra = f"，估计能量 {ee} kcal/mol，稳定性评分 {sc}"
            lines.append(
                f"- **Stem #{sid}**: 5' [{s.get('start_5p')}–{s.get('end_5p')}] / "
                f"3' [{s.get('start_3p')}–{s.get('end_3p')}]，长度 **{le}** bp{extra}"
            )
        if len(stems) > 30:
            lines.append("")
            lines.append(f"*（其余 {len(stems) - 30} 条茎区见 JSON）*")
    else:
        lines.append("*（未检测到连续茎区）*")
    lines.append("")

    loops = data.get("loops") or []
    lines.append(f"#### 环区（Loops）— 共 **{len(loops)}** 个")
    lines.append("")
    if loops:
        for lp in loops[:20]:
            lines.append(
                f"- **Loop #{lp.get('id')}** ({lp.get('loop_type')}): "
                f"位置 [{lp.get('start')}–{lp.get('end')}]，大小 **{lp.get('size')}** nt"
            )
        if len(loops) > 20:
            lines.append("")
            lines.append(f"*（其余 {len(loops) - 20} 个环区见 JSON）*")
    else:
        lines.append("*（未识别独立环区）*")
    lines.append("")

    lines.append("#### 数据下载")
    lines.append("")
    if json_url:
        lines.append(f"- [完整 JSON 结果]({json_url})（含 stems/loops 全量字段）")
    else:
        lines.append("- （本次未生成可下载 JSON 链接）")
    lines.append("")
    lines.append(
        "> **说明**: 自由能与环区罚分基于仓库内简化参数；发表级预测请使用 ViennaRNA RNAfold / UNAFold 等专业工具交叉验证。"
    )
    return "\n".join(lines).strip()


def _md_embed_png_base64(b64: str, title: str) -> str:
    b64 = (b64 or "").strip()
    if not b64:
        return ""
    return f"\n#### {title}\n\n![{title}](data:image/png;base64,{b64})\n\n"


def _build_itc_markdown(data: Dict[str, Any]) -> str:
    lines: list[str] = ["### ITC 一结合位点模型（热力学拟合）", ""]
    if not data.get("converged"):
        lines.append(f"- **状态**: 未收敛 — {data.get('error', '未知错误')}")
        return "\n".join(lines).strip()
    lines.append(f"- **Kd** (M): `{data.get('Kd')}`（σ≈`{data.get('Kd_err')}`）")
    lines.append(f"- **ΔH** (J/mol): `{data.get('dH')}`（σ≈`{data.get('dH_err')}`）")
    lines.append(f"- **n**（位点化学计量）: `{data.get('n')}`（σ≈`{data.get('n_err')}`）")
    lines.append(f"- **ΔG** (J/mol): `{data.get('dG')}`")
    lines.append(f"- **ΔS** (J/(mol·K)): `{data.get('dS')}`")
    lines.append(f"- **R²**: `{data.get('r_squared')}`，**RMSE**: `{data.get('rmse')}`")
    lines.append("")
    lines.append(_md_embed_png_base64(str(data.get("fit_plot_png_base64") or ""), "差分热：观测 vs 拟合曲线"))
    return "\n".join(lines).strip()


def _build_cell_cycle_markdown(data: Dict[str, Any]) -> str:
    lines: list[str] = ["### 细胞周期各阶段时长估计（EdU/BrdU 双脉冲流式）", ""]
    if data.get("success") is False:
        lines.append(f"- **优化状态**: ⚠️ `{data.get('message', '')}`")
    lines.append(f"- **G1** (h): `{data.get('g1_duration')}`")
    lines.append(f"- **S** (h): `{data.get('s_duration')}`")
    lines.append(f"- **G2/M** (h): `{data.get('g2m_duration')}`")
    lines.append(f"- **总周期 Tc** (h): `{data.get('total_cycle_time')}`")
    lines.append(f"- **死亡率** (/h): `{data.get('death_rate')}`（≈ `{data.get('death_rate_percent')}` %/h）")
    pf = data.get("phase_fractions") or {}
    if pf:
        lines.append(
            f"- **相占比**: G1 `{pf.get('g1')}`，S `{pf.get('s')}`，G2/M `{pf.get('g2m')}`（分数，和为 1）"
        )
    fs = data.get("fit_statistics") or {}
    if fs:
        lines.append(
            f"- **R²**: EdU `{fs.get('r2_edu')}`，BrdU `{fs.get('r2_brdu')}`，双阳 `{fs.get('r2_double')}`，平均 `{fs.get('r2_average')}`"
        )
    lines.append(f"- **最终损失**: `{data.get('final_loss')}`，迭代: `{data.get('iterations')}`")
    lines.append("")
    lines.append("> 注：完整观测与预测曲线见脚本 JSON 字段 `observed_data` / `predictions`。")
    return "\n".join(lines).strip()


def _build_protease_markdown(data: Dict[str, Any]) -> str:
    lines: list[str] = ["### 蛋白酶动力学（Michaelis-Menten 拟合）", ""]
    if data.get("error"):
        lines.append(f"- **错误**: {data.get('error')}")
        return "\n".join(lines).strip()
    lines.append(f"- **Vmax** (a.u./s): `{data.get('Vmax')}`")
    lines.append(f"- **Km** (μM): `{data.get('Km')}`")
    lines.append(f"- **kcat** (s⁻¹): `{data.get('kcat')}`")
    lines.append(f"- **kcat/Km** (μM⁻¹·s⁻¹): `{data.get('kcat_over_Km')}`")
    lines.append(f"- **R²**: `{data.get('r_squared')}`")
    lines.append("")
    lines.append(
        _md_embed_png_base64(str(data.get("michaelis_menten_fit_plot_png_base64") or ""), "MM 拟合曲线（v–[S]）")
    )
    log = (data.get("research_log") or "").strip()
    if log:
        lines.append("#### 分析日志（节选）\n\n```text\n" + log[:8000] + "\n```\n")
    return "\n".join(lines).strip()


def _build_enzyme_markdown(data: Dict[str, Any]) -> str:
    lines: list[str] = ["### 酶动力学（时间历程 + MM 拟合）", ""]
    if data.get("error"):
        lines.append(f"- **错误**: {data.get('error')}")
        return "\n".join(lines).strip()
    lines.append(f"- **酶**: `{data.get('enzyme_name')}`")
    fp = data.get("fit_params") or {}
    if fp.get("Vmax") is not None:
        lines.append(f"- **Vmax** (a.u./min): `{fp.get('Vmax')}`；**Km** (μM): `{fp.get('Km')}`；**R²**: `{fp.get('r_squared')}`")
    lines.append("")
    lines.append(
        _md_embed_png_base64(
            str(data.get("substrate_kinetics_fit_plot_png_base64") or ""),
            "底物动力学与 MM 拟合曲线",
        )
    )
    log = (data.get("research_log") or "").strip()
    if log:
        lines.append("#### 分析日志（节选）\n\n```text\n" + log[:8000] + "\n```\n")
    return "\n".join(lines).strip()


def _build_sequence_conservation_markdown(data: Dict[str, Any], *, rel: str) -> str:
    lines: list[str] = ["### 蛋白质序列保守性（MSA / 距离 / 保守性剖面）", ""]
    if data.get("error"):
        lines.append(f"- **错误**: {data.get('error')}")
        return "\n".join(lines).strip()
    log = (data.get("research_log") or "").strip()
    if log:
        lines.append("#### 分析日志（节选）\n\n```text\n" + log[:6000] + "\n```\n")
    extras: list[str] = []
    for label, key in (
        ("多序列比对 FASTA", "alignment_file_path"),
        ("系统发育树 Newick", "tree_file_path"),
        ("保守性文本报告", "conservation_file_path"),
        ("保守性 JSON", "conservation_json_path"),
    ):
        p = data.get(key)
        if p and Path(str(p)).is_file():
            extras.append(f"- [{label}]({_public_path_url('/results/' + rel + '/' + Path(str(p)).name)})")
    if extras:
        lines.append("#### 下载\n\n" + "\n".join(extras))
    return "\n".join(lines).strip()


def _build_protein_homology_markdown(data: Dict[str, Any]) -> str:
    lines: list[str] = ["### 蛋白同源与结构/PTM 综合评估", ""]
    if data.get("error"):
        lines.append(f"- **错误**: {data.get('error')}")
        return "\n".join(lines).strip()
    tgt = data.get("target") or {}
    lines.append(f"- **靶蛋白**: `{tgt.get('uniprot_id')}` — {tgt.get('name', '—')}")
    lines.append(
        f"- **序列长度**: {tgt.get('sequence_length')} aa；**PDB 条目**: {tgt.get('pdb_count')}；"
        f"**PTM 位点**: {tgt.get('ptm_count')}"
    )
    if tgt.get("alphafold_model"):
        lines.append(f"- **AlphaFold 模型**: {tgt.get('alphafold_model')}")
    summ = data.get("summary") or {}
    if summ:
        lines.append(
            f"- **BLAST 命中**: {summ.get('total_hits')}；**通过阈值**: {summ.get('passing_hits')}；"
            f"**平均同一性**: {summ.get('mean_identity', 0):.1f}%"
        )
    lines.append("")
    hits = data.get("high_confidence_hits") or data.get("hits") or []
    if not hits:
        lines.append(data.get("message") or "*（无显著同源命中）*")
        return "\n".join(lines).strip()
    lines.append("#### Top 同源命中（综合评分）")
    lines.append("")
    lines.append("| UniProt | 综合分 | 同一性% | E-value | 结构相似 | PTM 保守 |")
    lines.append("| --- | ---: | ---: | --- | ---: | ---: |")
    for h in hits[:15]:
        bid = h.get("uniprot_id") or "—"
        comp = (h.get("composite_score") or {}).get("total", "—")
        bl = h.get("blast") or {}
        st = (h.get("structural_similarity") or {}).get("structural_similarity_score", "—")
        pt = (h.get("ptm_analysis") or {}).get("ptm_conservation_score", "—")
        lines.append(
            f"| `{bid}` | {comp} | {bl.get('identity', '—')} | {bl.get('evalue', '—')} | {st} | {pt} |"
        )
    if data.get("test_mode"):
        lines.append("")
        lines.append("> **说明**: 本次为 **test_mode** 模拟 BLAST 命中，未调用 NCBI BLAST 远程 API。")
    return "\n".join(lines).strip()


def _build_oasis_humanness_markdown(data: Dict[str, Any]) -> str:
    lines: list[str] = ["### OASis 抗体人源性评估", ""]
    lines.append(f"- **抗体**: {data.get('antibody_name', '—')}")
    lines.append(
        f"- **编号方案**: `{data.get('scheme')}`；**CDR 定义**: `{data.get('cdr_definition')}`；"
        f"**阈值**: `{data.get('threshold')}`（{data.get('threshold_description', '')}）"
    )
    if data.get("warning"):
        lines.append(f"- **注意**: {data.get('warning')}")
    hc = data.get("heavy_chain") or {}
    lc = data.get("light_chain") or {}
    ov = data.get("overall") or {}
    lines.append("")
    lines.append(f"- **重链 OASis**: {float(hc.get('oasis_identity', 0)):.2%}（{hc.get('length')} aa）")
    lines.append(f"- **轻链 OASis**: {float(lc.get('oasis_identity', 0)):.2%}（{lc.get('length')} aa）")
    if ov:
        lines.append(
            f"- **整体平均 (H+L)**: {float(ov.get('average_identity', 0)):.2%}；"
            f"合并序列: {float(ov.get('oasis_identity', 0)):.2%}"
        )
        if ov.get("cdr_identity") is not None:
            lines.append(f"- **CDR 区**: {float(ov['cdr_identity']):.2%}")
        if ov.get("framework_identity") is not None:
            lines.append(f"- **Framework**: {float(ov['framework_identity']):.2%}")
    avg = float(ov.get("average_identity", 0))
    lines.append("")
    if avg >= 0.85:
        lines.append("**解读**: 人源性较高，预测免疫原性风险相对较低。")
    elif avg >= 0.65:
        lines.append("**解读**: 人源性中等，可考虑人源化优化。")
    else:
        lines.append("**解读**: 人源性偏低，预测免疫原性风险较高。")
    return "\n".join(lines).strip()


@registry.register(
    name="immune_cell_isolation_simulation",
    description="模拟组织来源免疫细胞的酶解、分离、纯化与质控指标（文本/JSON 报告 + CSV）。",
    category="Biomedicine",
)
@safe_tool_execution
async def immune_cell_isolation_simulation(
    simulation_json: str,
    output_format: str = "text",
    random_seed: Optional[int] = None,
) -> Dict[str, Any]:
    """
    simulation_json: JSON，需含 tissue_type、target_cell_type；可选 enzyme_type、digestion_time_min 等。
    output_format: text 或 json。
    """
    script = _asset_script("immune_cell_isolation.py")
    if not script.is_file():
        return {"status": "error", "message": f"未找到脚本: {script}"}

    run_id = uuid.uuid4().hex[:16]
    out_dir = _results_dir() / "immune_cell" / run_id
    try:
        out_dir.mkdir(parents=True, exist_ok=False)
    except OSError as exc:
        return {"status": "error", "message": f"无法创建输出目录: {exc}"}

    csv_path = out_dir / "immune_cell_isolation.csv"
    fmt = (output_format or "text").strip().lower()
    if fmt not in ("text", "json"):
        fmt = "text"
    argv = [
        "--data",
        simulation_json.strip(),
        "--output-csv",
        str(csv_path),
        "--output-format",
        fmt,
    ]
    if random_seed is not None:
        argv.extend(["--seed", str(int(random_seed))])

    timeout = float(os.getenv("IMMUNE_CELL_SUBPROCESS_TIMEOUT", "120"))
    try:
        code, stdout, stderr = await _run_script(script, argv, timeout=timeout)
    except TimeoutError as exc:
        return {"status": "error", "message": str(exc)}
    if code != 0:
        return {"status": "error", "message": (stderr or stdout or f"进程退出码 {code}")[:8000]}

    rel = f"immune_cell/{run_id}"
    csv_url = _public_path_url(f"/results/{rel}/{csv_path.name}") if csv_path.is_file() else ""
    md = stdout or "（无标准输出）"
    if csv_url:
        md += f"\n\n#### 下载\n\n- [CSV 报告]({csv_url})"
    return {
        "status": "success",
        "message": "免疫细胞分离与纯化仿真已完成。",
        "markdown": md,
        "csv_url": csv_url,
        "output_dir": str(out_dir),
    }


@registry.register(
    name="rna_secondary_structure_analysis",
    description="基于点括号表示与简化 Turner 参数的 RNA 二级结构解析（碱基对、茎环、估计自由能等）。",
    category="Biomedicine",
)
@safe_tool_execution
async def rna_secondary_structure_analysis(
    structure_json: str = "",
    file_path: str = "",
    result_format: str = "json",
) -> Dict[str, Any]:
    """
    structure_json: 可选，JSON 字符串，含 structure 与可选 sequence。
    file_path: 可选，已上传的 .fa / 点括号文本文件绝对路径。
    二者至少填其一。
    result_format: json（默认）| txt | csv | all。
    """
    script = _asset_script("rna_structure.py")
    if not script.is_file():
        return {"status": "error", "message": f"未找到脚本: {script}"}

    sj = (structure_json or "").strip()
    fp = (file_path or "").strip()
    if not sj and not fp:
        return {"status": "error", "message": "请提供 structure_json（内联 JSON）或 file_path（上传文件路径）。"}

    run_id = uuid.uuid4().hex[:16]
    out_dir = _results_dir() / "rna_structure" / run_id
    try:
        out_dir.mkdir(parents=True, exist_ok=False)
    except OSError as exc:
        return {"status": "error", "message": f"无法创建输出目录: {exc}"}

    fmt = (result_format or "json").strip().lower()
    if fmt not in ("json", "txt", "csv", "all"):
        fmt = "json"

    out_txt = out_dir / "report.txt"
    out_csv = out_dir / "base_pairs.csv"
    argv: list[str] = [
        "--format",
        fmt,
        "--output-txt",
        str(out_txt),
        "--output-csv",
        str(out_csv),
    ]
    if fp:
        argv.extend(["--file", _resolve_uploaded_path(fp)])
    else:
        argv.extend(["--input", sj])

    timeout = float(os.getenv("RNA_STRUCTURE_SUBPROCESS_TIMEOUT", "120"))
    try:
        code, stdout, stderr = await _run_script(script, argv, timeout=timeout)
    except TimeoutError as exc:
        return {"status": "error", "message": str(exc)}
    if code != 0:
        return {"status": "error", "message": (stderr or stdout or f"进程退出码 {code}")[:8000]}

    rel = f"rna_structure/{run_id}"
    json_url = ""
    if fmt == "json" and stdout:
        jp = out_dir / "summary.json"
        try:
            jp.write_text(stdout, encoding="utf-8")
            json_url = _public_path_url(f"/results/{rel}/{jp.name}")
        except OSError:
            pass

    md = ""
    if fmt == "json" and stdout.strip():
        try:
            summary = json.loads(stdout)
            md = _build_rna_expert_markdown(summary, json_url=json_url)
        except json.JSONDecodeError:
            md = "### RNA 二级结构分析\n\n```json\n" + stdout[:120000] + "\n```"
    elif stdout.strip():
        md = "### RNA 二级结构分析\n\n```\n" + stdout[:120000] + "\n```"
    else:
        md = "### RNA 二级结构分析\n\n（无标准输出）"

    if out_txt.is_file():
        u = _public_path_url(f"/results/{rel}/{out_txt.name}")
        md += f"\n\n#### 附加下载\n\n- [文本报告（机器生成）]({u})"
    if out_csv.is_file():
        u2 = _public_path_url(f"/results/{rel}/{out_csv.name}")
        md += f"\n- [茎区汇总 CSV]({u2})"

    out: Dict[str, Any] = {
        "status": "success",
        "message": "RNA 二级结构分析已完成。",
        "markdown": md.strip(),
        "output_dir": str(out_dir),
    }
    if json_url:
        out["json_url"] = json_url
    if out_csv.is_file():
        out["csv_url"] = _public_path_url(f"/results/{rel}/{out_csv.name}")
    return out


@registry.register(
    name="circular_dichroism_analysis",
    description="圆二色谱（CD）数据分析：二级结构组分估计、谱特征与可选热变性曲线；输出报告与 CSV。",
    category="Biomedicine",
)
@safe_tool_execution
async def circular_dichroism_analysis(
    spectrum_json: str,
    output_format: str = "text",
) -> Dict[str, Any]:
    """
    spectrum_json: JSON，需含 sample_name、sample_type(protein|nucleic_acid)、wavelength_data、cd_signal_data 数组；
    可选 temperature_data + thermal_cd_data。
    output_format: text 或 json。
    """
    script = _asset_script("cd_analysis.py")
    if not script.is_file():
        return {"status": "error", "message": f"未找到脚本: {script}"}

    run_id = uuid.uuid4().hex[:16]
    out_dir = _results_dir() / "cd_analysis" / run_id
    try:
        out_dir.mkdir(parents=True, exist_ok=False)
    except OSError as exc:
        return {"status": "error", "message": f"无法创建输出目录: {exc}"}

    fmt = (output_format or "text").strip().lower()
    if fmt not in ("text", "json"):
        fmt = "text"
    argv = [
        "--data",
        spectrum_json.strip(),
        "--output-dir",
        str(out_dir),
        "--output-format",
        fmt,
    ]
    timeout = float(os.getenv("CD_ANALYSIS_SUBPROCESS_TIMEOUT", "120"))
    try:
        code, stdout, stderr = await _run_script(script, argv, timeout=timeout)
    except TimeoutError as exc:
        return {"status": "error", "message": str(exc)}
    if code != 0:
        return {"status": "error", "message": (stderr or stdout or f"进程退出码 {code}")[:8000]}

    rel = f"cd_analysis/{run_id}"
    md = stdout or "（无标准输出）"
    extras: list[str] = []
    cd_csv = out_dir / "cd_spectrum.csv"
    th_csv = out_dir / "thermal_denaturation.csv"
    if cd_csv.is_file():
        extras.append(f"- [CD 光谱 CSV]({_public_path_url(f'/results/{rel}/{cd_csv.name}')})")
    if th_csv.is_file():
        extras.append(f"- [热变性 CSV]({_public_path_url(f'/results/{rel}/{th_csv.name}')})")
    if extras:
        md += "\n\n#### 下载\n\n" + "\n".join(extras)

    out: Dict[str, Any] = {
        "status": "success",
        "message": "圆二色谱分析已完成。",
        "markdown": md,
        "output_dir": str(out_dir),
    }
    if cd_csv.is_file():
        out["csv_url"] = _public_path_url(f"/results/{rel}/{cd_csv.name}")
    return out


_DEFAULT_CELL_CYCLE_INITIAL_JSON = (
    '{"g1_duration":6,"s_duration":8,"g2m_duration":4,"death_rate":0.02}'
)


@registry.register(
    name="itc_binding_thermodynamic_analysis",
    description="ITC 一结合位点模型拟合：Kd、ΔH、ΔS、ΔG 等；输出拟合曲线 PNG（Base64 嵌入 Markdown）。",
    category="Biomedicine",
)
@safe_tool_execution
async def itc_binding_thermodynamic_analysis(
    protein_conc: float,
    ligand_conc: float,
    itc_injections_json: str = "",
    file_path: str = "",
    cell_volume: float = 1.4e-3,
    temperature: float = 25.0,
    volume_unit: str = "L",
) -> Dict[str, Any]:
    """
    itc_injections_json: JSON 数组 [[injection#, volume_L, heat], ...]（与 file_path 二选一）。
    file_path: 已上传 CSV/TSV（三列：injection, volume, heat）。
    """
    script = _asset_script("itc_analysis.py")
    if not script.is_file():
        return {"status": "error", "message": f"未找到脚本: {script}"}

    ij = (itc_injections_json or "").strip()
    fp = (file_path or "").strip()
    if not ij and not fp:
        return {
            "status": "error",
            "message": "请提供 itc_injections_json（内联滴定数据 JSON 数组）或 file_path（上传文件路径）。",
        }

    run_id = uuid.uuid4().hex[:16]
    out_dir = _results_dir() / "itc_analysis" / run_id
    try:
        out_dir.mkdir(parents=True, exist_ok=False)
    except OSError as exc:
        return {"status": "error", "message": f"无法创建输出目录: {exc}"}

    argv: list[str] = [
        "--output-format",
        "json",
        "--protein-conc",
        str(float(protein_conc)),
        "--ligand-conc",
        str(float(ligand_conc)),
        "--cell-volume",
        str(float(cell_volume)),
        "--temperature",
        str(float(temperature)),
        "--volume-unit",
        (volume_unit or "L").strip(),
    ]
    if fp:
        argv.extend(["--file", _resolve_uploaded_path(fp)])
    else:
        argv.extend(["--data", ij])

    timeout = float(os.getenv("ITC_ANALYSIS_SUBPROCESS_TIMEOUT", "180"))
    try:
        code, stdout, stderr = await _run_script(script, argv, timeout=timeout)
    except TimeoutError as exc:
        return {"status": "error", "message": str(exc)}
    if code != 0:
        return {"status": "error", "message": (stderr or stdout or f"进程退出码 {code}")[:8000]}

    rel = f"itc_analysis/{run_id}"
    md = ""
    try:
        payload = json.loads(stdout or "{}")
        md = _build_itc_markdown(payload)
    except json.JSONDecodeError:
        md = "### ITC 分析\n\n```\n" + (stdout or "")[:120000] + "\n```"

    return {
        "status": "success",
        "message": "ITC 热力学拟合已完成。",
        "markdown": md.strip(),
        "output_dir": str(out_dir),
        "results_rel": rel,
    }


@registry.register(
    name="cell_cycle_phase_duration_estimation",
    description="基于 EdU/BrdU 双脉冲流式数据估计 G1/S/G2M 时长与死亡率（scipy 优化）。",
    category="Biomedicine",
)
@safe_tool_execution
async def cell_cycle_phase_duration_estimation(
    flow_cytometry_json: str,
    initial_estimates_json: str = "",
) -> Dict[str, Any]:
    """
    flow_cytometry_json: 含 time_points、edu_positive、brdu_positive、double_positive（百分比序列，等长）。
    initial_estimates_json: 可选，含 g1_duration、s_duration、g2m_duration、death_rate；缺省使用仓库内默认初值 JSON。
    """
    script = _asset_script("cell_cycle_analysis.py")
    if not script.is_file():
        return {"status": "error", "message": f"未找到脚本: {script}"}

    data_s = (flow_cytometry_json or "").strip()
    if not data_s:
        return {"status": "error", "message": "flow_cytometry_json 不能为空。"}
    init_s = (initial_estimates_json or "").strip() or _DEFAULT_CELL_CYCLE_INITIAL_JSON

    run_id = uuid.uuid4().hex[:16]
    out_dir = _results_dir() / "cell_cycle" / run_id
    try:
        out_dir.mkdir(parents=True, exist_ok=False)
    except OSError as exc:
        return {"status": "error", "message": f"无法创建输出目录: {exc}"}

    argv = ["--data", data_s, "--initial", init_s, "--output-format", "json"]
    timeout = float(os.getenv("CELL_CYCLE_SUBPROCESS_TIMEOUT", "180"))
    try:
        code, stdout, stderr = await _run_script(script, argv, timeout=timeout)
    except TimeoutError as exc:
        return {"status": "error", "message": str(exc)}
    if code != 0:
        return {"status": "error", "message": (stderr or stdout or f"进程退出码 {code}")[:8000]}

    rel = f"cell_cycle/{run_id}"
    try:
        payload = json.loads(stdout or "{}")
        md = _build_cell_cycle_markdown(payload)
    except json.JSONDecodeError:
        md = "### 细胞周期估计\n\n```\n" + (stdout or "")[:120000] + "\n```"

    return {
        "status": "success",
        "message": "细胞周期相时长估计已完成。",
        "markdown": md.strip(),
        "output_dir": str(out_dir),
        "results_rel": rel,
    }


@registry.register(
    name="protease_kinetics_analysis",
    description="蛋白酶荧光动力学：由时间曲线求初速率并拟合 Michaelis-Menten；返回 MM 拟合图 Base64。",
    category="Biomedicine",
)
@safe_tool_execution
async def protease_kinetics_analysis(kinetics_input_json: str) -> Dict[str, Any]:
    """
    kinetics_input_json: 单行 JSON，需含 time_points、substrate_concentrations、fluorescence_data（二维）、enzyme_concentration（μM）。
    """
    script = _asset_script("protease_kinetics.py")
    if not script.is_file():
        return {"status": "error", "message": f"未找到脚本: {script}"}
    js = (kinetics_input_json or "").strip()
    if not js:
        return {"status": "error", "message": "kinetics_input_json 不能为空。"}

    run_id = uuid.uuid4().hex[:16]
    out_dir = _results_dir() / "protease_kinetics" / run_id
    try:
        out_dir.mkdir(parents=True, exist_ok=False)
    except OSError as exc:
        return {"status": "error", "message": f"无法创建输出目录: {exc}"}

    argv = ["--json-input", js, "--output-dir", str(out_dir)]
    timeout = float(os.getenv("PROTEASE_KINETICS_SUBPROCESS_TIMEOUT", "180"))
    try:
        code, stdout, stderr = await _run_script(script, argv, timeout=timeout)
    except TimeoutError as exc:
        return {"status": "error", "message": str(exc)}
    if code != 0:
        return {"status": "error", "message": (stderr or stdout or f"进程退出码 {code}")[:8000]}

    rel = f"protease_kinetics/{run_id}"
    try:
        payload = json.loads(stdout or "{}")
        md = _build_protease_markdown(payload)
    except json.JSONDecodeError:
        md = "### 蛋白酶动力学\n\n```\n" + (stdout or "")[:120000] + "\n```"

    extras: list[str] = []
    try:
        pobj = json.loads(stdout or "{}")
        for label, key in (
            ("MM 图 PNG", "mm_plot_path"),
            ("时间历程 PNG", "timecourse_plot_path"),
            ("文本报告", "report_path"),
            ("结果 JSON", "json_path"),
        ):
            p = pobj.get(key)
            if p and Path(str(p)).is_file():
                extras.append(f"- [{label}]({_public_path_url('/results/' + rel + '/' + Path(str(p)).name)})")
    except Exception:
        pass
    if extras:
        md += "\n\n#### 下载\n\n" + "\n".join(extras)

    return {
        "status": "success",
        "message": "蛋白酶动力学分析已完成。",
        "markdown": md.strip(),
        "output_dir": str(out_dir),
        "results_rel": rel,
    }


@registry.register(
    name="enzyme_kinetics_analysis",
    description="酶动力学仿真：时间历程、底物饱和曲线与 MM 参数拟合；返回底物动力学拟合图 Base64。",
    category="Biomedicine",
)
@safe_tool_execution
async def enzyme_kinetics_analysis(kinetics_input_json: str) -> Dict[str, Any]:
    """
    kinetics_input_json: 单行 JSON，需含 enzyme_name、substrate_concentrations（μM 列表）、enzyme_concentration（nM）；可选 time_points、modulators 等。
    """
    script = _asset_script("enzyme_kinetics.py")
    if not script.is_file():
        return {"status": "error", "message": f"未找到脚本: {script}"}
    js = (kinetics_input_json or "").strip()
    if not js:
        return {"status": "error", "message": "kinetics_input_json 不能为空。"}

    run_id = uuid.uuid4().hex[:16]
    out_dir = _results_dir() / "enzyme_kinetics" / run_id
    try:
        out_dir.mkdir(parents=True, exist_ok=False)
    except OSError as exc:
        return {"status": "error", "message": f"无法创建输出目录: {exc}"}

    argv = ["--json-input", js, "--output-dir", str(out_dir)]
    timeout = float(os.getenv("ENZYME_KINETICS_SUBPROCESS_TIMEOUT", "180"))
    try:
        code, stdout, stderr = await _run_script(script, argv, timeout=timeout)
    except TimeoutError as exc:
        return {"status": "error", "message": str(exc)}
    if code != 0:
        return {"status": "error", "message": (stderr or stdout or f"进程退出码 {code}")[:8000]}

    rel = f"enzyme_kinetics/{run_id}"
    try:
        payload = json.loads(stdout or "{}")
        md = _build_enzyme_markdown(payload)
    except json.JSONDecodeError:
        md = "### 酶动力学\n\n```\n" + (stdout or "")[:120000] + "\n```"

    extras: list[str] = []
    try:
        pobj = json.loads(stdout or "{}")
        ofs = pobj.get("output_files") or {}
        for label, key in (
            ("底物动力学 PNG", "substrate_kinetics_plot"),
            ("时间历程 PNG", "time_course_plot"),
            ("时间历程 CSV", "time_course_csv"),
            ("底物动力学 CSV", "substrate_kinetics_csv"),
            ("文本报告", "report"),
        ):
            p = ofs.get(key)
            if p and Path(str(p)).is_file():
                extras.append(f"- [{label}]({_public_path_url('/results/' + rel + '/' + Path(str(p)).name)})")
    except Exception:
        pass
    if extras:
        md += "\n\n#### 下载\n\n" + "\n".join(extras)

    return {
        "status": "success",
        "message": "酶动力学分析已完成。",
        "markdown": md.strip(),
        "output_dir": str(out_dir),
        "results_rel": rel,
    }


@registry.register(
    name="protein_sequence_conservation_analysis",
    description="蛋白质多序列比对、距离矩阵、系统发育树与位点保守性（Biopython + scipy）。",
    category="Biomedicine",
)
@safe_tool_execution
async def protein_sequence_conservation_analysis(alignment_json: str) -> Dict[str, Any]:
    """
    alignment_json: 单行 JSON，需含 protein_sequences（字符串数组或带 > 头的 FASTA 串列表）。
    """
    script = _asset_script("sequence_conservation.py")
    if not script.is_file():
        return {"status": "error", "message": f"未找到脚本: {script}"}
    js = (alignment_json or "").strip()
    if not js:
        return {"status": "error", "message": "alignment_json 不能为空。"}

    run_id = uuid.uuid4().hex[:16]
    out_dir = _results_dir() / "sequence_conservation" / run_id
    try:
        out_dir.mkdir(parents=True, exist_ok=False)
    except OSError as exc:
        return {"status": "error", "message": f"无法创建输出目录: {exc}"}

    argv = ["--json-input", js, "--output-dir", str(out_dir)]
    timeout = float(os.getenv("SEQUENCE_CONSERVATION_SUBPROCESS_TIMEOUT", "300"))
    try:
        code, stdout, stderr = await _run_script(script, argv, timeout=timeout)
    except TimeoutError as exc:
        return {"status": "error", "message": str(exc)}
    if code != 0:
        return {"status": "error", "message": (stderr or stdout or f"进程退出码 {code}")[:8000]}

    rel = f"sequence_conservation/{run_id}"
    try:
        payload = json.loads(stdout or "{}")
        md = _build_sequence_conservation_markdown(payload, rel=rel)
    except json.JSONDecodeError:
        md = "### 序列保守性\n\n```\n" + (stdout or "")[:120000] + "\n```"

    return {
        "status": "success",
        "message": "蛋白质序列保守性分析已完成。",
        "markdown": md.strip(),
        "output_dir": str(out_dir),
        "results_rel": rel,
    }


@registry.register(
    name="protein_homology_structure_assessment",
    description="UniProt 靶蛋白 BLAST 同源检索 + 结构相似性与 PTM 保守性综合评分（可选 test_mode 离线演示）。",
    category="Biomedicine",
)
@safe_tool_execution
async def protein_homology_structure_assessment(
    uniprot_id: str,
    blast_evalue: float = 1e-5,
    struct_threshold: float = 0.8,
    ptm_threshold: float = 0.6,
    top_n: int = 10,
    test_mode: bool = False,
) -> Dict[str, Any]:
    """
    uniprot_id: 靶蛋白 UniProt 登录号（如 O95292）。
    test_mode: True 时使用脚本内模拟 BLAST 命中（无 NCBI 远程调用，适合演示/无网环境）。
    """
    script = _asset_script("protein_assess.py")
    if not script.is_file():
        return {"status": "error", "message": f"未找到脚本: {script}"}
    uid = (uniprot_id or "").strip()
    if not uid:
        return {"status": "error", "message": "uniprot_id 不能为空。"}

    run_id = uuid.uuid4().hex[:16]
    out_dir = _results_dir() / "protein_assess" / run_id
    try:
        out_dir.mkdir(parents=True, exist_ok=False)
    except OSError as exc:
        return {"status": "error", "message": f"无法创建输出目录: {exc}"}

    out_json = out_dir / "assessment.json"
    argv = [
        "--uniprot-id",
        uid,
        "--blast-evalue",
        str(float(blast_evalue)),
        "--struct-threshold",
        str(float(struct_threshold)),
        "--ptm-threshold",
        str(float(ptm_threshold)),
        "--top-n",
        str(int(top_n)),
        "--output",
        str(out_json),
    ]
    if test_mode:
        argv.append("--test-mode")

    timeout = float(os.getenv("PROTEIN_ASSESS_SUBPROCESS_TIMEOUT", "600"))
    try:
        code, stdout, stderr = await _run_script(script, argv, timeout=timeout)
    except TimeoutError as exc:
        return {"status": "error", "message": str(exc)}
    if code != 0:
        return {"status": "error", "message": (stderr or stdout or f"进程退出码 {code}")[:8000]}

    rel = f"protein_assess/{run_id}"
    raw = ""
    if out_json.is_file():
        raw = out_json.read_text(encoding="utf-8", errors="replace")
    elif stdout.strip():
        raw = stdout
    try:
        payload = json.loads(raw or "{}")
        if test_mode and isinstance(payload, dict):
            payload["test_mode"] = True
        md = _build_protein_homology_markdown(payload)
    except json.JSONDecodeError:
        md = "### 蛋白同源评估\n\n```\n" + (raw or stdout or "")[:120000] + "\n```"

    json_url = ""
    if out_json.is_file():
        json_url = _public_path_url(f"/results/{rel}/{out_json.name}")
        md += f"\n\n#### 下载\n\n- [完整 JSON]({json_url})"

    return {
        "status": "success",
        "message": "蛋白同源结构评估已完成。",
        "markdown": md.strip(),
        "output_dir": str(out_dir),
        "results_rel": rel,
        "json_url": json_url,
    }


def _resolve_antibody_chains(
    heavy_chain: str,
    light_chain: str,
    antibody_sequences_json: str,
) -> tuple[str, str, Optional[str]]:
    """
    解析 H/L 链序列。返回 (hc, lc, error_message)。
    antibody_sequences_json 优先于分字段（便于 LLM 一次输出 JSON）。
    """
    hc = (heavy_chain or "").strip().replace(" ", "").upper()
    lc = (light_chain or "").strip().replace(" ", "").upper()
    js = (antibody_sequences_json or "").strip()
    if js:
        try:
            obj = json.loads(js)
            if not hc:
                hc = str(
                    obj.get("heavy_chain") or obj.get("heavy") or obj.get("H") or ""
                ).strip().replace(" ", "").upper()
            if not lc:
                lc = str(
                    obj.get("light_chain") or obj.get("light") or obj.get("L") or ""
                ).strip().replace(" ", "").upper()
        except json.JSONDecodeError as exc:
            return "", "", f"antibody_sequences_json 不是合法 JSON: {exc}"
    if not hc or not lc:
        return (
            "",
            "",
            "heavy_chain 与 light_chain 均不能为空。"
            "请提供两条氨基酸序列，或使用 antibody_sequences_json："
            '{"heavy_chain":"EVQLV...","light_chain":"DIQMT..."}',
        )
    return hc, lc, None


@registry.register(
    name="antibody_humanness_oasis_evaluation",
    description="基于 OAS（promb human-oas）与可选 ANARCI 的抗体链人源性（OASis）评分。",
    category="Biomedicine",
)
@safe_tool_execution
async def antibody_humanness_oasis_evaluation(
    heavy_chain: str = "",
    light_chain: str = "",
    antibody_sequences_json: str = "",
    antibody_name: str = "Antibody",
    scheme: str = "kabat",
    cdr_definition: str = "kabat",
    threshold: str = "relaxed",
) -> Dict[str, Any]:
    """
    heavy_chain / light_chain: 抗体可变区氨基酸序列（单字母）；可与 antibody_sequences_json 二选一或互补。
    antibody_sequences_json: 单行 JSON，含 heavy_chain、light_chain（或 heavy/light、H/L）。
    scheme: kabat|chothia|imgt|aho|martin；cdr_definition: kabat|chothia|imgt|north；threshold: loose|relaxed|medium|strict。
    依赖：pip 包 promb；可选 anarci + 系统 hmmer（用于 CDR/框架分区）。
    """
    script = _asset_script("oasis_analysis.py")
    if not script.is_file():
        return {"status": "error", "message": f"未找到脚本: {script}"}
    hc, lc, err = _resolve_antibody_chains(heavy_chain, light_chain, antibody_sequences_json)
    if err:
        return {"status": "error", "message": err}

    run_id = uuid.uuid4().hex[:16]
    out_dir = _results_dir() / "oasis_humanness" / run_id
    try:
        out_dir.mkdir(parents=True, exist_ok=False)
    except OSError as exc:
        return {"status": "error", "message": f"无法创建输出目录: {exc}"}

    argv = [
        "--heavy-chain",
        hc,
        "--light-chain",
        lc,
        "--name",
        (antibody_name or "Antibody").strip(),
        "--scheme",
        (scheme or "kabat").strip(),
        "--cdr-definition",
        (cdr_definition or "kabat").strip(),
        "--threshold",
        (threshold or "relaxed").strip(),
        "--output-format",
        "json",
    ]
    timeout = float(os.getenv("OASIS_HUMANNESS_SUBPROCESS_TIMEOUT", "300"))
    try:
        code, stdout, stderr = await _run_script(script, argv, timeout=timeout)
    except TimeoutError as exc:
        return {"status": "error", "message": str(exc)}
    if code != 0:
        err = (stderr or stdout or f"进程退出码 {code}")[:8000]
        if "promb" in err.lower() or "No module named" in err:
            err += (
                "\n\n请在 API 容器内安装：pip install promb；"
                "可选 anarci（需 apt 安装 hmmer）。或重建镜像：docker compose build --no-cache api-server"
            )
        return {"status": "error", "message": err}

    rel = f"oasis_humanness/{run_id}"
    try:
        payload = json.loads(stdout or "{}")
        md = _build_oasis_humanness_markdown(payload)
    except json.JSONDecodeError:
        md = "### OASis 人源性\n\n```\n" + (stdout or "")[:120000] + "\n```"

    jp = out_dir / "oasis_result.json"
    try:
        jp.write_text(stdout or "{}", encoding="utf-8")
        json_url = _public_path_url(f"/results/{rel}/{jp.name}")
        md += f"\n\n#### 下载\n\n- [完整 JSON]({json_url})"
    except OSError:
        json_url = ""

    return {
        "status": "success",
        "message": "抗体人源性（OASis）评估已完成。",
        "markdown": md.strip(),
        "output_dir": str(out_dir),
        "results_rel": rel,
        "json_url": json_url,
    }
