# -*- coding: utf-8 -*-
"""
药物化学 / RDKit 原子工具：子进程调用 ``gibh_agent/assets/chem_rdkit/run_chem_tool.py``（或暂存区
``uploads/skills_staging/.../rdkit/run_chem_tool.py`` / ``CHEM_RDKIT_RUNNER``），
成功返回 ``markdown`` + ``image_urls`` + ``json_url``（与 DRT 右栏标准一致）。
"""
from __future__ import annotations

import json
import logging
import os
import subprocess
import sys
import uuid
from html import escape
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

from gibh_agent.core.tool_registry import registry
from gibh_agent.core.utils import safe_tool_execution, sanitize_for_json

logger = logging.getLogger(__name__)


def _repo_root() -> Path:
    return Path(__file__).resolve().parents[2]


def _resolve_chem_runner_script() -> Optional[Path]:
    env = (os.getenv("CHEM_RDKIT_RUNNER") or "").strip()
    if env:
        p = Path(env).expanduser().resolve()
        if p.is_file():
            return p
        logger.warning("CHEM_RDKIT_RUNNER 已设置但不是文件: %s", p)
    root = _repo_root()
    staging = (
        root
        / "uploads"
        / "skills_staging"
        / "2026.04.20 skills"
        / "rdkit"
        / "run_chem_tool.py"
    )
    if staging.is_file():
        return staging.resolve()
    bundled = root / "gibh_agent" / "assets" / "chem_rdkit" / "run_chem_tool.py"
    if bundled.is_file():
        return bundled.resolve()
    return None


def _pick_chem_python() -> str:
    py = (os.getenv("CHEM_RDKIT_PYTHON") or "").strip()
    return py if py else sys.executable


def _prepare_out_dir(tag: str) -> Path:
    results = Path(os.getenv("RESULTS_DIR", "/app/results")).expanduser().resolve()
    rid = uuid.uuid4().hex[:12]
    out = results / "chem_rdkit" / f"{tag}_{rid}"
    out.mkdir(parents=True, exist_ok=True)
    return out


def _collect_artifacts(out_dir: Path) -> Dict[str, List[str]]:
    json_paths: List[str] = []
    image_paths: List[str] = []
    if not out_dir.is_dir():
        return {"json_paths": json_paths, "image_paths": image_paths}
    for p in sorted(out_dir.rglob("*")):
        if not p.is_file():
            continue
        suf = p.suffix.lower()
        if suf == ".json":
            json_paths.append(str(p.resolve()))
        elif suf in (".png", ".jpg", ".jpeg", ".svg", ".pdf"):
            image_paths.append(str(p.resolve()))
    return {"json_paths": json_paths, "image_paths": image_paths}


def _artifact_abs_path_to_browser_url(abs_path: str) -> str:
    p = Path(abs_path).resolve()
    s = str(p)
    try:
        rdir = Path(os.getenv("RESULTS_DIR", "/app/results")).expanduser().resolve()
        rel = p.relative_to(rdir)
        return "/results/" + rel.as_posix()
    except ValueError:
        pass
    try:
        udir = Path(os.getenv("UPLOAD_DIR", "/app/uploads")).expanduser().resolve()
        rel = p.relative_to(udir)
        return "/uploads/" + rel.as_posix()
    except ValueError:
        pass
    if "/app/results/" in s:
        return "/results/" + s.split("/app/results/", 1)[-1].lstrip("/")
    if "/app/uploads/" in s:
        return "/uploads/" + s.split("/app/uploads/", 1)[-1].lstrip("/")
    return s


def _public_urls_from_artifacts(artifacts: Dict[str, List[str]]) -> Tuple[List[str], Optional[str]]:
    urls: List[str] = []
    for ip in artifacts.get("image_paths") or []:
        u = _artifact_abs_path_to_browser_url(ip)
        if u not in urls:
            urls.append(u)
    jp = (artifacts.get("json_paths") or [])[:1]
    json_u = _artifact_abs_path_to_browser_url(jp[0]) if jp else None
    return urls, json_u


def _safe_load_chem_summary_json(out_dir: Path, max_bytes: int = 2_000_000) -> Optional[Dict[str, Any]]:
    """
    仅从本工具创建的输出目录读取 chem_summary.json，限制体积，避免异常 JSON 拖垮进程。
    """
    p = out_dir / "chem_summary.json"
    if not p.is_file():
        return None
    try:
        raw = p.read_bytes()
        if len(raw) > max_bytes:
            logger.warning("chem_summary.json 过大 (%s bytes)，跳过摘要嵌入", len(raw))
            return None
        obj = json.loads(raw.decode("utf-8"))
        if not isinstance(obj, dict):
            return None
        return obj
    except (OSError, UnicodeDecodeError, json.JSONDecodeError, MemoryError) as e:
        logger.warning("读取 chem_summary.json 失败: %s", e)
        return None


def _html_chem_summary_for_tool(cli_tool: str, summary: Dict[str, Any]) -> str:
    """根据子进程 JSON 生成置于 2D 图之上的内联样式 HTML（Markdown 中可混排）。"""
    res = summary.get("result")
    if not isinstance(res, dict):
        return ""

    if cli_tool == "pains":
        n = int(res.get("hit_count") or 0)
        if n == 0:
            return (
                '<p style="margin:12px 0 8px 0;">'
                '<span style="background:#def7ec; color:#03543f; padding:2px 8px; border-radius:12px; '
                "font-weight:bold; font-size:12px;"
                ">✅ PAINS 检测通过 (0 警报)</span></p>"
            )
        parts: List[str] = []
        for h in (res.get("hits") or [])[:30]:
            if isinstance(h, dict):
                name = (h.get("name") or h.get("description") or "").strip()
            else:
                name = str(h).strip()
            if name:
                parts.append(escape(name))
        frag = "、".join(parts) if parts else escape("（未解析到片段名称）")
        return (
            '<p style="margin:12px 0 8px 0;">'
            '<span style="background:#fde8e8; color:#9b1c1c; padding:2px 8px; border-radius:12px; '
            "font-weight:bold; font-size:12px;"
            f">⚠️ PAINS 未通过 ({n} 条警报)</span></p>"
            '<p style="margin:4px 0 12px 0; font-size:13px; color:#374151; line-height:1.6;">'
            f"<strong>匹配结构片段：</strong>{frag}</p>"
        )

    if cli_tool == "brenk":
        n = int(res.get("hit_count") or 0)
        if n == 0:
            return (
                '<p style="margin:12px 0 8px 0;">'
                '<span style="background:#def7ec; color:#03543f; padding:2px 8px; border-radius:12px; '
                "font-weight:bold; font-size:12px;"
                ">✅ Brenk 筛查通过 (0 毒性/不良片段)</span></p>"
            )
        parts = []
        for h in (res.get("hits") or [])[:30]:
            if isinstance(h, dict):
                d = (h.get("description") or h.get("name") or "").strip()
            else:
                d = str(h).strip()
            if d:
                parts.append(escape(d))
        frag = "、".join(parts) if parts else escape("（未解析到片段）")
        return (
            '<p style="margin:12px 0 8px 0;">'
            '<span style="background:#fde8e8; color:#9b1c1c; padding:2px 8px; border-radius:12px; '
            "font-weight:bold; font-size:12px;"
            f">⚠️ Brenk 未通过 ({n} 条警报)</span></p>"
            '<p style="margin:4px 0 12px 0; font-size:13px; color:#374151; line-height:1.6;">'
            f"<strong>Toxic Fragments：</strong>{frag}</p>"
        )

    if cli_tool == "bbb":
        mw = res.get("MW")
        logp = res.get("MolLogP")
        tpsa = res.get("TPSA")
        label = str(res.get("summary_label") or "").strip()
        score = res.get("bbb_heuristic_score")
        try:
            sc = float(score) if score is not None else -1.0
        except (TypeError, ValueError):
            sc = -1.0
        # 与 run_chem_tool 启发式一致：高分或文案含「利于 CNS」视为偏正面
        positive = sc >= 0.65 or ("利于" in label and "CNS" in label)
        lbl_color = "#03543f" if positive else "#b45309"
        lbl_html = escape(label) if label else escape("—")
        def _fmt_num(x: Any, nd: int = 2) -> str:
            try:
                if x is None:
                    return "—"
                return escape(f"{float(x):.{nd}f}")
            except (TypeError, ValueError):
                return escape(str(x))

        mw_s = _fmt_num(mw, 2)
        lp_s = _fmt_num(logp, 2)
        tp_s = _fmt_num(tpsa, 1)
        return (
            '<div style="margin:12px 0 14px 0; padding:12px 14px; background:#f8fafc; '
            'border:1px solid #e2e8f0; border-radius:10px; max-width:520px;">'
            '<p style="margin:0 0 10px 0; font-size:13px; color:#64748b;">'
            "<strong>BBB（启发式）核心理化参数</strong></p>"
            '<table style="width:100%; border-collapse:collapse; font-size:13px; color:#1e293b;">'
            "<thead><tr>"
            '<th style="text-align:left; padding:6px 8px; border-bottom:1px solid #e2e8f0;">MW</th>'
            '<th style="text-align:left; padding:6px 8px; border-bottom:1px solid #e2e8f0;">LogP</th>'
            '<th style="text-align:left; padding:6px 8px; border-bottom:1px solid #e2e8f0;">'
            "TPSA (Å²)</th>"
            "</tr></thead><tbody><tr>"
            f'<td style="padding:8px; border-bottom:1px solid #f1f5f9;">{mw_s}</td>'
            f'<td style="padding:8px; border-bottom:1px solid #f1f5f9;">{lp_s}</td>'
            f'<td style="padding:8px; border-bottom:1px solid #f1f5f9;">{tp_s}</td>'
            "</tr></tbody></table>"
            f'<p style="margin:12px 0 0 0; font-size:14px; line-height:1.5;">'
            f"<strong>预测结论：</strong>"
            f'<span style="color:{lbl_color}; font-weight:700;">{lbl_html}</span>'
            "</p>"
            "</div>"
        )

    if cli_tool == "morgan_fp":
        r = int(res.get("radius") or 2)
        nb = int(res.get("n_bits") or 2048)
        nob = int(res.get("num_on_bits") or 0)
        return (
            f'<p style="margin:12px 0 10px 0; font-size:14px; color:#1f2937; line-height:1.55;">'
            f"🧬 Morgan 指纹 (Radius={r}, {nb}-bit)：共提取到 "
            f"<strong>{nob}</strong> 个激活的分子特征片段。</p>"
        )

    return ""


def _markdown_core_metrics_table(cli_tool: str, summary: Optional[Dict[str, Any]]) -> str:
    """
    从 chem_summary.json 结构生成 **Markdown 表格**（便于打印时可见，不仅依赖下载链接）。
    """
    if not summary or not isinstance(summary, dict):
        return ""
    res = summary.get("result")
    if not isinstance(res, dict):
        return ""
    tool = (summary.get("tool") or cli_tool or "").strip().lower()

    def _cell(v: Any) -> str:
        if v is None:
            return "—"
        s = str(v).strip()
        if not s:
            return "—"
        return s.replace("|", "\\|")

    rows: List[Tuple[str, str]] = []

    if tool == "pains":
        rows.append(("PAINS 警报数", _cell(res.get("hit_count"))))
    elif tool == "brenk":
        rows.append(("Brenk 警报数", _cell(res.get("hit_count"))))
    elif tool == "bbb":
        rows.extend(
            [
                ("MW (Da)", _cell(res.get("MW"))),
                ("MolLogP", _cell(res.get("MolLogP"))),
                ("TPSA (Å²)", _cell(res.get("TPSA"))),
                ("HBD", _cell(res.get("HBD"))),
                ("可旋转键数", _cell(res.get("NumRotatableBonds"))),
                ("BBB 启发式得分 (0–1)", _cell(res.get("bbb_heuristic_score"))),
                ("预测结论", _cell(res.get("summary_label"))),
            ]
        )
    elif tool == "morgan_fp":
        rows.extend(
            [
                ("Morgan 半径", _cell(res.get("radius"))),
                ("指纹位数", _cell(res.get("n_bits"))),
                ("激活位数", _cell(res.get("num_on_bits"))),
                ("位密度", _cell(res.get("density"))),
            ]
        )
    elif tool == "mol_weight":
        rows.extend(
            [
                ("精确分子量", _cell(res.get("exact_molecular_weight"))),
                ("平均分子量", _cell(res.get("average_molecular_weight"))),
            ]
        )
    elif tool == "tanimoto":
        rows.extend(
            [
                ("Tanimoto 相似度", _cell(res.get("tanimoto_similarity"))),
                ("距离 (1−Sim)", _cell(res.get("tanimoto_distance"))),
                ("指纹", _cell(res.get("fingerprint"))),
            ]
        )
    else:
        return ""

    canon = summary.get("canonical_smiles")
    if isinstance(canon, str) and canon.strip():
        rows.insert(0, ("Canonical SMILES（截断展示）", _cell(canon[:400] + ("…" if len(canon) > 400 else ""))))

    lines = ["| 指标 | 数值 |", "| --- | --- |"]
    for k, v in rows:
        lines.append(f"| {k} | {v} |")
    return "\n".join(["#### 核心数值摘要", "", "\n".join(lines), ""])


def _markdown_chem_report(
    title: str,
    image_urls: List[str],
    json_url: Optional[str],
    image_labels: Optional[Tuple[str, ...]] = None,
    html_summary: Optional[str] = None,
    *,
    core_metrics_md: Optional[str] = None,
) -> str:
    chunks: List[str] = [f"### {title}", "", "以下为 **RDKit** 计算结果（同源 `/results/` 静态资源）。", ""]
    link_parts: List[str] = []
    if json_url:
        link_parts.append(f"数值摘要 JSON：[下载/打开]({json_url})（`chem_summary.json`）")
    for idx, u in enumerate(image_urls):
        lab = image_labels[idx] if image_labels and idx < len(image_labels) else f"结果图 {idx + 1}"
        link_parts.append(f"{lab}：[下载]({u})")
    if link_parts:
        chunks.append("- " + " · ".join(link_parts))
        chunks.append("")
    if core_metrics_md and core_metrics_md.strip():
        chunks.append(core_metrics_md.strip())
        chunks.append("")
    if html_summary:
        chunks.append(html_summary.strip())
        chunks.append("")
    for image_url in image_urls:
        src = escape(image_url, quote=True)
        chunks.append(
            '<div style="max-height: 450px; width: 100%; overflow: hidden; '
            'border: 1px solid #e5e7eb; border-radius: 8px; display: flex; '
            "justify-content: center; align-items: center; background-color: #f8fafc; "
            f'margin-top: 16px;">\n'
            f'    <img src="{src}" alt="RDKit Structure" '
            'style="max-width: 100%; max-height: 450px; object-fit: contain; cursor: zoom-in;" '
            'onclick="window.open(this.src, \'_blank\')" />\n'
            "</div>\n"
            '<p style="text-align: center; font-size: 12px; color: #6b7280; margin-top: 8px;">'
            "点击图片可放大查看原图</p>"
        )
        chunks.append("")
    return "\n".join(chunks).strip()


def _validate_single_smiles_input(file_path: str, smiles_text: str) -> Optional[str]:
    fp = (file_path or "").strip()
    if fp and Path(fp).expanduser().is_file():
        return None
    if (smiles_text or "").strip():
        return None
    return "请提供 **smiles_text**（单行 SMILES）或已上传文件的 **file_path**。"


def _run_chem_tool(
    cli_tool: str,
    out_tag: str,
    md_title: str,
    image_labels: Tuple[str, ...],
    *,
    file_path: str = "",
    smiles_text: str = "",
    smiles2_text: str = "",
    timeout_seconds: int = 240,
) -> Dict[str, Any]:
    script = _resolve_chem_runner_script()
    if script is None:
        root = _repo_root()
        return {
            "status": "error",
            "message": (
                "未找到 run_chem_tool.py。请确认已包含 "
                f"{root / 'gibh_agent/assets/chem_rdkit/run_chem_tool.py'} "
                "或设置 CHEM_RDKIT_RUNNER 指向脚本绝对路径。"
            ),
        }

    if cli_tool != "tanimoto":
        ve = _validate_single_smiles_input(file_path, smiles_text)
        if ve:
            return {"status": "error", "message": ve}
    else:
        has_file = (file_path or "").strip() and Path((file_path or "").strip()).expanduser().is_file()
        st = (smiles_text or "").strip()
        s2 = (smiles2_text or "").strip()
        if not has_file and ((not st or not s2) and "|" not in st):
            return {
                "status": "error",
                "message": (
                    "Tanimoto 需要：**file_path**（两行 SMILES）、或 **smiles_text + smiles2_text**、"
                    "或在 **smiles_text** 内用 `|` 分隔两条 SMILES。"
                ),
            }

    out_dir = _prepare_out_dir(out_tag)
    py = _pick_chem_python()
    cmd: List[str] = [
        py,
        str(script),
        "--tool",
        cli_tool,
        "--output-dir",
        str(out_dir),
    ]
    if (file_path or "").strip():
        cmd += ["--file", str(Path(file_path.strip()).expanduser().resolve())]
    if (smiles_text or "").strip():
        cmd += ["--smiles", smiles_text.strip()]
    if cli_tool == "tanimoto" and (smiles2_text or "").strip():
        cmd += ["--smiles2", smiles2_text.strip()]

    logger.info("chem_rdkit subprocess: %s", " ".join(cmd[:10]))

    try:
        proc = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=max(45, int(timeout_seconds)),
            env={**os.environ, "PYTHONUNBUFFERED": "1"},
            cwd=str(out_dir),
        )
    except subprocess.TimeoutExpired:
        return {"status": "error", "message": f"RDKit 子进程超时（>{timeout_seconds}s）。"}
    except OSError as e:
        return {"status": "error", "message": f"无法启动子进程: {e}"}

    stdout_tail = (proc.stdout or "")[-6000:]
    stderr_tail = (proc.stderr or "")[-6000:]
    artifacts = _collect_artifacts(out_dir)

    if proc.returncode != 0:
        hint = (stderr_tail or stdout_tail or "").strip()[-800:]
        return {
            "status": "error",
            "message": f"RDKit 脚本失败（exit {proc.returncode}）。{hint}",
            "data": sanitize_for_json(
                {
                    "script_path": str(script),
                    "output_dir": str(out_dir.resolve()),
                    "stdout_tail": stdout_tail,
                    "stderr_tail": stderr_tail,
                    **artifacts,
                }
            ),
        }

    image_urls, json_url = _public_urls_from_artifacts(artifacts)
    labs = image_labels
    if cli_tool == "tanimoto" and len(image_urls) == 2:
        labs = ("分子 A（2D）", "分子 B（2D）")

    summary_obj = _safe_load_chem_summary_json(out_dir)
    html_summary = ""
    if summary_obj:
        try:
            html_summary = _html_chem_summary_for_tool(cli_tool, summary_obj)
        except Exception as embed_err:
            logger.warning("chem HTML 摘要生成失败（忽略）: %s", embed_err)

    msg = (
        f"计算完成。输出目录: {out_dir.resolve()}；"
        f"JSON {len(artifacts['json_paths'])}，图 {len(artifacts['image_paths'])}。"
    )
    core_md = _markdown_core_metrics_table(cli_tool, summary_obj) if summary_obj else ""
    return {
        "status": "success",
        "message": msg,
        "markdown": _markdown_chem_report(
            md_title,
            image_urls,
            json_url,
            labs,
            html_summary or None,
            core_metrics_md=core_md or None,
        ),
        "image_urls": image_urls,
        "json_url": json_url,
        "data": sanitize_for_json(
            {
                "script_path": str(script),
                "output_dir": str(out_dir.resolve()),
                "stdout_tail": stdout_tail,
                "stderr_tail": stderr_tail,
                **artifacts,
            }
        ),
    }


_LABEL_2D = ("2D 结构",)


@registry.register(
    name="chem_bbb_assessment",
    description=(
        "BBB 血脑屏障渗透 **启发式** 评估（RDKit 理化描述符：MW/logP/TPSA/HBD 等综合）。"
        "输入：**file_path**（含 SMILES 文本）或 **smiles_text**。子进程隔离执行。"
    ),
    category="MedicinalChemistry",
    output_type="mixed",
)
@safe_tool_execution
def chem_bbb_assessment(
    file_path: str = "",
    smiles_text: str = "",
    timeout_seconds: int = 240,
) -> Dict[str, Any]:
    return _run_chem_tool(
        "bbb",
        "bbb",
        "BBB 血脑屏障（启发式）筛选",
        _LABEL_2D,
        file_path=file_path,
        smiles_text=smiles_text,
        timeout_seconds=timeout_seconds,
    )


@registry.register(
    name="chem_pains_filter",
    description=(
        "PAINS（泛测定干扰结构）过滤：RDKit FilterCatalog PAINS A/B/C。"
        "输入：**file_path** 或 **smiles_text**。"
    ),
    category="MedicinalChemistry",
    output_type="mixed",
)
@safe_tool_execution
def chem_pains_filter(
    file_path: str = "",
    smiles_text: str = "",
    timeout_seconds: int = 240,
) -> Dict[str, Any]:
    return _run_chem_tool(
        "pains",
        "pains",
        "PAINS 泛测定干扰过滤",
        _LABEL_2D,
        file_path=file_path,
        smiles_text=smiles_text,
        timeout_seconds=timeout_seconds,
    )


@registry.register(
    name="chem_brenk_filter",
    description=(
        "Brenk 毒性/不良片段筛查：优先使用 RDKit FilterCatalog.BRENK（若可用），否则 SMARTS 兜底。"
        "输入：**file_path** 或 **smiles_text**。"
    ),
    category="MedicinalChemistry",
    output_type="mixed",
)
@safe_tool_execution
def chem_brenk_filter(
    file_path: str = "",
    smiles_text: str = "",
    timeout_seconds: int = 240,
) -> Dict[str, Any]:
    return _run_chem_tool(
        "brenk",
        "brenk",
        "Brenk 毒性/不良片段筛查",
        _LABEL_2D,
        file_path=file_path,
        smiles_text=smiles_text,
        timeout_seconds=timeout_seconds,
    )


@registry.register(
    name="chem_morgan_fingerprint",
    description=(
        "生成分子 Morgan（ECFP 类）指纹比特向量摘要（SHA256 预览 + 开位数）。"
        "输入：**file_path** 或 **smiles_text**。"
    ),
    category="MedicinalChemistry",
    output_type="mixed",
)
@safe_tool_execution
def chem_morgan_fingerprint(
    file_path: str = "",
    smiles_text: str = "",
    timeout_seconds: int = 240,
) -> Dict[str, Any]:
    return _run_chem_tool(
        "morgan_fp",
        "morgan",
        "Morgan 指纹（比特向量摘要）",
        _LABEL_2D,
        file_path=file_path,
        smiles_text=smiles_text,
        timeout_seconds=timeout_seconds,
    )


@registry.register(
    name="chem_molecular_weight",
    description=(
        "分子量（平均分子量 / 精确分子量）。输入：**file_path** 或 **smiles_text**。"
    ),
    category="MedicinalChemistry",
    output_type="mixed",
)
@safe_tool_execution
def chem_molecular_weight(
    file_path: str = "",
    smiles_text: str = "",
    timeout_seconds: int = 240,
) -> Dict[str, Any]:
    return _run_chem_tool(
        "mol_weight",
        "molwt",
        "分子量计算",
        _LABEL_2D,
        file_path=file_path,
        smiles_text=smiles_text,
        timeout_seconds=timeout_seconds,
    )


@registry.register(
    name="chem_tanimoto_similarity",
    description=(
        "基于 Morgan 指纹（2048-bit）的双分子 **Tanimoto 相似度**。"
        "输入：**file_path**（两行 SMILES）、或 **smiles_text + smiles2_text**、或 smiles_text 内含 `|` 分隔两条。"
    ),
    category="MedicinalChemistry",
    output_type="mixed",
)
@safe_tool_execution
def chem_tanimoto_similarity(
    file_path: str = "",
    smiles_text: str = "",
    smiles2_text: str = "",
    timeout_seconds: int = 240,
) -> Dict[str, Any]:
    return _run_chem_tool(
        "tanimoto",
        "tanimoto",
        "Tanimoto 相似度（双分子）",
        ("分子 A（2D）", "分子 B（2D）"),
        file_path=file_path,
        smiles_text=smiles_text,
        smiles2_text=smiles2_text,
        timeout_seconds=timeout_seconds,
    )
