# -*- coding: utf-8 -*-
"""
化学扩展原子工具（molmass / 元素查询 / 分子 2D 图 / Open Babel）：子进程调用 ``gibh_agent/assets/...`` 脚本，
返回 ``markdown`` + ``image_urls`` + ``json_url``（与 GI / RDKit 右栏一致）。
"""
from __future__ import annotations

import importlib.util
import json
import logging
import os
import subprocess
import sys
import uuid
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

from pydantic import BaseModel, Field, model_validator

from gibh_agent.core.tool_registry import registry
from gibh_agent.core.utils import safe_tool_execution, sanitize_for_json
from gibh_agent.tools.chem_rdkit_tools import (
    _collect_artifacts,
    _markdown_chem_report,
    _public_urls_from_artifacts,
    _safe_load_chem_summary_json,
)

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Pydantic 输入契约（与各 @registry.register 签名字段对齐）
# ---------------------------------------------------------------------------


class ChemMolmassInput(BaseModel):
    """化学式分子量（molmass 库）：输入为 **化学式**（如 C6H12O6），非 SMILES。"""

    formula_text: str = Field(
        ...,
        min_length=1,
        description="化学式字符串（Hill / 常规写法），如 C6H12O6、H2SO4",
    )
    include_spectrum: bool = Field(False, description="是否输出同位素分布谱")
    timeout_seconds: int = Field(240, ge=15, le=3600)


class ChemElementQueryInput(BaseModel):
    """元素周期表查询：符号 / 中英文名称 / 原子序数。"""

    query_text: str = Field(..., min_length=1, description="如 Au、金、79、Gold")
    timeout_seconds: int = Field(120, ge=15, le=3600)


class ChemMoleculeImageInput(BaseModel):
    """由 SMILES 生成 2D PNG（RDKit Draw）。"""

    smiles_text: str = Field(..., min_length=2, description="SMILES")
    image_size: int = Field(400, ge=64, le=2048)
    dpi: int = Field(150, ge=72, le=600)
    timeout_seconds: int = Field(240, ge=15, le=3600)


class ChemOpenbabelInput(BaseModel):
    """Open Babel（系统 obabel）：文件路径 **或** 内联 SMILES/InChI 二选一。"""

    file_path: str = Field(default="", description="输入分子文件路径")
    smiles_text: str = Field(default="", description="内联 SMILES（或其它文本格式）")
    in_format: str = Field("smi", description="输入格式标识，如 smi mol sdf pdb inchi")
    out_format: str = Field(
        "mol",
        description="输出格式：mol pdb xyz sdf png svg smiles inchi 等",
    )
    compute_properties: bool = Field(
        False,
        description="为 True 时优先计算 molwt/logP/TPSA 等（仍可按需转换）",
    )
    add_hydrogens: bool = False
    remove_hydrogens: bool = False
    gen3d: bool = False
    optimize: bool = False
    timeout_seconds: int = Field(180, ge=15, le=3600)

    @model_validator(mode="after")
    def _one_input(self) -> ChemOpenbabelInput:
        fp = (self.file_path or "").strip()
        st = (self.smiles_text or "").strip()
        if bool(fp) == bool(st):
            raise ValueError("必须且仅能指定 **file_path** 或 **smiles_text** 之一")
        return self


def _repo_root() -> Path:
    return Path(__file__).resolve().parents[2]


def _pick_python(*, need_rdkit: bool) -> str:
    if need_rdkit:
        py = (os.getenv("CHEM_RDKIT_PYTHON") or "").strip()
        return py if py else sys.executable
    return sys.executable


def _prepare_misc_out_dir(tag: str) -> Path:
    results = Path(os.getenv("RESULTS_DIR", "/app/results")).expanduser().resolve()
    rid = uuid.uuid4().hex[:12]
    out = results / "chem_misc" / f"{tag}_{rid}"
    out.mkdir(parents=True, exist_ok=True)
    return out


def _resolve_script(env_key: str, relative: Path) -> Optional[Path]:
    env = (os.getenv(env_key) or "").strip()
    if env:
        p = Path(env).expanduser().resolve()
        if p.is_file():
            return p
        logger.warning("%s 已设置但不是文件: %s", env_key, p)
    bundled = _repo_root() / "gibh_agent" / "assets" / relative
    if bundled.is_file():
        return bundled.resolve()
    return None


def _load_format_markdown(script: Path, func_name: str):
    spec = importlib.util.spec_from_file_location("_misc_fmt", script)
    if spec is None or spec.loader is None:
        return None
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return getattr(mod, func_name, None)


def _write_chem_summary(
    out_dir: Path,
    tool: str,
    flat: Dict[str, Any],
    *,
    label_key: str = "formula",
) -> None:
    inner = {k: v for k, v in flat.items() if k not in ("Status", "Message")}
    lab = flat.get(label_key) or flat.get("smiles") or flat.get("input") or ""
    summary = {
        "tool": tool,
        "canonical_smiles": flat.get("smiles") if tool != "molmass" else None,
        "canonical_formula": flat.get("formula") if tool == "molmass" else None,
        "label": lab,
        "result": inner,
        "status_line": flat.get("Status"),
        "message_line": flat.get("Message"),
    }
    (out_dir / "chem_summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2),
        encoding="utf-8",
    )


def _md_metrics_table(tool: str, summary: Optional[Dict[str, Any]]) -> str:
    if not summary or not isinstance(summary.get("result"), dict):
        return ""
    res = summary["result"]
    rows: List[Tuple[str, str]] = []
    if tool == "molmass":
        rows = [
            ("化学式", str(res.get("formula") or summary.get("canonical_formula") or "")),
            ("平均分子量 (g/mol)", str(res.get("average_mass") or "—")),
            ("单同位素质量 (Da)", str(res.get("monoisotopic_mass") or "—")),
        ]
    elif tool == "element_query":
        rows = [
            ("原子序数", str(res.get("number") or "—")),
            ("符号", str(res.get("symbol") or "—")),
            ("名称", f"{res.get('name', '')} ({res.get('name_cn', '')})"),
            ("原子量", str(res.get("mass") or "—")),
        ]
    elif tool == "molecule_image":
        rows = [
            ("SMILES", str(res.get("smiles") or "—")),
            ("尺寸", str(res.get("image_size") or "—")),
        ]
    elif tool == "openbabel":
        props = res.get("properties") if isinstance(res.get("properties"), dict) else {}
        rows = [
            ("输入", str(res.get("input") or "—")),
            ("输出格式", str(res.get("output_format") or "—")),
        ]
        for k, v in list(props.items())[:8]:
            rows.append((k, str(v)))
    else:
        return ""
    lines = ["| 指标 | 数值 |", "| --- | --- |"]
    for k, v in rows:
        vv = str(v).replace("|", "\\|")
        lines.append(f"| {k} | {vv} |")
    return "\n".join(["#### 核心摘要", "", "\n".join(lines), ""])


def _run_stdout_json_tool(
    *,
    tag: str,
    tool_key: str,
    md_title: str,
    cmd: List[str],
    timeout_seconds: int,
    script_for_fmt: Optional[Path],
    fmt_name: str,
    summary_label_key: str = "formula",
) -> Dict[str, Any]:
    out_dir = _prepare_misc_out_dir(tag)
    proc = subprocess.run(
        cmd,
        capture_output=True,
        text=True,
        timeout=max(30, int(timeout_seconds)),
        env={**os.environ, "PYTHONUNBUFFERED": "1"},
        cwd=str(out_dir),
    )
    stdout_tail = (proc.stdout or "")[-8000:]
    stderr_tail = (proc.stderr or "")[-8000:]
    flat: Optional[Dict[str, Any]] = None
    if proc.stdout and proc.stdout.strip():
        try:
            flat = json.loads(proc.stdout.strip())
        except json.JSONDecodeError:
            flat = None

    if proc.returncode != 0 or not flat or flat.get("Status") == "error":
        msg = (flat or {}).get("Message") if isinstance(flat, dict) else ""
        hint = (msg or stderr_tail or stdout_tail or "").strip()[-800:]
        return {
            "status": "error",
            "message": f"子进程失败：{hint or '未知错误'}",
            "data": sanitize_for_json({"stdout_tail": stdout_tail, "stderr_tail": stderr_tail}),
        }

    _write_chem_summary(out_dir, tool_key, flat, label_key=summary_label_key)
    artifacts = _collect_artifacts(out_dir)
    image_urls, json_url = _public_urls_from_artifacts(artifacts)
    summary_obj = _safe_load_chem_summary_json(out_dir)
    core_md = _md_metrics_table(tool_key, summary_obj)
    body_md = ""
    if script_for_fmt and script_for_fmt.is_file():
        fmt = _load_format_markdown(script_for_fmt, fmt_name)
        if callable(fmt):
            try:
                body_md = fmt(flat)
            except Exception as e:  # noqa: BLE001
                logger.warning("format_markdown 失败: %s", e)

    md = _markdown_chem_report(
        md_title,
        image_urls,
        json_url,
        ("输出图",),
        None,
        core_metrics_md=(core_md + "\n\n" + body_md).strip() if (core_md or body_md) else body_md or None,
    )
    return {
        "status": "success",
        "message": f"{md_title} 完成。输出目录: {out_dir}",
        "markdown": md,
        "image_urls": image_urls,
        "json_url": json_url,
        "data": sanitize_for_json(
            {
                "result": flat,
                "output_dir": str(out_dir.resolve()),
                "stdout_tail": stdout_tail,
                "stderr_tail": stderr_tail,
                **artifacts,
            }
        ),
    }


@registry.register(
    name="chem_molmass",
    description=(
        "基于 **molmass** 库：由 **化学式**（如 C6H12O6）计算平均/名义/单同位素质量、元素组成，"
        "可选同位素分布谱。**非 SMILES**；SMILES 请用「分子量计算工具」RDKit 路径。"
    ),
    category="MedicinalChemistry",
    output_type="mixed",
)
@safe_tool_execution
def chem_molmass(
    formula_text: str = "",
    include_spectrum: bool = False,
    timeout_seconds: int = 240,
) -> Dict[str, Any]:
    _ = ChemMolmassInput(
        formula_text=formula_text,
        include_spectrum=include_spectrum,
        timeout_seconds=timeout_seconds,
    )
    script = _resolve_script("MOLMASS_CLI_SCRIPT", Path("molmass") / "molmass_cli.py")
    if script is None:
        return {
            "status": "error",
            "message": f"未找到 molmass_cli.py：{_repo_root() / 'gibh_agent/assets/molmass/molmass_cli.py'}",
        }
    py = _pick_python(need_rdkit=False)
    cmd: List[str] = [
        py,
        str(script),
        "--formula",
        formula_text.strip(),
        "--format",
        "json",
    ]
    if include_spectrum:
        cmd.append("--spectrum")
    return _run_stdout_json_tool(
        tag="molmass",
        tool_key="molmass",
        md_title="化学式分子量与同位素（Molmass）",
        cmd=cmd,
        timeout_seconds=timeout_seconds,
        script_for_fmt=script,
        fmt_name="format_markdown",
        summary_label_key="formula",
    )


@registry.register(
    name="chem_element_query",
    description="查询元素周期表：**符号 / 中英文名称 / 原子序数** → 原子量、电子排布、简介等。",
    category="MedicinalChemistry",
    output_type="mixed",
)
@safe_tool_execution
def chem_element_query(
    query_text: str = "",
    timeout_seconds: int = 120,
) -> Dict[str, Any]:
    _ = ChemElementQueryInput(query_text=query_text, timeout_seconds=timeout_seconds)
    script = _resolve_script(
        "ELEMENT_QUERY_SCRIPT",
        Path("chemical_element_query") / "element_query.py",
    )
    if script is None:
        return {
            "status": "error",
            "message": "未找到 element_query.py",
        }
    py = sys.executable
    cmd = [py, str(script), "--query", query_text.strip(), "--format", "json"]
    return _run_stdout_json_tool(
        tag="element",
        tool_key="element_query",
        md_title="化学元素查询（周期表）",
        cmd=cmd,
        timeout_seconds=timeout_seconds,
        script_for_fmt=script,
        fmt_name="format_markdown",
        summary_label_key="symbol",
    )


@registry.register(
    name="chem_molecule_image",
    description=(
        "由 **SMILES** 生成 **2D 分子结构 PNG**（RDKit Draw），输出挂载 `/results/`。"
    ),
    category="MedicinalChemistry",
    output_type="mixed",
)
@safe_tool_execution
def chem_molecule_image(
    smiles_text: str = "",
    image_size: int = 400,
    dpi: int = 150,
    timeout_seconds: int = 240,
) -> Dict[str, Any]:
    _ = ChemMoleculeImageInput(
        smiles_text=smiles_text,
        image_size=image_size,
        dpi=dpi,
        timeout_seconds=timeout_seconds,
    )
    script = _resolve_script(
        "MOLECULE_IMAGE_SCRIPT",
        Path("molecular_image_generation") / "generate_molecule_image.py",
    )
    if script is None:
        return {"status": "error", "message": "未找到 generate_molecule_image.py"}
    py = _pick_python(need_rdkit=True)
    out_dir = _prepare_misc_out_dir("molimg")
    png_path = out_dir / "molecule_2d.png"
    cmd = [
        py,
        str(script),
        "--smiles",
        smiles_text.strip(),
        "--output",
        str(png_path),
        "--size",
        str(int(image_size)),
        "--dpi",
        str(int(dpi)),
        "--format",
        "json",
    ]
    proc = subprocess.run(
        cmd,
        capture_output=True,
        text=True,
        timeout=max(30, int(timeout_seconds)),
        env={**os.environ, "PYTHONUNBUFFERED": "1"},
        cwd=str(out_dir),
    )
    stdout_tail = (proc.stdout or "")[-8000:]
    stderr_tail = (proc.stderr or "")[-8000:]
    flat: Optional[Dict[str, Any]] = None
    if proc.stdout and proc.stdout.strip():
        try:
            flat = json.loads(proc.stdout.strip())
        except json.JSONDecodeError:
            flat = None

    if proc.returncode != 0 or not flat or flat.get("Status") == "error":
        msg = (flat or {}).get("Message") if isinstance(flat, dict) else ""
        hint = (msg or stderr_tail or stdout_tail or "").strip()[-800:]
        return {
            "status": "error",
            "message": f"图像生成失败：{hint}",
            "data": sanitize_for_json({"stdout_tail": stdout_tail, "stderr_tail": stderr_tail}),
        }

    _write_chem_summary(out_dir, "molecule_image", flat, label_key="smiles")
    artifacts = _collect_artifacts(out_dir)
    image_urls, json_url = _public_urls_from_artifacts(artifacts)
    summary_obj = _safe_load_chem_summary_json(out_dir)
    core_md = _md_metrics_table("molecule_image", summary_obj)
    md = _markdown_chem_report(
        "分子 2D 结构图（SMILES → PNG）",
        image_urls,
        json_url,
        ("2D 结构",),
        None,
        core_metrics_md=core_md or None,
    )
    return {
        "status": "success",
        "message": "分子图像已生成。",
        "markdown": md,
        "image_urls": image_urls,
        "json_url": json_url,
        "data": sanitize_for_json(
            {
                "result": flat,
                "output_dir": str(out_dir.resolve()),
                **artifacts,
            }
        ),
    }


@registry.register(
    name="chem_openbabel",
    description=(
        "Open Babel（**系统 obabel**）：格式互转、加氢、--gen3d、性质计算等。"
        "输入二选一：**file_path**（上传文件，缺省 **in_format** 时按扩展名猜测）或 **smiles_text**（内联）；"
        "**compute_properties** 为真时走性质表；否则须 **out_format** 并写出结果文件。"
    ),
    category="MedicinalChemistry",
    output_type="mixed",
)
@safe_tool_execution
def chem_openbabel(
    file_path: str = "",
    smiles_text: str = "",
    in_format: str = "smi",
    out_format: str = "mol",
    compute_properties: bool = False,
    add_hydrogens: bool = False,
    remove_hydrogens: bool = False,
    gen3d: bool = False,
    optimize: bool = False,
    timeout_seconds: int = 180,
) -> Dict[str, Any]:
    _ = ChemOpenbabelInput(
        file_path=file_path or "",
        smiles_text=smiles_text or "",
        in_format=in_format,
        out_format=out_format,
        compute_properties=compute_properties,
        add_hydrogens=add_hydrogens,
        remove_hydrogens=remove_hydrogens,
        gen3d=gen3d,
        optimize=optimize,
        timeout_seconds=timeout_seconds,
    )
    script = _resolve_script(
        "OPENBABEL_CONVERTER_SCRIPT",
        Path("openbabel_converter") / "openbabel_converter.py",
    )
    if script is None:
        return {"status": "error", "message": "未找到 openbabel_converter.py"}
    py = sys.executable
    out_dir = _prepare_misc_out_dir("obabel")
    fp = (file_path or "").strip()
    st = (smiles_text or "").strip()

    def _guess_in_format(pth: str) -> str:
        ext = Path(pth).suffix.lower()
        m = {
            ".mol": "mol",
            ".sdf": "sdf",
            ".smi": "smi",
            ".smiles": "smi",
            ".pdb": "pdb",
            ".xyz": "xyz",
            ".inchi": "inchi",
            ".cdx": "cdx",
        }
        return m.get(ext, "mol")

    out_ext = (out_format or "mol").lower().strip()
    if out_ext in ("smiles", "smi"):
        out_file = out_dir / "converted.smi"
    elif out_ext == "png":
        out_file = out_dir / "out.png"
    elif out_ext == "svg":
        out_file = out_dir / "out.svg"
    elif out_ext == "pdb":
        out_file = out_dir / "out.pdb"
    else:
        out_file = out_dir / f"out.{out_ext}"

    inf = (in_format or "").strip()
    cmd: List[str] = [py, str(script), "--format", "json"]
    if fp:
        pth = Path(fp).expanduser().resolve()
        if not pth.is_file():
            return {"status": "error", "message": f"输入文件不存在: {pth}"}
        inf_use = inf or _guess_in_format(str(pth))
        cmd += ["--input", str(pth), "--in-format", inf_use]
    else:
        inf_use = inf or "smi"
        cmd += ["--input", st, "--in-format", inf_use]

    if compute_properties:
        cmd.append("--properties")
    else:
        cmd += ["--out-format", out_ext, "--output", str(out_file)]

    if add_hydrogens:
        cmd.append("--add-hydrogens")
    if remove_hydrogens:
        cmd.append("--remove-hydrogens")
    if gen3d:
        cmd.append("--gen3d")
    if optimize:
        cmd.append("--optimize")

    proc = subprocess.run(
        cmd,
        capture_output=True,
        text=True,
        timeout=max(30, int(timeout_seconds)),
        env={**os.environ, "PYTHONUNBUFFERED": "1"},
        cwd=str(out_dir),
    )
    stdout_tail = (proc.stdout or "")[-12000:]
    stderr_tail = (proc.stderr or "")[-8000:]
    flat: Optional[Dict[str, Any]] = None
    if proc.stdout and proc.stdout.strip():
        try:
            flat = json.loads(proc.stdout.strip())
        except json.JSONDecodeError:
            flat = None

    if proc.returncode != 0 or not flat or flat.get("Status") == "error":
        msg = (flat or {}).get("Message") if isinstance(flat, dict) else ""
        hint = (msg or stderr_tail or stdout_tail or "").strip()[-800:]
        return {
            "status": "error",
            "message": f"Open Babel 失败：{hint}",
            "data": sanitize_for_json({"stdout_tail": stdout_tail, "stderr_tail": stderr_tail}),
        }

    inner = {k: v for k, v in flat.items() if k not in ("Status", "Message")}
    summary = {
        "tool": "openbabel",
        "canonical_smiles": st or None,
        "label": flat.get("input"),
        "result": inner,
        "status_line": flat.get("Status"),
        "message_line": flat.get("Message"),
    }
    (out_dir / "chem_summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2),
        encoding="utf-8",
    )
    artifacts = _collect_artifacts(out_dir)
    image_urls, json_url = _public_urls_from_artifacts(artifacts)
    summary_obj = _safe_load_chem_summary_json(out_dir)
    core_md = _md_metrics_table("openbabel", summary_obj)
    md = _markdown_chem_report(
        "Open Babel 转换 / 性质",
        image_urls,
        json_url,
        ("输出文件",),
        None,
        core_metrics_md=core_md or None,
    )
    return {
        "status": "success",
        "message": "Open Babel 执行完成。",
        "markdown": md,
        "image_urls": image_urls,
        "json_url": json_url,
        "data": sanitize_for_json(
            {
                "result": flat,
                "output_dir": str(out_dir.resolve()),
                **artifacts,
            }
        ),
    }
