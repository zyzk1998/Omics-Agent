# -*- coding: utf-8 -*-
"""
第二批 RDKit 技能（芳香性感知 / Lipinski 五规则 / 官能团 / Kekulization / Pattern 指纹 / Tanimoto 矩阵）：
子进程调用 ``gibh_agent/assets/chem_rdkit_batch2/*.py``，返回与 GI / chem_misc 一致的右栏契约。
"""
from __future__ import annotations

import json
import logging
import os
import subprocess
import uuid
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

from pydantic import BaseModel, Field, model_validator

from gibh_agent.core.tool_registry import registry
from gibh_agent.core.utils import safe_tool_execution, sanitize_for_json
from gibh_agent.tools.chem_misc_tools import _write_chem_summary
from gibh_agent.tools.chem_rdkit_tools import (
    _collect_artifacts,
    _markdown_chem_report,
    _pick_chem_python,
    _public_urls_from_artifacts,
    _safe_load_chem_summary_json,
)

logger = logging.getLogger(__name__)

_TOOL_SLUG = {
    "chem_aromaticity_perception": "aromaticity_perception",
    "chem_lipinski_five": "lipinski_five",
    "chem_functional_groups": "functional_groups",
    "chem_kekulization": "kekulization",
    "chem_pattern_fingerprint": "pattern_fingerprint",
    "chem_tanimoto_matrix": "tanimoto_matrix",
}


def _repo_root() -> Path:
    return Path(__file__).resolve().parents[2]


def _prepare_batch2_out_dir(tag: str) -> Path:
    results = Path(os.getenv("RESULTS_DIR", "/app/results")).expanduser().resolve()
    rid = uuid.uuid4().hex[:12]
    out = results / "chem_rdkit_batch2" / f"{tag}_{rid}"
    out.mkdir(parents=True, exist_ok=True)
    return out


def _resolve_batch2_script(filename: str) -> Optional[Path]:
    env_key = f"BATCH2_{filename.upper().replace('.', '_')}_SCRIPT"
    env = (os.getenv(env_key) or "").strip()
    if env:
        p = Path(env).expanduser().resolve()
        if p.is_file():
            return p
        logger.warning("%s 已设置但不是文件: %s", env_key, p)
    bundled = _repo_root() / "gibh_agent" / "assets" / "chem_rdkit_batch2" / filename
    if bundled.is_file():
        return bundled.resolve()
    return None


def _metrics_md(tool_id: str, summary: Optional[Dict[str, Any]]) -> str:
    if not summary or not isinstance(summary.get("result"), dict):
        return ""
    res = summary["result"]
    rows: List[Tuple[str, str]] = []
    if tool_id == "chem_aromaticity_perception":
        rows = [
            ("输入 SMILES", str(res.get("InputSMILES") or "—")),
            ("芳香式 SMILES", str(res.get("AromaticSMILES") or "—")),
        ]
    elif tool_id == "chem_lipinski_five":
        rows = [
            ("MW (Da)", str(res.get("MW") if res.get("MW") is not None else "—")),
            ("logP", str(res.get("logP") if res.get("logP") is not None else "—")),
            ("类药性 (≤1 条违规)", str(res.get("is_druglike"))),
            ("违规条数", str(res.get("violation_count") if res.get("violation_count") is not None else "—")),
        ]
    elif tool_id == "chem_functional_groups":
        groups = res.get("functional_groups")
        if isinstance(groups, list):
            preview = ", ".join(str(x) for x in groups[:12])
            if len(groups) > 12:
                preview += " …"
        else:
            preview = "—"
        rows = [
            ("检出官能团数", str(len(groups)) if isinstance(groups, list) else "—"),
            ("预览", preview or "—"),
        ]
    elif tool_id == "chem_kekulization":
        rows = [
            ("输入 SMILES", str(res.get("InputSMILES") or "—")),
            ("Kekulé 式 SMILES", str(res.get("KekulizedSMILES") or "—")),
        ]
    elif tool_id == "chem_pattern_fingerprint":
        rows = [
            ("SMILES", str(res.get("SMILES") or res.get("smiles") or "—")),
            ("位数", str(res.get("NumBits") or "—")),
            ("置位位数", str(res.get("NumOnBits") or "—")),
        ]
    elif tool_id == "chem_tanimoto_matrix":
        sl = res.get("smiles_list")
        dm = res.get("distance_matrix")
        n = len(sl) if isinstance(sl, list) else 0
        rows = [
            ("分子数", str(n)),
            ("矩阵阶数", f"{len(dm)}×{len(dm[0])}" if isinstance(dm, list) and dm and isinstance(dm[0], list) else "—"),
        ]
    else:
        return ""
    lines = ["| 指标 | 数值 |", "| --- | --- |"]
    for k, v in rows:
        vv = str(v).replace("|", "\\|")
        lines.append(f"| {k} | {vv} |")
    return "\n".join(["#### 核心摘要", "", "\n".join(lines), ""])


def _finalize_batch2_run(
    *,
    out_dir: Path,
    tool_id: str,
    flat: Dict[str, Any],
    md_title: str,
    caption: Tuple[str, ...],
    summary_label_key: str,
) -> Dict[str, Any]:
    inner = {k: v for k, v in flat.items() if k not in ("Status", "Message")}
    summary_payload = {**flat}
    _write_chem_summary(
        out_dir,
        _TOOL_SLUG.get(tool_id, tool_id.replace("chem_", "", 1)),
        summary_payload,
        label_key=summary_label_key,
    )
    artifacts = _collect_artifacts(out_dir)
    image_urls, json_url = _public_urls_from_artifacts(artifacts)
    summary_obj = _safe_load_chem_summary_json(out_dir)
    core_md = _metrics_md(tool_id, summary_obj)
    md = _markdown_chem_report(
        md_title,
        image_urls,
        json_url,
        caption,
        None,
        core_metrics_md=core_md or None,
    )
    return {
        "status": "success",
        "message": f"{md_title} 已完成。",
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


class ChemSmilesOnlyInput(BaseModel):
    smiles_text: str = Field(..., min_length=1, description="单条 SMILES")
    timeout_seconds: int = Field(240, ge=15, le=3600)


class ChemPatternFpInput(BaseModel):
    smiles_text: str = Field(..., min_length=1, description="SMILES")
    n_bits: int = Field(2048, ge=64, le=8192)
    timeout_seconds: int = Field(240, ge=15, le=3600)


class ChemTanimotoMatrixInput(BaseModel):
    """多分子：上传 ``file_path``（每行一条 SMILES）或内联 ``smiles_text``（多行 / 竖线分隔）。"""

    file_path: str = Field(default="", description="每行一条 SMILES 的文本文件")
    smiles_text: str = Field(
        default="",
        description="多条 SMILES：多行，或一行内用 | 分隔",
    )
    radius: int = Field(2, ge=0, le=5)
    n_bits: int = Field(2048, ge=64, le=8192)
    timeout_seconds: int = Field(360, ge=30, le=7200)

    @model_validator(mode="after")
    def _src(self) -> ChemTanimotoMatrixInput:
        fp = (self.file_path or "").strip()
        st = (self.smiles_text or "").strip()
        if bool(fp) == bool(st):
            raise ValueError("必须且仅能指定 **file_path** 或 **smiles_text** 之一")
        return self


def _split_multi_smiles(raw: str) -> List[str]:
    t = (raw or "").strip()
    if "|" in t and "\n" not in t:
        parts = [x.strip() for x in t.split("|")]
    else:
        parts = [ln.strip() for ln in t.splitlines() if ln.strip()]
    out = [p for p in parts if p and not p.startswith("#")]
    return out


def _err_missing_smiles_line() -> Dict[str, Any]:
    return {
        "status": "error",
        "message": (
            "请提供一条有效的 **分子结构线性表示**（对话正文或表格中的结构串均可）。"
            "若已上传单行结构的文本附件，请在消息中说明使用附件，以便助手填入路径。"
        ),
    }


def _err_matrix_input() -> Dict[str, Any]:
    return {
        "status": "error",
        "message": (
            "请任选其一：**①** 在对话中写出 **至少两条** 结构（多行书写，或同一行用 **|** 分隔）；"
            "**②** 上传每行一条结构的 `.txt`/`.smi` 并在消息中确认。"
            "勿同时留空「文件路径」与「多行结构正文」。"
        ),
    }


def _err_matrix_both_sources() -> Dict[str, Any]:
    return {
        "status": "error",
        "message": "「上传文件」与「正文多行结构」请 **二选一**，请勿同时填写两类输入。",
    }


# ---------------------------------------------------------------------------
# 工具实现
# ---------------------------------------------------------------------------


@registry.register(
    name="chem_aromaticity_perception",
    description="将凯库勒式 SMILES 规范为 **芳香表示**（RDKit 感知），便于检索与标准化。",
    category="MedicinalChemistry",
    output_type="mixed",
)
@safe_tool_execution
def chem_aromaticity_perception(
    smiles_text: str = "",
    timeout_seconds: int = 240,
) -> Dict[str, Any]:
    if not (smiles_text or "").strip():
        return _err_missing_smiles_line()
    try:
        _ = ChemSmilesOnlyInput(smiles_text=smiles_text, timeout_seconds=timeout_seconds)
    except Exception as e:
        return {"status": "error", "message": str(e)}
    script = _resolve_batch2_script("aromaticity_perception.py")
    if script is None:
        return {"status": "error", "message": "未找到 aromaticity_perception.py"}
    py = _pick_chem_python()
    out_dir = _prepare_batch2_out_dir("aromaticity")
    cmd = [py, str(script), "--smiles", smiles_text.strip(), "--format", "json"]
    proc = subprocess.run(
        cmd,
        capture_output=True,
        text=True,
        timeout=max(30, int(timeout_seconds)),
        env={**os.environ, "PYTHONUNBUFFERED": "1"},
        cwd=str(out_dir),
    )
    tail = (proc.stdout or "")[-8000:]
    try:
        flat = json.loads(proc.stdout.strip()) if proc.stdout.strip() else None
    except json.JSONDecodeError:
        flat = None
    if proc.returncode != 0 or not flat or flat.get("Status") == "error":
        msg = (flat or {}).get("Message") if isinstance(flat, dict) else ""
        hint = (msg or proc.stderr or tail or "").strip()[-800:]
        return {"status": "error", "message": f"芳香性感知失败：{hint}"}
    return _finalize_batch2_run(
        out_dir=out_dir,
        tool_id="chem_aromaticity_perception",
        flat=flat,
        md_title="芳香性规范化（Kekulé → Aromatic）",
        caption=("结构表示",),
        summary_label_key="InputSMILES",
    )


@registry.register(
    name="chem_lipinski_five",
    description="**Lipinski 五规则** 快速评估（MW、logP、氢键供体/受体及违规条数），用于早期口服类药性粗筛。",
    category="MedicinalChemistry",
    output_type="mixed",
)
@safe_tool_execution
def chem_lipinski_five(
    smiles_text: str = "",
    timeout_seconds: int = 240,
) -> Dict[str, Any]:
    if not (smiles_text or "").strip():
        return _err_missing_smiles_line()
    try:
        _ = ChemSmilesOnlyInput(smiles_text=smiles_text, timeout_seconds=timeout_seconds)
    except Exception as e:
        return {"status": "error", "message": str(e)}
    script = _resolve_batch2_script("druglikeness.py")
    if script is None:
        return {"status": "error", "message": "未找到 druglikeness.py"}
    py = _pick_chem_python()
    out_dir = _prepare_batch2_out_dir("lipinski")
    cmd = [py, str(script), "--smiles", smiles_text.strip(), "--format", "json"]
    proc = subprocess.run(
        cmd,
        capture_output=True,
        text=True,
        timeout=max(30, int(timeout_seconds)),
        env={**os.environ, "PYTHONUNBUFFERED": "1"},
        cwd=str(out_dir),
    )
    try:
        flat = json.loads(proc.stdout.strip()) if proc.stdout.strip() else None
    except json.JSONDecodeError:
        flat = None
    if proc.returncode != 0 or not flat or flat.get("Status") == "error":
        msg = (flat or {}).get("Message") if isinstance(flat, dict) else ""
        hint = (msg or proc.stderr or proc.stdout or "").strip()[-800:]
        return {"status": "error", "message": f"Lipinski 评估失败：{hint}"}
    return _finalize_batch2_run(
        out_dir=out_dir,
        tool_id="chem_lipinski_five",
        flat=flat,
        md_title="Lipinski 五规则类药性",
        caption=("理化阈值",),
        summary_label_key="smiles",
    )


@registry.register(
    name="chem_functional_groups",
    description="基于 SMARTS 识别分子中 **常见官能团**（醇/酚、羰基、胺、芳香杂环等），用于结构解析与合成路线讨论。",
    category="MedicinalChemistry",
    output_type="mixed",
)
@safe_tool_execution
def chem_functional_groups(
    smiles_text: str = "",
    timeout_seconds: int = 240,
) -> Dict[str, Any]:
    if not (smiles_text or "").strip():
        return _err_missing_smiles_line()
    try:
        _ = ChemSmilesOnlyInput(smiles_text=smiles_text, timeout_seconds=timeout_seconds)
    except Exception as e:
        return {"status": "error", "message": str(e)}
    script = _resolve_batch2_script("functional_groups.py")
    if script is None:
        return {"status": "error", "message": "未找到 functional_groups.py"}
    py = _pick_chem_python()
    out_dir = _prepare_batch2_out_dir("fgroups")
    cmd = [py, str(script), "--smiles", smiles_text.strip(), "--format", "json"]
    proc = subprocess.run(
        cmd,
        capture_output=True,
        text=True,
        timeout=max(30, int(timeout_seconds)),
        env={**os.environ, "PYTHONUNBUFFERED": "1"},
        cwd=str(out_dir),
    )
    try:
        flat = json.loads(proc.stdout.strip()) if proc.stdout.strip() else None
    except json.JSONDecodeError:
        flat = None
    if proc.returncode != 0 or not flat or flat.get("Status") == "error":
        msg = (flat or {}).get("Message") if isinstance(flat, dict) else ""
        hint = (msg or proc.stderr or proc.stdout or "").strip()[-800:]
        return {"status": "error", "message": f"官能团识别失败：{hint}"}
    return _finalize_batch2_run(
        out_dir=out_dir,
        tool_id="chem_functional_groups",
        flat=flat,
        md_title="官能团识别",
        caption=("子结构标签",),
        summary_label_key="smiles",
    )


@registry.register(
    name="chem_kekulization",
    description="**Kekulization**：将芳香 SMILES 展开为 **显式单双键交替** 表示，便于部分绘图或导出流程。",
    category="MedicinalChemistry",
    output_type="mixed",
)
@safe_tool_execution
def chem_kekulization(
    smiles_text: str = "",
    timeout_seconds: int = 240,
) -> Dict[str, Any]:
    if not (smiles_text or "").strip():
        return _err_missing_smiles_line()
    try:
        _ = ChemSmilesOnlyInput(smiles_text=smiles_text, timeout_seconds=timeout_seconds)
    except Exception as e:
        return {"status": "error", "message": str(e)}
    script = _resolve_batch2_script("kekulization.py")
    if script is None:
        return {"status": "error", "message": "未找到 kekulization.py"}
    py = _pick_chem_python()
    out_dir = _prepare_batch2_out_dir("kekule")
    cmd = [py, str(script), "--smiles", smiles_text.strip(), "--format", "json"]
    proc = subprocess.run(
        cmd,
        capture_output=True,
        text=True,
        timeout=max(30, int(timeout_seconds)),
        env={**os.environ, "PYTHONUNBUFFERED": "1"},
        cwd=str(out_dir),
    )
    try:
        flat = json.loads(proc.stdout.strip()) if proc.stdout.strip() else None
    except json.JSONDecodeError:
        flat = None
    if proc.returncode != 0 or not flat or flat.get("Status") == "error":
        msg = (flat or {}).get("Message") if isinstance(flat, dict) else ""
        hint = (msg or proc.stderr or proc.stdout or "").strip()[-800:]
        return {"status": "error", "message": f"Kekulization 失败：{hint}"}
    return _finalize_batch2_run(
        out_dir=out_dir,
        tool_id="chem_kekulization",
        flat=flat,
        md_title="Kekulé 式展开",
        caption=("键级显式化",),
        summary_label_key="InputSMILES",
    )


@registry.register(
    name="chem_pattern_fingerprint",
    description="生成 **Pattern 指纹**（基于预定义子结构的固定长度比特向量），用于子结构特征与机器学习输入摘要。",
    category="MedicinalChemistry",
    output_type="mixed",
)
@safe_tool_execution
def chem_pattern_fingerprint(
    smiles_text: str = "",
    n_bits: int = 2048,
    timeout_seconds: int = 240,
) -> Dict[str, Any]:
    if not (smiles_text or "").strip():
        return _err_missing_smiles_line()
    try:
        _ = ChemPatternFpInput(
            smiles_text=smiles_text,
            n_bits=n_bits,
            timeout_seconds=timeout_seconds,
        )
    except Exception as e:
        return {"status": "error", "message": str(e)}
    script = _resolve_batch2_script("generate_fingerprint.py")
    if script is None:
        return {"status": "error", "message": "未找到 generate_fingerprint.py"}
    py = _pick_chem_python()
    out_dir = _prepare_batch2_out_dir("patternfp")
    cmd = [
        py,
        str(script),
        "--smiles",
        smiles_text.strip(),
        "--format",
        "json",
        "--bits",
        str(int(n_bits)),
    ]
    proc = subprocess.run(
        cmd,
        capture_output=True,
        text=True,
        timeout=max(30, int(timeout_seconds)),
        env={**os.environ, "PYTHONUNBUFFERED": "1"},
        cwd=str(out_dir),
    )
    try:
        flat = json.loads(proc.stdout.strip()) if proc.stdout.strip() else None
    except json.JSONDecodeError:
        flat = None
    if proc.returncode != 0 or not flat or flat.get("Status") == "error":
        msg = (flat or {}).get("Message") if isinstance(flat, dict) else ""
        hint = (msg or proc.stderr or proc.stdout or "").strip()[-800:]
        return {"status": "error", "message": f"Pattern 指纹失败：{hint}"}
    return _finalize_batch2_run(
        out_dir=out_dir,
        tool_id="chem_pattern_fingerprint",
        flat=flat,
        md_title="Pattern 分子指纹",
        caption=("比特摘要",),
        summary_label_key="SMILES",
    )


@registry.register(
    name="chem_tanimoto_matrix",
    description="基于 Morgan 指纹计算 **多分子两两 Tanimoto 距离矩阵**（距离 = 1 − 相似度），用于小型化合物集的多样性粗测。",
    category="MedicinalChemistry",
    output_type="mixed",
)
@safe_tool_execution
def chem_tanimoto_matrix(
    file_path: str = "",
    smiles_text: str = "",
    radius: int = 2,
    n_bits: int = 2048,
    timeout_seconds: int = 360,
) -> Dict[str, Any]:
    fp = (file_path or "").strip()
    st = (smiles_text or "").strip()
    if not fp and not st:
        return _err_matrix_input()
    if fp and st:
        return _err_matrix_both_sources()
    try:
        _ = ChemTanimotoMatrixInput(
            file_path=file_path or "",
            smiles_text=smiles_text or "",
            radius=radius,
            n_bits=n_bits,
            timeout_seconds=timeout_seconds,
        )
    except Exception as e:
        return {"status": "error", "message": str(e)}
    script = _resolve_batch2_script("tanimoto_matrix.py")
    if script is None:
        return {"status": "error", "message": "未找到 tanimoto_matrix.py"}
    py = _pick_chem_python()
    out_dir = _prepare_batch2_out_dir("tanimat")
    cmd: List[str] = [
        py,
        str(script),
        "--format",
        "json",
        "--radius",
        str(int(radius)),
        "--n-bits",
        str(int(n_bits)),
    ]
    if fp:
        pth = Path(fp).expanduser().resolve()
        if not pth.is_file():
            return {"status": "error", "message": f"输入文件不存在：{pth}"}
        cmd += ["--file", str(pth)]
    else:
        parts = _split_multi_smiles(st)
        if len(parts) < 2:
            return {
                "status": "error",
                "message": "距离矩阵至少需要 **2 条**有效 SMILES（多行书写或用 | 分隔）。",
            }
        cmd += ["--smiles", *parts]
    proc = subprocess.run(
        cmd,
        capture_output=True,
        text=True,
        timeout=max(45, int(timeout_seconds)),
        env={**os.environ, "PYTHONUNBUFFERED": "1"},
        cwd=str(out_dir),
    )
    raw_out = (proc.stdout or "").strip()
    try:
        parsed = json.loads(raw_out) if raw_out else None
    except json.JSONDecodeError:
        parsed = None
    if proc.returncode != 0 or not parsed:
        hint = (proc.stderr or raw_out or "").strip()[-800:]
        return {"status": "error", "message": f"Tanimoto 矩阵失败：{hint}"}
    flat = {
        "Status": "success",
        "Message": "Tanimoto distance matrix completed",
        "smiles_list": parsed.get("smiles_list"),
        "distance_matrix": parsed.get("distance_matrix"),
        "parameters": parsed.get("parameters"),
    }
    return _finalize_batch2_run(
        out_dir=out_dir,
        tool_id="chem_tanimoto_matrix",
        flat=flat,
        md_title="Tanimoto 距离矩阵",
        caption=("结构多样性",),
        summary_label_key="Message",
    )
