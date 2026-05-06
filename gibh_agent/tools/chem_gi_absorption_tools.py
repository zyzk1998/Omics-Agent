# -*- coding: utf-8 -*-
"""
胃肠道被动吸收启发式评估（gi-absorption 技能包）：子进程调用 ``gibh_agent/assets/gi_absorption/gi_absorption.py``。
成功返回 ``markdown`` + ``image_urls`` + ``json_url``（与化学 RDKit 原子工具右栏一致）。
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
from typing import Any, Dict, List, Optional

from pydantic import BaseModel, Field

from gibh_agent.core.tool_registry import registry
from gibh_agent.core.utils import safe_tool_execution, sanitize_for_json
from gibh_agent.tools.chem_rdkit_tools import (
    _collect_artifacts,
    _markdown_chem_report,
    _public_urls_from_artifacts,
    _safe_load_chem_summary_json,
)

logger = logging.getLogger(__name__)


class GIAbsorptionInput(BaseModel):
    """
    GI 吸收评估参数契约：与 ``chem_gi_absorption`` 注册签名一致；
    路径类字段命名遵循架构宪法词汇表。
    """

    file_path: str = Field(
        default="",
        description="含首行 SMILES 的文本/.csv 路径（与会话附件解析路径一致）",
    )
    smiles_text: str = Field(default="", description="单行 SMILES 字符串")
    timeout_seconds: int = Field(default=240, ge=15, le=3600, description="子进程超时（秒）")


def _repo_root() -> Path:
    return Path(__file__).resolve().parents[2]


def _resolve_gi_script() -> Optional[Path]:
    env = (os.getenv("GI_ABSORPTION_SCRIPT") or "").strip()
    if env:
        p = Path(env).expanduser().resolve()
        if p.is_file():
            return p
        logger.warning("GI_ABSORPTION_SCRIPT 已设置但不是文件: %s", p)
    bundled = _repo_root() / "gibh_agent" / "assets" / "gi_absorption" / "gi_absorption.py"
    if bundled.is_file():
        return bundled.resolve()
    return None


def _pick_chem_python() -> str:
    py = (os.getenv("CHEM_RDKIT_PYTHON") or "").strip()
    return py if py else sys.executable


def _prepare_gi_out_dir() -> Path:
    results = Path(os.getenv("RESULTS_DIR", "/app/results")).expanduser().resolve()
    rid = uuid.uuid4().hex[:12]
    out = results / "gi_absorption" / rid
    out.mkdir(parents=True, exist_ok=True)
    return out


def _load_format_markdown_from_script(script: Path):
    spec = importlib.util.spec_from_file_location("_gi_absorption_fmt", script)
    if spec is None or spec.loader is None:
        return None
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return getattr(mod, "format_markdown", None)


def _validate_input(file_path: str, smiles_text: str) -> Optional[str]:
    fp = (file_path or "").strip()
    if fp and Path(fp).expanduser().is_file():
        return None
    if (smiles_text or "").strip():
        return None
    return "请提供 **smiles_text**（单行 SMILES）或已上传文件的 **file_path**。"


def _markdown_gi_metrics(summary: Optional[Dict[str, Any]]) -> str:
    if not summary or not isinstance(summary.get("result"), dict):
        return ""
    res = summary["result"]
    rows = [
        (
            "Canonical SMILES（截断）",
            _md_cell((str(summary.get("canonical_smiles") or res.get("smiles") or ""))[:400]),
        ),
        ("TPSA (Å²)", _md_cell(res.get("TPSA"))),
        ("logP", _md_cell(res.get("logP"))),
        ("高 GI 被动吸收（TPSA≤140 且 logP≤5）", _md_cell(res.get("is_high_gi_absorption"))),
    ]
    lines = ["| 指标 | 数值 |", "| --- | --- |"]
    for k, v in rows:
        lines.append(f"| {k} | {v} |")
    return "\n".join(["#### 核心数值摘要", "", "\n".join(lines), ""])


def _md_cell(v: Any) -> str:
    if v is None:
        return "—"
    if isinstance(v, bool):
        return "是" if v else "否"
    s = str(v).strip()
    if not s:
        return "—"
    return s.replace("|", "\\|")


def _html_gi_banner(summary: Dict[str, Any]) -> str:
    from html import escape as html_esc

    res = summary.get("result")
    if not isinstance(res, dict):
        return ""
    tpsa = res.get("TPSA")
    logp = res.get("logP")
    high = res.get("is_high_gi_absorption")
    ok = bool(high) if high is not None else False
    bg = "#def7ec" if ok else "#fef3c7"
    fg = "#03543f" if ok else "#92400e"
    lbl = "高胃肠道被动吸收可能性（启发式）" if ok else "口服被动吸收可能偏低（启发式）"
    tp_s = html_esc(str(tpsa)) if tpsa is not None else "—"
    lp_s = html_esc(str(logp)) if logp is not None else "—"
    return (
        f'<p style="margin:12px 0 8px 0;">'
        f'<span style="background:{bg}; color:{fg}; padding:2px 8px; border-radius:12px; '
        f"font-weight:bold; font-size:12px;"
        f">{html_esc(lbl)}</span></p>"
        f'<p style="margin:4px 0 12px 0; font-size:13px; color:#374151;">'
        f"<strong>TPSA</strong> {tp_s} Å² · <strong>logP</strong> {lp_s}"
        f"</p>"
    )


def _run_gi_absorption_subprocess(
    *,
    file_path: str = "",
    smiles_text: str = "",
    timeout_seconds: int = 240,
) -> Dict[str, Any]:
    script = _resolve_gi_script()
    if script is None:
        root = _repo_root()
        return {
            "status": "error",
            "message": (
                "未找到 gi_absorption.py。请确认已包含 "
                f"{root / 'gibh_agent/assets/gi_absorption/gi_absorption.py'} "
                "或设置 GI_ABSORPTION_SCRIPT。"
            ),
        }

    ve = _validate_input(file_path, smiles_text)
    if ve:
        return {"status": "error", "message": ve}

    out_dir = _prepare_gi_out_dir()
    py = _pick_chem_python()
    cmd: List[str] = [py, str(script), "--output-dir", str(out_dir), "--format", "json"]
    if (file_path or "").strip():
        cmd += ["--file", str(Path(file_path.strip()).expanduser().resolve())]
    else:
        cmd += ["--smiles", (smiles_text or "").strip()]

    logger.info("gi_absorption subprocess: %s", " ".join(cmd[:12]))

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
        return {"status": "error", "message": f"GI 吸收脚本超时（>{timeout_seconds}s）。"}
    except OSError as e:
        return {"status": "error", "message": f"无法启动子进程: {e}"}

    stdout_tail = (proc.stdout or "")[-6000:]
    stderr_tail = (proc.stderr or "")[-6000:]
    artifacts = _collect_artifacts(out_dir)

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
            "message": f"GI 吸收评估失败：{hint or '未知错误'}",
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

    summary_obj = _safe_load_chem_summary_json(out_dir)
    html_summary = ""
    if summary_obj:
        try:
            html_summary = _html_gi_banner(summary_obj)
        except Exception as embed_err:  # noqa: BLE001
            logger.warning("GI HTML 摘要生成失败（忽略）: %s", embed_err)

    fmt = _load_format_markdown_from_script(script)
    body_md = ""
    if callable(fmt):
        try:
            body_md = fmt(flat)
        except Exception as e:  # noqa: BLE001
            logger.warning("format_markdown 调用失败: %s", e)

    image_urls, json_url = _public_urls_from_artifacts(artifacts)
    core_md = _markdown_gi_metrics(summary_obj)
    title = "胃肠道吸收（GI）评估（被动扩散启发式）"
    md = _markdown_chem_report(
        title,
        image_urls,
        json_url,
        ("2D 结构",),
        html_summary or None,
        core_metrics_md=(core_md + "\n\n" + body_md).strip() if (core_md or body_md) else body_md or None,
    )

    msg = (
        f"GI 吸收评估完成。输出目录: {out_dir.resolve()}；"
        f"JSON {len(artifacts['json_paths'])}，图 {len(artifacts['image_paths'])}。"
    )
    return {
        "status": "success",
        "message": msg,
        "markdown": md,
        "image_urls": image_urls,
        "json_url": json_url,
        "data": sanitize_for_json(
            {
                "gi_result": flat,
                "script_path": str(script),
                "output_dir": str(out_dir.resolve()),
                "stdout_tail": stdout_tail,
                "stderr_tail": stderr_tail,
                **artifacts,
            }
        ),
    }


@registry.register(
    name="chem_gi_absorption",
    description=(
        "胃肠道（GI）被动吸收倾向启发式评估：基于 RDKit **TPSA** 与 **logP**，"
        "采用经验阈值 **TPSA≤140 Å² 且 logP≤5** 判定高吸收可能性（非 PBPK）。"
        "输入：**file_path**（首行 SMILES）或 **smiles_text**。"
    ),
    category="MedicinalChemistry",
    output_type="mixed",
)
@safe_tool_execution
def chem_gi_absorption(
    file_path: str = "",
    smiles_text: str = "",
    timeout_seconds: int = 240,
) -> Dict[str, Any]:
    # 参数模型见 ``GIAbsorptionInput``（与注册器生成的 schema 对齐）
    _ = GIAbsorptionInput(
        file_path=file_path or "",
        smiles_text=smiles_text or "",
        timeout_seconds=timeout_seconds,
    )
    return _run_gi_absorption_subprocess(
        file_path=file_path,
        smiles_text=smiles_text,
        timeout_seconds=timeout_seconds,
    )
