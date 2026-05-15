# -*- coding: utf-8 -*-
"""
CRISPR-Cas9 编辑流程仿真（第三批技能包）：子进程调用 `assets/bioml_batch3/crispr_cas9.py`，仅标准库 + 本地随机模拟。
"""
from __future__ import annotations

import asyncio
import logging
import os
import sys
import uuid
from pathlib import Path
from typing import Any, Dict, Optional

from ..core.tool_registry import registry
from ..core.utils import safe_tool_execution

logger = logging.getLogger(__name__)


def _package_dir() -> Path:
    return Path(__file__).resolve().parent.parent


def _crispr_script() -> Path:
    return _package_dir() / "assets" / "bioml_batch3" / "crispr_cas9.py"


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


@registry.register(
    name="crispr_cas9_simulation",
    description=(
        "CRISPR-Cas9 基因组编辑仿真：验证 gRNA、扫描 NGG PAM、估计递送效率并模拟 "
        "DSB 修复（NHEJ/HDR），输出 Markdown 摘要及原始/编辑序列 FASTA 下载链接。"
    ),
    category="Biomedicine",
)
@safe_tool_execution
async def crispr_cas9_simulation(
    guides_text: str,
    target_sequence: str,
    cell_line: str = "HEK293",
    result_format: str = "markdown",
    random_seed: Optional[int] = None,
) -> Dict[str, Any]:
    """
    guides_text: 一条或多条 20nt gRNA（DNA 字母），多条英文逗号分隔。
    target_sequence: 靶标基因组 DNA 片段。
    cell_line: 细胞系关键字（与脚本内置表一致，默认 HEK293）。
    result_format: markdown 或 json（默认 markdown，便于右侧工作台渲染）。
    random_seed: 可选，固定随机修复细节。
    """
    script = _crispr_script()
    if not script.is_file():
        return {
            "status": "error",
            "message": f"未找到仿真脚本: {script}（请确认镜像已包含 gibh_agent/assets/bioml_batch3）。",
        }

    guides = (guides_text or "").strip()
    target = (target_sequence or "").strip()
    if not guides or not target:
        return {"status": "error", "message": "guides_text 与 target_sequence 不能为空。"}

    fmt = (result_format or "markdown").strip().lower()
    if fmt not in ("markdown", "json"):
        fmt = "markdown"

    run_id = uuid.uuid4().hex[:16]
    out_dir = _results_dir() / "crispr_cas9" / run_id
    try:
        out_dir.mkdir(parents=True, exist_ok=False)
    except OSError as exc:
        return {"status": "error", "message": f"无法创建输出目录: {exc}"}

    py = (sys.executable or "python3").strip()
    cmd: list[str] = [
        py,
        str(script),
        "--guides",
        guides,
        "--target",
        target,
        "--cell",
        (cell_line or "HEK293").strip() or "HEK293",
        "--output",
        str(out_dir),
        "--format",
        fmt,
    ]
    if random_seed is not None:
        cmd.extend(["--seed", str(int(random_seed))])

    timeout = float(os.getenv("CRISPR_CAS9_SUBPROCESS_TIMEOUT", "120"))
    env = {**os.environ, "PYTHONUNBUFFERED": "1"}
    logger.info("CRISPR-Cas9 子进程: %s", " ".join(cmd[:6]) + " ...")

    proc = await asyncio.create_subprocess_exec(
        *cmd,
        cwd=str(script.parent),
        env=env,
        stdout=asyncio.subprocess.PIPE,
        stderr=asyncio.subprocess.PIPE,
    )
    try:
        stdout_b, stderr_b = await asyncio.wait_for(proc.communicate(), timeout=timeout)
    except asyncio.TimeoutError:
        proc.kill()
        await proc.communicate()
        return {"status": "error", "message": f"CRISPR-Cas9 子进程超时（>{timeout}s）"}

    stderr_text = (stderr_b or b"").decode("utf-8", errors="replace").strip()
    if stderr_text:
        logger.warning("CRISPR-Cas9 stderr: %s", stderr_text[:4000])

    stdout_text = (stdout_b or b"").decode("utf-8", errors="replace").strip()
    if proc.returncode != 0:
        return {
            "status": "error",
            "message": (
                f"仿真进程退出码 {proc.returncode}。"
                + (f" stderr: {stderr_text[:2000]}" if stderr_text else "")
                + (f" stdout: {stdout_text[:2000]}" if stdout_text else "")
            ),
        }

    rel_base = f"crispr_cas9/{run_id}"
    extras: Dict[str, str] = {}
    for fn in ("original.txt", "modified.txt"):
        if (out_dir / fn).is_file():
            extras[f"{fn.replace('.', '_')}_url"] = _public_path_url(f"/results/{rel_base}/{fn}")

    md_body = stdout_text
    if fmt == "markdown" and extras:
        lines = [md_body, "", "#### 序列文件下载", ""]
        if extras.get("original_txt_url"):
            lines.append(f"- [原始靶序列 FASTA]({extras['original_txt_url']})")
        if extras.get("modified_txt_url"):
            lines.append(f"- [编辑后序列 FASTA]({extras['modified_txt_url']})")
        md_body = "\n".join(lines)

    out: Dict[str, Any] = {
        "status": "success",
        "message": "CRISPR-Cas9 仿真已完成。",
        "markdown": md_body,
        "output_dir": str(out_dir),
    }
    out.update(extras)
    return out
