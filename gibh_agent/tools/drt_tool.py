# -*- coding: utf-8 -*-
"""
电化学 EIS → DRT（弛豫时间分布）分析 — 薄封装工具。

解析顺序：``DRT_ANALYSIS_SCRIPT`` 环境变量 → 随仓库分发的
``gibh_agent/assets/drt/drt_analysis.py`` → 历史暂存区
``uploads/skills_staging/.../drt-tools/scripts/drt_analysis.py``。

通过 **子进程** 调用底层脚本，避免在主 API 进程内 ``import`` 庞杂依赖。

环境变量：
- ``DRT_ANALYSIS_SCRIPT``：覆盖底层脚本绝对路径（可选；未设则使用仓库内置脚本）。
- ``DRT_PYTHON``：子进程 Python 解释器（默认 ``sys.executable``；若缺 numpy/scipy 可指向独立 venv）。
"""
from __future__ import annotations

import logging
import os
import subprocess
import sys
import uuid
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

from gibh_agent.core.tool_registry import registry
from gibh_agent.core.utils import safe_tool_execution, sanitize_for_json

logger = logging.getLogger(__name__)


def _repo_root() -> Path:
    # gibh_agent/tools/drt_tool.py -> parents[2] == 仓库根
    return Path(__file__).resolve().parents[2]


def resolve_drt_analysis_script() -> Optional[Path]:
    """解析底层 ``drt_analysis.py``；存在则返回绝对路径。"""
    env = (os.getenv("DRT_ANALYSIS_SCRIPT") or "").strip()
    if env:
        p = Path(env).expanduser().resolve()
        if p.is_file():
            return p
        logger.warning("DRT_ANALYSIS_SCRIPT 已设置但不是文件: %s", p)
    root = _repo_root()
    candidates = [
        # 随镜像 / 挂载仓库发布（推荐，免手工解压技能 ZIP）
        root / "gibh_agent" / "assets" / "drt" / "drt_analysis.py",
        root
        / "uploads"
        / "skills_staging"
        / "2026.04.20 skills"
        / "drt-tools"
        / "scripts"
        / "drt_analysis.py",
    ]
    for c in candidates:
        if c.is_file() and c.suffix == ".py":
            return c.resolve()
    # 常见：仅解压了子目录到 staging
    alt = (
        root
        / "uploads"
        / "skills_staging"
        / "drt-tools"
        / "scripts"
        / "drt_analysis.py"
    )
    if alt.is_file():
        return alt.resolve()
    return None


def _pick_python_executable() -> str:
    py = (os.getenv("DRT_PYTHON") or "").strip()
    if py:
        return py
    return sys.executable


def _prepare_run_output_dir(input_file: Path) -> Path:
    """在 RESULTS_DIR 下创建本次运行的输出目录。"""
    results = Path(os.getenv("RESULTS_DIR", "/app/results")).expanduser().resolve()
    stem = input_file.stem or "eis"
    rid = uuid.uuid4().hex[:12]
    out = results / "drt_analysis" / f"{stem}_{rid}"
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
    """将容器内绝对路径映射为 nginx 同源 URL（/results/ 或 /uploads/）。"""
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


def _markdown_for_drt(image_urls: List[str], json_url: Optional[str]) -> str:
    lines = [
        "### DRT（弛豫时间分布）分析结果",
        "",
        "- 以下为 **DRT 分布图** 与 **Nyquist 拟合对比**（同源 `/results/` 静态资源）。",
    ]
    if json_url:
        lines.append(f"- 数值摘要 JSON：[下载/打开]({json_url})（`drt_summary.json`）。")
        lines.append("")
    for u in image_urls:
        lines.append(f"![]({u})")
        lines.append("")
    return "\n".join(lines).strip()


@registry.register(
    name="eis_drt_analysis",
    description=(
        "对电化学阻抗谱(EIS) CSV 数据执行弛豫时间分布(DRT)分析：输入文件需包含频率与阻抗实部/虚部列"
        "（常见列名 frequency、Z_real、Z_imag；具体以底层脚本解析为准）。"
        "通过独立子进程调用 DRTtools 脚本，在 RESULTS_DIR 下生成输出目录，返回 JSON 与图像路径列表。"
        "触发场景：用户上传 EIS 数据、需要 DRT 分解、弛豫时间分布、Nyquist/Bode 解读、电化学阻抗谱建模等。"
    ),
    category="Electrochemistry",
    output_type="mixed",
)
@safe_tool_execution
def eis_drt_analysis(
    file_path: str,
    regularization_lambda: float = 0.1,
    method: str = "tikhonov",
    timeout_seconds: int = 900,
) -> Dict[str, Any]:
    """
    Args:
        file_path: 已上传到服务可访问路径的 EIS 数据文件（通常为 CSV）。
        regularization_lambda: Tikhonov 等正则强度（对应底层脚本 ``--lambda``）。
        method: 反演/正则方法名（对应底层脚本 ``--method``，默认 tikhonov）。
        timeout_seconds: 子进程超时（秒）。
    """
    fp = Path((file_path or "").strip()).expanduser()
    if not fp.is_file():
        return {"status": "error", "message": f"数据文件不存在或不可读: {fp}"}

    script = resolve_drt_analysis_script()
    if script is None:
        root = _repo_root()
        return {
            "status": "error",
            "message": (
                "未找到 drt_analysis.py。请确认镜像/挂载中含 "
                f"{root / 'gibh_agent/assets/drt/drt_analysis.py'}；"
                "或将旧版技能包解压至 "
                f"{root / 'uploads/skills_staging/2026.04.20 skills/drt-tools/scripts/'}"
                "；也可设置环境变量 DRT_ANALYSIS_SCRIPT 指向脚本绝对路径。"
            ),
        }

    out_dir = _prepare_run_output_dir(fp.resolve())
    py_exe = _pick_python_executable()

    cmd: List[str] = [
        py_exe,
        str(script),
        "--input",
        str(fp.resolve()),
        "--output",
        str(out_dir),
        "--lambda",
        str(float(regularization_lambda)),
        "--method",
        (method or "tikhonov").strip(),
    ]

    logger.info("eis_drt_analysis subprocess: %s", " ".join(cmd[:6]) + " ...")

    try:
        proc = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=max(30, int(timeout_seconds)),
            env={**os.environ, "PYTHONUNBUFFERED": "1"},
            cwd=str(out_dir),
        )
    except subprocess.TimeoutExpired:
        return {
            "status": "error",
            "message": f"DRT 子进程超时（>{timeout_seconds}s）。可适当提高 timeout_seconds 或检查数据规模。",
        }
    except OSError as e:
        return {"status": "error", "message": f"无法启动子进程: {e}"}

    stdout_tail = (proc.stdout or "")[-8000:]
    stderr_tail = (proc.stderr or "")[-8000:]
    artifacts = _collect_artifacts(out_dir)

    if proc.returncode != 0:
        err_hint = (stderr_tail or stdout_tail or "").strip()
        if len(err_hint) > 600:
            err_hint = err_hint[-600:]
        extra = f" 详情: {err_hint}" if err_hint else ""
        return {
            "status": "error",
            "message": (
                "DRT 脚本执行失败（非零退出码）。请检查 CSV 列名/编码/数值，"
                "并确认环境含 numpy、scipy、matplotlib、 pandas（或设置 DRT_PYTHON）。"
                f"{extra}"
            ),
            "data": sanitize_for_json(
                {
                    "script_path": str(script),
                    "exit_code": proc.returncode,
                    "output_dir": str(out_dir.resolve()),
                    "stdout_tail": stdout_tail,
                    "stderr_tail": stderr_tail,
                    **artifacts,
                }
            ),
        }

    image_urls, json_url = _public_urls_from_artifacts(artifacts)
    msg = (
        f"DRT 分析完成。输出目录: {out_dir.resolve()}；"
        f"JSON {len(artifacts['json_paths'])} 个，图像 {len(artifacts['image_paths'])} 个。"
    )
    return {
        "status": "success",
        "message": msg,
        "markdown": _markdown_for_drt(image_urls, json_url),
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
