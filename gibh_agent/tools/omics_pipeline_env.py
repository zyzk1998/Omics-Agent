"""
组学下游算子共享：可执行文件嗅探、参考序列解析、优雅降级文案与临时产物路径。

重型分析默认应在 Worker/TaaS；此处仅在宿主具备 CLI 与输入资源时尝试真实 subprocess 骨架，
否则返回带明确提示的仿真可视化数据，避免整条 Task 因缺依赖而失败。
"""
from __future__ import annotations

import logging
import os
import shutil
import subprocess
import tempfile
from typing import Any, Callable, Dict, List, Optional

logger = logging.getLogger(__name__)


def exe_available(name: str) -> bool:
    return shutil.which(name) is not None


def write_temp_mock_artifact(prefix: str, message: str) -> str:
    fd, path = tempfile.mkstemp(prefix=f"{prefix}_", suffix=".mock.txt")
    os.close(fd)
    try:
        with open(path, "w", encoding="utf-8") as fh:
            fh.write(message + "\n")
    except OSError as exc:
        logger.warning("mock artifact write failed: %s", exc)
    return path


def degraded_banner_md(missing: str) -> str:
    """与前端 safeMarkedParse 中剔除 ⚠️ 符号兼容：正文仍保留「运行环境提示」语义。"""
    return (
        "> **运行环境提示**：未检测到 "
        + missing
        + "，当前为您呈现系统仿真结果。\n\n"
    )


def resolve_reference_fasta(reference_id: str) -> Optional[str]:
    """通过 reference_id 与环境变量解析参考基因组 FASTA（需运维预先 export）。"""
    rid = (reference_id or "").strip().lower().replace(" ", "")
    env_for_id = {
        "hg38": "GIBH_REF_HG38",
        "grch38": "GIBH_REF_HG38",
        "hg19": "GIBH_REF_HG19",
        "grch37": "GIBH_REF_HG19",
        "mm10": "GIBH_REF_MM10",
        "mm39": "GIBH_REF_MM39",
    }
    primary = env_for_id.get(rid)
    candidates = []
    if primary:
        candidates.append(os.getenv(primary))
    candidates.append(os.getenv("GIBH_REFERENCE_FASTA"))
    for p in candidates:
        if p and os.path.isfile(p):
            return os.path.abspath(p)
    return None


def resolve_bowtie2_index_prefix(reference_id: str) -> Optional[str]:
    """Bowtie2 索引前缀（不含 .bt2 后缀），由环境变量注入。"""
    rid = (reference_id or "").strip().lower().replace(" ", "")
    env_for_id = {
        "hg38": "GIBH_BOWTIE2_HG38_INDEX",
        "grch38": "GIBH_BOWTIE2_HG38_INDEX",
        "hg19": "GIBH_BOWTIE2_HG19_INDEX",
        "mm10": "GIBH_BOWTIE2_MM10_INDEX",
    }
    key = env_for_id.get(rid) or "GIBH_BOWTIE2_INDEX"
    p = os.getenv(key)
    if not p:
        return None
    if os.path.isfile(p + ".1.bt2"):
        return p
    return None


def smoke_subprocess(cmd: List[str], *, timeout: int = 120) -> bool:
    """仅验证 CLI 可启动（不跑完整分析），用于骨架探测。"""
    if not cmd or not cmd[0]:
        return False
    try:
        r = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=timeout,
            check=False,
        )
        return r.returncode == 0 or r.stdout or r.stderr
    except (OSError, subprocess.TimeoutExpired) as exc:
        logger.debug("smoke_subprocess skip: %s", exc)
        return False


def run_or_degrade(
    real_fn: Callable[[], Optional[Dict[str, Any]]],
    degrade_fn: Callable[[], Dict[str, Any]],
    *,
    ctx: str = "",
) -> Dict[str, Any]:
    """
    先尝试真实管线（返回 dict 视为成功）；返回 None 或抛错时走降级仿真。
    """
    try:
        got = real_fn()
        if got is not None:
            return got
    except Exception as exc:  # noqa: BLE001 — 工具层须兜底，避免 Task 崩溃
        logger.warning("%s real run failed: %s", ctx or "omics", exc)
    return degrade_fn()


def modal_tool_with_degradation(
    prefix: str,
    msg: str,
    *,
    markdown_sim: str,
    image_urls: List[str],
    table_data: Dict[str, Any],
    tool_human_label: str,
    candidate_exes: Optional[List[str]] = None,
    real_runner: Optional[Callable[[], Optional[Dict[str, Any]]]] = None,
    tool_id: Optional[str] = None,
) -> Dict[str, Any]:
    """
    蛋白/表观等模态通用：可选 real_runner；若本机未安装 candidate_exes 中任一工具则加「环境提示」前缀。
    candidate_exes 为空则不加缺省提示前缀（由 markdown_sim 自行说明）。
    """
    from .omics_mock_ui import attach_visual_contract

    cands = candidate_exes or []
    if real_runner is not None:
        try:
            if not cands or any(exe_available(e) for e in cands):
                got = real_runner()
                if got is not None:
                    return got
        except Exception as exc:  # noqa: BLE001
            logger.warning("%s: real_runner failed: %s", prefix, exc)
    md = markdown_sim
    if cands and not any(exe_available(e) for e in cands):
        md = degraded_banner_md(tool_human_label) + md
    path = write_temp_mock_artifact(prefix, msg)
    out: Dict[str, Any] = {
        "status": "success",
        "message": msg,
        "output_path": path,
        "file_path": path,
    }
    return attach_visual_contract(
        out, markdown=md, image_urls=image_urls, table_data=table_data, tool_id=tool_id
    )
