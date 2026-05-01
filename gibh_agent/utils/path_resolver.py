"""
三态路径解析：Web（容器 uploads）/ Local（Windows 客户端绝对路径）/ HPC（集群 POSIX 前缀）。

根因修复：在 Linux 容器内 pathlib.Path('C:/Users/...').is_absolute() == False，
会被误拼接 UPLOAD_DIR，产生诸如 /app/uploads/C:/Users/.../matrix.mtx/matrix.mtx 的畸形路径。
"""
from __future__ import annotations

import logging
import os
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, Literal, Optional, Tuple

logger = logging.getLogger(__name__)

Kind = Literal["web", "local", "hpc"]


@dataclass(frozen=True)
class ResolvedMount:
    """resolve_real_path 的返回值。"""

    kind: Kind
    original: str
    skip_local_fs_check: bool
    """True：不要在 Docker 内对磁盘做存在性假设（Local/HPC 语义）。"""
    web_absolute_path: Optional[str] = None
    """仅 kind=='web'：容器内可读绝对路径（已拼接 UPLOAD_DIR）。"""


def _hpc_prefixes() -> Tuple[str, ...]:
    raw = os.getenv("OMICS_HPC_PATH_PREFIXES", "/public,/scratch,/work,/gpfs,/lustre,/mnt/hpc")
    parts = tuple(p.strip() for p in raw.split(",") if p.strip())
    return parts if parts else ("/public",)


def is_windows_abs_path(path_str: str) -> bool:
    """兼容 Linux 容器：不依赖 pathlib.is_absolute()。"""
    s = (path_str or "").strip()
    if len(s) >= 3 and s[0].isalpha() and s[1] == ":" and s[2] in "\\/":
        return True
    if s.startswith("\\\\"):
        return True
    return False


def is_hpc_style_path(path_str: str) -> bool:
    """POSIX 集群路径：不以客户端盘符形式出现。"""
    u = path_str.replace("\\", "/").strip()
    if not u.startswith("/"):
        return False
    if is_windows_abs_path(path_str):
        return False
    for pref in _hpc_prefixes():
        if u.startswith(pref.rstrip("/") + "/") or u == pref.rstrip("/"):
            return True
    return False


def normalize_duplicate_tail_filename(path_str: str) -> str:
    """
    修复错误 join 导致的尾部重复，例如：
    .../matrix.mtx/matrix.mtx 或 ...\\matrix.mtx\\matrix.mtx
    """
    if not path_str:
        return path_str
    s = path_str.strip()
    use_win_backslash = "\\" in s and not s.startswith("\\\\")
    n = s.replace("\\", "/")
    # 反复折叠末尾重复段：parent/name/name -> parent/name
    dup_pat = re.compile(r"(.+)/([^/]+)/\2$")
    guard = 0
    while guard < 16:
        guard += 1
        n2 = dup_pat.sub(r"\1/\2", n)
        if n2 == n:
            break
        n = n2
    if use_win_backslash and len(n) >= 3 and n[1] == ":" and n[0].isalpha():
        return n[0] + ":" + n[2:].replace("/", "\\")
    return n


def resolve_real_path(path_str: str, upload_dir: str = "/app/uploads") -> ResolvedMount:
    """
    - Web：相对路径 → 拼接到 upload_dir；已在容器内的 POSIX 绝对路径 → 直接使用。
    - Local：Windows 绝对路径 → 禁止与 UPLOAD_DIR 拼接。
    - HPC：符合 OMICS_HPC_PATH_PREFIXES 的 POSIX 路径 → 禁止与 UPLOAD_DIR 拼接。
    """
    original = normalize_duplicate_tail_filename((path_str or "").strip())
    if not original:
        raise ValueError("empty path")

    if is_windows_abs_path(original):
        return ResolvedMount(
            kind="local",
            original=original,
            skip_local_fs_check=True,
            web_absolute_path=None,
        )

    if is_hpc_style_path(original):
        return ResolvedMount(
            kind="hpc",
            original=original,
            skip_local_fs_check=True,
            web_absolute_path=None,
        )

    up = Path(upload_dir).expanduser()
    p = Path(original)

    try:
        if p.is_absolute():
            return ResolvedMount(
                kind="web",
                original=original,
                skip_local_fs_check=False,
                web_absolute_path=str(p.expanduser().resolve()),
            )
    except (OSError, ValueError) as e:
        logger.debug("resolve_real_path absolute resolve fallback: %s", e)

    rel = original.lstrip("/\\")
    joined = (up / rel).resolve()
    return ResolvedMount(
        kind="web",
        original=original,
        skip_local_fs_check=False,
        web_absolute_path=str(joined),
    )


def resolve_upload_path_for_container(raw: str, upload_dir: str) -> str:
    """编排器 / Agent：统一得到用于容器内访问的路径字符串；Local/HPC 返回原始串。"""
    r = resolve_real_path(normalize_duplicate_tail_filename((raw or "").strip()), upload_dir)
    if r.kind == "web":
        return r.web_absolute_path or ""
    return r.original


# --- 格式门禁（仅路径字符串，不读盘）---

_FORBIDDEN_EXT = {".txt"}
_ALLOWED_EXT = {
    ".h5ad",
    ".h5",
    ".mtx",
    ".csv",
    ".tsv",
    ".xlsx",
    ".parquet",
    ".vcf",
    ".bam",
    ".sam",
    ".bed",
    ".bw",
    ".bigwig",
    ".gz",
    ".nii",
    ".dcm",
    ".tif",
    ".tiff",
    ".png",
    ".jpg",
    ".jpeg",
    ".pdf",
    ".zip",
    ".tar",
    ".tar.gz",
}

_COMPUND_OK = (
    ".nii.gz",
    ".mtx.gz",
    ".tsv.gz",
    ".csv.gz",
    ".fq.gz",
    ".fastq.gz",
    ".tar.gz",
)


def _suffix_for_inspector_gate(path_str: str) -> str:
    """用于白名单判断的后缀（尽量贴近 pathlib 的多后缀语义）。"""
    name = Path(path_str.replace("\\", "/")).name.lower()
    for compound in _COMPUND_OK:
        if name.endswith(compound):
            return compound
    return Path(name).suffix.lower()


def validate_inspector_file_format(path_str: str) -> Optional[Dict[str, Any]]:
    """
    第一步：仅解析后缀，不访问磁盘。
    返回 None 表示放行；否则返回 FileInspector 风格的 error dict。
    """
    if not path_str or not str(path_str).strip():
        return None
    raw = str(path_str).strip()
    # 无前缀点号的「目录风格」名称（如 my_run_folder）放行给 10x 目录
    base = Path(raw.replace("\\", "/")).name
    if base and "." not in base:
        return None

    suf = _suffix_for_inspector_gate(raw)
    if suf in _FORBIDDEN_EXT:
        return {
            "status": "error",
            "success": False,
            "error": (
                "不支持该文件扩展名（已禁止 .txt 作为数据入口）。\n\n"
                "支持格式：单细胞 RNA 请上传 .h5ad 或解压后的 10x 目录（含 matrix.mtx、barcodes.tsv、features.tsv）；"
                "代谢组学请上传 CSV；影像组学请上传 .nii / .nii.gz / .dcm。"
            ),
            "file_type": "unknown",
            "file_path": raw,
        }
    if not suf:
        return None
    if suf in _ALLOWED_EXT or suf in _COMPUND_OK:
        return None
    # 未知后缀：保持硬编码提示风格，与 orchestrator SSE 文案对齐
    return {
        "status": "error",
        "success": False,
        "error": (
            f"不支持的文件扩展名: {suf}\n\n"
            "支持格式：单细胞 RNA 请上传 .h5ad 或解压后的 10x 目录（含 matrix.mtx、barcodes.tsv、features.tsv）；"
            "代谢组学请上传 CSV；影像组学请上传 .nii / .nii.gz / .dcm。"
        ),
        "file_type": "unknown",
        "file_path": raw,
    }


def _local_sidecar_base() -> str:
    return os.getenv("OMICS_LOCAL_SIDECAR_URL", "http://127.0.0.1:8019").rstrip("/")


def verify_path_exists_after_resolve(resolved: ResolvedMount, timeout: float = 8.0) -> Tuple[bool, str]:
    """
    第二步：存在性校验。
    Web：容器内 os.path.exists。
    Local：HTTP GET 本地 Sidecar /api/tools/check_file。
    HPC：容器内 exists（若集群挂载可见）；否则失败提示。
    """
    if resolved.kind == "web":
        wp = resolved.web_absolute_path
        if not wp:
            return False, "内部错误：web 路径未解析"
        ok = Path(wp).exists()
        return (
            ok,
            ""
            if ok
            else (
                f"File or directory not found: '{resolved.original}'\n\n"
                f"**Upload directory (configured):** {os.getenv('UPLOAD_DIR', '/app/uploads')}"
            ),
        )

    if resolved.kind == "local":
        try:
            import requests
        except ImportError:
            return False, "requests 未安装，无法呼叫本地侧车校验路径"

        url = f"{_local_sidecar_base()}/api/tools/check_file"
        try:
            r = requests.get(url, params={"path": resolved.original}, timeout=timeout)
            if r.status_code == 404:
                return False, (
                    f"本地文件不存在（侧车校验）: '{resolved.original}'"
                )
            if r.status_code >= 400:
                return False, (
                    f"本地侧车校验失败 (HTTP {r.status_code}): {resolved.original}"
                )
            data = r.json()
            if data.get("exists"):
                return True, ""
            return False, f"本地文件不存在: '{resolved.original}'"
        except Exception as e:
            logger.warning("本地侧车 check_file 失败: %s", e)
            return (
                False,
                "无法连接本地侧车进行路径校验（请确认桌面客户端已启动并已挂载工作区）。"
                f" 详情: {e}",
            )

    # hpc
    try:
        ok = Path(resolved.original).exists()
        return (
            ok,
            ""
            if ok
            else (
                f"路径在当前服务端不可访问（请确认集群挂载或路径是否正确）: '{resolved.original}'"
            ),
        )
    except Exception as e:
        return False, str(e)


def infer_loose_file_type_from_path_string(path_str: str) -> str:
    """无法在容器内 stat 时，仅凭路径字符串推断类型（供远程挂载占位体检）。"""
    lower = path_str.lower().replace("\\", "/")
    if lower.endswith(".h5ad") or lower.endswith(".h5"):
        return "h5ad"
    if "matrix.mtx" in lower:
        return "10x_mtx"
    if lower.endswith(".csv"):
        return "tabular"
    if ".nii" in lower or lower.endswith(".dcm"):
        return "radiomics"
    if lower.endswith(".mtx"):
        return "mtx"
    return "unknown"
