"""
全组学 I/O 注册表：扩展名字典、格式提示、Basename 模糊路径自愈。

供 path_resolver / asset_locator / tool_input_validator / file_inspector 统一引用，
避免各模块硬编码狭隘白名单。
"""
from __future__ import annotations

import logging
import os
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Set, Tuple

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# 全组学合法扩展名（复合后缀优先匹配）
# ---------------------------------------------------------------------------

OMICS_COMPOUND_EXTENSIONS: frozenset[str] = frozenset(
    {
        # 测序原始/比对
        ".fastq.gz",
        ".fq.gz",
        ".sam.gz",
        ".bam.bai",  # 索引，仅识别不阻断
        ".vcf.gz",
        ".bed.gz",
        # 单细胞 / 空间 / 表格
        ".mtx.gz",
        ".tsv.gz",
        ".csv.gz",
        ".tar.gz",
        # 影像 / 代谢 / 蛋白
        ".nii.gz",
        ".mzxml.gz",
        ".mzml.gz",
    }
)

OMICS_SIMPLE_EXTENSIONS: frozenset[str] = frozenset(
    {
        # 测序原始/比对
        ".fastq",
        ".fq",
        ".sam",
        ".bam",
        ".cram",
        ".crai",
        ".bai",
        ".vcf",
        ".bed",
        ".bw",
        ".bigwig",
        ".gff",
        ".gff3",
        ".gtf",
        # 单细胞 / 空间
        ".h5ad",
        ".h5",
        ".mtx",
        ".loom",
        ".tsv",
        ".csv",
        ".xlsx",
        ".parquet",
        ".rds",
        ".rdata",
        # 影像
        ".nii",
        ".dcm",
        ".dicom",
        ".tif",
        ".tiff",
        ".png",
        ".jpg",
        ".jpeg",
        ".jpe",
        ".jfif",
        ".gif",
        ".bmp",
        ".webp",
        ".svg",
        ".ico",
        ".heic",
        ".heif",
        ".avif",
        ".ppm",
        ".pgm",
        ".pbm",
        # 代谢 / 蛋白
        ".mzml",
        ".mzxml",
        ".raw",
        ".mgf",
        ".mgf.gz",
        # 结构 / 序列
        ".pdb",
        ".cif",
        ".mmcif",
        ".fasta",
        ".fa",
        ".faa",
        # 归档
        ".zip",
        ".tar",
        ".gz",
        ".tgz",
        ".7z",
        # 文档（只警告不阻断）
        ".pdf",
        ".json",
        ".txt",
        ".md",
    }
)

# 明确禁止作为「主数据入口」的扩展（仍允许 warning 后继续，不 hard error）
DISCOURAGED_DATA_EXTENSIONS: frozenset[str] = frozenset({".txt", ".log", ".md"})

# 弹性跳过工具标记（与 Executor 对齐，供 ToolRegistry 软校验）
ELASTIC_SKIP_TOOL_MARKERS: Tuple[str, ...] = (
    "_data_validation",
    "data_validation",
    "preview",
    "nifti_preview",
    "clustering_comparison",
    "model_comparison",
)

# HITL 工具：image_path 可为 http(s) URL，不做容器内路径存在性校验
PATH_VALIDATION_EXEMPT_TOOLS: frozenset[str] = frozenset(
    {
        "Trigger_Expert_Annotation",
    }
)


def is_elastic_skip_tool(tool_id: Optional[str]) -> bool:
    tid = (tool_id or "").strip().lower()
    if not tid:
        return False
    return any(m in tid for m in ELASTIC_SKIP_TOOL_MARKERS)


def get_path_suffix(path_str: str) -> str:
    """提取用于分类的后缀（支持 .fastq.gz / .nii.gz 等复合后缀）。"""
    name = Path(str(path_str or "").replace("\\", "/")).name.lower()
    for compound in sorted(OMICS_COMPOUND_EXTENSIONS, key=len, reverse=True):
        if name.endswith(compound):
            return compound
    return Path(name).suffix.lower()


def is_known_omics_extension(path_str: str) -> bool:
    suf = get_path_suffix(path_str)
    if not suf:
        return False
    return suf in OMICS_COMPOUND_EXTENSIONS or suf in OMICS_SIMPLE_EXTENSIONS


def is_upload_allowed_filename(filename: str) -> bool:
    """HTTP 上传白名单（server.py）：组学数据 + 全量常见图片 + 语料 json/txt/md。"""
    return is_known_omics_extension(filename)


def upload_rejection_hint(filename: str) -> str:
    """上传被拒时返回给前端的简短说明。"""
    suf = get_path_suffix(filename) or Path(str(filename or "")).suffix.lower() or "(无后缀)"
    return (
        f"不允许的文件类型: {suf}。"
        "支持：组学矩阵/测序（h5ad、mtx、fastq.gz 等）、医学影像（nii、dcm、tif）、"
        "图片（png、jpg、jpeg、gif、webp、bmp、svg、tif 等）及语料 json/txt/md。"
    )


def path_looks_like_data_file(path_str: str) -> bool:
    """路径是否像数据文件（含扩展名或 uploads 语义）。"""
    raw = str(path_str or "").strip()
    if not raw:
        return False
    if raw.startswith("/app/uploads") or raw.startswith("uploads/"):
        return True
    base = Path(raw.replace("\\", "/")).name
    if base and "." not in base:
        return True  # 10x 目录等
    return bool(get_path_suffix(raw))


def get_extension_warning(path_str: str) -> Optional[str]:
    """
    未知后缀 → 仅返回 warning 文案；已知或目录风格 → None。
    绝不返回 blocking error。
    """
    raw = str(path_str or "").strip()
    if not raw:
        return None
    base = Path(raw.replace("\\", "/")).name
    if base and "." not in base:
        return None
    suf = get_path_suffix(raw)
    if not suf:
        return None
    if is_known_omics_extension(raw):
        if suf in DISCOURAGED_DATA_EXTENSIONS:
            return f"扩展名 {suf} 通常不是主数据格式，已继续处理。"
        return None
    return (
        f"未在组学扩展名字典中识别到后缀 {suf!r}，已放行并尝试读取；"
        f"若后续步骤失败请检查文件格式。{get_omics_format_hint(short=True)}"
    )


def get_omics_format_hint(*, short: bool = False) -> str:
    """动态格式提示，替代 orchestrator / path_resolver 硬编码长串。"""
    if short:
        return "常见格式：FASTQ/BAM/VCF、.h5ad/10x、.nii/.dcm、.mzML/.csv 等。"
    return (
        "常见组学格式：测序原始/比对 (.fastq/.fastq.gz/.fq/.bam/.sam/.vcf)；"
        "单细胞/空间 (.h5ad/.h5/.mtx/.loom/.tsv/.csv 或 10x 目录)；"
        "影像 (.nii/.nii.gz/.dcm)；代谢/蛋白 (.mzML/.mzXML/.raw)；"
        "未识别后缀将仅警告，不阻断流程。"
    )


def _basename_key(path_str: str) -> str:
    return Path(str(path_str or "").replace("\\", "/")).name.lower()


def _path_usable(p: Path, *, allow_dir: bool) -> bool:
    try:
        if not p.exists():
            return False
        if p.is_file():
            return True
        return allow_dir and p.is_dir()
    except OSError:
        return False


def _collect_basename_index(
    context_paths: Optional[Sequence[str]],
    uploaded_entries: Optional[Sequence[Dict[str, Any]]],
) -> Dict[str, List[str]]:
    """basename(lower) → 候选绝对路径列表（保序）。"""
    index: Dict[str, List[str]] = {}
    seen: Set[str] = set()

    def _add(raw: object) -> None:
        if raw is None:
            return
        s = str(raw).strip()
        if not s:
            return
        key = s.rstrip("/")
        if key in seen:
            return
        seen.add(key)
        bn = _basename_key(s)
        if not bn:
            return
        index.setdefault(bn, []).append(s)

    for p in context_paths or []:
        _add(p)
    for entry in uploaded_entries or []:
        if not isinstance(entry, dict):
            continue
        _add(entry.get("path") or entry.get("file_path"))
        _add(entry.get("name") or entry.get("file_name"))
        if entry.get("group_dir"):
            _add(entry["group_dir"])

    return index


def fuzzy_heal_path_by_basename(
    raw: str,
    *,
    upload_dir: Optional[Path] = None,
    owner_id: Optional[str] = None,
    context_paths: Optional[Sequence[str]] = None,
    uploaded_entries: Optional[Sequence[Dict[str, Any]]] = None,
    allow_dir: bool = False,
) -> Tuple[Optional[str], List[str]]:
    """
    Basename 降维模糊自愈：绝对路径不存在时，按文件名在上下文与 owner 目录内找回真实路径。

    Returns:
        (healed_absolute_path_or_none, attempted_paths)
    """
    attempted: List[str] = []
    s = str(raw or "").strip()
    if not s:
        return None, attempted

    p = Path(s)
    # 1) 直查原路径
    try:
        rp = p.resolve()
        attempted.append(str(rp))
        if _path_usable(rp, allow_dir=allow_dir):
            return str(rp), attempted
    except OSError:
        attempted.append(s)

    basename = p.name or Path(s.replace("\\", "/")).name
    if not basename or len(basename) < 2:
        return None, attempted

    bn_lower = basename.lower()

    # 2) 上下文 uploaded_files / omics_resolved 等同名匹配
    index = _collect_basename_index(context_paths, uploaded_entries)
    for cand in index.get(bn_lower, []):
        attempted.append(cand)
        cp = Path(cand)
        try:
            cr = cp.resolve()
            if _path_usable(cr, allow_dir=allow_dir):
                logger.info("[FuzzyHeal] 上下文 basename 命中: %s -> %s", s, cr)
                return str(cr), attempted
        except OSError:
            if _path_usable(cp, allow_dir=allow_dir):
                logger.info("[FuzzyHeal] 上下文 basename 命中: %s -> %s", s, cp)
                return str(cp.resolve() if cp.exists() else cp), attempted

    # 3) owner 目录下 rglob（严禁全局 uploads 盲搜）
    if upload_dir is not None and owner_id and str(owner_id).strip():
        owner_root = (upload_dir / str(owner_id).strip()).resolve()
        attempted.append(str(owner_root / "**" / basename))
        try:
            if owner_root.exists():
                for hit in owner_root.rglob(basename):
                    if _path_usable(hit, allow_dir=allow_dir):
                        hit_s = str(hit.resolve())
                        attempted.append(hit_s)
                        logger.info("[FuzzyHeal] owner rglob 命中: %s -> %s", s, hit_s)
                        return hit_s, attempted
        except OSError as exc:
            logger.debug("[FuzzyHeal] owner rglob 失败: %s", exc)

    # 4) 无 owner 时：仅在 upload_dir 一级子目录下按 basename 扫描（不深扫全树）
    if upload_dir is not None:
        ud = upload_dir.resolve()
        attempted.append(str(ud / basename))
        flat = ud / basename
        if _path_usable(flat, allow_dir=allow_dir):
            return str(flat.resolve()), attempted
        try:
            if ud.is_dir():
                for child in ud.iterdir():
                    if not child.is_dir():
                        continue
                    candidate = child / basename
                    attempted.append(str(candidate))
                    if _path_usable(candidate, allow_dir=allow_dir):
                        logger.info("[FuzzyHeal] upload 子目录命中: %s -> %s", s, candidate)
                        return str(candidate.resolve()), attempted
        except OSError:
            pass

    return None, attempted


def heal_web_path_if_missing(
    raw: str,
    upload_dir: str,
    *,
    owner_id: Optional[str] = None,
    context_paths: Optional[Sequence[str]] = None,
    uploaded_entries: Optional[Sequence[Dict[str, Any]]] = None,
    allow_dir: bool = True,
) -> Tuple[Optional[str], List[str]]:
    """FileInspector / verify 链路专用：容器 Web 路径自愈包装。"""
    ud = Path(os.path.expanduser(upload_dir))
    return fuzzy_heal_path_by_basename(
        raw,
        upload_dir=ud,
        owner_id=owner_id,
        context_paths=context_paths,
        uploaded_entries=uploaded_entries,
        allow_dir=allow_dir,
    )
