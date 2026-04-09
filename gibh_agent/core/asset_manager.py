"""
多组学资产总线 — OmicsAssetManager

将零散路径规范为 DataAsset 列表（捆绑包优先于单文件），并通过 to_legacy_resolved_dict
与 resolve_omics_paths 的六键字典结构兼容，并额外包含 ``fasta`` 桶供蛋白序列下游使用。
"""
from __future__ import annotations

import logging
import os
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, FrozenSet, List, Optional, Set, Tuple

logger = logging.getLogger(__name__)

# Level-1：未显式传 lineage 时整次 classify_assets 共用一个虚拟血缘，允许跨目录配对
REQUEST_BATCH_LINEAGE_KEY = "__request_batch__"


LEGACY_BUCKET_KEYS = (
    "tables",
    "h5ad",
    "10x_mtx",
    "images",
    "masks",
    "unknown",
)

# 额外 legacy 键（resolve_omics_paths 六键之外，向下兼容扩展）
LEGACY_EXTRA_KEYS = ("fasta",)


@dataclass
class DataAsset:
    """标准化组学数据资产。"""

    asset_type: str
    primary_path: str
    associated_files: List[str] = field(default_factory=list)
    metadata: Dict[str, Any] = field(default_factory=dict)


def _suffix_lower(path: Path) -> str:
    suf = path.suffix.lower()
    if suf == ".gz" and len(path.suffixes) >= 2:
        suf = "".join(path.suffixes[-2:]).lower()
    return suf


def _name_lower(path: Path) -> str:
    return path.name.lower()


def _is_tenx_candidate(path: Path) -> bool:
    nl = _name_lower(path)
    return (
        "matrix.mtx" in nl
        or "barcodes" in nl
        or "features" in nl
        or "genes" in nl
    )


def _try_tenx_bundle(
    paths: List[Path],
    allowed_indices: Optional[Set[int]] = None,
    assigned: Optional[List[bool]] = None,
) -> Optional[Tuple[DataAsset, Set[int]]]:
    """
    若同批凑齐 matrix + (barcodes|features|genes)，返回 bundle 资产及应移除的下标集合。
    规则与 TenXPathResolver 一致。

    allowed_indices:
        非空时仅在指定全局下标内凑 TenX（血缘隔离）；None 表示整批路径均可参与。
    assigned:
        已占用下标跳过（同血缘分组内多轮抽 bundle）。
    """
    cand_idx = [
        i
        for i, p in enumerate(paths)
        if (allowed_indices is None or i in allowed_indices)
        and (assigned is None or not assigned[i])
        and _is_tenx_candidate(p)
    ]
    if not cand_idx:
        return None
    tenx_strs = [str(paths[i]) for i in cand_idx]
    has_matrix = any("matrix.mtx" in Path(s).name.lower() for s in tenx_strs)
    has_bf = any(
        "barcodes" in Path(s).name.lower()
        or "features" in Path(s).name.lower()
        or "genes" in Path(s).name.lower()
        for s in tenx_strs
    )
    if not (has_matrix and has_bf):
        return None
    dirs = [os.path.dirname(s) for s in tenx_strs]
    try:
        common = os.path.commonpath(dirs) if len(dirs) > 1 else dirs[0]
    except (ValueError, OSError):
        common = os.path.dirname(tenx_strs[0]) if tenx_strs else ""
    asset = DataAsset(
        asset_type="10x_bundle",
        primary_path=common,
        associated_files=list(tenx_strs),
        metadata={},
    )
    return asset, set(cand_idx)


RADIOMICS_SUFFIXES = (".nii", ".nii.gz", ".dcm", ".tiff", ".tif", ".nrrd")
MASK_KEYWORDS = ("mask", "roi", "seg", "label")
IMAGE_KEYWORDS = ("image", "img", "mri", "ct")


def _is_radiomics_file(path: Path) -> bool:
    nl = _name_lower(path)
    sl = _suffix_lower(path)
    return (
        sl in RADIOMICS_SUFFIXES
        or ".nii" in nl
        or nl.endswith(".nii.gz")
        or nl.endswith(".nrrd")
    )


def _is_radiomics_mask_file(path: Path) -> bool:
    nl = _name_lower(path)
    if not _is_radiomics_file(path):
        return False
    if any(kw in nl for kw in MASK_KEYWORDS) and not any(
        k in nl for k in IMAGE_KEYWORDS
    ):
        return True
    return False


_RADIOMICS_STRIP_TOKENS: Tuple[str, ...] = (
    "_mask",
    "-mask",
    "_seg",
    "-seg",
    "_segmentation",
    "-segmentation",
    "_label",
    "-label",
    "_roi",
    "-roi",
    "_mri",
    "-mri",
    "_ct",
    "-ct",
    "_pet",
    "-pet",
    "_t1",
    "-t1",
    "_t2",
    "-t2",
    "_t1w",
    "-t1w",
    "_t2w",
    "-t2w",
    "mask",
    "seg",
)


def _radiomics_strip_extensions(name: str) -> str:
    for ext in (".nii.gz", ".nrrd", ".nii", ".dcm", ".tiff", ".tif"):
        if name.endswith(ext):
            return name[: -len(ext)]
    return name


def _radiomics_identity_base(path: Path) -> str:
    """
    Level-2 强语义：迭代剥除影像/掩膜常见后缀 token，得到用于「同源同病例」比对的 identity。
    不依赖物理目录；两文件 identity 相等即候选 Bundle。
    """
    name = _radiomics_strip_extensions(_name_lower(path))
    changed = True
    while changed and name:
        changed = False
        for token in _RADIOMICS_STRIP_TOKENS:
            if name.endswith(token):
                name = name[: -len(token)].rstrip("_-")
                changed = True
                break
    return name


def _radiomics_stem_core(path: Path) -> str:
    """用于配对的弱化 stem（去掉扩展名与常见 mask/seg 后缀）；兼容旧子串匹配。"""
    name = _name_lower(path)
    for ext in (".nii.gz", ".nrrd", ".nii", ".dcm", ".tiff", ".tif"):
        if name.endswith(ext):
            name = name[: -len(ext)]
            break
    for token in _RADIOMICS_STRIP_TOKENS:
        if name.endswith(token):
            name = name[: -len(token)].rstrip("_-")
            break
    return name


def _probe_radiomics_volume_alignment(
    image_path: Path, mask_path: Path
) -> Tuple[bool, Optional[str]]:
    """
    Level-3：体数据维度探针。一致 → 可确认 radiomics_pair；不一致或不可读 → 降级为孤立单文件。
    未安装 nibabel 时跳过探针（允许配对），避免无头 CI/最小镜像全量拆对。
    返回 (ok, reason_if_not_ok)。
    """
    try:
        import nibabel as nib  # type: ignore
    except ImportError:
        logger.debug(
            "[OmicsBus] nibabel not available; skip radiomics spatial probe (pair allowed)"
        )
        return True, None
    try:
        img = nib.load(str(image_path))
        msk = nib.load(str(mask_path))
        ishape = getattr(img, "shape", None)
        mshape = getattr(msk, "shape", None)
        if ishape is None or mshape is None:
            return False, "missing shape on nibabel image"
        # 对齐空间维（3D）；允许通道维差异
        ituple = tuple(int(x) for x in ishape[:3])
        mtuple = tuple(int(x) for x in mshape[:3])
        if ituple == mtuple:
            return True, None
        return False, f"spatial_shape mismatch img{ituple} vs mask{mtuple}"
    except Exception as exc:
        logger.warning(
            "⚠️ [OmicsBus] radiomics probe failed (degrade to singletons): %s",
            exc,
            exc_info=False,
        )
        return False, str(exc)


def _radiomics_abs(p: Path) -> str:
    try:
        return str(p.resolve()) if p.exists() else str(p)
    except OSError:
        return str(p)


def _commit_radiomics_pair(
    assets: List[DataAsset],
    consumed: Set[int],
    ii: int,
    mi: int,
    paths: List[Path],
    used_images: Set[int],
    used_masks: Set[int],
    *,
    tier: str,
) -> bool:
    """探针通过则提交 radiomics_pair 并标记消费；失败则返回 False（两文件保持未配对）。"""
    img_p, mask_p = paths[ii], paths[mi]
    ok, reason = _probe_radiomics_volume_alignment(img_p, mask_p)
    if not ok:
        logger.warning(
            "⚠️ [OmicsBus] radiomics pair rejected (tier=%s): %s | img=%s mask=%s",
            tier,
            reason,
            img_p.name,
            mask_p.name,
        )
        return False
    ia, ma = _radiomics_abs(img_p), _radiomics_abs(mask_p)
    assets.append(
        DataAsset(
            asset_type="radiomics_pair",
            primary_path=ia,
            associated_files=[ia, ma],
            metadata={
                "mask_path": ma,
                "image_path": ia,
                "radiomics_probe_ok": True,
                "pairing_tier": tier,
            },
        )
    )
    consumed.add(mi)
    consumed.add(ii)
    used_masks.add(mi)
    used_images.add(ii)
    return True


def _pair_radiomics(
    indices: List[int], paths: List[Path]
) -> Tuple[List[DataAsset], Set[int]]:
    """
    同源分组内的影像配对（Level 2 + 3）。
    不依赖 path.parent 一致：跨 UUID 目录仅靠语义 + 探针确认。
    """
    assets: List[DataAsset] = []
    consumed: Set[int] = set()
    mask_idxs = [i for i in indices if _is_radiomics_mask_file(paths[i])]
    img_idxs = [i for i in indices if i not in mask_idxs and _is_radiomics_file(paths[i])]

    used_masks: Set[int] = set()
    used_images: Set[int] = set()

    # Tier A：同目录 1 图 + 1 mask（快路径，仍走探针）
    by_dir: Dict[str, List[int]] = {}
    for i in mask_idxs + img_idxs:
        d = str(paths[i].parent.resolve()) if paths[i].exists() else str(paths[i].parent)
        by_dir.setdefault(d, []).append(i)
    for _d, idxs in by_dir.items():
        ms = [i for i in idxs if i in mask_idxs]
        ims = [i for i in idxs if i in img_idxs]
        if len(ms) == 1 and len(ims) == 1:
            mi, ii = ms[0], ims[0]
            _commit_radiomics_pair(
                assets, consumed, ii, mi, paths, used_images, used_masks, tier="same_directory"
            )

    # Tier B：强语义 identity base 完全一致（跨目录）
    id_groups: Dict[str, List[int]] = defaultdict(list)
    for i in mask_idxs + img_idxs:
        if i in consumed:
            continue
        bid = _radiomics_identity_base(paths[i])
        if bid:
            id_groups[bid].append(i)
    for _bid, idxs in id_groups.items():
        if len(idxs) < 2:
            continue
        ms = [i for i in idxs if i in mask_idxs and i not in used_masks]
        ims = [i for i in idxs if i in img_idxs and i not in used_images]
        if len(ms) == 1 and len(ims) == 1:
            mi, ii = ms[0], ims[0]
            _commit_radiomics_pair(
                assets,
                consumed,
                ii,
                mi,
                paths,
                used_images,
                used_masks,
                tier="identity_base",
            )

    # Tier C：弱化 stem 子串匹配（兼容历史命名）
    for mi in mask_idxs:
        if mi in used_masks:
            continue
        core_m = _radiomics_stem_core(paths[mi])
        best_ii: Optional[int] = None
        for ii in img_idxs:
            if ii in used_images:
                continue
            core_i = _radiomics_stem_core(paths[ii])
            if core_m and core_i and (core_m == core_i or core_m in core_i or core_i in core_m):
                best_ii = ii
                break
        if best_ii is not None:
            _commit_radiomics_pair(
                assets,
                consumed,
                best_ii,
                mi,
                paths,
                used_images,
                used_masks,
                tier="stem_fuzzy",
            )

    # Tier D：分组内恰好 2 个 radiomics 且能区分 1 图 + 1 mask
    pending = [i for i in indices if i not in consumed and _is_radiomics_file(paths[i])]
    if len(pending) == 2:
        a, b = pending[0], pending[1]
        a_mask = _is_radiomics_mask_file(paths[a])
        b_mask = _is_radiomics_mask_file(paths[b])
        if a_mask ^ b_mask:
            mi, ii = (a, b) if a_mask else (b, a)
            _commit_radiomics_pair(
                assets,
                consumed,
                ii,
                mi,
                paths,
                used_images,
                used_masks,
                tier="binary_complement",
            )

    # 剩余 → 单文件资产
    for ii in img_idxs:
        if ii in used_images:
            continue
        assets.append(
            DataAsset(
                asset_type="radiomics_image",
                primary_path=_radiomics_abs(paths[ii]),
                associated_files=[_radiomics_abs(paths[ii])],
                metadata={"image_path": _radiomics_abs(paths[ii])},
            )
        )
        consumed.add(ii)
    for mi in mask_idxs:
        if mi in used_masks:
            continue
        assets.append(
            DataAsset(
                asset_type="radiomics_mask",
                primary_path=_radiomics_abs(paths[mi]),
                associated_files=[_radiomics_abs(paths[mi])],
                metadata={"mask_path": _radiomics_abs(paths[mi])},
            )
        )
        consumed.add(mi)

    return assets, consumed


def _is_table_file(path: Path) -> bool:
    nl = _name_lower(path)
    sl = _suffix_lower(path)
    return sl in (".csv", ".tsv", ".xlsx", ".txt") or nl.endswith(".tsv.gz")


def _is_h5ad_file(path: Path) -> bool:
    """仅 AnnData .h5ad；原生 10x HDF5 矩阵 .h5 见 _is_visium_matrix_h5_file / 10x_h5_matrix。"""
    return _suffix_lower(path) == ".h5ad"


def _is_visium_matrix_h5_file(path: Path) -> bool:
    """10x 官方 Visium / HD 风格 HDF5 矩阵（单文件）。"""
    if not path.is_file():
        return False
    if _suffix_lower(path) != ".h5":
        return False
    nl = path.name.lower()
    return (
        nl == "filtered_feature_bc_matrix.h5"
        or nl == "raw_feature_bc_matrix.h5"
        or nl.endswith("_feature_bc_matrix.h5")
    )


def _is_10x_mtx_directory(path: Path) -> bool:
    """目录内含 matrix.mtx / matrix.mtx.gz 等，视为 10x 表达矩阵目录（spatial_bundle 矩阵端）。"""
    if not path.is_dir():
        return False
    try:
        if not path.exists():
            return False
        for c in path.iterdir():
            if not c.is_file():
                continue
            n = c.name.lower()
            if "matrix.mtx" in n:
                return True
    except OSError:
        return False
    return False


_SPATIAL_DIR_MARKERS = (
    "tissue_positions",
    "tissue_hires_image",
    "tissue_lowres_image",
    "scalefactors",
)


def _is_spatial_image_directory(path: Path) -> bool:
    """
    空间成像侧候选：目录名含 spatial，或目录下存在典型 Visium spatial 文件。
    """
    if not path.is_dir():
        return False
    if "spatial" in path.name.lower():
        return True
    try:
        if not path.exists():
            return False
        for c in path.iterdir():
            if not c.is_file():
                continue
            cl = c.name.lower()
            for mk in _SPATIAL_DIR_MARKERS:
                if mk in cl:
                    return True
    except OSError:
        return False
    return False


def _spatial_matrix_compatible(matrix_path: Path, spatial_dir: Path) -> bool:
    """矩阵与 spatial 目录是否属于同一次上传会话（启发式）。"""
    if not spatial_dir.is_dir():
        return False
    try:
        m = matrix_path.resolve()
        s = spatial_dir.resolve()
    except OSError:
        m, s = matrix_path, spatial_dir
    if matrix_path.is_file():
        mpar = m.parent
    else:
        mpar = m
    spar = s.parent
    if mpar == spar or spar == mpar:
        return True
    if spar.parent == mpar:
        return True
    if mpar.parent == spar.parent:
        return True
    try:
        s.relative_to(mpar)
        return True
    except ValueError:
        pass
    return False


def _spatial_pair_score(matrix_path: Path, spatial_dir: Path) -> int:
    """越小越优（在已通过 compatible 的前提下）。"""
    try:
        m = matrix_path.resolve()
        s = spatial_dir.resolve()
    except OSError:
        m, s = matrix_path, spatial_dir
    if matrix_path.is_file():
        mpar = m.parent
    else:
        mpar = m
    spar = s.parent
    if mpar == spar:
        return 0
    if spar.parent == mpar:
        return 1
    if mpar.parent == spar.parent:
        return 2
    try:
        s.relative_to(mpar)
        return 3
    except ValueError:
        return 10


def _spatial_bundle_visium_root(matrix_path: Path, spatial_dir: Path) -> Path:
    """供体检：SpatialVisiumHandler 需要含 spatial/ 与矩阵的目录根。"""
    try:
        m_base = matrix_path.resolve() if matrix_path.is_dir() else matrix_path.parent.resolve()
        s = spatial_dir.resolve()
        return Path(os.path.commonpath([str(m_base), str(s)]))
    except (ValueError, OSError):
        return matrix_path.parent if matrix_path.is_file() else matrix_path


def _try_spatial_bundle(
    paths: List[Path],
    assigned: List[bool],
    allowed_indices: Optional[Set[int]] = None,
) -> Optional[Tuple[DataAsset, Set[int]]]:
    """
    矩阵（.h5 Visium 或 10x_mtx 目录）+ 空间成像目录 → spatial_bundle。
    allowed_indices: 非空时仅在同血缘分组内配对（禁止跨血缘交叉绑定）。
    """
    n = len(paths)
    matrix_idxs: List[int] = []
    spatial_idxs: List[int] = []
    for i in range(n):
        if assigned[i]:
            continue
        if allowed_indices is not None and i not in allowed_indices:
            continue
        p = paths[i]
        if p.is_file() and _is_visium_matrix_h5_file(p):
            matrix_idxs.append(i)
        elif p.is_dir() and _is_10x_mtx_directory(p):
            matrix_idxs.append(i)
        elif p.is_dir() and _is_spatial_image_directory(p):
            spatial_idxs.append(i)

    if not matrix_idxs or not spatial_idxs:
        return None

    best: Optional[Tuple[int, int, int]] = None  # (score, mi, si)
    for mi in matrix_idxs:
        mp = paths[mi]
        for si in spatial_idxs:
            sd = paths[si]
            if not _spatial_matrix_compatible(mp, sd):
                continue
            sc = _spatial_pair_score(mp, sd)
            cand = (sc, mi, si)
            if best is None or cand < best:
                best = cand

    if best is None:
        if len(matrix_idxs) == 1 and len(spatial_idxs) == 1:
            mi, si = matrix_idxs[0], spatial_idxs[0]
            if not _spatial_matrix_compatible(paths[mi], paths[si]):
                logger.warning(
                    "⚠️ [OmicsBus] spatial_bundle: 1 matrix + 1 spatial but "
                    "directory-heuristic incompatible; still bundling (cross-UUID / lineage)."
                )
        else:
            return None
    else:
        _sc, mi, si = best

    mpath = paths[mi]
    spath = paths[si]
    m_str = str(mpath.resolve()) if mpath.exists() else str(mpath)
    s_str = str(spath.resolve()) if spath.exists() else str(spath)
    vroot = _spatial_bundle_visium_root(mpath, spath)
    v_str = str(vroot.resolve()) if vroot.exists() else str(vroot)
    asset = DataAsset(
        asset_type="spatial_bundle",
        primary_path=m_str,
        associated_files=[m_str, s_str],
        metadata={
            "spatial_dir": s_str,
            "visium_root": v_str,
        },
    )
    return asset, {mi, si}


PROTEIN_FASTA_SUFFIXES = frozenset({".fasta", ".fa", ".faa", ".ffn"})


def _is_protein_fasta_file(path: Path) -> bool:
    return _suffix_lower(path) in PROTEIN_FASTA_SUFFIXES


# ---------------------------------------------------------------------------
# 动态嗅探：后缀 → 资产类型（配置化，扩展时只改此字典）
# 说明：.txt 仍由上方 _is_table_file 优先归为 metabolomics_csv，不在此映射中。
# ---------------------------------------------------------------------------
_DYNAMIC_ASSET_SUFFIXES_BY_TYPE: Dict[str, Tuple[str, ...]] = {
    "protein_structure": (
        ".pdb",
        ".cif",
        ".mmcif",
        ".ent",
        ".pdb.gz",
        ".cif.gz",
        ".mmcif.gz",
    ),
    "document": (
        ".pdf",
        ".docx",
        ".doc",
        ".pptx",
        ".ppt",
        ".odt",
        ".rtf",
    ),
    "plain_text": (
        ".json",
        ".jsonl",
        ".yaml",
        ".yml",
        ".md",
        ".log",
        ".xml",
        ".toml",
        ".ini",
        ".cfg",
    ),
}


def _build_suffix_to_asset_type_map() -> Dict[str, str]:
    m: Dict[str, str] = {}
    for atype, suffixes in _DYNAMIC_ASSET_SUFFIXES_BY_TYPE.items():
        for suf in suffixes:
            m[suf] = atype
    return m


SUFFIX_TO_ASSET_TYPE: Dict[str, str] = _build_suffix_to_asset_type_map()


def sniff_asset_type_from_path(path: Path) -> str:
    """按扩展名嗅探资产类型；无法识别则 generic_unknown（不再使用单一 unknown_file）。"""
    sl = _suffix_lower(path)
    if sl in SUFFIX_TO_ASSET_TYPE:
        return SUFFIX_TO_ASSET_TYPE[sl]
    return "generic_unknown"


# 工作流 DAG 原生资产（用于全局路由：与技能/对话工具区分）
ROUTING_WORKFLOW_NATIVE_ASSET_TYPES: FrozenSet[str] = frozenset(
    {
        "10x_bundle",
        "spatial_bundle",
        "10x_h5_matrix",
        "h5ad_single",
        "metabolomics_csv",
        "radiomics_pair",
        "radiomics_image",
        "radiomics_mask",
    }
)

# 非工作流管线典型资产（易误判进 RNA/影像等 DAG 时，供 orchestrator / 提示词参考）
ROUTING_NON_WORKFLOW_ASSET_TYPES: FrozenSet[str] = frozenset(
    {
        "protein_structure",
        "protein_fasta",
        "document",
        "plain_text",
        "generic_unknown",
    }
)


def routing_asset_inventory(assets: List[DataAsset]) -> List[Dict[str, Any]]:
    """供 LLM / JSON 注入：精确资产类型 + 文件名（不含硬编码业务分支）。"""
    out: List[Dict[str, Any]] = []
    for a in assets:
        out.append(
            {
                "asset_type": a.asset_type,
                "file_name": Path(a.primary_path).name,
                "primary_path": a.primary_path,
            }
        )
    return out


def format_routing_asset_digest(assets: List[DataAsset]) -> str:
    """人类可读摘要，用于 system/user 提示。"""
    if not assets:
        return "(本次路径列表未解析出任何资产)"
    lines: List[str] = []
    for a in assets:
        lines.append(
            f"- asset_type={a.asset_type} | file={Path(a.primary_path).name!r}"
        )
    return "\n".join(lines)


def _normalize_lineage_keys(
    n: int, lineage_group_ids: Optional[List[Optional[str]]]
) -> List[str]:
    """
    Level-1 血缘分组：与 raw_file_paths 等长。
    - 未传或长度不符 → 全部为 REQUEST_BATCH_LINEAGE_KEY（单次请求单一虚拟池，可跨物理目录配对）。
    - 显式传入时：同非空字符串同组；空串/None → ``__default__``。
    """
    if lineage_group_ids is None or len(lineage_group_ids) != n:
        return [REQUEST_BATCH_LINEAGE_KEY] * n
    out: List[str] = []
    for x in lineage_group_ids:
        s = (x or "").strip()
        out.append(s if s else "__default__")
    return out


def _classify_lineage_group(
    paths: List[Path],
    gidx: List[int],
    assigned: List[bool],
    assets: List[DataAsset],
) -> None:
    """对单血缘组内下标执行 TenX → Spatial → Radiomics → 单文件阶梯。"""
    allowed: Set[int] = set(gidx)

    while True:
        tenx = _try_tenx_bundle(paths, allowed_indices=allowed, assigned=assigned)
        if not tenx:
            break
        bundle, idxs = tenx
        if not idxs <= allowed:
            break
        assets.append(bundle)
        for i in idxs:
            assigned[i] = True

    while True:
        sb = _try_spatial_bundle(paths, assigned, allowed_indices=allowed)
        if not sb:
            break
        bundle, idxs = sb
        if not idxs <= allowed:
            break
        assets.append(bundle)
        for i in idxs:
            assigned[i] = True

    pending_rad = [
        i for i in gidx if not assigned[i] and _is_radiomics_file(paths[i])
    ]
    if pending_rad:
        rad_assets, rad_consumed = _pair_radiomics(pending_rad, paths)
        assets.extend(rad_assets)
        for i in rad_consumed:
            assigned[i] = True

    for i in gidx:
        if assigned[i]:
            continue
        if _is_table_file(paths[i]):
            assets.append(
                DataAsset(
                    asset_type="metabolomics_csv",
                    primary_path=str(paths[i]),
                    associated_files=[str(paths[i])],
                    metadata={},
                )
            )
            assigned[i] = True

    for i in gidx:
        if assigned[i]:
            continue
        if _is_h5ad_file(paths[i]):
            assets.append(
                DataAsset(
                    asset_type="h5ad_single",
                    primary_path=str(paths[i]),
                    associated_files=[str(paths[i])],
                    metadata={},
                )
            )
            assigned[i] = True

    for i in gidx:
        if assigned[i]:
            continue
        if _is_visium_matrix_h5_file(paths[i]):
            assets.append(
                DataAsset(
                    asset_type="10x_h5_matrix",
                    primary_path=str(paths[i]),
                    associated_files=[str(paths[i])],
                    metadata={},
                )
            )
            assigned[i] = True

    for i in gidx:
        if assigned[i]:
            continue
        if _is_protein_fasta_file(paths[i]):
            assets.append(
                DataAsset(
                    asset_type="protein_fasta",
                    primary_path=str(paths[i]),
                    associated_files=[str(paths[i])],
                    metadata={},
                )
            )
            assigned[i] = True

    for i in gidx:
        if assigned[i]:
            continue
        sniffed = sniff_asset_type_from_path(paths[i])
        assets.append(
            DataAsset(
                asset_type=sniffed,
                primary_path=str(paths[i]),
                associated_files=[str(paths[i])],
                metadata={"sniff_rule": "suffix_map"},
            )
        )
        assigned[i] = True


def _empty_legacy() -> Dict[str, List[str]]:
    out = {k: [] for k in LEGACY_BUCKET_KEYS}
    for k in LEGACY_EXTRA_KEYS:
        out[k] = []
    return out


def to_legacy_resolved_dict(assets: List[DataAsset]) -> Dict[str, List[str]]:
    """
    将 DataAsset 列表转为与 resolve_omics_paths 对齐的桶字典，并含 ``fasta`` 列表。

    六键部分与 path_resolvers.resolve_omics_paths 取值语义一致；``fasta`` 为扩展桶。
    """
    out = _empty_legacy()
    for a in assets:
        t = a.asset_type
        if t == "10x_bundle":
            out["10x_mtx"].append(a.primary_path)
        elif t == "spatial_bundle":
            primary = Path(a.primary_path)
            if primary.is_file() and _suffix_lower(primary) == ".h5":
                out["h5ad"].append(a.primary_path)
            else:
                out["10x_mtx"].append(a.primary_path)
        elif t == "10x_h5_matrix":
            out["h5ad"].append(a.primary_path)
        elif t == "metabolomics_csv":
            out["tables"].append(a.primary_path)
        elif t == "h5ad_single":
            out["h5ad"].append(a.primary_path)
        elif t == "protein_fasta":
            out.setdefault("fasta", [])
            out["fasta"].append(a.primary_path)
        elif t == "radiomics_pair":
            out["images"].append(a.primary_path)
            mp = a.metadata.get("mask_path")
            if mp:
                out["masks"].append(str(mp))
        elif t == "radiomics_image":
            out["images"].append(a.primary_path)
        elif t == "radiomics_mask":
            out["masks"].append(a.primary_path)
        elif t in (
            "generic_unknown",
            "unknown_file",
            "protein_structure",
            "document",
            "plain_text",
        ):
            out["unknown"].append(a.primary_path)
        else:
            out["unknown"].append(a.primary_path)
    return out


class OmicsAssetManager:
    """
    接收规范化路径列表，产出 DataAsset 列表（捆绑包优先）。
    不在此阶段拼接逗号字符串；仅使用 List[str] 入参。

    lineage_group_ids（可选）与路径等长时启用多血缘隔离；省略时整次请求为单一虚拟池
    （BATCH_LINEAGE_KEY），允许跨 UUID 物理目录的 radiomics / 空间配对。
    """

    BATCH_LINEAGE_KEY = REQUEST_BATCH_LINEAGE_KEY

    def classify_assets(
        self,
        raw_file_paths: List[str],
        user_intent_domain: Optional[str] = None,
        lineage_group_ids: Optional[List[Optional[str]]] = None,
    ) -> List[DataAsset]:
        paths = [Path(p.strip()) for p in raw_file_paths if p and str(p).strip()]
        if not paths:
            return []

        n = len(paths)
        lineage_keys = _normalize_lineage_keys(n, lineage_group_ids)
        groups: Dict[str, List[int]] = defaultdict(list)
        for i, lk in enumerate(lineage_keys):
            groups[lk].append(i)

        assigned = [False] * n
        assets: List[DataAsset] = []

        for _gk in sorted(groups.keys()):
            gidx = sorted(groups[_gk])
            _classify_lineage_group(paths, gidx, assigned, assets)

        if user_intent_domain:
            assets = self._sort_by_intent(assets, user_intent_domain)
        return assets

    def _sort_by_intent(self, assets: List[DataAsset], domain: str) -> List[DataAsset]:
        d = domain.strip().lower()
        priority = {
            "radiomics": ("radiomics_pair", "radiomics_image", "radiomics_mask"),
            "metabolomics": ("metabolomics_csv",),
            "spatial": ("spatial_bundle", "10x_bundle", "10x_h5_matrix", "h5ad_single"),
            "scrna": ("spatial_bundle", "10x_bundle", "10x_h5_matrix", "h5ad_single"),
            "sc_rna": ("spatial_bundle", "10x_bundle", "10x_h5_matrix", "h5ad_single"),
            "10x": ("spatial_bundle", "10x_bundle", "10x_h5_matrix", "h5ad_single"),
            "protein": ("protein_structure", "protein_fasta"),
        }
        pref = priority.get(d, ())

        def key(a: DataAsset) -> Tuple[int, str]:
            try:
                rank = pref.index(a.asset_type)
            except ValueError:
                rank = len(pref) + 1
            return (rank, a.asset_type)

        return sorted(assets, key=key)

    def to_legacy_resolved_dict(self, assets: List[DataAsset]) -> Dict[str, List[str]]:
        return to_legacy_resolved_dict(assets)
