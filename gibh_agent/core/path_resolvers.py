"""
组学路径解析：策略模式（替代 resolve_omics_paths 巨型 if-else）。

执行顺序（捆绑包优先于单文件语义）：
1. TenXPathResolver — 同批凑齐 matrix + barcodes/features/genes 时写入 10x_mtx 并消费路径；否则不占用 .tsv，留给表格策略。
2. RadiomicsPathResolver
3. MetabolomicsPathResolver — .csv / .tsv / .xlsx / .txt 等
4. AnnDataPathResolver — .h5ad / .h5
5. 门面将仍未分配的路径写入 unknown

门面：OmicsPathResolver、resolve_omics_paths（返回六键结构不变）。
"""
from __future__ import annotations

import os
from abc import ABC, abstractmethod
from pathlib import Path
from typing import Dict, List, Optional, Union


class OmicsPathContext:
    """单次解析的可变上下文：路径列表 + 分配标记 + 六类输出桶。"""

    __slots__ = ("paths", "n", "assigned", "out")

    def __init__(self, resolved_paths: List[Path]) -> None:
        self.paths = resolved_paths
        self.n = len(resolved_paths)
        self.assigned = [False] * self.n
        self.out: Dict[str, List[str]] = {
            "tables": [],
            "h5ad": [],
            "10x_mtx": [],
            "images": [],
            "masks": [],
            "unknown": [],
        }

    def path_str(self, i: int) -> str:
        return str(self.paths[i])

    def name_lower(self, i: int) -> str:
        return self.paths[i].name.lower()

    def suffix_lower(self, i: int) -> str:
        po = self.paths[i]
        suf = po.suffix.lower()
        if suf == ".gz" and len(po.suffixes) >= 2:
            suf = "".join(po.suffixes[-2:]).lower()
        return suf


class BasePathResolver(ABC):
    """路径解析策略基类：仅处理 assigned[i] 为 False 的项，匹配后必须置 assigned[i]=True。"""

    @abstractmethod
    def resolve(self, ctx: OmicsPathContext) -> None:
        """根据规则将未分配路径写入 ctx.out 并标记 ctx.assigned。"""
        raise NotImplementedError


class MetabolomicsPathResolver(BasePathResolver):
    """代谢组学 / 单表：.csv, .tsv, .xlsx, .txt, .tsv.gz"""

    def resolve(self, ctx: OmicsPathContext) -> None:
        for i in range(ctx.n):
            if ctx.assigned[i]:
                continue
            nl = ctx.name_lower(i)
            sl = ctx.suffix_lower(i)
            if sl in (".csv", ".tsv", ".xlsx", ".txt") or nl.endswith(".tsv.gz"):
                ctx.out["tables"].append(ctx.path_str(i))
                ctx.assigned[i] = True


class AnnDataPathResolver(BasePathResolver):
    """单细胞矩阵文件：.h5ad, .h5（与旧版 h5ad 桶一致；独立于代谢表）。"""

    def resolve(self, ctx: OmicsPathContext) -> None:
        for i in range(ctx.n):
            if ctx.assigned[i]:
                continue
            if ctx.suffix_lower(i) in (".h5ad", ".h5"):
                ctx.out["h5ad"].append(ctx.path_str(i))
                ctx.assigned[i] = True


class TenXPathResolver(BasePathResolver):
    """
    10x Genomics：收集文件名含 matrix.mtx / barcodes / features / genes 的候选路径。
    仅当同批凑齐 matrix + (barcodes|features|genes) 时：写入 10x_mtx 公共父目录，并标记这些路径已消费。
    若不凑齐：不标记、不写入 unknown，交由 Metabolomics 等后续策略（例如独立 barcodes.tsv → tables）。
    """

    def resolve(self, ctx: OmicsPathContext) -> None:
        cand_idx: List[int] = []
        for i in range(ctx.n):
            if ctx.assigned[i]:
                continue
            nl = ctx.name_lower(i)
            if (
                "matrix.mtx" in nl
                or "barcodes" in nl
                or "features" in nl
                or "genes" in nl
            ):
                cand_idx.append(i)

        if not cand_idx:
            return

        tenx_candidates = [ctx.path_str(i) for i in cand_idx]
        has_matrix = any("matrix.mtx" in Path(s).name.lower() for s in tenx_candidates)
        has_barcodes_or_features = any(
            "barcodes" in Path(s).name.lower()
            or "features" in Path(s).name.lower()
            or "genes" in Path(s).name.lower()
            for s in tenx_candidates
        )
        if not (has_matrix and has_barcodes_or_features):
            return

        for i in cand_idx:
            ctx.assigned[i] = True
        try:
            dirs = [os.path.dirname(s) for s in tenx_candidates]
            common = os.path.commonpath(dirs) if len(dirs) > 1 else dirs[0]
            ctx.out["10x_mtx"] = [common]
        except (ValueError, OSError):
            ctx.out["10x_mtx"] = (
                [os.path.dirname(tenx_candidates[0])] if tenx_candidates else []
            )


class RadiomicsPathResolver(BasePathResolver):
    """影像组学：.nii / .nii.gz / .dcm / tiff / nrrd 等；mask/seg 与 image 关键词分流。"""

    RADIOMICS_SUFFIXES = (".nii", ".nii.gz", ".dcm", ".tiff", ".tif", ".nrrd")
    MASK_KEYWORDS = ("mask", "roi", "seg", "label")
    IMAGE_KEYWORDS = ("image", "img", "mri", "ct")

    def resolve(self, ctx: OmicsPathContext) -> None:
        for i in range(ctx.n):
            if ctx.assigned[i]:
                continue
            p_str = ctx.path_str(i)
            nl = ctx.name_lower(i)
            sl = ctx.suffix_lower(i)
            is_radiomics = (
                sl in self.RADIOMICS_SUFFIXES
                or ".nii" in nl
                or nl.endswith(".nii.gz")
                or nl.endswith(".nrrd")
            )
            if not is_radiomics:
                continue
            if any(kw in nl for kw in self.MASK_KEYWORDS) and not any(
                k in nl for k in self.IMAGE_KEYWORDS
            ):
                ctx.out["masks"].append(p_str)
            else:
                ctx.out["images"].append(p_str)
            ctx.assigned[i] = True


def _resolve_path_objects(paths: List[str]) -> List[Path]:
    resolved_paths: List[Path] = []
    for p in paths:
        try:
            path_obj = Path(p)
            resolved_paths.append(path_obj.resolve() if path_obj.exists() else path_obj)
        except Exception:
            resolved_paths.append(Path(p))
    return resolved_paths


# 捆绑包优先：10x → 影像 → 代谢表 → AnnData；门面末尾写入 unknown
_DEFAULT_PIPELINE: List[BasePathResolver] = [
    TenXPathResolver(),
    RadiomicsPathResolver(),
    MetabolomicsPathResolver(),
    AnnDataPathResolver(),
]

def _fresh_buckets() -> Dict[str, List[str]]:
    return {
        "tables": [],
        "h5ad": [],
        "10x_mtx": [],
        "images": [],
        "masks": [],
        "unknown": [],
    }


class OmicsPathResolver:
    """默认策略链门面；可注入自定义 pipeline 供测试。"""

    def __init__(self, pipeline: Optional[List[BasePathResolver]] = None) -> None:
        self._pipeline = pipeline if pipeline is not None else list(_DEFAULT_PIPELINE)

    def resolve(self, raw_path_input: Union[str, List[str]]) -> Dict[str, List[str]]:
        return resolve_omics_paths_with_pipeline(raw_path_input, self._pipeline)


def resolve_omics_paths_with_pipeline(
    raw_path_input: Union[str, List[str]],
    pipeline: List[BasePathResolver],
) -> Dict[str, List[str]]:
    paths: List[str] = []
    if isinstance(raw_path_input, str):
        paths = [p.strip() for p in raw_path_input.replace(";", ",").split(",") if p.strip()]
    elif isinstance(raw_path_input, list):
        paths = [str(p).strip() for p in raw_path_input if str(p).strip()]
    else:
        return _fresh_buckets()

    if not paths:
        return _fresh_buckets()

    resolved = _resolve_path_objects(paths)
    ctx = OmicsPathContext(resolved)
    for resolver in pipeline:
        resolver.resolve(ctx)

    for i in range(ctx.n):
        if not ctx.assigned[i]:
            ctx.out["unknown"].append(ctx.path_str(i))

    return ctx.out


def resolve_omics_paths(raw_path_input: Union[str, List[str]]) -> Dict[str, List[str]]:
    """
    多文件解析总线（门面）：拆解入参后依次执行策略，合并为固定六键字典。

    Returns:
        {
            "tables": [],
            "h5ad": [],
            "10x_mtx": [],
            "images": [],
            "masks": [],
            "unknown": [],
        }
    """
    return resolve_omics_paths_with_pipeline(raw_path_input, list(_DEFAULT_PIPELINE))
