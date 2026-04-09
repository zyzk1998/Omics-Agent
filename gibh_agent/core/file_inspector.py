"""
文件检测器 - Universal Eyes of the System
多模态文件检查器，支持表格数据、单细胞数据、图像等

架构：策略模式（Strategy Pattern）
- InspectorRegistry: 单例注册表
- BaseFileHandler: 抽象基类
- 具体策略：TenXDirectoryHandler, AnnDataHandler, TabularHandler
- FileInspector: 门面类
"""
import os
import json
import gzip
import logging
import tarfile
import zipfile
from pathlib import Path
from typing import Dict, Optional, Any, List, Tuple, Set, Union
from abc import ABC, abstractmethod
import numpy as np

from .path_resolvers import resolve_omics_paths
from .asset_manager import (
    DataAsset,
    OmicsAssetManager,
    to_legacy_resolved_dict,
)

logger = logging.getLogger(__name__)


def parse_multiple_files(file_paths_string: str) -> Dict[str, List[str]]:
    """
    将前端逗号/分号拼接的路径字符串拆分为列表，并按后缀/文件名特征归类为结构化字典。
    供影像组学、单细胞 10x、空间组学等统一使用；无法识别的后缀归入 unknown，不丢弃。

    Returns:
        {
            "images": [...],   # 原图（nii/dcm/tiff 等，且文件名非 mask）
            "masks": [...],    # 掩膜（文件名含 mask/roi/seg 等）
            "matrix": [...],   # 矩阵/主数据（h5ad, mtx, 或 10x 目录路径）
            "barcodes": [...], # 10x barcodes.tsv
            "features": [...], # 10x features.tsv / genes.tsv
            "tables": [...],   # CSV 等表格
            "unknown": [...]  # 未识别后缀一律保留，防止误删用户数据
        }
    """
    if not file_paths_string or not isinstance(file_paths_string, str):
        return {
            "images": [], "masks": [], "matrix": [], "barcodes": [], "features": [],
            "tables": [], "unknown": [],
        }
    raw = file_paths_string.replace(";", ",").strip()
    parts = [x.strip() for x in raw.split(",") if x.strip()]
    paths: List[str] = []
    for p in parts:
        try:
            path_obj = Path(p)
            paths.append(str(path_obj.resolve()) if path_obj.exists() else p)
        except Exception:
            paths.append(p)
    paths = [p for p in paths if p]

    mask_keywords = ("mask", "roi", "seg", "segmentation", "label")
    image_extensions = (".nii", ".nii.gz", ".dcm", ".tiff", ".tif", ".png", ".jpg", ".jpeg", ".bmp", ".webp")
    out: Dict[str, List[str]] = {
        "images": [],
        "masks": [],
        "matrix": [],
        "barcodes": [],
        "features": [],
        "tables": [],
        "unknown": [],
    }

    for p in paths:
        path_obj = Path(p)
        name_lower = path_obj.name.lower()
        suffix_lower = path_obj.suffix.lower()
        if path_obj.suffix.lower() == ".gz" and len(path_obj.suffixes) >= 2:
            suffix_lower = "".join(path_obj.suffixes[-2:]).lower()

        # 10x
        if "barcodes" in name_lower and (name_lower.endswith(".tsv") or name_lower.endswith(".tsv.gz")):
            out["barcodes"].append(p)
            continue
        if ("features" in name_lower or "genes" in name_lower) and (name_lower.endswith(".tsv") or name_lower.endswith(".tsv.gz")):
            out["features"].append(p)
            continue
        if "matrix.mtx" in name_lower or name_lower == "matrix.mtx.gz":
            out["matrix"].append(p)
            continue

        # 单文件数据
        if suffix_lower in (".h5ad", ".h5"):
            out["matrix"].append(p)
            continue
        if suffix_lower == ".csv":
            out["tables"].append(p)
            continue

        # 影像：按扩展名 + 文件名区分 image / mask（.nii.gz 的 Path.suffix 为 .gz，用 name 判断）
        if (
            suffix_lower in image_extensions
            or ".nii" in name_lower
            or name_lower.endswith(".nii.gz")
            or suffix_lower == ".dcm"
        ):
            if any(kw in name_lower for kw in mask_keywords) and "image" not in name_lower and "img" not in name_lower and "mri" not in name_lower and "ct" not in name_lower:
                out["masks"].append(p)
            else:
                out["images"].append(p)
            continue

        # 未识别：保留到 unknown
        out["unknown"].append(p)

    return out


def path_looks_like_medical_imaging(path: Path) -> bool:
    """
    基于文件名识别影像组学常见后缀（不读取文件内容）。
    正确处理 .nii.gz（Path.suffix 仅为 .gz 的情况）。
    """
    if not path.is_file():
        return False
    name = path.name.lower()
    if name.endswith(".nii.gz"):
        return True
    if name.endswith(".nii"):
        return True
    if name.endswith(".dcm"):
        return True
    return False


def build_medical_imaging_inspection_result(path: Path) -> Dict[str, Any]:
    """
    影像文件体检：优先 nibabel 读取维度；ImportError 或其它失败时仅返回大小等元数据。
    绝不使用 pandas / scanpy。
    """
    resolved = path.resolve()
    file_size = resolved.stat().st_size if resolved.is_file() else 0
    shape_list: Optional[List[int]] = None
    try:
        import nibabel as nib  # type: ignore
    except ImportError:
        nib = None  # type: ignore
    if nib is not None:
        try:
            img = nib.load(str(resolved))
            shape_list = [int(x) for x in img.shape]
        except Exception as e:
            logger.debug("nibabel 读取影像失败（仍将返回基础元数据）: %s", e)
            shape_list = None

    abs_image = str(resolved)
    mask_path_resolved: Optional[str] = None
    try:
        from .file_handlers.structure_normalizer import pair_radiomics_files
        if pair_radiomics_files is not None:
            img_cand, mask_cand = pair_radiomics_files(resolved.parent)
            if img_cand is not None:
                abs_image = str(Path(img_cand).resolve())
            if mask_cand is not None:
                mask_path_resolved = str(Path(mask_cand).resolve())
    except Exception as e:
        logger.debug("pair_radiomics_files 跳过: %s", e)

    # Level-3：同目录配对后做体维度探针；不一致则丢弃 mask（降级为单图），避免下游致命假设
    if mask_path_resolved and abs_image:
        from .asset_manager import _probe_radiomics_volume_alignment

        _ok, _reason = _probe_radiomics_volume_alignment(
            Path(abs_image), Path(mask_path_resolved)
        )
        if not _ok:
            logger.warning(
                "⚠️ [FileInspector] radiomics mask dropped after probe: %s",
                _reason,
            )
            mask_path_resolved = None

    if shape_list is not None:
        message = f"影像文件校验通过，空间维度: {shape_list}"
    else:
        message = (
            "影像文件校验通过（未安装 nibabel 或无法解析体数据维度，已返回文件名与大小）"
        )

    shape_info: Dict[str, Any] = {
        "spatial_shape": shape_list,
        "file_size_bytes": int(file_size),
    }

    return {
        "status": "success",
        "success": True,
        "file_path": abs_image,
        "mask_path": mask_path_resolved,
        "file_type": "medical_imaging",
        "domain": "Radiomics",
        "modality": "Radiomics",
        "message": message,
        "shape": shape_info,
        "head": {
            "markdown": f"医学影像\n- 文件: {resolved.name}\n- {message}",
            "json": {"file_type": "medical_imaging", "shape": shape_info},
        },
    }


# ============================================
# Part 1: Registry & Base Class
# ============================================

class InspectorRegistry:
    """
    检查器注册表（单例模式）
    管理所有文件检查策略，按优先级排序
    """
    _instance = None
    _handlers: List[Tuple[int, type]] = []  # [(priority, handler_class), ...]
    
    def __new__(cls):
        if cls._instance is None:
            cls._instance = super().__new__(cls)
        return cls._instance
    
    def register(self, priority: int, handler_class: type):
        """注册检查器"""
        self._handlers.append((priority, handler_class))
        # 按优先级降序排序（优先级高的先执行）
        self._handlers.sort(key=lambda x: x[0], reverse=True)
        logger.debug(f"📝 Registered handler: {handler_class.__name__} (priority={priority})")
    
    def get_handlers(self) -> List[type]:
        """获取所有已注册的检查器（按优先级排序）"""
        return [handler_class for _, handler_class in self._handlers]
    
    def clear(self):
        """清空注册表（用于测试）"""
        self._handlers.clear()


# 全局注册表实例
_registry = InspectorRegistry()


def register_inspector(priority: int):
    """
    装饰器：注册文件检查器
    
    Args:
        priority: 优先级（数字越大，优先级越高）
    """
    def decorator(cls):
        _registry.register(priority, cls)
        return cls
    return decorator


class BaseFileHandler(ABC):
    """
    文件检查器基类
    
    所有具体的文件检查器都应该继承此类并实现：
    - can_handle(path) -> bool: 检查是否可以处理该路径
    - inspect(path) -> dict: 执行检查并返回结果
    """
    
    @abstractmethod
    def can_handle(self, path: Path) -> bool:
        """
        检查是否可以处理该路径
        
        Args:
            path: 文件或目录路径
            
        Returns:
            True 如果可以处理，False 否则
            
        Note:
            MUST check file content/structure, not just extension
        """
        pass
    
    @abstractmethod
    def inspect(self, path: Path) -> Dict[str, Any]:
        """
        执行文件检查
        
        Args:
            path: 文件或目录路径
            
        Returns:
            包含检查结果的字典，必须包含：
            - status: "success" 或 "error"
            - shape: {"rows": int, "cols": int} 或 {"n_obs": int, "n_vars": int}
            - columns: List[str] 或 None
            - file_type: str
            - preview: Optional[Dict] 或 str
            - error: Optional[str] (如果 status == "error")
        """
        pass


# ============================================
# Part 2: Concrete Strategies
# ============================================

@register_inspector(priority=10)
class TenXDirectoryHandler(BaseFileHandler):
    """
    10x Genomics 目录检查器（优先级：10）
    
    检查包含 matrix.mtx, barcodes.tsv, features.tsv 的目录
    支持递归搜索嵌套目录
    """
    
    def _find_10x_data_path(self, root_path: Path) -> Optional[Path]:
        """
        递归搜索 10x 数据目录
        
        🔥 CRITICAL: 必须找到包含 matrix.mtx 的实际数据目录
        支持嵌套目录结构（如 MyData/Sample1/outs/）
        
        Args:
            root_path: 根目录路径
            
        Returns:
            包含 matrix.mtx 的子目录路径，如果未找到则返回 None
        """
        print(f"DEBUG: Searching for 10x data in: {root_path}")  # 🔥 Loud logging
        logger.debug(f"DEBUG: Inspecting 10x path: {root_path}")
        
        # Step 1: 首先检查根目录
        try:
            dir_contents = os.listdir(root_path)
            has_matrix = any(f in dir_contents for f in ['matrix.mtx', 'matrix.mtx.gz'])
            has_barcodes = any(f in dir_contents for f in ['barcodes.tsv', 'barcodes.tsv.gz'])
            has_features = any(f in dir_contents for f in ['features.tsv', 'features.tsv.gz', 'genes.tsv', 'genes.tsv.gz'])
            
            if has_matrix and (has_barcodes or has_features):
                print(f"DEBUG: Found 10x data at root: {root_path}")  # 🔥 Loud logging
                logger.debug(f"DEBUG: Found 10x data at root: {root_path}")
                return root_path
        except Exception as e:
            logger.warning(f"⚠️ Failed to list root directory: {e}")
        
        # Step 2: 递归搜索子目录（使用 os.walk）
        print(f"DEBUG: Recursively searching subdirectories...")  # 🔥 Loud logging
        for root, dirs, files in os.walk(root_path):
            root_path_obj = Path(root)
            
            # 🔥 CRITICAL: 检查当前目录是否包含 matrix.mtx（未压缩或压缩）
            has_matrix_local = 'matrix.mtx' in files or 'matrix.mtx.gz' in files
            
            if has_matrix_local:
                # 🔥 健壮的文件名匹配：检查 barcodes 和 features（支持多种命名）
                has_barcodes_local = any(f in files for f in ['barcodes.tsv', 'barcodes.tsv.gz'])
                has_features_local = any(f in files for f in ['features.tsv', 'features.tsv.gz', 'genes.tsv', 'genes.tsv.gz'])
                
                # 🔥 如果找到 matrix.mtx，即使没有 barcodes/features 也返回（可能在其他位置）
                # 但优先返回同时包含所有文件的目录
                if has_barcodes_local or has_features_local:
                    print(f"DEBUG: Found data root at: {root_path_obj}")  # 🔥 Loud logging
                    logger.debug(f"DEBUG: Found matrix at: {root_path_obj}")
                    return root_path_obj
                else:
                    # 如果只有 matrix.mtx，也返回（barcodes/features 可能在父目录或子目录）
                    print(f"DEBUG: Found matrix.mtx at: {root_path_obj} (checking for barcodes/features nearby)")  # 🔥 Loud logging
                    # 检查父目录和当前目录
                    parent_dir = root_path_obj.parent
                    if parent_dir.exists():
                        try:
                            parent_files = os.listdir(parent_dir)
                            if any(f in parent_files for f in ['barcodes.tsv', 'barcodes.tsv.gz', 'features.tsv', 'features.tsv.gz', 'genes.tsv', 'genes.tsv.gz']):
                                print(f"DEBUG: Found barcodes/features in parent: {parent_dir}")  # 🔥 Loud logging
                                return root_path_obj
                        except Exception:
                            pass
                    
                    # 如果当前目录有 matrix.mtx，返回它（让后续逻辑处理文件查找）
                    print(f"DEBUG: Using directory with matrix.mtx: {root_path_obj}")  # 🔥 Loud logging
                    return root_path_obj
        
        print(f"DEBUG: No 10x data found in {root_path}")  # 🔥 Loud logging
        logger.debug(f"DEBUG: No 10x data found in {root_path}")
        return None
    
    def can_handle(self, path: Path) -> bool:
        """检查是否为 10x Genomics 目录（支持嵌套目录）"""
        if not path.is_dir():
            return False
        
        try:
            real_data_path = self._find_10x_data_path(path)
            return real_data_path is not None
        except Exception as e:
            logger.debug(f"DEBUG: can_handle failed: {e}")
            return False
    
    def _smart_open(self, filepath: Path):
        """
        智能打开文件，自动处理 gzip 压缩
        
        Args:
            filepath: 文件路径（Path 对象）
            
        Returns:
            文件对象（上下文管理器）
        """
        filepath_str = str(filepath)
        if filepath_str.endswith('.gz'):
            # 🔥 CRITICAL: 'rt' 模式 + encoding='utf-8' 对于 gzip 文件至关重要
            return gzip.open(filepath_str, 'rt', encoding='utf-8')
        else:
            return open(filepath_str, 'r', encoding='utf-8')
    
    def _count_lines_skip_comments(self, file_path: Path) -> int:
        """
        统计文件行数，跳过注释行（以 % 或 # 开头）
        
        支持 gzip 压缩文件（透明处理）
        """
        count = 0
        try:
            with self._smart_open(file_path) as f:
                for line in f:
                    line = line.strip()
                    if line and not line.startswith('%') and not line.startswith('#'):
                        count += 1
            logger.debug(f"DEBUG: Counted {count} lines in {file_path}")
        except Exception as e:
            logger.error(f"❌ Failed to count lines in {file_path}: {e}", exc_info=True)
            # 返回 0 而不是抛出异常，让调用者决定如何处理
        return count
    
    def inspect(self, path: Path) -> Dict[str, Any]:
        """
        检查 10x Genomics 目录（支持嵌套目录）
        
        返回包含 debug_trace 字段的字典，用于调试
        """
        debug_logs = []  # 🔥 创建调试日志列表
        debug_logs.append(f"Input path: {path}")
        
        try:
            logger.debug(f"DEBUG: Inspecting 10x path: {path}")
            debug_logs.append(f"Starting inspection of: {path}")
            
            # 记录根目录内容
            try:
                dir_contents = os.listdir(path)
                debug_logs.append(f"Directory contents: {dir_contents}")
                logger.debug(f"DEBUG: Directory contents: {dir_contents}")
            except Exception as e:
                debug_logs.append(f"Failed to list directory: {str(e)}")
                logger.warning(f"⚠️ Failed to list directory: {e}")
            
            # 🔥 递归搜索找到实际的数据目录
            debug_logs.append(f"Starting recursive search for matrix.mtx...")
            real_data_path = self._find_10x_data_path(path)
            
            if real_data_path is None:
                debug_logs.append(f"ERROR: Could not find matrix.mtx in directory tree")
                debug_trace = "\n".join(debug_logs)
                return {
                    "status": "error",
                    "error": f"在目录 {path} 中未找到 10x 数据文件（matrix.mtx, barcodes.tsv, features.tsv）",
                    "file_type": "directory",
                    "debug_trace": debug_trace  # 🔥 返回调试跟踪
                }
            
            debug_logs.append(f"Found data root at: {real_data_path}")
            logger.debug(f"DEBUG: Found 10x root at: {real_data_path}")
            print(f"DEBUG: Found data root at: {real_data_path}")  # 🔥 Loud logging
            
            logger.debug(f"DEBUG: Found 10x root at: {real_data_path}")
            print(f"DEBUG: Found data root at: {real_data_path}")  # 🔥 Loud logging
            
            # Step B: 在找到的目录中查找必需文件（支持压缩和未压缩）
            # 🔥 健壮的文件名匹配：支持多种命名约定
            matrix_path = None
            barcodes_path = None
            features_path = None
            
            # 查找 matrix.mtx（压缩或未压缩）
            debug_logs.append(f"Looking for matrix.mtx in: {real_data_path}")
            print(f"DEBUG: Looking for matrix.mtx in: {real_data_path}")  # 🔥 Loud logging
            try:
                data_dir_contents = os.listdir(real_data_path)
                debug_logs.append(f"Data directory contents: {data_dir_contents}")
            except Exception as e:
                debug_logs.append(f"Failed to list data directory: {str(e)}")
            
            for f in ['matrix.mtx', 'matrix.mtx.gz']:
                candidate = real_data_path / f
                if candidate.exists():
                    matrix_path = candidate
                    debug_logs.append(f"Matrix found at: {matrix_path}")
                    print(f"DEBUG: Found matrix at: {matrix_path}")  # 🔥 Loud logging
                    break
            
            if not matrix_path:
                debug_logs.append(f"ERROR: matrix.mtx not found in {real_data_path}")
            
            # 查找 barcodes.tsv（压缩或未压缩）
            # 🔥 如果不在当前目录，也检查父目录
            debug_logs.append(f"Looking for barcodes.tsv in: {real_data_path}")
            print(f"DEBUG: Looking for barcodes.tsv in: {real_data_path}")  # 🔥 Loud logging
            search_dirs = [real_data_path]
            if real_data_path.parent.exists():
                search_dirs.append(real_data_path.parent)
                debug_logs.append(f"Also checking parent directory: {real_data_path.parent}")
            
            for search_dir in search_dirs:
                debug_logs.append(f"Scanning subdir: {search_dir}")
                for f in ['barcodes.tsv', 'barcodes.tsv.gz']:
                    candidate = search_dir / f
                    if candidate.exists():
                        barcodes_path = candidate
                        debug_logs.append(f"Barcodes found at: {barcodes_path}")
                        print(f"DEBUG: Found barcodes at: {barcodes_path}")  # 🔥 Loud logging
                        break
                if barcodes_path:
                    break
            
            if not barcodes_path:
                debug_logs.append(f"ERROR: barcodes.tsv not found in search directories")
            
            # 查找 features.tsv 或 genes.tsv（压缩或未压缩）
            # 🔥 健壮匹配：支持 features.tsv, genes.tsv（标准 10x 格式）
            debug_logs.append(f"Looking for features.tsv/genes.tsv in: {real_data_path}")
            print(f"DEBUG: Looking for features.tsv/genes.tsv in: {real_data_path}")  # 🔥 Loud logging
            for search_dir in search_dirs:
                debug_logs.append(f"Scanning subdir for features: {search_dir}")
                for f in ['features.tsv', 'features.tsv.gz', 'genes.tsv', 'genes.tsv.gz']:
                    candidate = search_dir / f
                    if candidate.exists():
                        features_path = candidate
                        debug_logs.append(f"Features found at: {features_path}")
                        print(f"DEBUG: Found features at: {features_path}")  # 🔥 Loud logging
                        break
                if features_path:
                    break
            
            if not features_path:
                debug_logs.append(f"ERROR: features.tsv/genes.tsv not found in search directories")
            
            if not all([matrix_path, barcodes_path, features_path]):
                missing = []
                if not matrix_path:
                    missing.append("matrix.mtx")
                if not barcodes_path:
                    missing.append("barcodes.tsv")
                if not features_path:
                    missing.append("features.tsv/genes.tsv")
                
                debug_logs.append(f"ERROR: Missing required files: {', '.join(missing)}")
                debug_trace = "\n".join(debug_logs)
                return {
                    "status": "error",
                    "error": f"10x目录缺少必需文件: {', '.join(missing)}。搜索路径: {real_data_path}",
                    "file_type": "directory",
                    "debug_trace": debug_trace  # 🔥 返回调试跟踪
                }
            
            # Step C: 统计细胞数（barcodes.tsv 的行数，跳过注释）
            debug_logs.append(f"Attempting to read barcodes file: {barcodes_path}")
            logger.debug(f"DEBUG: Counting cells from: {barcodes_path}")
            print(f"DEBUG: Counting cells from: {barcodes_path}")  # 🔥 Loud logging
            try:
                n_cells = self._count_lines_skip_comments(barcodes_path)
                debug_logs.append(f"Successfully counted {n_cells} cells from barcodes.tsv")
            except Exception as e:
                debug_logs.append(f"Error reading barcodes file: {str(e)}")
                logger.error(f"❌ Error reading barcodes: {e}", exc_info=True)
                n_cells = 0
            
            # 🔥 CRITICAL FIX: 如果行数统计失败（返回0），回退到 scanpy.read_10x_mtx
            if n_cells == 0:
                debug_logs.append(f"WARNING: Line counting returned 0 cells, falling back to scanpy.read_10x_mtx")
                logger.warning(f"⚠️ Line counting returned 0, falling back to scanpy.read_10x_mtx")
                try:
                    import scanpy as sc
                    from .rna_utils import read_10x_data
                    debug_logs.append(f"Attempting scanpy.read_10x_mtx fallback...")
                    logger.info(f"🔄 Using scanpy.read_10x_mtx as fallback...")
                    adata = read_10x_data(str(real_data_path), var_names='gene_symbols', cache=False)
                    n_cells = adata.n_obs
                    n_genes = adata.n_vars
                    debug_logs.append(f"Fallback successful: {n_cells} cells, {n_genes} genes")
                    logger.info(f"✅ Fallback successful: {n_cells} cells, {n_genes} genes")
                    print(f"DEBUG: Detected {n_cells} cells (via scanpy fallback)")  # 🔥 Loud logging
                    
                    # 使用 scanpy 的结果，跳过手动统计
                    feature_types = ['Gene Expression']  # 默认值
                    columns = list(adata.var_names[:20]) if hasattr(adata, 'var_names') else []
                    
                    if n_cells == 0 or n_genes == 0:
                        debug_logs.append(f"ERROR: Fallback also returned 0 cells or genes")
                        debug_trace = "\n".join(debug_logs)
                        return {
                            "status": "error",
                            "error": f"10x数据为空: {n_cells} 个细胞, {n_genes} 个基因",
                            "file_type": "directory",
                            "n_obs": n_cells,
                            "n_vars": n_genes,
                            "n_samples": n_cells,
                            "n_features": n_genes,
                            "debug_trace": debug_trace  # 🔥 返回调试跟踪
                        }
                    
                    # 构建预览
                    preview = {
                        "n_cells": n_cells,
                        "n_genes": n_genes,
                        "feature_types": feature_types,
                        "sample_genes": columns[:10]
                    }
                    
                    debug_trace = "\n".join(debug_logs)  # 🔥 生成调试跟踪
                    
                    # 直接返回结果（跳过下面的手动统计逻辑）
                    result = {
                        "status": "success",
                        "file_path": str(path.resolve()),
                        "real_data_path": str(real_data_path.resolve()),
                        "file_type": "anndata",
                        "shape": {
                            "rows": n_cells,
                            "cols": n_genes
                        },
                        "n_obs": n_cells,
                        "n_vars": n_genes,
                        "n_samples": n_cells,  # 🔥 关键：确保存在
                        "n_features": n_genes,  # 🔥 关键：确保存在
                        "columns": columns[:20] if columns else None,
                        "head": {
                            "markdown": f"10x Genomics 数据\n- 细胞数: {n_cells}\n- 基因数: {n_genes}\n- 数据路径: {real_data_path}",
                            "json": preview
                        },
                        "feature_types": feature_types,
                        "debug_trace": debug_trace,  # 🔥 返回调试跟踪
                        "data": {
                            "summary": {
                                "n_samples": n_cells,
                                "n_features": n_genes,
                                "feature_types": feature_types
                            }
                        }
                    }
                    
                    logger.info(f"✅ [TenXDirectoryHandler] Success (fallback): {n_cells} cells, {n_genes} genes")
                    return result
                    
                except Exception as fallback_error:
                    debug_logs.append(f"EXCEPTION in fallback: {str(fallback_error)}")
                    debug_trace = "\n".join(debug_logs)
                    logger.error(f"❌ Fallback to scanpy also failed: {fallback_error}", exc_info=True)
                    return {
                        "status": "error",
                        "error": f"无法读取10x数据: {str(fallback_error)}",
                        "file_type": "directory",
                        "debug_trace": debug_trace  # 🔥 返回调试跟踪
                    }
            
            print(f"DEBUG: Detected {n_cells} cells.")  # 🔥 Loud logging
            
            # 读取 features.tsv 获取基因信息（使用 _smart_open）
            debug_logs.append(f"Attempting to read features file: {features_path}")
            n_genes = 0
            feature_types = []
            columns = []
            
            try:
                logger.debug(f"DEBUG: Reading features from: {features_path}")
                with self._smart_open(features_path) as f:
                    for line in f:
                        line = line.strip()
                        if line and not line.startswith('%') and not line.startswith('#'):
                            n_genes += 1
                            parts = line.split('\t')
                            if len(parts) >= 1:
                                gene_id = parts[0]
                                gene_symbol = parts[1] if len(parts) > 1 else gene_id
                                feature_type = parts[2] if len(parts) > 2 else 'Gene Expression'
                                
                                if feature_type not in feature_types:
                                    feature_types.append(feature_type)
                                
                                # 只保存前100个基因名（避免返回过多数据）
                                if len(columns) < 100:
                                    columns.append(gene_symbol or gene_id)
                debug_logs.append(f"Successfully read {n_genes} genes from features.tsv")
            except Exception as e:
                debug_logs.append(f"Error reading features file: {str(e)}")
                logger.warning(f"⚠️ Failed to read features.tsv: {e}")
                # 如果读取失败，至少统计行数
                try:
                    n_genes = self._count_lines_skip_comments(features_path)
                    debug_logs.append(f"Fallback: counted {n_genes} genes using line counting")
                except Exception as count_err:
                    debug_logs.append(f"Error counting lines in features file: {str(count_err)}")
                    n_genes = 0
            
            print(f"DEBUG: Detected {n_genes} genes.")  # 🔥 Loud logging
            debug_logs.append(f"Final counts: {n_cells} cells, {n_genes} genes")
            
            if n_cells == 0 or n_genes == 0:
                return {
                    "status": "error",
                    "error": f"10x数据为空: {n_cells} 个细胞, {n_genes} 个基因",
                    "file_type": "directory",
                    "n_obs": n_cells,
                    "n_vars": n_genes
                }
            
            # 构建预览（前5个细胞和基因）
            preview = {
                "n_cells": n_cells,
                "n_genes": n_genes,
                "feature_types": feature_types,
                "sample_genes": columns[:10]
            }
            
            # 🔥 CRITICAL FIX: 确保返回所有必需的键（兼容 DataDiagnostician）
            # 同时设置 n_obs/n_vars 和 n_samples/n_features
            debug_trace = "\n".join(debug_logs)  # 🔥 生成调试跟踪字符串
            
            result = {
                "status": "success",
                "file_path": str(path.resolve()),  # 返回原始路径（用户上传的路径）
                "real_data_path": str(real_data_path.resolve()),  # 实际数据所在路径（用于调试）
                "file_type": "anndata",
                "shape": {
                    "rows": n_cells,
                    "cols": n_genes
                },
                # scRNA-seq 格式（ScanpyTool 期望的）
                "n_obs": n_cells,  # 🔥 关键：确保存在
                "n_vars": n_genes,  # 🔥 关键：确保存在
                # 通用格式（DataDiagnostician 期望的）
                "n_samples": n_cells,  # 🔥 关键：确保存在
                "n_features": n_genes,  # 🔥 关键：确保存在
                "columns": columns[:20] if columns else None,  # 前20个基因名
                "head": {
                    "markdown": f"10x Genomics 数据\n- 细胞数: {n_cells}\n- 基因数: {n_genes}\n- 特征类型: {', '.join(feature_types)}\n- 数据路径: {real_data_path}",
                    "json": preview
                },
                "feature_types": feature_types,
                "debug_trace": debug_trace,  # 🔥 返回调试跟踪
                "data": {
                    "summary": {
                        "n_samples": n_cells,  # 🔥 确保嵌套字典中也存在
                        "n_features": n_genes,  # 🔥 确保嵌套字典中也存在
                        "feature_types": feature_types
                    }
                }
            }
            
            logger.info(f"✅ [TenXDirectoryHandler] Success: {n_cells} cells, {n_genes} genes (data_path: {real_data_path})")
            print(f"DEBUG: Final result - {n_cells} cells, {n_genes} genes")  # 🔥 Loud logging
            return result
            
        except Exception as e:
            debug_logs.append(f"EXCEPTION: {str(e)}")
            debug_trace = "\n".join(debug_logs)
            logger.error(f"❌ [TenXDirectoryHandler] Failed: {e}", exc_info=True)
            return {
                "status": "error",
                "error": f"无法读取10x数据: {str(e)}",
                "file_type": "directory",
                "debug_trace": debug_trace  # 🔥 返回调试跟踪
            }


def _archive_member_basenames(archive_path: Path) -> Optional[Set[str]]:
    """
    Peek into archive and return the set of member basenames (lowercase).
    Returns None on open/list error.
    """
    names: Set[str] = set()
    name_lower = archive_path.name.lower()
    try:
        if name_lower.endswith(".zip"):
            with zipfile.ZipFile(archive_path, "r") as z:
                for n in z.namelist():
                    base = n.replace("\\", "/").split("/")[-1].lower()
                    if base:
                        names.add(base)
            return names
        if name_lower.endswith((".tar.gz", ".tgz", ".tar")):
            with tarfile.open(archive_path, "r:*") as t:
                for m in t.getmembers():
                    base = (m.name or "").replace("\\", "/").split("/")[-1].lower()
                    if base:
                        names.add(base)
            return names
    except Exception as e:
        logger.debug("Archive peek failed %s: %s", archive_path, e)
    return None


def _archive_contains_10x_scrna(archive_path: Path) -> bool:
    """
    Content sniffing: return True only if the archive contains standard 10x scRNA-seq
    files (matrix.mtx[/.gz], barcodes.tsv[/.gz], features.tsv or genes.tsv[/.gz]).
    Returns False for spatial-only archives (images, scalefactors, no matrix.mtx).
    """
    basenames = _archive_member_basenames(archive_path)
    if not basenames:
        return False
    has_matrix = any(
        b in basenames for b in ("matrix.mtx", "matrix.mtx.gz")
    )
    has_barcodes = any(
        b in basenames for b in ("barcodes.tsv", "barcodes.tsv.gz")
    )
    has_features = any(
        b in basenames for b in (
            "features.tsv", "features.tsv.gz",
            "genes.tsv", "genes.tsv.gz",
        )
    )
    return bool(has_matrix and (has_barcodes or has_features))


@register_inspector(priority=9)
class Archive10xHandler(BaseFileHandler):
    """
    10x / Cell Ranger 压缩包检查器（优先级：9）
    
    委托 core.rna_utils.load_10x_from_tarball 解压并加载，与执行阶段共用同一套逻辑，
    保证 Inspector 能读则 Executor 一定能读。persist_h5ad=True 以便执行阶段直接使用 .h5ad，避免重复解压。
    can_handle 使用内容嗅探（peek）：仅当压缩包内含有 matrix.mtx + barcodes/features 时返回 True，
    不含 MTX 的 spatial 包由 StructureNormalizer 解压后由 SpatialVisiumHandler 识别。
    """
    
    def can_handle(self, path: Path) -> bool:
        """Only handle archives that contain standard 10x scRNA-seq files (content sniffing)."""
        if not path.is_file():
            return False
        name = path.name.lower()
        if not (name.endswith(".tar.gz") or name.endswith(".tgz") or name.endswith(".tar") or name.endswith(".zip")):
            return False
        return _archive_contains_10x_scrna(path)
    
    def inspect(self, path: Path) -> Dict[str, Any]:
        """
        使用共享的 load_10x_from_tarball 解压并加载 10x 数据，持久化为 .h5ad 供执行阶段使用。
        返回与 TenXDirectoryHandler 兼容的元数据结构。
        """
        try:
            from .rna_utils import load_10x_from_tarball
            adata, output_base_dir = load_10x_from_tarball(
                str(path.resolve()),
                var_names="gene_symbols",
                persist_h5ad=True,
            )
            n_cells = adata.n_obs
            n_genes = adata.n_vars
            name = path.name.lower()
            if name.endswith(".tar.gz"):
                stem = path.name[:-7]
            elif name.endswith(".zip"):
                stem = path.name[:-4]
            else:
                stem = path.stem
            persisted_h5ad = Path(output_base_dir) / f"{stem}.h5ad"
            file_path = str(persisted_h5ad.resolve())
            columns = list(adata.var_names[:20]) if hasattr(adata, "var_names") else []
            try:
                feature_types = list(adata.var["feature_types"].astype(str).unique()) if "feature_types" in adata.var.columns else ["Gene Expression"]
            except Exception:
                feature_types = ["Gene Expression"]
            preview = {
                "n_cells": n_cells,
                "n_genes": n_genes,
                "feature_types": feature_types,
                "sample_genes": columns[:10],
            }
            result = {
                "status": "success",
                "success": True,
                "file_path": file_path,
                "file_type": "anndata",
                "shape": {"rows": n_cells, "cols": n_genes},
                "n_obs": n_cells,
                "n_vars": n_genes,
                "n_samples": n_cells,
                "n_features": n_genes,
                "columns": columns[:20] if columns else None,
                "head": {
                    "markdown": f"10x Genomics 数据（来自压缩包）\n- 细胞数: {n_cells}\n- 基因数: {n_genes}\n- 数据路径: {file_path}",
                    "json": preview,
                },
                "feature_types": feature_types,
                "data": {
                    "summary": {
                        "n_samples": n_cells,
                        "n_features": n_genes,
                        "feature_types": feature_types,
                    }
                },
                "unpacked_from": str(path.resolve()),
            }
            logger.info("✅ [Archive10xHandler] 通过 load_10x_from_tarball 完成体检，持久化: %s", file_path)
            return result
        except Exception as e:
            logger.exception("❌ [Archive10xHandler] 解压或识别失败: %s", path)
            return {
                "status": "error",
                "success": False,
                "error": f"压缩包处理失败: {str(e)}",
                "file_type": "unknown",
                "file_path": str(path.resolve()),
            }


@register_inspector(priority=8)
class FastqDirectoryHandler(BaseFileHandler):
    """
    FASTQ 目录检查器（优先级：8，高于其他目录检查器）
    
    检查包含 FASTQ 文件的目录（.fastq, .fq, .fastq.gz, .fq.gz）
    提供轻量级诊断：文件数量、总大小、文件类型等
    """
    
    def can_handle(self, path: Path) -> bool:
        """检查是否为 FASTQ 目录"""
        if not path.is_dir():
            return False
        
        try:
            # 检查目录中是否包含 FASTQ 文件
            fastq_extensions = ('.fastq', '.fq', '.fastq.gz', '.fq.gz')
            for item in path.iterdir():
                if item.is_file() and item.name.lower().endswith(fastq_extensions):
                    return True
                # 也检查子目录（递归一层）
                elif item.is_dir():
                    try:
                        for subitem in item.iterdir():
                            if subitem.is_file() and subitem.name.lower().endswith(fastq_extensions):
                                return True
                    except (PermissionError, OSError):
                        continue
            return False
        except (PermissionError, OSError) as e:
            logger.warning(f"⚠️ 无法访问目录 {path}: {e}")
            return False
    
    def inspect(self, path: Path) -> Dict[str, Any]:
        """
        检查 FASTQ 目录
        
        返回轻量级诊断信息：
        - 文件数量
        - 总大小
        - 文件类型（R1, R2, I1等）
        - 是否包含完整的配对端数据
        """
        try:
            fastq_files = []
            total_size = 0
            file_types = set()  # R1, R2, I1等
            
            fastq_extensions = ('.fastq', '.fq', '.fastq.gz', '.fq.gz')
            
            # 收集所有FASTQ文件
            for item in path.rglob('*'):  # 递归搜索
                if item.is_file() and item.name.lower().endswith(fastq_extensions):
                    fastq_files.append(item)
                    total_size += item.stat().st_size
                    
                    # 识别文件类型（R1, R2, I1等）
                    name_lower = item.name.lower()
                    if '_r1_' in name_lower or '_r1.' in name_lower or name_lower.endswith('_r1.fastq') or name_lower.endswith('_r1.fq'):
                        file_types.add('R1')
                    elif '_r2_' in name_lower or '_r2.' in name_lower or name_lower.endswith('_r2.fastq') or name_lower.endswith('_r2.fq'):
                        file_types.add('R2')
                    elif '_i1_' in name_lower or '_i1.' in name_lower or name_lower.endswith('_i1.fastq') or name_lower.endswith('_i1.fq'):
                        file_types.add('I1')
                    elif '_i2_' in name_lower or '_i2.' in name_lower:
                        file_types.add('I2')
            
            if not fastq_files:
                return {
                    "status": "error",
                    "error": f"目录 {path} 中未找到 FASTQ 文件",
                    "file_type": "directory"
                }
            
            # 检查是否包含配对端数据
            has_paired_end = 'R1' in file_types and 'R2' in file_types
            has_index = 'I1' in file_types or 'I2' in file_types
            
            # 构建诊断信息
            diagnosis_info = {
                "file_count": len(fastq_files),
                "total_size_bytes": total_size,
                "total_size_mb": round(total_size / (1024 * 1024), 2),
                "total_size_gb": round(total_size / (1024 * 1024 * 1024), 2),
                "file_types": sorted(list(file_types)),
                "has_paired_end": has_paired_end,
                "has_index": has_index,
                "is_10x_format": has_paired_end and has_index,  # 10x格式通常包含R1, R2, I1
                "sample_names": sorted(list(set(f.name.split('_')[0] for f in fastq_files)))[:10]  # 前10个样本名
            }
            
            return {
                "status": "success",
                "success": True,
                "file_type": "fastq",
                "file_path": str(path.resolve()),
                "shape": {"n_files": len(fastq_files)},  # 兼容性
                "columns": None,
                "preview": None,
                "diagnosis_info": diagnosis_info,
                "message": f"检测到 {len(fastq_files)} 个 FASTQ 文件，总大小 {diagnosis_info['total_size_gb']:.2f} GB"
            }
            
        except Exception as e:
            logger.error(f"❌ [FastqDirectoryHandler] 检查失败: {e}", exc_info=True)
            return {
                "status": "error",
                "error": f"无法检查 FASTQ 目录: {str(e)}",
                "file_type": "directory"
            }


@register_inspector(priority=5)
class AnnDataHandler(BaseFileHandler):
    """
    AnnData (H5AD) 文件检查器（优先级：5）
    
    检查 .h5ad 文件，使用 backed='r' 模式避免加载全部数据
    """
    
    def can_handle(self, path: Path) -> bool:
        """检查是否为 H5AD 文件"""
        if not path.is_file():
            return False
        
        if not path.suffix.lower() == '.h5ad':
            return False
        
        # 内容检查：验证是否为有效的 HDF5 文件
        try:
            import h5py
            return h5py.is_hdf5(str(path))
        except ImportError:
            # 如果没有 h5py，至少检查扩展名
            logger.debug("⚠️ h5py not available, using extension check only")
            return True
        except Exception:
            return False
    
    def inspect(self, path: Path) -> Dict[str, Any]:
        """检查 H5AD 文件"""
        try:
            import scanpy as sc
            
            logger.info(f"🔍 [AnnDataHandler] Loading H5AD: {path}")
            
            # 使用 backed='r' 模式（只读，不加载全部数据）
            try:
                adata = sc.read_h5ad(str(path), backed='r')
            except Exception as e:
                logger.warning(f"⚠️ Backed mode failed: {e}, using normal mode")
                adata = sc.read_h5ad(str(path))
            
            if adata is None:
                return {
                    "status": "error",
                    "error": "数据读取失败：返回了空对象",
                    "file_type": "anndata"
                }
            
            if adata.n_obs == 0 or adata.n_vars == 0:
                return {
                    "status": "error",
                    "error": f"数据为空: {adata.n_obs} 个细胞, {adata.n_vars} 个基因",
                    "file_type": "anndata",
                    "n_obs": adata.n_obs,
                    "n_vars": adata.n_vars
                }
            
            # 提取基本信息
            n_obs = adata.n_obs
            n_vars = adata.n_vars
            obs_keys = list(adata.obs.columns) if hasattr(adata.obs, 'columns') else []
            var_keys = list(adata.var.columns) if hasattr(adata.var, 'columns') else []
            
            # 获取基因名（前20个）
            columns = list(adata.var_names[:20]) if hasattr(adata, 'var_names') else []
            
            # 计算稀疏度
            sparsity = 0
            if hasattr(adata.X, 'nnz'):
                total_cells = n_obs * n_vars
                sparsity = (1 - adata.X.nnz / total_cells) * 100 if total_cells > 0 else 0
            
            # 检查数据值范围（采样）
            data_range = {}
            try:
                sample_size = min(1000, n_obs * n_vars)
                if sample_size > 0:
                    if hasattr(adata.X, 'toarray'):
                        sample_data = adata.X[:min(100, n_obs), :min(100, n_vars)].toarray()
                    else:
                        sample_data = adata.X[:min(100, n_obs), :min(100, n_vars)]
                    
                    data_range = {
                        "min": float(np.min(sample_data)),
                        "max": float(np.max(sample_data)),
                        "mean": float(np.mean(sample_data)),
                        "median": float(np.median(sample_data))
                    }
            except Exception as e:
                logger.warning(f"⚠️ Could not calculate data range: {e}")
            
            # 检查是否有聚类结果
            has_clusters = 'leiden' in obs_keys or 'louvain' in obs_keys
            has_umap = 'X_umap' in adata.obsm_keys() if hasattr(adata, 'obsm_keys') else False
            
            # 构建预览
            preview = {
                "n_obs": n_obs,
                "n_vars": n_vars,
                "obs_keys": obs_keys[:10],
                "var_keys": var_keys[:10],
                "sparsity": round(sparsity, 2),
                "has_clusters": has_clusters,
                "has_umap": has_umap
            }
            
            return {
                "status": "success",
                "file_path": str(path.resolve()),
                "file_type": "anndata",
                "shape": {
                    "rows": n_obs,
                    "cols": n_vars
                },
                "n_obs": n_obs,
                "n_vars": n_vars,
                "n_samples": n_obs,  # 兼容性
                "n_features": n_vars,  # 兼容性
                "columns": columns,
                "head": {
                    "markdown": f"H5AD 文件\n- 细胞数: {n_obs}\n- 基因数: {n_vars}\n- 稀疏度: {sparsity:.2f}%",
                    "json": preview
                },
                "obs_keys": obs_keys,
                "var_keys": var_keys,
                "sparsity": sparsity,
                "has_clusters": has_clusters,
                "has_umap": has_umap,
                "data_range": data_range,
                "data": {
                    "summary": {
                        "n_samples": n_obs,
                        "n_features": n_vars,
                        "sparsity": round(sparsity, 2)
                    }
                }
            }
            
        except ImportError:
            return {
                "status": "error",
                "error": "scanpy not installed. Please install: pip install scanpy",
                "file_type": "anndata"
            }
        except Exception as e:
            logger.error(f"❌ [AnnDataHandler] Failed: {e}", exc_info=True)
            return {
                "status": "error",
                "error": f"无法读取H5AD文件: {str(e)}",
                "file_type": "anndata"
            }


@register_inspector(priority=1)
class TabularHandler(BaseFileHandler):
    """
    表格文件检查器（优先级：1）
    
    支持 CSV, TSV, TXT, XLSX 文件
    """
    
    # 支持的文件扩展名
    SUPPORTED_EXTENSIONS = {'.csv', '.tsv', '.txt', '.xlsx', '.xls'}
    
    def can_handle(self, path: Path) -> bool:
        """检查是否为支持的表格文件"""
        if not path.is_file():
            return False
        if path_looks_like_medical_imaging(path):
            return False
        return path.suffix.lower() in self.SUPPORTED_EXTENSIONS
    
    def _count_csv_lines(self, file_path: Path, separator: str = ',') -> Optional[int]:
        """统计 CSV/TSV 文件行数（避免加载全部数据）"""
        try:
            count = 0
            with open(file_path, 'r', encoding='utf-8', errors='ignore') as f:
                # 跳过第一行（表头）
                next(f, None)
                for _ in f:
                    count += 1
            return count
        except Exception:
            return None
    
    def inspect(self, path: Path) -> Dict[str, Any]:
        """检查表格文件"""
        if path_looks_like_medical_imaging(path):
            return build_medical_imaging_inspection_result(path)
        try:
            import pandas as pd
            
            logger.info(f"🔍 [TabularHandler] Loading: {path}")
            
            absolute_path = str(path.resolve())
            file_size_mb = path.stat().st_size / (1024 * 1024)
            
            # 检测分隔符
            separator = ','
            try:
                with open(absolute_path, 'r', encoding='utf-8', errors='ignore') as f:
                    first_line = f.readline()
                    if '\t' in first_line:
                        separator = '\t'
                    elif ',' in first_line:
                        separator = ','
                    elif ';' in first_line:
                        separator = ';'
            except Exception:
                pass
            
            # 读取预览（最多 100 行，避免规划阶段阻塞）
            LARGE_FILE_THRESHOLD_MB = 200
            PREVIEW_MAX_ROWS = 100
            preview_rows = 10

            if file_size_mb < LARGE_FILE_THRESHOLD_MB:
                # 小文件：仅读前 100 行用于元数据，避免大 CSV 阻塞
                df = pd.read_csv(absolute_path, sep=separator, nrows=PREVIEW_MAX_ROWS)
                is_sampled = True
                total_rows = self._count_csv_lines(path, separator)
                if total_rows is None:
                    total_rows = len(df)
            else:
                # 大文件：采样读取
                df = pd.read_csv(absolute_path, sep=separator, nrows=preview_rows)
                is_sampled = True
                total_rows = self._count_csv_lines(path, separator)
                if total_rows is None:
                    total_rows = len(df)
            
            # 识别列类型
            metadata_cols = []
            numeric_cols = []
            
            for col in df.columns:
                if pd.api.types.is_numeric_dtype(df[col]):
                    numeric_cols.append(col)
                else:
                    metadata_cols.append(col)
            
            n_samples = total_rows if not is_sampled else (total_rows if total_rows else len(df))
            n_features = len(numeric_cols)
            
            # 缺失值统计
            missing_rate = 0
            if len(numeric_cols) > 0:
                numeric_data = df[numeric_cols]
                total_cells = numeric_data.size
                missing_cells = numeric_data.isnull().sum().sum()
                missing_rate = (missing_cells / total_cells * 100) if total_cells > 0 else 0
            
            # 数据范围
            data_range = {}
            if len(numeric_cols) > 0:
                numeric_data = df[numeric_cols]
                data_range = {
                    "min": float(numeric_data.min().min()),
                    "max": float(numeric_data.max().max()),
                    "mean": float(numeric_data.mean().mean()),
                    "median": float(numeric_data.median().median())
                }
            
            # 构建预览
            try:
                from tabulate import tabulate
                head_markdown = tabulate(df.head(preview_rows), headers='keys', tablefmt='grid', showindex=False)
            except ImportError:
                # 回退到 CSV 格式
                head_markdown = df.head(preview_rows).to_csv(index=False)
            
            head_json = df.head(preview_rows).to_dict(orient='records')
            
            return {
                "status": "success",
                "file_path": absolute_path,
                "file_type": "tabular",
                "file_size_mb": round(file_size_mb, 2),
                "is_sampled": is_sampled,
                "separator": separator,
                "columns": list(df.columns),
                "shape": {
                    "rows": n_samples,
                    "cols": len(df.columns)
                },
                "n_samples": n_samples,
                "n_features": n_features,
                "head": {
                    "markdown": head_markdown,
                    "json": head_json
                },
                "metadata_columns": metadata_cols,
                "feature_columns": numeric_cols[:20],
                "total_feature_columns": len(numeric_cols),
                "missing_rate": round(missing_rate, 2),
                "data_range": data_range,
                "data": {
                    "summary": {
                        "n_samples": n_samples,
                        "n_features": n_features,
                        "missing_rate": round(missing_rate, 2),
                        "data_range": data_range,
                        "is_sampled": is_sampled
                    }
                }
            }
            
        except ImportError:
            return {
                "status": "error",
                "error": "pandas not installed. Please install: pip install pandas",
                "file_type": "tabular"
            }
        except Exception as e:
            logger.error(f"❌ [TabularHandler] Failed: {e}", exc_info=True)
            return {
                "status": "error",
                "error": str(e),
                "error_type": type(e).__name__,
                "file_type": "tabular"
            }


# ============================================
# Part 3: The Facade (FileInspector)
# ============================================

class FileInspector:
    """
    文件检测器 - 系统的"通用眼睛"（门面模式）
    
    使用策略模式，自动选择合适的检查器处理文件
    """
    
    def __init__(self, upload_dir: str):
        """
        初始化文件检测器
        
        Args:
            upload_dir: 上传文件目录
        """
        self.upload_dir = Path(upload_dir)
        try:
            self.upload_dir.mkdir(parents=True, exist_ok=True)
        except (PermissionError, OSError) as e:
            logger.warning(f"⚠️ 无法创建上传目录 {self.upload_dir}: {e}")
        
        # 常见 Docker 挂载路径列表
        self.common_mount_paths = [
            "/app/uploads",
            "/app/data/uploads",
            "/app/data",
            "/workspace/uploads",
            "./uploads",
            "./data"
        ]
        
        # 获取所有已注册的检查器
        self.handlers = [handler_class() for handler_class in _registry.get_handlers()]
        logger.info(f"✅ [FileInspector] Loaded {len(self.handlers)} handlers")
    
    def _resolve_actual_path(self, file_path: str) -> Tuple[Optional[str], List[str]]:
        """
        智能路径解析：尝试在多个常见路径中查找文件或目录
        
        Args:
            file_path: 原始文件路径
            
        Returns:
            (actual_path, searched_paths)
        """
        searched_paths = []
        original_path = Path(file_path)
        
        # Step 1: 检查原始路径
        if original_path.exists():
            resolved_path = str(original_path.resolve())
            searched_paths.append(resolved_path)
            logger.info(f"✅ [Smart Path Resolution] Found: {resolved_path}")
            return resolved_path, searched_paths
        
        if original_path.is_absolute():
            searched_paths.append(str(original_path))
        else:
            cwd_path = Path(os.getcwd()) / original_path
            searched_paths.append(str(cwd_path.resolve()))
        
        # Step 1.5: 相对路径时优先在 upload_dir 下按完整路径查找（修复 guest/session/10x_data_xxx 无法识别）
        if not original_path.is_absolute() and file_path and self.upload_dir:
            upload_base = Path(self.upload_dir)
            try:
                candidate = (upload_base / file_path).resolve()
                candidate = candidate.resolve()
                if candidate.exists():
                    searched_paths.append(str(candidate))
                    logger.info(f"✅ [Smart Path Resolution] Found under upload_dir: {candidate}")
                    return str(candidate), searched_paths
                searched_paths.append(str(candidate))
            except (OSError, ValueError) as e:
                logger.debug(f"⚠️ upload_dir + file_path resolve failed: {e}")
        
        # Step 2: 提取文件名或目录名
        path_name = original_path.name if original_path.name else str(original_path)
        if not path_name or path_name == '.':
            return None, searched_paths
        
        # Step 3: 在常见挂载路径中搜索
        for mount_path in self.common_mount_paths:
            mount_path_obj = Path(mount_path)
            if not mount_path_obj.is_absolute():
                mount_path_obj = Path(os.getcwd()) / mount_path_obj
            
            try:
                resolved_mount = mount_path_obj.resolve()
                candidate_path = resolved_mount / path_name
                searched_paths.append(str(candidate_path))
                
                if candidate_path.exists():
                    logger.info(f"✅ [Smart Path Resolution] Found: {candidate_path}")
                    return str(candidate_path), searched_paths
            except (OSError, ValueError) as e:
                logger.debug(f"⚠️ Invalid path {mount_path}: {e}")
                continue
        
        logger.warning(f"❌ [Smart Path Resolution] Not found: {path_name}")
        return None, searched_paths
    
    def inspect_file(self, file_path: str) -> Dict[str, Any]:
        """
        多模态文件检查主入口（分发器）
        
        Args:
            file_path: 文件路径（相对或绝对）
        
        Returns:
            包含检查结果的字典
        """
        # Step 1: 使用智能路径解析
        actual_path, searched_paths = self._resolve_actual_path(file_path)
        
        if actual_path is None:
            current_cwd = os.getcwd()
            error_msg = (
                f"File or directory not found: '{file_path}'\n\n"
                f"**Searched locations ({len(searched_paths)}):**\n"
                + "\n".join(f"  - {path}" for path in searched_paths[:10])
                + f"\n\n**Current working directory:** {current_cwd}\n"
                f"**Upload directory (configured):** {self.upload_dir}\n"
                f"**Environment UPLOAD_DIR:** {os.getenv('UPLOAD_DIR', 'Not set')}"
            )
            
            logger.error(f"❌ [FileInspector] {error_msg}")
            
            return {
                "status": "error",
                "success": False,
                "error": error_msg,
                "file_type": "unknown",
                "file_path": file_path,
                "searched_paths": searched_paths,
                "current_cwd": current_cwd
            }
        
        # Step 2: 转换为 Path 对象
        path = Path(actual_path)

        # Step 2b: 医学影像（.nii / .nii.gz / .dcm）严禁走表格类 pandas 逻辑，优先旁路体检
        if path.is_file() and path_looks_like_medical_imaging(path):
            result = build_medical_imaging_inspection_result(path)
            if result.get("status") == "success":
                result["file_path"] = str(path.resolve())
                result["success"] = True
            return result
        
        # Step 3: 遍历所有检查器，找到第一个可以处理的
        for handler in self.handlers:
            try:
                if handler.can_handle(path):
                    logger.info(f"✅ [FileInspector] Using handler: {handler.__class__.__name__}")
                    result = handler.inspect(path)
                    # 确保返回绝对路径
                    if result.get("status") == "success":
                        result["file_path"] = str(path.resolve())
                        result["success"] = True
                    return result
            except Exception as e:
                logger.warning(f"⚠️ Handler {handler.__class__.__name__} failed: {e}")
                # 继续尝试下一个检查器
                continue
        
        # Step 4: 特殊处理：单独的 .mtx 文件
        if path.suffix == '.mtx' or path.name.lower() == 'matrix.mtx':
            error_msg = (
                "⚠️ **单独的 matrix.mtx 文件无法单独使用**\n\n"
                "10x Genomics 数据需要三个文件：\n"
                "1. `matrix.mtx` - 表达矩阵\n"
                "2. `barcodes.tsv` - 细胞条形码\n"
                "3. `features.tsv` 或 `genes.tsv` - 基因信息\n\n"
                "**解决方案：**\n"
                "- 请同时上传这三个文件（或包含这三个文件的目录）\n"
                "- 如果文件已解压，请将它们放在同一个目录中上传\n"
                "- 如果文件已压缩，请上传压缩后的文件（.gz 格式）\n\n"
                f"**当前文件：** {path.name}"
            )
            logger.warning(f"⚠️ {error_msg}")
            return {
                "status": "error",
                "success": False,
                "error": error_msg,
                "file_type": "mtx",
                "file_path": str(path.resolve())
            }
        
        # Step 5: 如果没有检查器可以处理，返回错误
        if not path.suffix and path.is_file():
            error_msg = (
                f"无法识别文件类型: {path.name}（无扩展名）。"
                "请确认：若为数据文件请使用正确扩展名（如 .h5ad、.csv、.nii）；"
                "若为 10x 目录请上传包含 matrix.mtx、barcodes.tsv、features.tsv 的目录或压缩包。"
            )
        else:
            error_msg = f"无法识别文件类型: {path.name}"
        logger.error(f"❌ [FileInspector] {error_msg}")
        
        return {
            "status": "error",
            "success": False,
            "error": error_msg,
            "file_type": "unknown",
            "file_path": str(path.resolve())
        }


# ============================================
# Extended handlers (Spatial Visium before generic 10x)
# ============================================
try:
    from .file_handlers.extended_handlers import SpatialVisiumHandler
    _registry.register(20, SpatialVisiumHandler)
    logger.info("✅ [FileInspector] Registered SpatialVisiumHandler (priority=20)")
except ImportError as e:
    logger.debug("SpatialVisiumHandler not loaded: %s", e)

try:
    from .file_handlers.extended_handlers import RadiomicsHandler
    _registry.register(10, RadiomicsHandler)
    logger.info("✅ [FileInspector] Registered RadiomicsHandler (priority=10)")
except ImportError as e:
    logger.debug("RadiomicsHandler not loaded: %s", e)
