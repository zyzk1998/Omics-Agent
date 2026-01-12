"""
文件检测器 - Universal Eyes of the System
多模态文件检查器，支持表格数据、单细胞数据、图像等
"""
import os
import json
import gzip
import logging
from pathlib import Path
from typing import Dict, Optional, Any, List, Tuple
import numpy as np

logger = logging.getLogger(__name__)


class FileInspector:
    """
    文件检测器 - 系统的"通用眼睛"
    
    支持多模态文件检查：
    - Tabular: CSV, TSV, TXT, XLSX
    - Single-Cell: H5AD, MTX, 10x Genomics
    - Image: JPG, PNG, TIFF
    """
    
    # 文件大小阈值（MB）
    LARGE_FILE_THRESHOLD_MB = 200
    SAMPLE_SIZE_LARGE_FILE = 10000
    
    def __init__(self, upload_dir: str):
        """
        初始化文件检测器
        
        Args:
            upload_dir: 上传文件目录
        """
        self.upload_dir = Path(upload_dir)
        # 尝试创建目录，如果失败则记录警告（可能在容器外测试）
        try:
            self.upload_dir.mkdir(parents=True, exist_ok=True)
        except (PermissionError, OSError) as e:
            logger.warning(f"⚠️ 无法创建上传目录 {self.upload_dir}: {e}。将在运行时尝试查找文件。")
        
        # 常见 Docker 挂载路径列表（用于智能路径解析）
        self.common_mount_paths = [
            "/app/uploads",
            "/app/data/uploads",
            "/app/data",
            "/workspace/uploads",
            "./uploads",
            "./data"
        ]
    
    def _resolve_actual_path(self, file_path: str) -> Tuple[Optional[str], List[str]]:
        """
        智能路径解析：尝试在多个常见路径中查找文件
        
        Args:
            file_path: 原始文件路径（可能是相对路径、绝对路径或仅文件名）
        
        Returns:
            (actual_path, searched_paths): 
            - actual_path: 找到的实际路径，如果未找到则为 None
            - searched_paths: 已搜索的路径列表（用于错误报告）
        """
        import os
        from pathlib import Path
        
        searched_paths = []
        
        # Step 1: 检查原始路径是否存在
        original_path = Path(file_path)
        if original_path.exists():
            searched_paths.append(str(original_path.resolve()))
            return str(original_path.resolve()), searched_paths
        
        # 如果原始路径不存在，记录它
        if original_path.is_absolute():
            searched_paths.append(str(original_path))
        else:
            # 尝试相对于当前工作目录
            cwd_path = Path(os.getcwd()) / original_path
            searched_paths.append(str(cwd_path.resolve()))
        
        # Step 2: 提取文件名
        filename = original_path.name
        if not filename:
            # 如果路径是目录或无效，返回 None
            return None, searched_paths
        
        # Step 3: 在常见挂载路径中搜索
        for mount_path in self.common_mount_paths:
            mount_path_obj = Path(mount_path)
            
            # 如果是相对路径，转换为绝对路径
            if not mount_path_obj.is_absolute():
                mount_path_obj = Path(os.getcwd()) / mount_path_obj
            
            # 尝试解析路径
            try:
                resolved_mount = mount_path_obj.resolve()
                candidate_path = resolved_mount / filename
                
                searched_paths.append(str(candidate_path))
                
                if candidate_path.exists() and candidate_path.is_file():
                    logger.info(f"✅ [Smart Path Resolution] Found file at: {candidate_path}")
                    return str(candidate_path), searched_paths
            except (OSError, ValueError) as e:
                # 路径无效，跳过
                logger.debug(f"⚠️ [Smart Path Resolution] Invalid path {mount_path}: {e}")
                continue
        
        # Step 4: 如果仍未找到，返回 None 和已搜索的路径列表
        logger.warning(f"❌ [Smart Path Resolution] File not found: {filename}")
        logger.warning(f"   Searched in {len(searched_paths)} locations")
        return None, searched_paths
    
    def inspect_file(self, file_path: str) -> Dict[str, Any]:
        """
        多模态文件检查主入口（分发器）
        
        🔥 升级：使用智能路径解析，自动在多个常见路径中查找文件
        
        Args:
            file_path: 文件路径（相对或绝对）
        
        Returns:
            包含检查结果的字典，包含绝对路径（file_path字段）
        """
        import os
        
        # 🔥 Step 1: 使用智能路径解析
        actual_path, searched_paths = self._resolve_actual_path(file_path)
        
        if actual_path is None:
            # 🔥 CRITICAL: 返回详细的错误信息，列出所有搜索过的路径
            current_cwd = os.getcwd()
            error_msg = (
                f"File not found: '{file_path}'\n\n"
                f"**Searched locations ({len(searched_paths)}):**\n"
                + "\n".join(f"  - {path}" for path in searched_paths[:10])  # 最多显示10个
                + f"\n\n**Current working directory:** {current_cwd}\n"
                f"**Upload directory (configured):** {self.upload_dir}\n"
                f"**Environment UPLOAD_DIR:** {os.getenv('UPLOAD_DIR', 'Not set')}"
            )
            
            logger.error(f"❌ [FileInspector] {error_msg}")
            
            return {
                "status": "error",
                "success": False,  # 🔥 添加 success 字段用于前端检查
                "error": error_msg,
                "file_type": "unknown",
                "file_path": file_path,  # 返回原始路径
                "searched_paths": searched_paths,  # 调试信息
                "current_cwd": current_cwd
            }
        
        # Step 2: 使用找到的实际路径继续处理
        filepath = Path(actual_path)
        absolute_path = str(filepath.resolve())
        
        logger.info(f"✅ [FileInspector] Using resolved path: {absolute_path}")
        
        # 检查文件类型并分发到相应的检查器（传递绝对路径）
        if filepath.is_dir():
            # 10x Genomics 目录
            return self._inspect_anndata(absolute_path)
        
        file_ext = filepath.suffix.lower()
        
        # Tabular 数据
        if file_ext in ['.csv', '.tsv', '.txt', '.xlsx']:
            return self._inspect_tabular(absolute_path)
        
        # Single-Cell 数据
        elif file_ext in ['.h5ad', '.mtx'] or file_ext.endswith('.mtx.gz'):
            return self._inspect_anndata(absolute_path)
        
        # 图像数据
        elif file_ext in ['.jpg', '.jpeg', '.png', '.tiff', '.tif']:
            return self._inspect_image(absolute_path)
        
        else:
            return {
                "status": "error",
                "error": f"Unsupported file type: {file_ext}",
                "file_type": "unknown",
                "file_path": absolute_path  # 返回绝对路径
            }
    
    def _inspect_tabular(
        self,
        file_path: str
    ) -> Dict[str, Any]:
        """
        检查表格数据（CSV, TSV, TXT, XLSX）
        
        策略：
        - 小文件（<200MB）：完整读取以获得准确统计
        - 大文件（>200MB）：采样10000行防止OOM
        
        Args:
            file_path: 文件路径（相对或绝对）
        
        Returns:
            包含统计信息的字典，包括：
            - file_path: 绝对路径
            - columns: 列名列表
            - shape: (rows, cols) 元组
            - head: 前10行（markdown格式）
            - separator: 分隔符（',' 或 '\t'）
        """
        import pandas as pd
        
        try:
            # 🔥 Step 1: 转换为绝对路径
            filepath = Path(file_path)
            if not filepath.is_absolute():
                # 如果是相对路径，尝试相对于 upload_dir
                filepath = self.upload_dir / filepath
            absolute_path = str(filepath.resolve())
            
            if not filepath.exists():
                return {
                    "status": "error",
                    "error": f"File not found: {absolute_path}",
                    "file_type": "tabular"
                }
            
            file_size_mb = filepath.stat().st_size / (1024 * 1024)
            
            logger.info(f"🔍 [Tabular Inspector] File: {absolute_path}, Size: {file_size_mb:.2f} MB")
            
            # 🔥 Step 2: 检测分隔符（先读第一行）
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
            except:
                pass  # 使用默认分隔符
            
            # 决定读取策略
            if file_size_mb < self.LARGE_FILE_THRESHOLD_MB:
                # 小文件：完整读取
                logger.info(f"   📖 Reading full file (size < {self.LARGE_FILE_THRESHOLD_MB}MB)")
                df = pd.read_csv(absolute_path, sep=separator)
                is_sampled = False
            else:
                # 大文件：采样读取
                logger.info(f"   📖 Sampling {self.SAMPLE_SIZE_LARGE_FILE} rows (size >= {self.LARGE_FILE_THRESHOLD_MB}MB)")
                df = pd.read_csv(absolute_path, sep=separator, nrows=self.SAMPLE_SIZE_LARGE_FILE)
                is_sampled = True
                # 估算总行数
                try:
                    with open(absolute_path, 'r', encoding='utf-8', errors='ignore') as f:
                        total_lines = sum(1 for _ in f) - 1  # 减去表头
                except:
                    total_lines = None
            
            # 识别列类型
            metadata_cols = []
            numeric_cols = []
            
            for col in df.columns:
                if pd.api.types.is_numeric_dtype(df[col]):
                    numeric_cols.append(col)
                else:
                    metadata_cols.append(col)
            
            # 计算统计信息
            n_samples = len(df) if not is_sampled else (total_lines if total_lines else len(df))
            n_features = len(numeric_cols)
            
            # 缺失值统计（关键：需要完整数据或足够大的采样）
            if len(numeric_cols) > 0:
                numeric_data = df[numeric_cols]
                total_cells = numeric_data.size
                missing_cells = numeric_data.isnull().sum().sum()
                missing_rate = (missing_cells / total_cells * 100) if total_cells > 0 else 0
            else:
                missing_rate = 0
            
            # 数据范围（关键：用于判断是否需要 Log2 变换）
            data_range = {}
            if len(numeric_cols) > 0:
                numeric_data = df[numeric_cols]
                data_range = {
                    "min": float(numeric_data.min().min()),
                    "max": float(numeric_data.max().max()),
                    "mean": float(numeric_data.mean().mean()),
                    "median": float(numeric_data.median().median())
                }
            
            # 识别潜在的分组列（唯一值较少的列）
            potential_groups = {}
            for col in metadata_cols:
                unique_count = df[col].nunique()
                if unique_count > 1 and unique_count <= min(20, len(df) // 2):
                    potential_groups[col] = {
                        "n_unique": int(unique_count),
                        "values": df[col].unique().tolist()[:10]  # 只显示前10个
                    }
            
            # 🔥 Step 3: 提取前10行作为 head（markdown格式）
            head_df = df.head(10)
            head_json = head_df.to_dict(orient='records')
            
            # 🔥 修复：处理 tabulate 依赖缺失问题
            try:
                head_markdown = head_df.to_markdown(index=False)
            except ImportError:
                # 如果 tabulate 不可用，生成简单的文本表格
                logger.warning("⚠️ tabulate 不可用，使用文本格式替代 markdown")
                try:
                    # 使用 pandas 的 to_string 方法生成文本表格
                    head_markdown = head_df.to_string(index=False, max_rows=10)
                except Exception as e:
                    # 如果 to_string 也失败，使用简单的 CSV 格式
                    logger.warning(f"⚠️ to_string 失败，使用 CSV 格式: {e}")
                    head_markdown = head_df.to_csv(index=False)
            
            # 构建结果（包含完整元数据，作为单一数据源）
            result = {
                "status": "success",
                "file_path": absolute_path,  # 🔥 绝对路径
                "file_type": "tabular",
                "file_size_mb": round(file_size_mb, 2),
                "is_sampled": is_sampled,
                "separator": separator,  # 🔥 分隔符
                "columns": list(df.columns),  # 🔥 所有列名
                "shape": {
                    "rows": n_samples if isinstance(n_samples, int) else (total_lines if total_lines else len(df)),
                    "cols": len(df.columns)
                },
                "head": {  # 🔥 前10行数据
                    "markdown": head_markdown,
                    "json": head_json
                },
                "n_samples": n_samples if isinstance(n_samples, int) else f"~{n_samples}",
                "n_features": n_features,
                "n_metadata_cols": len(metadata_cols),
                "metadata_columns": metadata_cols,
                "feature_columns": numeric_cols[:20],  # 只显示前20个特征列名
                "total_feature_columns": len(numeric_cols),
                "missing_rate": round(missing_rate, 2),
                "data_range": data_range,
                "potential_groups": potential_groups,
                # 前端可用的摘要数据
                "data": {
                    "summary": {
                        "n_samples": n_samples if isinstance(n_samples, int) else f"~{n_samples}",
                        "n_features": n_features,
                        "missing_rate": round(missing_rate, 2),
                        "data_range": data_range,
                        "is_sampled": is_sampled
                    }
                }
            }
            
            logger.info(f"✅ [Tabular Inspector] Success: {n_samples} samples, {n_features} features, missing_rate={missing_rate:.2f}%")
            return result
            
        except Exception as e:
            logger.error(f"❌ [Tabular Inspector] Failed: {e}", exc_info=True)
            return {
                "status": "error",
                "error": str(e),
                "error_type": type(e).__name__,
                "file_type": "tabular"
            }
    
    def _inspect_anndata(
        self,
        file_path: str
    ) -> Dict[str, Any]:
        """
        检查单细胞数据（H5AD, MTX, 10x Genomics）
        
        Args:
            file_path: 文件路径或目录路径
        
        Returns:
            包含检查结果的字典
        """
        try:
            import scanpy as sc
            
            filepath = Path(file_path)
            
            logger.info(f"🔍 [Anndata Inspector] File: {file_path}")
            
            # 处理目录（10x Genomics）
            if filepath.is_dir():
                dir_contents = os.listdir(filepath)
                if any(f in dir_contents for f in ['matrix.mtx', 'matrix.mtx.gz']):
                    # 10x 格式需要完整加载
                    logger.info("   📖 Loading 10x Genomics directory...")
                    adata = sc.read_10x_mtx(filepath)
                else:
                    return {
                        "status": "error",
                        "error": f"Unknown directory format: {file_path}",
                        "file_type": "directory"
                    }
            # 处理 H5AD 文件
            elif filepath.suffix == '.h5ad':
                # 使用 backed='r' 模式（只读，不加载全部数据到内存）
                logger.info("   📖 Loading H5AD file (backed mode)...")
                try:
                    adata = sc.read_h5ad(file_path, backed='r')
                except:
                    # 如果 backed 模式失败，使用普通模式
                    logger.warning("   ⚠️ Backed mode failed, using normal mode")
                    adata = sc.read_h5ad(file_path)
            # 处理 MTX 文件
            else:
                logger.info("   📖 Loading MTX file...")
                adata = sc.read(file_path)
            
            # 提取基本信息
            n_obs = adata.n_obs
            n_vars = adata.n_vars
            obs_keys = list(adata.obs.columns) if hasattr(adata.obs, 'columns') else []
            var_keys = list(adata.var.columns) if hasattr(adata.var, 'columns') else []
            
            # 计算稀疏度
            if hasattr(adata.X, 'nnz'):
                # 稀疏矩阵
                total_cells = n_obs * n_vars
                sparsity = (1 - adata.X.nnz / total_cells) * 100 if total_cells > 0 else 0
            else:
                # 密集矩阵
                sparsity = 0
            
            # 检查数据值范围（采样检查，不加载全部数据）
            data_range = {}
            try:
                sample_size = min(1000, n_obs * n_vars)
                if sample_size > 0:
                    if hasattr(adata.X, 'toarray'):
                        # 稀疏矩阵：采样检查
                        sample_data = adata.X[:min(100, n_obs), :min(100, n_vars)].toarray()
                    else:
                        # 密集矩阵：采样检查
                        sample_data = adata.X[:min(100, n_obs), :min(100, n_vars)]
                    
                    data_range = {
                        "min": float(np.min(sample_data)),
                        "max": float(np.max(sample_data)),
                        "mean": float(np.mean(sample_data)),
                        "median": float(np.median(sample_data))
                    }
            except Exception as e:
                logger.warning(f"   ⚠️ Could not calculate data range: {e}")
            
            # 检查是否有聚类结果
            has_clusters = 'leiden' in obs_keys or 'louvain' in obs_keys
            has_umap = 'X_umap' in adata.obsm_keys() if hasattr(adata, 'obsm_keys') else False
            
            result = {
                "status": "success",
                "file_path": file_path,
                "file_type": "anndata",
                "n_obs": n_obs,
                "n_vars": n_vars,
                "obs_keys": obs_keys,
                "var_keys": var_keys,
                "sparsity": round(sparsity, 2),
                "data_range": data_range,
                "has_clusters": has_clusters,
                "has_umap": has_umap,
                # 前端可用的摘要数据
                "data": {
                    "summary": {
                        "n_obs": n_obs,
                        "n_vars": n_vars,
                        "sparsity": round(sparsity, 2),
                        "data_range": data_range,
                        "has_clusters": has_clusters,
                        "has_umap": has_umap
                    }
                }
            }
            
            logger.info(f"✅ [Anndata Inspector] Success: {n_obs} cells, {n_vars} genes, sparsity={sparsity:.2f}%")
            return result
            
        except ImportError:
            logger.error("❌ [Anndata Inspector] scanpy not available")
            return {
                "status": "error",
                "error": "scanpy not installed",
                "file_type": "anndata"
            }
        except Exception as e:
            logger.error(f"❌ [Anndata Inspector] Failed: {e}", exc_info=True)
            return {
                "status": "error",
                "error": str(e),
                "error_type": type(e).__name__,
                "file_type": "anndata"
            }
    
    def _inspect_image(
        self,
        file_path: str
    ) -> Dict[str, Any]:
        """
        检查图像文件（JPG, PNG, TIFF）
        
        Args:
            file_path: 文件路径
        
        Returns:
            包含检查结果的字典
        """
        try:
            from PIL import Image
            
            logger.info(f"🔍 [Image Inspector] File: {file_path}")
            
            with Image.open(file_path) as img:
                width, height = img.size
                mode = img.mode
                format_name = img.format
                
                # 获取文件大小
                filepath = Path(file_path)
                file_size_mb = filepath.stat().st_size / (1024 * 1024)
            
            result = {
                "status": "success",
                "file_path": file_path,
                "file_type": "image",
                "file_size_mb": round(file_size_mb, 2),
                "width": width,
                "height": height,
                "mode": mode,
                "format": format_name,
                # 前端可用的摘要数据
                "data": {
                    "summary": {
                        "width": width,
                        "height": height,
                        "mode": mode,
                        "format": format_name,
                        "file_size_mb": round(file_size_mb, 2)
                    }
                }
            }
            
            logger.info(f"✅ [Image Inspector] Success: {width}x{height}, mode={mode}, format={format_name}")
            return result
            
        except ImportError:
            logger.error("❌ [Image Inspector] PIL not available")
            return {
                "status": "error",
                "error": "PIL (Pillow) not installed",
                "file_type": "image"
            }
        except Exception as e:
            logger.error(f"❌ [Image Inspector] Failed: {e}", exc_info=True)
            return {
                "status": "error",
                "error": str(e),
                "error_type": type(e).__name__,
                "file_type": "image"
            }
    
    # ================= 保留旧方法以保持兼容性 =================
    
    def _is_gzipped(self, filepath: Path) -> bool:
        """检查文件是否为 gzip 压缩"""
        if filepath.is_dir():
            return False
        try:
            with open(filepath, 'rb') as f:
                return f.read(2) == b'\x1f\x8b'
        except:
            return False
    
    def _read_head(self, filepath: Path, lines: int = 5) -> list:
        """安全读取文件前几行"""
        if filepath.is_dir():
            return []
        try:
            if self._is_gzipped(filepath):
                with gzip.open(filepath, 'rt') as f:
                    return [next(f).strip() for _ in range(lines) if True]
            else:
                with open(filepath, 'r', encoding='utf-8', errors='ignore') as f:
                    return [next(f).strip() for _ in range(lines) if True]
        except Exception:
            return []
    
    def generate_metadata(self, filename: str) -> Optional[Dict]:
        """
        主动检测：分析单个文件并保存 .meta.json（保留兼容性）
        
        Args:
            filename: 文件名或目录名
        
        Returns:
            元数据字典，如果失败返回 None
        """
        filepath = self.upload_dir / filename
        if not filepath.exists():
            return None
        
        # 使用新的 inspect_file 方法
        result = self.inspect_file(str(filepath))
        
        if result.get("status") == "success":
            # 保存元数据
            meta_path = filepath.with_suffix(filepath.suffix + '.meta.json')
            try:
                with open(meta_path, 'w', encoding='utf-8') as f:
                    json.dump(result, f, indent=2, ensure_ascii=False)
            except Exception as e:
                logger.warning(f"⚠️ 保存元数据失败: {e}")
        
        return result
    
    def get_metadata(self, filename: str) -> Optional[Dict]:
        """
        获取文件的元数据（如果已生成）
        
        Args:
            filename: 文件名
        
        Returns:
            元数据字典，如果不存在返回 None
        """
        filepath = self.upload_dir / filename
        meta_path = filepath.with_suffix(filepath.suffix + '.meta.json')
        
        if meta_path.exists():
            try:
                with open(meta_path, 'r', encoding='utf-8') as f:
                    return json.load(f)
            except Exception:
                return None
        return None
