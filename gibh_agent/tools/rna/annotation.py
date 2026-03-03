"""
单细胞数据注释工具 - 细胞类型注释
"""
import os
import time
import logging
from typing import Dict, Any, Optional, List
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from ...core.tool_registry import registry

logger = logging.getLogger(__name__)


def _get_fallback_markers() -> Dict[str, str]:
    """
    内置 PBMC/免疫细胞 marker 字典：基因名 -> 细胞类型。
    用于无 .pkl 模型时基于表达的规则注释，输出具生物学意义的标签（如 T cells, B cells）。
    """
    return {
        "CD3D": "T cells",
        "CD3E": "T cells",
        "CD4": "CD4+ T cells",
        "CD8A": "CD8+ T cells",
        "MS4A1": "B cells",   # CD20
        "CD79A": "B cells",
        "GNLY": "NK cells",
        "NKG7": "NK cells",
        "CD14": "Monocytes",
        "LYZ": "Monocytes",
        "FCGR3A": "CD16+ Monocytes",
        "S100A8": "Monocytes",
        "S100A9": "Monocytes",
        "PPBP": "Platelets",
        "PF4": "Platelets",
        "FCER1A": "pDC",
        "CST3": "pDC",
        "IL3RA": "pDC",
    }


def _fallback_marker_annotation(
    adata,
    adata_path: str,
    output_dir: Optional[str],
    sc,
) -> Dict[str, Any]:
    """
    当 CellTypist 参考文件不可用时，使用内置 marker 字典按簇打分，
    为每个簇分配生物学细胞类型（T cells, B cells, NK cells 等），而非 Cluster_0/1。
    """
    markers = _get_fallback_markers()
    var_names = set(adata.var_names)
    cluster_key = None
    for key in ("leiden", "louvain", "cluster"):
        if key in adata.obs.columns:
            cluster_key = key
            break
    if cluster_key is None:
        adata.obs["cell_type"] = "Unknown"
        n_cell_types = 1
        label_counts = {"Unknown": len(adata)}
    else:
        # 每个基因可能对应同一类型，合并为 基因 -> 类型，按类型聚合表达
        gene_to_type: Dict[str, str] = {g: t for g, t in markers.items() if g in var_names}
        if not gene_to_type:
            adata.obs["cell_type"] = "Cluster_" + adata.obs[cluster_key].astype(str)
            label_counts = adata.obs["cell_type"].value_counts()
            n_cell_types = len(label_counts)
            logger.warning("No fallback marker genes found in adata.var_names, using cluster labels.")
        else:
            # 按类型取该类型所有 marker 的平均表达作为类型得分
            type_genes: Dict[str, List[str]] = {}
            for g, t in gene_to_type.items():
                type_genes.setdefault(t, []).append(g)
            clusters = adata.obs[cluster_key].astype(str).unique()
            cluster_to_label = {}
            import numpy as np
            for cl in clusters:
                mask = adata.obs[cluster_key].astype(str) == cl
                scores = {}
                for cell_type, genes in type_genes.items():
                    try:
                        sub = adata[mask, genes].X
                        if hasattr(sub, "toarray"):
                            sub = sub.toarray()
                        scores[cell_type] = float(np.mean(sub))
                    except Exception:
                        scores[cell_type] = 0.0
                best = max(scores, key=scores.get)
                cluster_to_label[cl] = best if scores[best] > 0 else f"Cluster_{cl}"
            adata.obs["cell_type"] = adata.obs[cluster_key].astype(str).map(cluster_to_label)
            label_counts = adata.obs["cell_type"].value_counts()
            n_cell_types = len(label_counts)
    plot_path = None
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)
        if "X_umap" in adata.obsm.keys():
            timestamp = int(time.time())
            plot_path = os.path.join(output_dir, f"umap_annotated_{timestamp}.png")
            fig, ax = plt.subplots(figsize=(10, 8))
            sc.pl.umap(
                adata,
                color="cell_type",
                ax=ax,
                show=False,
                title="UMAP: Cell Type (marker-based)",
                legend_loc="right margin",
                frameon=False,
                legend_fontsize=8,
            )
            plt.savefig(plot_path, bbox_inches="tight", dpi=300)
            plt.close()
        output_h5ad = os.path.join(output_dir, "annotated.h5ad")
        adata.write(output_h5ad)
    else:
        out_dir = os.path.dirname(adata_path) or "."
        output_h5ad = os.path.join(out_dir, "annotated.h5ad")
        adata.write(output_h5ad)
    return {
        "status": "success",
        "method": "marker_fallback",
        "model": None,
        "n_cell_types": n_cell_types,
        "cell_types": label_counts.to_dict() if hasattr(label_counts, "to_dict") else dict(label_counts),
        "plot_path": plot_path,
        "output_h5ad": output_h5ad,
        "summary": f"细胞类型注释完成（基于 marker 打分）: 识别到 {n_cell_types} 种类型",
    }


@registry.register(
    name="rna_cell_annotation",
    description="Annotates cell types in single-cell RNA-seq data using CellTypist or marker-based methods. Cell type annotation is crucial for interpreting cell populations.",
    category="scRNA-seq",
    output_type="mixed"
)
def run_cell_annotation(
    adata_path: str,
    method: str = "celltypist",
    model_name: str = "Immune_All_Low.pkl",
    cache_dir: Optional[str] = None,
    output_dir: Optional[str] = None
) -> Dict[str, Any]:
    """
    执行细胞类型注释
    
    Args:
        adata_path: AnnData 文件路径（.h5ad）
        method: 注释方法（"celltypist" 或 "marker"）
        model_name: CellTypist 模型名称（如果使用 celltypist）
        cache_dir: 模型缓存目录
        output_dir: 输出目录（可选）
    
    Returns:
        注释结果字典
    """
    try:
        import scanpy as sc
        
        # 🔥 TASK 3 FIX: 验证并解析输入文件路径
        adata_path_obj = Path(adata_path)
        if not adata_path_obj.is_absolute():
            # 如果是相对路径，尝试在多个位置查找
            potential_paths = [
                adata_path_obj,
                Path(os.getcwd()) / adata_path_obj,
                Path(os.getenv("RESULTS_DIR", "/app/results")) / adata_path_obj,
                Path(os.getenv("UPLOAD_DIR", "/app/uploads")) / adata_path_obj,
            ]
            for potential_path in potential_paths:
                if potential_path.exists():
                    adata_path = str(potential_path.resolve())
                    logger.info(f"✅ [Cell Annotation] 找到输入文件: {adata_path}")
                    break
            else:
                logger.error(f"❌ [Cell Annotation] 无法找到输入文件: {adata_path}")
                return {
                    "status": "error",
                    "error": f"File not found: {adata_path}",
                    "user_message": f"数据文件问题：细胞类型注释步骤无法找到或读取数据文件。",
                    "error_category": "data_issue",
                    "suggestion": f"请检查数据文件路径是否正确，文件是否存在且可读。尝试的路径: {[str(p) for p in potential_paths]}",
                    "can_skip": False
                }
        else:
            if not adata_path_obj.exists():
                logger.error(f"❌ [Cell Annotation] 输入文件不存在: {adata_path}")
                return {
                    "status": "error",
                    "error": f"File not found: {adata_path}",
                    "user_message": f"数据文件问题：细胞类型注释步骤无法找到或读取数据文件。",
                    "error_category": "data_issue",
                    "suggestion": f"请检查数据文件路径是否正确，文件是否存在且可读。路径: {adata_path}",
                    "can_skip": False
                }
        
        # 加载数据
        logger.info(f"📂 [Cell Annotation] 加载数据文件: {adata_path}")
        adata = sc.read_h5ad(adata_path)
        
        if method == "celltypist":
            try:
                import celltypist
                from celltypist import models
                
                # 1) 优先使用持久化参考目录（容器内 /app/data/references 或环境变量）
                ref_dirs = [
                    Path(os.getenv("CELLTYPIST_REFERENCE_DIR", "/app/data/references")),
                    Path(os.getenv("REFERENCE_DIR", "/app/data/references")),
                    Path(os.getcwd()) / "test_data" / "cache",
                ]
                if cache_dir:
                    ref_dirs.insert(0, Path(cache_dir))
                cache_path = None
                model_path = None
                for d in ref_dirs:
                    if d.exists():
                        p = d / model_name
                        if os.path.exists(str(p)):
                            model_path = p
                            cache_path = d
                            logger.info("✅ [Cell Annotation] 使用参考模型: %s", model_path)
                            break
                if model_path is None:
                    # 2) 使用默认缓存目录并确保存在
                    cache_path = Path(cache_dir or os.path.join(os.getcwd(), "test_data", "cache"))
                    cache_path.mkdir(parents=True, exist_ok=True)
                    model_path = cache_path / model_name
                
                # 下载或加载模型（仅当本地尚不存在时尝试下载）
                if not model_path.exists():
                    logger.info(f"📥 正在下载 CellTypist 模型: {model_name}")
                    # 🔥 FIX: celltypist 的 download_models 使用 model 参数，不使用 folder 参数
                    # 模型会自动下载到默认位置，然后我们需要移动到指定目录
                    try:
                        # 尝试使用新版本 API（如果支持）
                        models.download_models(model=model_name)
                        # 查找下载的模型文件并移动到指定目录
                        import shutil
                        # Path 已在文件顶部导入，无需重复导入
                        # 默认下载位置通常在用户目录下的 .celltypist 文件夹
                        default_cache = Path.home() / ".celltypist" / "models"
                        if default_cache.exists():
                            downloaded_model = default_cache / model_name
                            if downloaded_model.exists():
                                shutil.move(str(downloaded_model), str(model_path))
                                logger.info(f"✅ 模型已移动到: {model_path}")
                    except TypeError as e:
                        # 🔥 TASK 3 FIX: 如果新版本 API 不支持，尝试旧版本（可能使用不同的参数）
                        logger.debug(f"🔍 [CellAnnotation] 新版本API失败: {e}，尝试旧版本API")
                        try:
                            # 某些版本可能使用 path 而不是 folder
                            models.download_models(model=model_name, path=str(cache_path))
                            logger.info(f"✅ [CellAnnotation] 使用path参数下载模型成功")
                        except (TypeError, AttributeError) as e2:
                            # 🔥 TASK 3 FIX: 如果都不支持，尝试直接下载到默认位置，然后复制
                            logger.debug(f"🔍 [CellAnnotation] path参数也失败: {e2}，尝试默认位置下载")
                            try:
                                models.download_models(model=model_name)
                                # 尝试从默认位置复制到目标位置
                                default_cache = Path.home() / ".celltypist" / "models"
                                if default_cache.exists():
                                    downloaded_model = default_cache / model_name
                                    if downloaded_model.exists():
                                        import shutil
                                        shutil.copy2(str(downloaded_model), str(model_path))
                                        logger.info(f"✅ [CellAnnotation] 模型已从默认位置复制到: {model_path}")
                                    else:
                                        logger.warning(f"⚠️ [CellAnnotation] 模型已下载到默认位置，但未找到: {downloaded_model}")
                                else:
                                    logger.warning(f"⚠️ [CellAnnotation] 模型已下载到默认位置，请手动移动到: {cache_path}")
                            except Exception as e3:
                                logger.error(f"❌ [CellAnnotation] 所有下载方法都失败: {e3}")
                                # 不报错退出：使用 scanpy 聚类作为注释结果，保证步骤成功
                                logger.warning("CellTypist 模型不可用，改用基于聚类的注释。")
                                return _fallback_marker_annotation(adata, adata_path, output_dir, sc)
                
                # 若下载后仍无本地文件，使用 marker 规则注释
                if not model_path.exists():
                    logger.warning(
                        "Cell type annotation reference not found. Using marker-based fallback. path=%s",
                        model_path,
                    )
                    return _fallback_marker_annotation(adata, adata_path, output_dir, sc)
                
                # 加载模型（捕获 FileNotFoundError：路径存在检查后仍可能缺失，如权限/挂载问题）
                try:
                    model = celltypist.models.Model.load(str(model_path))
                except (FileNotFoundError, OSError) as e:
                    logger.warning("Reference model file not readable or missing: %s. Using marker-based fallback.", e)
                    return _fallback_marker_annotation(adata, adata_path, output_dir, sc)
                logger.info(f"✅ 模型加载成功: {model_name}")
                
                # 运行注释
                logger.info("🔬 正在运行 CellTypist 注释...")
                predictions = celltypist.annotate(
                    adata,
                    model=model,
                    majority_voting=True,
                    mode='probabilities'
                )
                
                # 保存预测结果
                adata.obs['predicted_labels'] = predictions.predicted_labels['majority_voting']
                if 'predicted_labels' in predictions.predicted_labels.columns:
                    adata.obs['predicted_labels_prob'] = predictions.predicted_labels['predicted_labels']
                
                # 统计注释结果
                label_counts = adata.obs['predicted_labels'].value_counts()
                n_cell_types = len(label_counts)
                
                # 生成 UMAP 图（按预测标签着色）
                plot_path = None
                if output_dir:
                    os.makedirs(output_dir, exist_ok=True)
                    
                    if 'X_umap' in adata.obsm.keys():
                        timestamp = int(time.time())
                        plot_path = os.path.join(output_dir, f"umap_annotated_{timestamp}.png")
                        
                        fig, ax = plt.subplots(figsize=(10, 8))
                        sc.pl.umap(
                            adata,
                            color='predicted_labels',
                            ax=ax,
                            show=False,
                            title="UMAP: Cell Type Annotation",
                            legend_loc='right margin',
                            frameon=False,
                            legend_fontsize=8
                        )
                        plt.savefig(plot_path, bbox_inches='tight', dpi=300)
                        plt.close()
                
                # 保存结果
                output_h5ad = None
                if output_dir:
                    output_h5ad = os.path.join(output_dir, "annotated.h5ad")
                    adata.write(output_h5ad)
                
                return {
                    "status": "success",
                    "method": "celltypist",
                    "model": model_name,
                    "n_cell_types": n_cell_types,
                    "cell_types": label_counts.to_dict(),
                    "plot_path": plot_path,
                    "output_h5ad": output_h5ad,
                    "summary": f"细胞类型注释完成: 识别到 {n_cell_types} 种细胞类型"
                }
            
            except ImportError:
                return {
                    "status": "error",
                    "error": "celltypist not installed. Please install: pip install celltypist"
                }
        
        elif method == "marker":
            # Marker-based 注释（需要用户提供 marker 基因列表）
            return {
                "status": "error",
                "error": "Marker-based annotation requires marker gene lists. This feature is not yet implemented."
            }
        
        else:
            return {
                "status": "error",
                "error": f"Unknown annotation method: {method}"
            }
    
    except ImportError:
        return {
            "status": "error",
            "error": "scanpy not installed. Please install: pip install scanpy"
        }
    except Exception as e:
        logger.error(f"❌ 细胞类型注释失败: {e}", exc_info=True)
        return {
            "status": "error",
            "error": str(e)
        }

