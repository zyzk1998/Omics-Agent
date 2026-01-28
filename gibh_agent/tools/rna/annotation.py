"""
单细胞数据注释工具 - 细胞类型注释
"""
import os
import time
import logging
from typing import Dict, Any, Optional
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from ...core.tool_registry import registry

logger = logging.getLogger(__name__)


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
        
        # 加载数据
        adata = sc.read_h5ad(adata_path)
        
        if method == "celltypist":
            try:
                import celltypist
                from celltypist import models
                
                # 设置缓存目录
                if cache_dir is None:
                    cache_dir = os.path.join(os.getcwd(), "test_data", "cache")
                cache_path = Path(cache_dir)
                cache_path.mkdir(parents=True, exist_ok=True)
                
                model_path = cache_path / model_name
                
                # 下载或加载模型
                if not model_path.exists():
                    logger.info(f"📥 正在下载 CellTypist 模型: {model_name}")
                    # 🔥 FIX: celltypist 的 download_models 使用 model 参数，不使用 folder 参数
                    # 模型会自动下载到默认位置，然后我们需要移动到指定目录
                    try:
                        # 尝试使用新版本 API（如果支持）
                        models.download_models(model=model_name)
                        # 查找下载的模型文件并移动到指定目录
                        import shutil
                        from pathlib import Path
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
                                return {
                                    "status": "error",
                                    "error": f"无法下载CellTypist模型: {str(e3)}。请检查网络连接或手动下载模型。"
                                }
                
                # 加载模型
                model = celltypist.models.Model.load(str(model_path))
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

