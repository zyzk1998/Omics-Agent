"""
单细胞数据导出工具
"""
import os
import zipfile
import logging
from typing import Dict, Any, Optional, List
from pathlib import Path

from ...core.tool_registry import registry

logger = logging.getLogger(__name__)


@registry.register(
    name="rna_export_results",
    description="Exports single-cell RNA-seq analysis results in multiple formats: H5AD file, CSV tables (cell metadata, gene metadata, marker genes), and a ZIP archive containing all figures.",
    category="scRNA-seq",
    output_type="mixed"
)
def export_results(
    adata_path: str,
    output_dir: str,
    export_h5ad: bool = True,
    export_csv: bool = True,
    export_figures: bool = True,
    figure_dir: Optional[str] = None
) -> Dict[str, Any]:
    """
    导出分析结果
    
    Args:
        adata_path: AnnData 文件路径（.h5ad）
        output_dir: 输出目录
        export_h5ad: 是否导出 H5AD 文件
        export_csv: 是否导出 CSV 表格
        export_figures: 是否导出图片（ZIP）
        figure_dir: 图片目录（如果为 None，则从 output_dir 推断）
    
    Returns:
        导出结果字典，包含所有导出文件的路径
    """
    try:
        import scanpy as sc
        import pandas as pd
        
        # 加载数据
        adata = sc.read_h5ad(adata_path)
        
        # 确保输出目录存在
        os.makedirs(output_dir, exist_ok=True)
        
        exported_files = []
        
        # 1. 导出 H5AD 文件
        if export_h5ad:
            h5ad_path = os.path.join(output_dir, "analysis_results.h5ad")
            adata.write(h5ad_path)
            exported_files.append(h5ad_path)
            logger.info(f"✅ 已导出 H5AD: {h5ad_path}")
        
        # 2. 导出 CSV 表格
        if export_csv:
            # 导出细胞元数据
            cell_metadata_path = os.path.join(output_dir, "cell_metadata.csv")
            adata.obs.to_csv(cell_metadata_path, index=True)
            exported_files.append(cell_metadata_path)
            logger.info(f"✅ 已导出细胞元数据: {cell_metadata_path}")
            
            # 导出基因元数据
            gene_metadata_path = os.path.join(output_dir, "gene_metadata.csv")
            adata.var.to_csv(gene_metadata_path, index=True)
            exported_files.append(gene_metadata_path)
            logger.info(f"✅ 已导出基因元数据: {gene_metadata_path}")
            
            # 导出 Marker 基因（如果存在）
            if 'rank_genes_groups' in adata.uns:
                marker_path = os.path.join(output_dir, "marker_genes.csv")
                marker_df = pd.DataFrame(adata.uns['rank_genes_groups'])
                marker_df.to_csv(marker_path, index=True)
                exported_files.append(marker_path)
                logger.info(f"✅ 已导出 Marker 基因: {marker_path}")
        
        # 3. 导出图片（ZIP）
        zip_path = None
        if export_figures:
            if figure_dir is None:
                # 尝试从 output_dir 推断图片目录
                figure_dir = os.path.join(output_dir, "figures")
            
            if os.path.exists(figure_dir):
                zip_path = os.path.join(output_dir, "figures.zip")
                with zipfile.ZipFile(zip_path, 'w', zipfile.ZIP_DEFLATED) as zipf:
                    for root, dirs, files in os.walk(figure_dir):
                        for file in files:
                            if file.endswith(('.png', '.jpg', '.jpeg', '.pdf', '.svg')):
                                file_path = os.path.join(root, file)
                                arcname = os.path.relpath(file_path, figure_dir)
                                zipf.write(file_path, arcname)
                
                exported_files.append(zip_path)
                logger.info(f"✅ 已导出图片 ZIP: {zip_path}")
            else:
                logger.warning(f"⚠️ 图片目录不存在: {figure_dir}")
        
        return {
            "status": "success",
            "output_dir": output_dir,
            "exported_files": exported_files,
            "h5ad_path": os.path.join(output_dir, "analysis_results.h5ad") if export_h5ad else None,
            "cell_metadata_path": os.path.join(output_dir, "cell_metadata.csv") if export_csv else None,
            "gene_metadata_path": os.path.join(output_dir, "gene_metadata.csv") if export_csv else None,
            "zip_path": zip_path,
            "summary": f"已导出 {len(exported_files)} 个文件到 {output_dir}"
        }
    
    except ImportError as e:
        missing_module = "scanpy" if "scanpy" in str(e) else "pandas"
        return {
            "status": "error",
            "error": f"{missing_module} not installed. Please install: pip install {missing_module}"
        }
    except Exception as e:
        logger.error(f"❌ 导出失败: {e}", exc_info=True)
        return {
            "status": "error",
            "error": str(e)
        }

