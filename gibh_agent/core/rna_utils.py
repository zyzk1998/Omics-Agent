"""
单细胞RNA-seq数据读取工具函数
提供统一的10x Genomics数据读取接口，支持压缩和未压缩格式
"""
import os
import logging
from pathlib import Path
from typing import Optional, Tuple

logger = logging.getLogger(__name__)


def read_10x_data(
    data_path: str,
    var_names: str = 'gene_symbols',
    cache: bool = False
):
    """
    读取10x Genomics数据，自动检测压缩/未压缩格式
    
    这个函数提供了统一的接口来读取10x数据，支持：
    - 压缩格式（.gz）
    - 未压缩格式
    - 自动检测和多重尝试策略
    - 手动读取作为后备方案
    
    Args:
        data_path: 10x数据目录路径
        var_names: 变量名类型，'gene_symbols' 或 'gene_ids'
        cache: 是否使用缓存
    
    Returns:
        AnnData对象
    
    Raises:
        FileNotFoundError: 如果找不到必需的文件
        Exception: 如果所有读取方法都失败
    """
    import scanpy as sc
    import pandas as pd
    
    data_path = Path(data_path)
    
    if not data_path.is_dir():
        raise ValueError(f"路径不是目录: {data_path}")
    
    # 获取目录内容
    dir_contents = os.listdir(data_path)
    
    # 检测文件是压缩还是未压缩的
    is_compressed = any(f.endswith('.gz') for f in dir_contents if 'matrix.mtx' in f)
    logger.info(f"📦 检测到10x文件格式: {'压缩' if is_compressed else '未压缩'}")
    
    # 策略1: 使用 scanpy 的标准方法，根据文件格式设置 compressed 参数
    try:
        adata = sc.read_10x_mtx(
            str(data_path),
            var_names=var_names,
            cache=cache,
            compressed=is_compressed
        )
        adata.var_names_make_unique()
        logger.info(f"✅ 成功使用 sc.read_10x_mtx (compressed={is_compressed}): {adata.n_obs} cells, {adata.n_vars} genes")
        return adata
    except Exception as e:
        logger.warning(f"⚠️ read_10x_mtx failed (compressed={is_compressed}): {e}, trying opposite compression setting...")
        last_error = e
    
    # 策略2: 尝试相反的压缩设置
    try:
        adata = sc.read_10x_mtx(
            str(data_path),
            var_names=var_names,
            cache=cache,
            compressed=not is_compressed
        )
        adata.var_names_make_unique()
        logger.info(f"✅ 成功使用相反的压缩设置 (compressed={not is_compressed}): {adata.n_obs} cells, {adata.n_vars} genes")
        return adata
    except Exception as e2:
        logger.warning(f"⚠️ 相反的压缩设置也失败: {e2}, trying without var_names...")
        last_error = e2
    
    # 策略3: 尝试不使用 var_names 参数
    try:
        adata = sc.read_10x_mtx(
            str(data_path),
            cache=cache,
            compressed=is_compressed
        )
        adata.var_names_make_unique()
        logger.info(f"✅ 成功不使用 var_names (compressed={is_compressed}): {adata.n_obs} cells, {adata.n_vars} genes")
        return adata
    except Exception as e3:
        logger.warning(f"⚠️ 所有 scanpy 方法都失败: {e3}, trying manual load...")
        last_error = e3
    
    # 策略4: 手动读取作为后备方案
    try:
        # 查找文件路径（支持压缩和非压缩）
        mtx_path = None
        for f in ['matrix.mtx', 'matrix.mtx.gz']:
            candidate = data_path / f
            if candidate.exists():
                mtx_path = candidate
                break
        
        if not mtx_path:
            raise FileNotFoundError("无法找到 matrix.mtx 文件")
        
        barcodes_path = None
        for f in ['barcodes.tsv', 'barcodes.tsv.gz']:
            candidate = data_path / f
            if candidate.exists():
                barcodes_path = candidate
                break
        
        if not barcodes_path:
            raise FileNotFoundError("无法找到 barcodes.tsv 文件")
        
        features_path = None
        for f in ['features.tsv', 'features.tsv.gz', 'genes.tsv', 'genes.tsv.gz']:
            candidate = data_path / f
            if candidate.exists():
                features_path = candidate
                break
        
        if not features_path:
            raise FileNotFoundError("无法找到 features.tsv 或 genes.tsv 文件")
        
        logger.info(f"📖 手动读取10x文件:")
        logger.info(f"   Matrix: {mtx_path}")
        logger.info(f"   Barcodes: {barcodes_path}")
        logger.info(f"   Features: {features_path}")
        
        # 读取矩阵文件
        adata = sc.read_mtx(str(mtx_path)).T
        
        # 读取基因/特征文件
        if str(features_path).endswith('.gz'):
            import gzip
            with gzip.open(features_path, 'rt') as f:
                features_df = pd.read_csv(f, header=None, sep='\t')
        else:
            features_df = pd.read_csv(features_path, header=None, sep='\t')
        
        # 设置基因名（features.tsv 格式：gene_id, gene_symbol, feature_type）
        # 或 genes.tsv 格式：gene_id, gene_symbol
        if len(features_df.columns) >= 2:
            if var_names == 'gene_symbols':
                adata.var_names = features_df[1].values  # 使用基因符号
            else:
                adata.var_names = features_df[0].values  # 使用基因ID
            if len(features_df.columns) >= 1:
                adata.var['gene_ids'] = features_df[0].values
        else:
            adata.var_names = features_df[0].values
        
        # 读取 barcodes 文件
        if str(barcodes_path).endswith('.gz'):
            import gzip
            with gzip.open(barcodes_path, 'rt') as f:
                barcodes_df = pd.read_csv(f, header=None, sep='\t')
        else:
            barcodes_df = pd.read_csv(barcodes_path, header=None, sep='\t')
        
        adata.obs_names = barcodes_df[0].values
        
        # 确保基因名唯一
        adata.var_names_make_unique()
        
        logger.info(f"✅ 手动读取成功: {adata.n_obs} cells, {adata.n_vars} genes")
        return adata
        
    except Exception as e4:
        logger.error(f"❌ 手动读取也失败: {e4}")
        raise Exception(f"所有10x数据读取方法都失败。最后错误: {str(e4)}。原始错误: {str(last_error)}")

