"""
Scanpy 分析工具
参考旧版本实现，直接执行单细胞转录组分析流程
支持十步标准流程
"""
# 🔧 修复：在容器环境中设置 Numba 缓存目录（避免权限问题）
import os
# 如果 NUMBA_CACHE_DIR 未设置，使用临时目录
if 'NUMBA_CACHE_DIR' not in os.environ:
    import tempfile
    cache_dir = tempfile.mkdtemp(prefix='numba_cache_')
    os.environ['NUMBA_CACHE_DIR'] = cache_dir
    os.makedirs(cache_dir, exist_ok=True)

# 在导入 scanpy 之前配置 Numba（避免缓存错误）
try:
    import numba
    # 设置缓存目录
    if 'NUMBA_CACHE_DIR' in os.environ:
        numba.config.CACHE_DIR = os.environ['NUMBA_CACHE_DIR']
    # 启用缓存调试（帮助诊断问题）
    numba.config.DEBUG_CACHE = 1
except (ImportError, AttributeError):
    # 如果 numba 未安装或没有该配置，忽略
    pass

import scanpy as sc
import os
import matplotlib
# 设置无头模式，防止服务器报错
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
import time
import warnings
from typing import Dict, Any, List, Optional
from pathlib import Path

# 忽略警告
warnings.filterwarnings("ignore")

# 全局绘图设置
sc.settings.verbosity = 3
sc.settings.set_figure_params(
    dpi=300,
    facecolor='white',
    frameon=True,
    vector_friendly=True
)


class ScanpyTool:
    """
    Scanpy 分析工具
    
    核心功能：直接执行单细胞转录组分析流程
    参考旧版本：/home/ubuntu/GIBH-AGENT/services/api/src/scrna_analysis.py
    """
    
    def __init__(self, config: Dict[str, Any] = None, cellranger_tool=None):
        """
        初始化 Scanpy 工具
        
        Args:
            config: 配置字典
            cellranger_tool: CellRangerTool 实例（可选，用于运行 Cell Ranger）
        """
        self.config = config or {}
        # 使用相对路径，避免权限问题
        default_output = os.path.join(os.getcwd(), "results")
        self.output_dir = self.config.get("output_dir", default_output)
        try:
            os.makedirs(self.output_dir, exist_ok=True)
        except PermissionError:
            # 如果权限不足，使用当前目录下的 results
            self.output_dir = os.path.join(os.getcwd(), "results")
            os.makedirs(self.output_dir, exist_ok=True)
        
        # Cell Ranger 工具（可选）
        self.cellranger_tool = cellranger_tool
        
        # 工具映射表：将 tool_id 映射到具体的处理函数
        self.tool_map = {
            "inspect_file": self.inspect_file,  # 数据检查工具
            "run_cellranger": self.run_cellranger,  # Cell Ranger 计数
            "convert_cellranger_to_h5ad": self.convert_cellranger_to_h5ad,  # 转换 Cell Ranger 输出
            "local_qc": self.step_qc,
            "local_normalize": self.step_normalize,
            "local_hvg": self.step_hvg,
            "local_scale": self.step_scale,
            "local_pca": self.step_pca,
            "local_neighbors": self.step_neighbors,
            "local_cluster": self.step_cluster,
            "local_umap": self.step_umap,
            "local_tsne": self.step_tsne,
            "local_markers": self.step_markers,
            "local_annotate": self.step_annotate
        }
    
    def _save_plot(self, name_prefix: str) -> str:
        """保存图片并返回文件路径（相对于 results 目录）"""
        timestamp = int(time.time())
        filename = f"{name_prefix}_{timestamp}.png"
        save_path = os.path.join(self.output_dir, filename)
        plt.savefig(save_path, bbox_inches='tight', dpi=300)
        plt.close()
        
        # 返回相对于 results 目录的路径
        # 如果 output_dir 是 results/run_xxx，返回 run_xxx/filename
        if "results" in self.output_dir:
            # 提取 run_xxx 部分
            parts = self.output_dir.split(os.sep)
            if "results" in parts:
                results_idx = parts.index("results")
                if results_idx + 1 < len(parts):
                    run_dir = parts[results_idx + 1]
                    return f"{run_dir}/{filename}".replace("\\", "/")
        
        # 如果无法提取，返回完整相对路径
        if os.path.isabs(save_path):
            # 尝试找到 results 目录
            current = save_path
            while current != os.path.dirname(current):
                if os.path.basename(current) == "results":
                    rel_path = os.path.relpath(save_path, current)
                    return rel_path.replace("\\", "/")
                current = os.path.dirname(current)
        
        # 最后返回文件名
        return filename
    
    # ================= 📦 数据加载 =================
    def load_data(self, data_input: str):
        """加载单细胞数据"""
        print(f"📂 Loading data from: {data_input}")
        
        # 检查是否是 FASTQ 目录（不应该直接加载）
        if os.path.isdir(data_input):
            fastq_files = [f for f in os.listdir(data_input) if f.endswith(('.fastq', '.fq', '.fastq.gz', '.fq.gz'))]
            if fastq_files:
                raise ValueError(
                    f"检测到 FASTQ 目录: {data_input}。"
                    f"FASTQ 文件需要先通过 Cell Ranger 处理，不能直接加载。"
                    f"请先运行 Cell Ranger count，然后转换输出为 .h5ad 格式。"
                )
        
        if os.path.isdir(data_input):
            try:
                adata = sc.read_10x_mtx(data_input, var_names='gene_symbols', cache=False)
            except FileNotFoundError:
                print("⚠️ read_10x_mtx failed, trying manual mtx load...")
                # 手动加载逻辑
                mtx_path = os.path.join(data_input, "matrix.mtx")
                if not os.path.exists(mtx_path):
                    mtx_path = os.path.join(data_input, "matrix.mtx.gz")
                if not os.path.exists(mtx_path):
                    raise FileNotFoundError(
                        f"无法在目录 {data_input} 中找到 matrix.mtx 文件。"
                        f"这可能是 FASTQ 目录，需要先运行 Cell Ranger。"
                    )
                adata = sc.read_mtx(mtx_path).T
                
                genes_path = os.path.join(data_input, "features.tsv")
                if not os.path.exists(genes_path):
                    genes_path = os.path.join(data_input, "genes.tsv")
                if not os.path.exists(genes_path):
                    raise FileNotFoundError(f"无法找到基因文件: {genes_path}")
                genes = pd.read_csv(genes_path, header=None, sep='\t')
                adata.var_names = genes[1].values
                adata.var['gene_ids'] = genes[0].values
                
                barcodes_path = os.path.join(data_input, "barcodes.tsv")
                if not os.path.exists(barcodes_path):
                    raise FileNotFoundError(f"无法找到 barcodes 文件: {barcodes_path}")
                barcodes = pd.read_csv(barcodes_path, header=None, sep='\t')
                adata.obs_names = barcodes[0].values
            
            adata.var_names_make_unique()
        elif data_input.endswith('.h5ad'):
            adata = sc.read_h5ad(data_input)
        else:
            adata = sc.read(data_input)
        return adata
    
    # ================= 🔍 数据检查工具 =================
    
    def inspect_file(self, file_path: str) -> Dict[str, Any]:
        """
        检查文件内容，返回数据摘要
        
        这是一个强制性的检查步骤，必须在执行任何分析之前调用。
        
        Args:
            file_path: 数据文件路径（.h5ad 文件或 10x 目录）
        
        Returns:
            包含数据摘要的字典：
            - n_obs: 细胞数量
            - n_vars: 基因数量
            - obs_keys: .obs 中的列名列表
            - var_keys: .var 中的列名列表
            - is_normalized: 是否已标准化（基于最大值猜测）
            - max_value: 数据最大值
            - min_value: 数据最小值
            - preview: .obs 的前5行预览
            - has_clusters: 是否已有聚类结果
            - has_umap: 是否已有 UMAP 坐标
        """
        try:
            # 高效加载：使用 backed='r' 模式只读取元数据，不加载全部数据到内存
            if file_path.endswith('.h5ad'):
                try:
                    # 尝试使用 backed 模式（只读模式，不加载全部数据）
                    adata = sc.read_h5ad(file_path, backed='r')
                except:
                    # 如果 backed 模式失败，使用普通模式
                    adata = sc.read_h5ad(file_path)
            elif os.path.isdir(file_path):
                # 10x 格式需要完整加载
                adata = self.load_data(file_path)
            else:
                adata = sc.read(file_path)
            
            # 提取基本信息
            n_obs = adata.n_obs
            n_vars = adata.n_vars
            obs_keys = list(adata.obs.columns) if hasattr(adata.obs, 'columns') else []
            var_keys = list(adata.var.columns) if hasattr(adata.var, 'columns') else []
            
            # 检查数据值范围（用于判断是否已标准化）
            # 只检查一个小样本以提高效率
            import numpy as np
            sample_size = min(1000, adata.n_obs * adata.n_vars)
            if sample_size > 0:
                # 随机采样检查
                if hasattr(adata.X, 'toarray'):
                    # 稀疏矩阵
                    sample_data = adata.X[:min(100, adata.n_obs), :min(100, adata.n_vars)]
                    if hasattr(sample_data, 'toarray'):
                        sample_data = sample_data.toarray()
                    else:
                        sample_data = np.array(sample_data)
                else:
                    sample_data = np.array(adata.X[:min(100, adata.n_obs), :min(100, adata.n_vars)])
                
                max_value = float(np.nanmax(sample_data)) if sample_data.size > 0 else 0.0
                min_value = float(np.nanmin(sample_data)) if sample_data.size > 0 else 0.0
            else:
                max_value = 0.0
                min_value = 0.0
            
            # 判断是否已标准化
            # 经验规则：如果最大值 < 20，可能是 log-transformed；如果最大值很大（>1000），可能是原始 counts
            is_normalized = max_value < 20 if max_value > 0 else False
            
            # 预览 .obs 的前5行
            preview = None
            if n_obs > 0:
                try:
                    preview_df = adata.obs.head(5)
                    preview = preview_df.to_dict('records') if hasattr(preview_df, 'to_dict') else str(preview_df)
                except:
                    preview = "无法生成预览"
            
            # 检查是否已有分析结果
            has_clusters = 'leiden' in adata.obs.columns or 'louvain' in adata.obs.columns
            has_umap = 'X_umap' in adata.obsm_keys() if hasattr(adata, 'obsm_keys') else False
            
            # 检查是否有 QC 指标
            has_qc_metrics = any(key in obs_keys for key in ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'])
            
            result = {
                "n_obs": n_obs,
                "n_vars": n_vars,
                "obs_keys": obs_keys,
                "var_keys": var_keys,
                "is_normalized": is_normalized,
                "max_value": max_value,
                "min_value": min_value,
                "preview": preview,
                "has_clusters": has_clusters,
                "has_umap": has_umap,
                "has_qc_metrics": has_qc_metrics,
                "file_path": file_path
            }
            
            return result
            
        except Exception as e:
            return {
                "error": str(e),
                "file_path": file_path
            }
    
    # ================= 🔧 原子化工具函数 =================
    
    def step_qc(self, adata, params: Dict[str, Any]):
        """步骤1: 质量控制"""
        adata.var['mt'] = adata.var_names.str.startswith(('MT-', 'mt-'))
        sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], inplace=True)
        
        # 绘图
        sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'], 
                     jitter=0.4, multi_panel=True, show=False)
        plot_path = self._save_plot("qc_violin")
        
        # 过滤
        min_genes = int(params.get('min_genes', 200))
        max_mt = float(params.get('max_mt', 20))
        sc.pp.filter_cells(adata, min_genes=min_genes)
        adata = adata[adata.obs.pct_counts_mt < max_mt, :]
        sc.pp.filter_genes(adata, min_cells=3)
        
        return {
            "summary": f"过滤后剩余 {adata.n_obs} 细胞",
            "plot": plot_path
        }
    
    def step_normalize(self, adata, params: Dict[str, Any]):
        """步骤2: 标准化"""
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
        return {"summary": "LogNormalize 完成"}
    
    def step_hvg(self, adata, params: Dict[str, Any]):
        """步骤3: 寻找高变基因"""
        n_top_genes = int(params.get('n_top_genes', 2000))
        sc.pp.highly_variable_genes(adata, n_top_genes=n_top_genes)
        sc.pl.highly_variable_genes(adata, show=False)
        plot_path = self._save_plot("hvg")
        # 过滤高变基因
        adata._inplace_subset_var(adata.var['highly_variable'])
        return {"summary": f"筛选 {n_top_genes} 高变基因", "plot": plot_path}
    
    def step_scale(self, adata, params: Dict[str, Any]):
        """步骤4: 数据缩放"""
        sc.pp.scale(adata, max_value=10)
        return {"summary": "数据缩放完成"}
    
    def step_pca(self, adata, params: Dict[str, Any]):
        """步骤5: PCA 降维"""
        sc.tl.pca(adata, svd_solver='arpack')
        sc.pl.pca_variance_ratio(adata, log=True, show=False)
        plot_path = self._save_plot("pca_variance")
        return {"summary": "PCA 降维完成", "plot": plot_path}
    
    def step_neighbors(self, adata, params: Dict[str, Any]):
        """步骤6: 计算邻居"""
        n_neighbors = int(params.get('n_neighbors', 10))
        n_pcs = int(params.get('n_pcs', 40))
        sc.pp.neighbors(adata, n_neighbors=n_neighbors, n_pcs=n_pcs)
        return {"summary": "邻接图构建完成"}
    
    def step_cluster(self, adata, params: Dict[str, Any]):
        """步骤7: Leiden 聚类"""
        resolution = float(params.get('resolution', 0.5))
        sc.tl.leiden(adata, resolution=resolution)
        n_clusters = len(adata.obs['leiden'].unique())
        return {"summary": f"Leiden 聚类 (Res={resolution}): {n_clusters} 簇"}
    
    def step_umap(self, adata, params: Dict[str, Any]):
        """步骤8: UMAP 可视化"""
        sc.tl.umap(adata)
        fig, ax = plt.subplots(figsize=(8, 6))
        sc.pl.umap(adata, color=['leiden'], ax=ax, show=False, 
                   title="UMAP", legend_loc='on data', frameon=False)
        plot_path = self._save_plot("final_umap")
        return {"summary": "UMAP 生成完毕", "plot": plot_path}
    
    def step_tsne(self, adata, params: Dict[str, Any]):
        """步骤9: t-SNE 可视化"""
        if adata.n_obs < 5000:
            sc.tl.tsne(adata)
            fig, ax = plt.subplots(figsize=(8, 6))
            sc.pl.tsne(adata, color=['leiden'], ax=ax, show=False, 
                       title="t-SNE", frameon=False)
            plot_path = self._save_plot("final_tsne")
            return {"summary": "t-SNE 生成完毕", "plot": plot_path}
        else:
            return {"summary": "细胞数过多，跳过 t-SNE"}
    
    def step_markers(self, adata, params: Dict[str, Any]):
        """步骤10: 寻找 Marker 基因"""
        method = params.get('method', 't-test')
        sc.tl.rank_genes_groups(adata, 'leiden', method=method)
        result = adata.uns['rank_genes_groups']
        groups = result['names'].dtype.names
        
        # 构建 Marker 基因表格
        markers_data = {}
        for group in groups:
            markers_data[f"{group}_names"] = result['names'][group][:5]
            markers_data[f"{group}_pvals"] = result['pvals'][group][:5]
        
        markers_df = pd.DataFrame(markers_data)
        return {
            "summary": "Marker 基因鉴定完成",
            "details": markers_df.to_html(classes="table table-sm", index=False)
        }
    
    def step_annotate(self, adata, params: Dict[str, Any]):
        """步骤11: 细胞类型注释 (CellTypist)"""
        try:
            import celltypist
            from celltypist import models
        except ImportError:
            return {
                "summary": "错误: 未安装 celltypist",
                "error": "请运行: pip install celltypist"
            }
        
        # 模型缓存目录
        cache_dir = Path(self.config.get("cache_dir", "test_data/cache"))
        cache_dir.mkdir(parents=True, exist_ok=True)
        model_name = "Immune_All_Low.pkl"
        model_path = cache_dir / model_name
        
        # 下载或加载模型
        try:
            if not model_path.exists():
                print(f"📥 正在下载 CellTypist 模型: {model_name}")
                models.download_models(model=model_name, folder=str(cache_dir))
            
            # 加载模型
            model = celltypist.models.Model.load(str(model_path))
            print(f"✅ 模型加载成功: {model_name}")
            
            # 运行注释
            print("🔬 正在运行 CellTypist 注释...")
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
            if 'X_umap' in adata.obsm.keys():
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
                plot_path = self._save_plot("umap_annotated")
                
                return {
                    "summary": f"细胞类型注释完成: 识别到 {n_cell_types} 种细胞类型",
                    "plot": plot_path,
                    "cell_types": label_counts.to_dict(),
                    "n_cell_types": n_cell_types
                }
            else:
                return {
                    "summary": f"细胞类型注释完成: 识别到 {n_cell_types} 种细胞类型（请先运行 UMAP）",
                    "cell_types": label_counts.to_dict(),
                    "n_cell_types": n_cell_types
                }
                
        except Exception as e:
            import traceback
            error_msg = f"CellTypist 注释失败: {str(e)}"
            print(f"❌ {error_msg}")
            print(traceback.format_exc())
            return {
                "summary": error_msg,
                "error": str(e)
            }
    
    # ================= 🔬 Cell Ranger 工具 =================
    
    def run_cellranger(
        self,
        fastq_dir: str,
        sample_id: str,
        output_dir: str,
        reference: Optional[str] = None,
        sample: Optional[str] = None,
        localcores: int = 8,
        localmem: int = 32,
        create_bam: bool = False,
        expect_cells: Optional[int] = None
    ) -> Dict[str, Any]:
        """
        运行 Cell Ranger count
        
        Args:
            fastq_dir: FASTQ 文件目录路径
            sample_id: 样本 ID
            output_dir: 输出目录路径
            reference: 参考基因组路径（可选）
            sample: 样本名称（可选）
            localcores: CPU 核心数
            localmem: 内存（GB）
            create_bam: 是否创建 BAM 文件
            expect_cells: 预期细胞数（可选）
        
        Returns:
            执行结果字典
        """
        if not self.cellranger_tool:
            return {
                "status": "error",
                "error": "CellRangerTool not initialized. Please provide cellranger_tool in ScanpyTool.__init__()",
                "output_dir": None,
                "matrix_dir": None
            }
        
        return self.cellranger_tool.run_count(
            fastq_dir=fastq_dir,
            sample_id=sample_id,
            output_dir=output_dir,
            reference=reference,
            sample=sample,
            localcores=localcores,
            localmem=localmem,
            create_bam=create_bam,
            expect_cells=expect_cells
        )
    
    def convert_cellranger_to_h5ad(
        self,
        cellranger_matrix_dir: str,
        output_h5ad_path: str
    ) -> Dict[str, Any]:
        """
        将 Cell Ranger 输出转换为 .h5ad 格式
        
        Args:
            cellranger_matrix_dir: Cell Ranger 矩阵目录路径（filtered_feature_bc_matrix）
            output_h5ad_path: 输出的 .h5ad 文件路径
        
        Returns:
            转换结果字典，包含：
            - status: "success" 或 "error"
            - output_path: 输出文件路径
            - n_obs: 细胞数
            - n_vars: 基因数
            - error: 错误信息（如果有）
        """
        try:
            print(f"📖 读取 Cell Ranger 输出: {cellranger_matrix_dir}")
            
            # 检查输入目录
            if not os.path.exists(cellranger_matrix_dir):
                return {
                    "status": "error",
                    "error": f"Cell Ranger matrix directory does not exist: {cellranger_matrix_dir}",
                    "output_path": None,
                    "n_obs": None,
                    "n_vars": None
                }
            
            # 读取 10x MTX 数据
            adata = sc.read_10x_mtx(
                cellranger_matrix_dir,
                var_names='gene_symbols',  # 使用基因符号作为变量名
                cache=True
            )
            
            # 确保基因名唯一
            adata.var_names_make_unique()
            
            # 保存为 .h5ad 格式
            print(f"💾 保存为 .h5ad 格式: {output_h5ad_path}")
            os.makedirs(os.path.dirname(output_h5ad_path), exist_ok=True)
            adata.write(output_h5ad_path)
            
            file_size_mb = os.path.getsize(output_h5ad_path) / (1024 * 1024)
            
            return {
                "status": "success",
                "output_path": output_h5ad_path,
                "n_obs": adata.n_obs,
                "n_vars": adata.n_vars,
                "matrix_type": type(adata.X).__name__,
                "file_size_mb": round(file_size_mb, 2)
            }
        except Exception as e:
            return {
                "status": "error",
                "error": f"Failed to convert Cell Ranger output: {str(e)}",
                "output_path": None,
                "n_obs": None,
                "n_vars": None
        }
    
    # ================= 🚀 主调度器 =================
    def run_pipeline(
        self,
        data_input: str,
        steps_config: Optional[List[Dict[str, Any]]] = None
    ) -> Dict[str, Any]:
        """
        执行完整的 Scanpy 工作流
        
        Args:
            data_input: 输入数据路径（.h5ad 文件或 10x 目录）
            steps_config: 步骤配置列表，每个步骤包含 tool_id 和 params
        
        Returns:
            分析报告字典
        """
        report = {
            "status": "running",
            "steps_details": [],
            "final_plot": None,
            "qc_metrics": {},
            "diagnosis": "",
            "error": None
        }
        
        try:
            # 1. 加载数据
            adata = self.load_data(data_input)
            report["qc_metrics"]["raw_cells"] = adata.n_obs
            report["qc_metrics"]["raw_genes"] = adata.n_vars
            
            if not steps_config:
                print("⚠️ No steps provided, returning raw data stats.")
                report["status"] = "success"
                return report
            
            # 2. 按顺序执行步骤
            print(f"📋 Pipeline Plan: {[s['tool_id'] for s in steps_config]}")
            
            for step in steps_config:
                tool_id = step['tool_id']
                params = step.get('params', {})
                
                print(f"▶️ Executing: {tool_id}")
                
                # 从映射表中找到函数并执行
                if tool_id in self.tool_map:
                    func = self.tool_map[tool_id]
                    result = func(adata, params)
                    
                    # 构造标准返回格式
                    step_report = {
                        "name": tool_id,
                        "status": "success",
                        "plot": result.get("plot"),
                        "details": result.get("details", ""),
                        "summary": result.get("summary", "完成")
                    }
                    
                    # 特殊处理：更新全局报告状态
                    if tool_id == "local_qc":
                        report["qc_metrics"]["filtered_cells"] = adata.n_obs
                        report["qc_metrics"]["filtered_genes"] = adata.n_vars
                    if tool_id == "local_umap":
                        report["final_plot"] = result.get("plot")
                    
                    report["steps_details"].append(step_report)
                else:
                    print(f"⚠️ Unknown tool_id: {tool_id}, skipping...")
            
            # 3. 保存处理后的数据
            output_file = os.path.join(self.output_dir, "processed.h5ad")
            adata.write(output_file)
            report["output_file"] = output_file
            
            # 4. 生成诊断
            report["diagnosis"] = f"""
            ### ✅ 分析完成
            - **执行步骤数**: {len(report['steps_details'])}
            - **剩余细胞**: {adata.n_obs}
            - **剩余基因**: {adata.n_vars}
            """
            
            report["status"] = "success"
            return report
        
        except Exception as e:
            print(f"❌ Pipeline Error: {e}")
            import traceback
            traceback.print_exc()
            report["status"] = "failed"
            report["error"] = str(e)
            return report
    
    def generate_workflow_script(
        self,
        input_path: str,
        output_dir: str,
        steps: List[Dict[str, Any]],
        **kwargs
    ) -> str:
        """
        生成完整的 Scanpy 工作流脚本（保留此方法以兼容旧接口）
        
        Args:
            input_path: 输入文件路径（.h5ad 或 10x 目录）
            output_dir: 输出目录
            steps: 步骤列表，每个步骤包含 tool_id 和 params
            **kwargs: 其他参数
        
        Returns:
            Python 脚本内容
        """
        script_lines = [
            "#!/usr/bin/env python3",
            "# Scanpy Workflow Script",
            "# Generated by GIBH-Agent",
            "",
            "import scanpy as sc",
            "import pandas as pd",
            "import numpy as np",
            "import matplotlib",
            "matplotlib.use('Agg')",
            "import matplotlib.pyplot as plt",
            "from pathlib import Path",
            "import json",
            "import warnings",
            "warnings.filterwarnings('ignore')",
            "",
            "# 配置",
            f"input_path = '{input_path}'",
            f"output_dir = Path('{output_dir}')",
            "output_dir.mkdir(parents=True, exist_ok=True)",
            "",
            "# 加载数据",
            "print('Loading data...')",
            "if input_path.endswith('.h5ad'):",
            "    adata = sc.read_h5ad(input_path)",
            "elif os.path.isdir(input_path):",
            "    adata = sc.read_10x_mtx(input_path, var_names='gene_symbols', cache=False)",
            "else:",
            "    adata = sc.read(input_path)",
            "adata.var_names_make_unique()",
            "",
        ]
        
        # 生成每个步骤的代码
        for step in steps:
            tool_id = step.get("tool_id", "")
            params = step.get("params", {})
            
            step_code = self._generate_step_code(tool_id, params)
            script_lines.extend([
                f"# {step.get('name', tool_id)}",
                step_code,
                ""
            ])
        
        # 保存结果
        script_lines.extend([
            "# 保存结果",
            "output_file = output_dir / 'processed.h5ad'",
            "adata.write(output_file)",
            "print(f'Results saved to: {output_file}')",
            "",
            "print('Workflow completed successfully!')"
        ])
        
        return "\n".join(script_lines)
    
    def _generate_step_code(self, tool_id: str, params: Dict[str, Any]) -> str:
        """生成单个步骤的代码（用于脚本生成）"""
        code_map = {
            "local_qc": self._qc_code(params),
            "local_normalize": self._normalize_code(),
            "local_hvg": self._hvg_code(params),
            "local_scale": self._scale_code(),
            "local_pca": self._pca_code(),
            "local_neighbors": self._neighbors_code(params),
            "local_cluster": self._cluster_code(params),
            "local_umap": self._umap_code(),
            "local_tsne": self._tsne_code(),
            "local_markers": self._markers_code(params)
        }
        
        return code_map.get(tool_id, f"# Unknown tool: {tool_id}")
    
    def _qc_code(self, params: Dict[str, Any]) -> str:
        """QC 步骤代码"""
        min_genes = params.get("min_genes", 200)
        max_mt = params.get("max_mt", 20)
        
        return f"""# 计算 QC 指标
adata.var['mt'] = adata.var_names.str.startswith(('MT-', 'mt-'))
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], inplace=True)

# 过滤
sc.pp.filter_cells(adata, min_genes={min_genes})
adata = adata[adata.obs.pct_counts_mt < {max_mt}, :]
sc.pp.filter_genes(adata, min_cells=3)"""
    
    def _normalize_code(self) -> str:
        """标准化步骤代码"""
        return """# 标准化
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)"""
    
    def _hvg_code(self, params: Dict[str, Any]) -> str:
        """高变基因步骤代码"""
        n_top_genes = params.get("n_top_genes", 2000)
        return f"""# 寻找高变基因
sc.pp.highly_variable_genes(adata, n_top_genes={n_top_genes})
adata._inplace_subset_var(adata.var['highly_variable'])"""
    
    def _scale_code(self) -> str:
        """缩放步骤代码"""
        return """# 缩放
sc.pp.scale(adata, max_value=10)"""
    
    def _pca_code(self) -> str:
        """PCA 步骤代码"""
        return """# PCA
sc.tl.pca(adata, svd_solver='arpack')"""
    
    def _neighbors_code(self, params: Dict[str, Any]) -> str:
        """Neighbors 步骤代码"""
        n_neighbors = params.get("n_neighbors", 10)
        n_pcs = params.get("n_pcs", 40)
        return f"""# 计算邻居
sc.pp.neighbors(adata, n_neighbors={n_neighbors}, n_pcs={n_pcs})"""
    
    def _cluster_code(self, params: Dict[str, Any]) -> str:
        """聚类步骤代码"""
        resolution = params.get("resolution", 0.5)
        return f"""# Leiden 聚类
sc.tl.leiden(adata, resolution={resolution})"""
    
    def _umap_code(self) -> str:
        """UMAP 步骤代码"""
        return """# UMAP
sc.tl.umap(adata)
fig, ax = plt.subplots(figsize=(8, 6))
sc.pl.umap(adata, color=['leiden'], ax=ax, show=False, title="UMAP", legend_loc='on data', frameon=False)
plt.savefig(output_dir / 'umap.png', bbox_inches='tight', dpi=300)
plt.close()"""
    
    def _tsne_code(self) -> str:
        """t-SNE 步骤代码"""
        return """# t-SNE
if adata.n_obs < 5000:
    sc.tl.tsne(adata)
    fig, ax = plt.subplots(figsize=(8, 6))
    sc.pl.tsne(adata, color=['leiden'], ax=ax, show=False, title="t-SNE", frameon=False)
    plt.savefig(output_dir / 'tsne.png', bbox_inches='tight', dpi=300)
    plt.close()"""
    
    def _markers_code(self, params: Dict[str, Any]) -> str:
        """Markers 步骤代码"""
        method = params.get("method", "t-test")
        return f"""# 寻找 Marker 基因
sc.tl.rank_genes_groups(adata, 'leiden', method='{method}')
result = adata.uns['rank_genes_groups']
groups = result['names'].dtype.names
markers_df = pd.DataFrame({{group + '_' + key: result[key][group][:5] 
                            for group in groups for key in ['names', 'pvals']}})
markers_df.to_csv(output_dir / 'markers.csv', index=False)"""
