"""
上游数据处理工具 - Cell Ranger
"""
import os
import sys
import subprocess
import shutil
import logging
from typing import Dict, Any, Optional
from pathlib import Path

from ...core.tool_registry import registry

logger = logging.getLogger(__name__)


@registry.register(
    name="rna_cellranger_count",
    description="Runs Cell Ranger count on FASTQ files to generate gene expression count matrices. This is the upstream processing step for 10x Genomics single-cell RNA-seq data. Executes asynchronously in the background to avoid blocking the UI.",
    category="scRNA-seq",
    output_type="json"
)
def run_cellranger_count(
    fastqs_path: str,
    sample_id: str,
    transcriptome_path: str,
    output_dir: str,
    localcores: int = 8,
    localmem: int = 32,
    create_bam: bool = False,
    expect_cells: Optional[int] = None,
    cellranger_path: str = "/opt/cellranger"
) -> Dict[str, Any]:
    """
    异步执行 Cell Ranger count（后台任务）
    
    Args:
        fastqs_path: FASTQ 文件目录路径（本地绝对路径）
        sample_id: 样本 ID（也是输出目录名）
        transcriptome_path: 参考转录组路径
        output_dir: 最终输出目录路径
        localcores: 本地核心数
        localmem: 本地内存（GB）
        create_bam: 是否创建 BAM 文件（Cell Ranger 10.0.0+ 必需参数）
        expect_cells: 预期细胞数（可选）
        cellranger_path: Cell Ranger 安装路径
    
    Returns:
        异步作业状态字典，包含：
        - status: "async_job_started"
        - job_id: 进程 ID（格式：PID_12345）
        - message: 用户友好的消息
        - output_dir: 输出目录路径
        - log_path: 日志文件路径
        - error: 错误信息（如果验证失败）
    """
    try:
        # Step 1: 验证路径存在
        if not os.path.exists(fastqs_path):
            return {
                "status": "error",
                "error": f"FASTQ 目录不存在: {fastqs_path}",
                "message": f"请检查 FASTQ 文件路径是否正确: {fastqs_path}"
            }
        
        if not os.path.exists(transcriptome_path):
            return {
                "status": "error",
                "error": f"参考转录组路径不存在: {transcriptome_path}",
                "message": f"请检查参考转录组路径是否正确: {transcriptome_path}"
            }
        
        if not os.path.exists(cellranger_path):
            return {
                "status": "error",
                "error": f"Cell Ranger 路径不存在: {cellranger_path}",
                "message": f"请检查 Cell Ranger 安装路径: {cellranger_path}"
            }
        
        # Step 2: 准备输出目录和日志文件
        os.makedirs(output_dir, exist_ok=True)
        log_path = os.path.join(output_dir, f"cellranger_{sample_id}.log")
        
        # Step 3: 构建 Cell Ranger 命令
        cellranger_cmd = os.path.join(cellranger_path, "cellranger")
        cmd = [
            cellranger_cmd,
            "count",
            "--id", sample_id,
            "--fastqs", fastqs_path,
            "--transcriptome", transcriptome_path,
            "--localcores", str(localcores),
            "--localmem", str(localmem)
        ]
        
        if create_bam:
            cmd.append("--include-introns")
        
        if expect_cells:
            cmd.extend(["--expect-cells", str(expect_cells)])
        
        # Step 4: 启动异步进程（Fire-and-Forget）
        logger.info(f"🚀 启动异步 Cell Ranger count 任务: {sample_id}")
        logger.info(f"   命令: {' '.join(cmd)}")
        logger.info(f"   输出目录: {output_dir}")
        logger.info(f"   日志文件: {log_path}")
        
        # 打开日志文件（追加模式）
        log_file = open(log_path, 'a')
        
        # 启动进程（不等待完成）
        process = subprocess.Popen(
            cmd,
            stdout=log_file,
            stderr=subprocess.STDOUT,  # 将 stderr 重定向到 stdout
            cwd=output_dir,  # 在输出目录中执行
            start_new_session=True  # 创建新的会话，使进程独立
        )
        
        # 获取进程 ID
        job_id = f"PID_{process.pid}"
        
        # 不关闭 log_file，让它在进程运行期间保持打开
        # 注意：在实际应用中，可能需要一个后台任务管理器来跟踪这些进程
        
        logger.info(f"✅ Cell Ranger 任务已启动: {job_id}")
        logger.info(f"   进程 PID: {process.pid}")
        logger.info(f"   日志文件: {log_path}")
        
        return {
            "status": "async_job_started",
            "job_id": job_id,
            "pid": process.pid,
            "message": f"Cell Ranger 任务已后台启动（{job_id}），请耐心等待。\n日志文件: {log_path}\n输出目录: {output_dir}",
            "output_dir": output_dir,
            "log_path": log_path,
            "sample_id": sample_id,
            "summary": f"Cell Ranger count 任务已启动，进程 ID: {process.pid}"
        }
    
    except Exception as e:
        logger.error(f"❌ 启动 Cell Ranger 任务失败: {e}", exc_info=True)
        return {
            "status": "error",
            "error": str(e),
            "message": f"启动 Cell Ranger 任务失败: {str(e)}"
        }


@registry.register(
    name="rna_cellranger_mkfastq",
    description="Runs Cell Ranger mkfastq to demultiplex FASTQ files from Illumina sequencer output. This is typically the first step in 10x Genomics workflows.",
    category="scRNA-seq",
    output_type="json"
)
def run_cellranger_mkfastq(
    run_dir: str,
    sample_sheet: str,
    output_dir: str,
    cellranger_path: str = "/opt/cellranger"
) -> Dict[str, Any]:
    """
    执行 Cell Ranger mkfastq（解复用）
    
    Args:
        run_dir: Illumina 运行目录
        sample_sheet: 样本表 CSV 文件路径
        output_dir: 输出目录路径
        cellranger_path: Cell Ranger 安装路径
    
    Returns:
        执行结果字典
    """
    # 注意：mkfastq 功能在原始 CellRangerTool 中未实现
    # 这里提供一个占位符实现
    return {
        "status": "error",
        "error": "Cell Ranger mkfastq is not yet implemented. Please use Cell Ranger count directly with FASTQ files."
    }


@registry.register(
    name="rna_convert_cellranger_to_h5ad",
    description="Converts Cell Ranger output (filtered_feature_bc_matrix) to Scanpy-compatible H5AD format for downstream analysis.",
    category="scRNA-seq",
    output_type="json"
)
def convert_cellranger_to_h5ad(
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
        import scanpy as sc
        
        logger.info(f"📖 读取 Cell Ranger 输出: {cellranger_matrix_dir}")
        
        # 检查输入目录
        if not os.path.exists(cellranger_matrix_dir):
            return {
                "status": "error",
                "error": f"Cell Ranger matrix directory does not exist: {cellranger_matrix_dir}",
                "output_path": None,
                "n_obs": None,
                "n_vars": None
            }
        
        # 🔥 使用统一的10x数据读取函数，支持压缩和未压缩格式
        from ...core.rna_utils import read_10x_data
        adata = read_10x_data(
            cellranger_matrix_dir,
            var_names='gene_symbols',  # 使用基因符号作为变量名
            cache=True
        )
        
        # 保存为 .h5ad 格式
        logger.info(f"💾 保存为 .h5ad 格式: {output_h5ad_path}")
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
    except ImportError:
        return {
            "status": "error",
            "error": "scanpy not installed. Please install: pip install scanpy",
            "output_path": None,
            "n_obs": None,
            "n_vars": None
        }
    except Exception as e:
        logger.error(f"❌ 转换失败: {e}", exc_info=True)
        return {
            "status": "error",
            "error": f"Failed to convert Cell Ranger output: {str(e)}",
            "output_path": None,
            "n_obs": None,
            "n_vars": None
        }

