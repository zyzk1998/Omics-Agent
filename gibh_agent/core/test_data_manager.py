"""
测试数据管理器
自动检测和管理可用的测试数据集
"""
import os
from pathlib import Path
from typing import Dict, List, Optional, Any
import json


class TestDataManager:
    """测试数据管理器"""
    
    def __init__(self, test_data_dir: Optional[str] = None):
        """
        初始化测试数据管理器
        
        Args:
            test_data_dir: 测试数据目录路径，默认为项目根目录下的 test_data
        """
        if test_data_dir is None:
            # 默认使用项目根目录下的 test_data
            current_file = Path(__file__)
            project_root = current_file.parent.parent.parent
            test_data_dir = str(project_root / "test_data")
        
        self.test_data_dir = Path(test_data_dir)
        self.test_data_dir.mkdir(parents=True, exist_ok=True)
    
    def scan_test_datasets(self) -> List[Dict[str, Any]]:
        """
        扫描可用的测试数据集
        
        Returns:
            测试数据集列表，每个数据集包含：
            - id: 数据集ID
            - name: 数据集名称
            - description: 描述
            - fastq_dir: FASTQ 目录路径（如果存在）
            - reference: 参考基因组路径（如果存在）
            - h5ad_file: .h5ad 文件路径（如果存在）
            - available: 是否可用
            - workflow_type: 适用工作流类型（"full_workflow" 或 "downstream"）
        """
        datasets = []
        
        # 扫描 test_data 目录
        if not self.test_data_dir.exists():
            return datasets
        
        # 查找 pbmc_1k_v3 数据集（FASTQ格式，适用于全流程）
        pbmc_fastq = self.test_data_dir / "pbmc_1k_v3_fastqs"
        pbmc_ref = self.test_data_dir / "refdata-gex-GRCh38-2024-A"
        pbmc_h5ad = self.test_data_dir / "pbmc_1k_v3_filtered.h5ad"
        pbmc_output = self.test_data_dir / "pbmc_1k_v3_output"
        
        # 1. FASTQ数据集（适用于Cell Ranger全流程）
        if pbmc_fastq.exists() and pbmc_ref.exists():
            datasets.append({
                "id": "pbmc_1k_v3_fastqs",
                "name": "PBMC 1k v3 FASTQ",
                "description": "10x Genomics PBMC 1k cells v3 chemistry dataset (FASTQ files)",
                "fastq_dir": str(pbmc_fastq),
                "reference": str(pbmc_ref),
                "h5ad_file": str(pbmc_h5ad) if pbmc_h5ad.exists() else None,
                "available": True,
                "cells": "~1,200",
                "chemistry": "v3",
                "workflow_type": "full_workflow",  # 适用于全流程（包含Cell Ranger）
                "data_type": "fastq"
            })
        
        # 2. h5ad数据集（适用于下游分析）
        if pbmc_h5ad.exists() and pbmc_h5ad.is_file():
            datasets.append({
                "id": "pbmc_1k_v3_h5ad",
                "name": "PBMC 1k v3 (h5ad)",
                "description": "10x Genomics PBMC 1k cells v3 chemistry dataset (processed h5ad file)",
                "fastq_dir": None,
                "reference": None,
                "h5ad_file": str(pbmc_h5ad),
                "available": True,
                "cells": "~1,200",
                "chemistry": "v3",
                "workflow_type": "downstream",  # 适用于下游分析（质控、标准化、聚类等）
                "data_type": "h5ad"
            })
        
        # 3. Cell Ranger输出（适用于10x格式转换）
        if pbmc_output.exists() and (pbmc_output / "outs").exists():
            outs_dir = pbmc_output / "outs"
            if (outs_dir / "filtered_feature_bc_matrix.h5").exists() or (outs_dir / "filtered_feature_bc_matrix").exists():
                datasets.append({
                    "id": "pbmc_1k_v3_output",
                    "name": "PBMC 1k v3 (Cell Ranger Output)",
                    "description": "10x Genomics PBMC 1k cells v3 chemistry dataset (Cell Ranger output)",
                    "fastq_dir": None,
                    "reference": None,
                    "h5ad_file": None,
                    "cellranger_output": str(pbmc_output),
                    "available": True,
                    "cells": "~1,200",
                    "chemistry": "v3",
                    "workflow_type": "downstream",  # 适用于10x格式转换和下游分析
                    "data_type": "cellranger_output"
                })
        
        return datasets
    
    def get_datasets_by_workflow_type(self, workflow_steps: List[str]) -> List[Dict[str, Any]]:
        """
        根据工作流步骤推荐对应的测试数据集
        
        Args:
            workflow_steps: 工作流步骤列表（工具ID或步骤名称）
        
        Returns:
            推荐的测试数据集列表
        """
        all_datasets = self.scan_test_datasets()
        if not all_datasets:
            return []
        
        # 检查是否需要Cell Ranger（全流程）
        needs_cellranger = any(
            step and ('cellranger' in step.lower() or 'cell_ranger' in step.lower())
            for step in workflow_steps
        )
        
        # 检查是否只需要下游分析（h5ad或10x格式）
        needs_downstream = any(
            step and (step.lower() in ['qc', 'quality_control', 'normalize', 'normalization', 
                                      'clustering', 'annotation', 'umap', 'tsne', 'pca'])
            for step in workflow_steps
        ) and not needs_cellranger
        
        # 根据需求筛选数据集
        recommended = []
        if needs_cellranger:
            # 全流程：推荐FASTQ数据集
            for dataset in all_datasets:
                if dataset.get('workflow_type') == 'full_workflow' or dataset.get('data_type') == 'fastq':
                    recommended.append(dataset)
        elif needs_downstream:
            # 下游分析：推荐h5ad或Cell Ranger输出
            for dataset in all_datasets:
                if dataset.get('workflow_type') == 'downstream' or dataset.get('data_type') in ['h5ad', 'cellranger_output']:
                    recommended.append(dataset)
        else:
            # 默认：返回所有可用数据集
            recommended = all_datasets
        
        return recommended if recommended else all_datasets
    
    def get_dataset_by_id(self, dataset_id: str) -> Optional[Dict[str, Any]]:
        """
        根据ID获取数据集信息
        
        Args:
            dataset_id: 数据集ID
        
        Returns:
            数据集信息字典，如果不存在则返回 None
        """
        datasets = self.scan_test_datasets()
        for dataset in datasets:
            if dataset["id"] == dataset_id:
                return dataset
        return None
    
    def format_datasets_for_display(self, datasets: List[Dict[str, Any]]) -> str:
        """
        格式化数据集列表用于显示
        
        Args:
            datasets: 数据集列表
        
        Returns:
            格式化的字符串
        """
        if not datasets:
            return "暂无可用的测试数据集。"
        
        lines = ["可用的测试数据集：\n"]
        for i, dataset in enumerate(datasets, 1):
            lines.append(f"{i}. **{dataset['name']}** (ID: {dataset['id']})")
            lines.append(f"   - 描述: {dataset['description']}")
            if dataset.get('cells'):
                lines.append(f"   - 细胞数: {dataset['cells']}")
            if dataset.get('chemistry'):
                lines.append(f"   - 化学版本: {dataset['chemistry']}")
            if dataset.get('fastq_dir'):
                lines.append(f"   - FASTQ: {dataset['fastq_dir']}")
            if dataset.get('reference'):
                lines.append(f"   - 参考基因组: {dataset['reference']}")
            if dataset.get('h5ad_file'):
                lines.append(f"   - 已转换: {dataset['h5ad_file']}")
            lines.append("")
        
        return "\n".join(lines)
    
    def format_datasets_for_selection(self, datasets: List[Dict[str, Any]]) -> str:
        """
        格式化数据集列表用于用户选择（JSON格式）
        
        Args:
            datasets: 数据集列表
        
        Returns:
            JSON格式的字符串
        """
        # 简化数据集信息，只保留必要的字段
        simplified = []
        for dataset in datasets:
            simplified.append({
                "id": dataset["id"],
                "name": dataset["name"],
                "description": dataset["description"],
                "available": dataset["available"]
            })
        return json.dumps(simplified, ensure_ascii=False, indent=2)


def create_test_data_manager(test_data_dir: Optional[str] = None) -> TestDataManager:
    """创建测试数据管理器实例"""
    return TestDataManager(test_data_dir)

