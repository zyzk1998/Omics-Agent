"""
工作流注册表和基础工作流类

提供：
- BaseWorkflow: 抽象工作流基类
- WorkflowRegistry: 工作流注册表（单例）
- MetabolomicsWorkflow: 代谢组学工作流
- RNAWorkflow: scRNA-seq 工作流
"""

from .base import BaseWorkflow
from .registry import WorkflowRegistry
from .metabolomics import MetabolomicsWorkflow
from .rna import RNAWorkflow

__all__ = [
    "BaseWorkflow",
    "WorkflowRegistry",
    "MetabolomicsWorkflow",
    "RNAWorkflow",
]

