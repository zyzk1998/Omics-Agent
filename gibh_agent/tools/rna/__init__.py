"""
单细胞转录组（scRNA-seq）工具模块

自动导入所有RNA工具以触发注册装饰器
"""
# 导入所有RNA工具模块以触发 @registry.register 装饰器
from . import quality_control
from . import analysis
from . import plotting
from . import annotation
from . import export
from . import upstream

__all__ = [
    "quality_control",
    "analysis",
    "plotting",
    "annotation",
    "export",
    "upstream"
]

