"""
Radiomics 子包：本地 PyRadiomics 原子工具（可选）。
重计算默认在 worker-pyskills（见顶层 radiomics_tools.py TaaS 适配器）；主 API 镜像可不安装 pyradiomics。
"""
import logging

logger = logging.getLogger(__name__)

try:
    import radiomics  # noqa: F401
except ImportError:
    logger.info("主进程未安装 PyRadiomics：跳过 radiomics 子包注册（影像 TaaS 见 radiomics_tools）。")
    __all__: list = []
else:
    from . import io
    from . import analysis
    from . import visualization
    from . import modeling

    __all__ = ["io", "analysis", "visualization", "modeling"]
