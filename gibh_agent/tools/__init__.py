"""
工具模块 - 模块化插件系统

自动发现和加载所有工具定义。
"""
import pkgutil
import importlib
import os
import logging

logger = logging.getLogger(__name__)


def load_all_tools():
    """
    自动导入所有模块以触发工具注册装饰器
    
    递归遍历 gibh_agent/tools 目录，导入所有 .py 文件（排除 __init__.py）
    """
    package_path = os.path.dirname(__file__)
    package_name = __name__
    
    logger.info("🔍 开始自动发现工具模块...")
    logger.info(f"   包路径: {package_path}")
    logger.info(f"   包名称: {package_name}")
    
    loaded_count = 0
    failed_count = 0
    
    # 使用 pkgutil.walk_packages 递归遍历所有子包和模块
    for importer, module_name, is_pkg in pkgutil.walk_packages(
        [package_path], 
        prefix=package_name + "."
    ):
        # 跳过 __init__.py 文件
        if module_name.endswith(".__init__"):
            continue
        
        try:
            # 导入模块（这会触发 @registry.register 装饰器）
            importlib.import_module(module_name)
            loaded_count += 1
            logger.info(f"✅ 已加载工具模块: {module_name}")
        except ImportError as e:
            # ImportError 可能是正常的（缺少可选依赖），只记录警告
            logger.warning(f"⚠️ 导入模块失败（可能是缺少依赖）: {module_name} - {e}")
            failed_count += 1
        except Exception as e:
            # 其他错误需要记录
            logger.error(f"❌ 加载模块失败: {module_name} - {e}", exc_info=True)
            failed_count += 1
    
    logger.info(f"📊 工具模块加载完成: 成功 {loaded_count} 个, 失败 {failed_count} 个")
    
    return {
        "loaded": loaded_count,
        "failed": failed_count
    }


# 显式导入 7 大组学原子工具模块（保证被注册并纳入 ChromaDB）
def _import_atomic_omics_tools():
    """Import atomic tool modules for 7 omics modalities (Genomics, Transcriptomics, Epigenomics, Proteomics, Metabolomics, Spatial, Radiomics)."""
    modules = [
        "genomics_tools",
        "transcriptomics_tools",
        "epigenomics_tools",
        "proteomics_tools",
        "protein_tools",
        "bio_tools",
        "metabolomics_tools",
        "radiomics_tools",
    ]
    for name in modules:
        try:
            importlib.import_module(f".{name}", package=__name__)
            logger.info("✅ 已加载原子工具模块: %s", name)
        except ImportError as e:
            logger.warning("⚠️ 可选原子工具模块未加载（可能缺少依赖）: %s - %s", name, e)


# 自动加载所有工具（当模块被导入时）
# 注意：这会在导入时立即执行，确保工具被注册
try:
    load_all_tools()
    _import_atomic_omics_tools()
except Exception as e:
    logger.error(f"❌ 自动加载工具失败: {e}", exc_info=True)

# Step 2 (Additive only): Explicitly ensure modality packages and flat modules are loaded.
# Genomics, Epigenomics, Proteomics, Spatial, Radiomics. Do not remove existing imports above.
for _mod in [
    "genomics_tools",
    "epigenomics_tools",
    "proteomics_tools",
    "protein_tools",
    "bio_tools",
    "radiomics_tools",
    "dynamic_proxy",
]:
    try:
        importlib.import_module(f".{_mod}", package=__name__)
    except ImportError:
        pass
# Radiomics package (atomic: io, analysis, visualization)
try:
    importlib.import_module(".radiomics", package=__name__)
    logger.info("✅ 已加载 Radiomics 工具包 (io, analysis, visualization)")
except ImportError as e:
    logger.warning("⚠️ Radiomics 工具包未加载（可能缺少依赖）: %s", e)

