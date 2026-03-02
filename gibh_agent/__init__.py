"""
GIBH-AGENT-V2 主模块
延迟导入 main，避免仅使用 tools（如 radiomics）时加载 scanpy 等重依赖。
"""
__all__ = ["GIBHAgent", "create_agent"]


def __getattr__(name: str):
    if name in ("GIBHAgent", "create_agent"):
        from .main import GIBHAgent, create_agent
        return GIBHAgent if name == "GIBHAgent" else create_agent
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")

