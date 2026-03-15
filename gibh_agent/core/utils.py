"""
通用工具函数
包含 JSON 序列化辅助函数、绘图路径矫正等。
"""
from datetime import date, datetime
from decimal import Decimal
from pathlib import Path
from typing import Union

import numpy as np
import pandas as pd
import math
import logging

logger = logging.getLogger(__name__)


def sanitize_for_json(obj):
    """
    递归清理对象，使其可以被 JSON 序列化
    
    处理：
    - Numpy 数据类型（int64, float32, ndarray 等）
    - NaN 和 Infinity 值
    - Pandas DataFrame 和 Series
    
    Args:
        obj: 需要清理的对象（可以是 dict, list, numpy 类型等）
    
    Returns:
        清理后的对象，可以被 JSON 序列化
    """
    # 处理字典
    if isinstance(obj, dict):
        return {k: sanitize_for_json(v) for k, v in obj.items()}
    
    # 处理列表
    elif isinstance(obj, list):
        return [sanitize_for_json(v) for v in obj]
    
    # 处理 Numpy 整数类型
    elif isinstance(obj, (np.integer, np.intc, np.intp, np.int8,
                         np.int16, np.int32, np.int64, np.uint8,
                         np.uint16, np.uint32, np.uint64)):
        return int(obj)
    
    # 处理 Numpy 浮点类型
    elif isinstance(obj, (np.floating, np.float16, np.float32, np.float64)):
        # JSON 不支持 NaN 和 Infinity，转换为 None
        if np.isnan(obj) or np.isinf(obj):
            logger.warning(f"⚠️ 检测到 NaN/Infinity 值，转换为 None: {obj}")
            return None
        return float(obj)
    
    # 🔥 处理 Numpy 布尔类型（numpy.bool_）
    elif isinstance(obj, np.bool_):
        return bool(obj)
    
    # 处理 Numpy 数组
    elif isinstance(obj, np.ndarray):
        return sanitize_for_json(obj.tolist())
    
    # 处理 Pandas DataFrame
    elif isinstance(obj, pd.DataFrame):
        # 转换为字典列表（records 格式）
        return sanitize_for_json(obj.to_dict(orient='records'))
    
    # 处理 Pandas Series
    elif isinstance(obj, pd.Series):
        return sanitize_for_json(obj.to_dict())
    
    # 处理 Python 内置的 float（可能包含 NaN/Inf）
    elif isinstance(obj, float):
        if math.isnan(obj) or math.isinf(obj):
            logger.warning(f"⚠️ 检测到 Python float NaN/Infinity 值，转换为 None: {obj}")
            return None
        return obj

    # 处理 datetime/date（API 返回 JSON 时需为字符串）
    elif isinstance(obj, (datetime, date)):
        return obj.isoformat() if hasattr(obj, "isoformat") else str(obj)

    # 处理 bytes（转为 UTF-8 字符串，避免 JSON 报错）
    elif isinstance(obj, bytes):
        try:
            return obj.decode("utf-8")
        except Exception:
            return f"<bytes len={len(obj)}>"

    # 处理 Decimal（MySQL/JSON 有时返回）
    elif isinstance(obj, Decimal):
        return float(obj)

    # 其他类型直接返回（str, int, bool, None 等）
    return obj


# 绘图工具允许的扩展名（matplotlib savefig 支持）；禁止 .csv 等导致 Format not supported
_PLOT_EXTENSIONS = (".png", ".pdf", ".svg", ".jpg", ".jpeg")


def sanitize_plot_path(path: Union[str, Path]) -> Path:
    """
    强制矫正绘图输出路径后缀，避免外部传入 .csv 等导致 matplotlib 报错。
    仅当后缀在允许列表内时保留；否则改为 .png。不影响合法 .pdf/.svg 等高精度输出。
    """
    p = Path(path) if isinstance(path, str) else path
    if not p.suffix or p.suffix.lower() not in _PLOT_EXTENSIONS:
        return (p.parent / (p.stem + ".png")) if p.stem else (p.parent / "plot.png")
    return p

