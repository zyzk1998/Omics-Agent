"""
通用工具函数
包含 JSON 序列化辅助函数、绘图路径矫正等。
"""
from datetime import date, datetime
from decimal import Decimal
from enum import Enum
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
    if isinstance(obj, dict):
        return {k: sanitize_for_json(v) for k, v in obj.items()}
    elif isinstance(obj, list):
        return [sanitize_for_json(v) for v in obj]
    elif isinstance(obj, (np.integer, np.intc, np.intp, np.int8,
                         np.int16, np.int32, np.int64, np.uint8,
                         np.uint16, np.uint32, np.uint64)):
        return int(obj)
    elif isinstance(obj, (np.floating, np.float16, np.float32, np.float64)):
        if np.isnan(obj) or np.isinf(obj):
            logger.warning("⚠️ 检测到 NaN/Infinity 值，转换为 None: %s", obj)
            return None
        return float(obj)
    elif isinstance(obj, np.bool_):
        return bool(obj)
    elif isinstance(obj, np.ndarray):
        return sanitize_for_json(obj.tolist())
    elif isinstance(obj, pd.DataFrame):
        return sanitize_for_json(obj.to_dict(orient='records'))
    elif isinstance(obj, pd.Series):
        return sanitize_for_json(obj.to_dict())
    elif isinstance(obj, float):
        if math.isnan(obj) or math.isinf(obj):
            logger.warning("⚠️ 检测到 Python float NaN/Infinity 值，转换为 None: %s", obj)
            return None
        return obj
    elif isinstance(obj, (datetime, date)):
        return obj.isoformat() if hasattr(obj, "isoformat") else str(obj)
    elif isinstance(obj, bytes):
        try:
            return obj.decode("utf-8")
        except Exception:
            return f"<bytes len={len(obj)}>"
    elif isinstance(obj, Decimal):
        return float(obj)
    elif isinstance(obj, (tuple, set, frozenset)):
        return [sanitize_for_json(v) for v in obj]
    elif isinstance(obj, Path):
        return str(obj)
    elif isinstance(obj, Enum):
        return sanitize_for_json(obj.value)
    elif isinstance(obj, (str, int, bool)) or obj is None:
        return obj
    try:
        s = str(obj)
        if len(s) > 8000:
            return s[:8000] + f"...<truncated, type={type(obj).__name__}>"
        return s
    except Exception:  # pragma: no cover
        return f"<non-serializable type={type(obj).__name__}>"


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
