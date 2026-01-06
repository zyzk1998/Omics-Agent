"""
通用工具函数
包含 JSON 序列化辅助函数
"""
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
    
    # 其他类型直接返回（str, int, bool, None 等）
    return obj

