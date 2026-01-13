"""
代谢组学数据预处理工具
"""
import logging
from typing import Dict, Any
import pandas as pd
import numpy as np
from sklearn.preprocessing import StandardScaler

from ...core.tool_registry import registry

logger = logging.getLogger(__name__)


@registry.register(
    name="metabolomics_preprocess",
    description="Preprocesses metabolite data: handles missing values, applies log2 transformation, and standardizes the data. Returns preprocessed DataFrame.",
    category="Metabolomics",
    output_type="json"
)
def preprocess_metabolite_data(
    file_path: str,
    missing_imputation: str = "min",
    log_transform: bool = True,
    standardize: bool = True
) -> Dict[str, Any]:
    """
    预处理代谢物数据
    
    Args:
        file_path: 输入数据文件路径（CSV）
        missing_imputation: 缺失值填充方法（"min", "median", "mean", "zero"）
        log_transform: 是否进行 log2 转换（默认 True）
        standardize: 是否标准化（默认 True）
    
    Returns:
        包含以下键的字典:
        - status: "success" 或 "error"
        - preprocessed_data: 预处理后的数据（JSON 格式）
        - output_path: 保存的文件路径（如果保存）
        - error: 错误信息（如果失败）
    """
    try:
        # 读取数据
        df = pd.read_csv(file_path, index_col=0)
        
        # 提取数值列
        numeric_cols = df.select_dtypes(include=[np.number]).columns
        data = df[numeric_cols].copy()
        
        # 1. 缺失值填充
        if missing_imputation == "min":
            data = data.fillna(data.min())
        elif missing_imputation == "median":
            data = data.fillna(data.median())
        elif missing_imputation == "mean":
            data = data.fillna(data.mean())
        else:  # "zero"
            data = data.fillna(0)
        
        # 2. Log2 转换
        if log_transform:
            data = data.apply(lambda x: np.log2(x + 1))
        
        # 3. 标准化
        if standardize:
            scaler = StandardScaler()
            data_scaled = scaler.fit_transform(data)
            data = pd.DataFrame(
                data_scaled,
                index=data.index,
                columns=data.columns
            )
        
        return {
            "status": "success",
            "preprocessed_data": data.to_dict(orient='index'),
            "shape": {
                "rows": len(data),
                "columns": len(data.columns)
            }
        }
    
    except Exception as e:
        logger.error(f"❌ 数据预处理失败: {e}", exc_info=True)
        return {
            "status": "error",
            "error": str(e)
        }

