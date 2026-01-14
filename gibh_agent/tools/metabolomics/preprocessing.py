"""
代谢组学数据预处理工具
"""
import logging
import os
from typing import Dict, Any, Optional
from pathlib import Path
import pandas as pd
import numpy as np
from sklearn.preprocessing import StandardScaler

from ...core.tool_registry import registry

logger = logging.getLogger(__name__)


@registry.register(
    name="preprocess_data",
    description="Preprocesses metabolite data: handles missing values, applies log2 transformation, and standardizes the data. Returns preprocessed DataFrame and saves to CSV file.",
    category="Metabolomics",
    output_type="json"
)
def preprocess_metabolite_data(
    file_path: str,
    missing_imputation: str = "min",
    log_transform: bool = True,
    standardize: bool = True,
    output_dir: Optional[str] = None
) -> Dict[str, Any]:
    """
    预处理代谢物数据
    
    Args:
        file_path: 输入数据文件路径（CSV）
        missing_imputation: 缺失值填充方法（"min", "median", "mean", "zero"）
        log_transform: 是否进行 log2 转换（默认 True）
        standardize: 是否标准化（默认 True）
        output_dir: 输出目录（如果提供，将保存预处理后的数据）
    
    Returns:
        包含以下键的字典:
        - status: "success" 或 "error"
        - preprocessed_data: 预处理后的数据（JSON 格式）
        - output_path: 保存的文件路径（如果保存）
        - output_file: 保存的文件路径（别名，用于数据流传递）
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
        
        # 4. 保存预处理后的数据到文件（用于数据流传递）
        output_path = None
        if output_dir:
            # 确保输出目录存在
            output_path_obj = Path(output_dir)
            output_path_obj.mkdir(parents=True, exist_ok=True)
            
            # 生成输出文件路径
            input_filename = Path(file_path).stem
            output_path = str(output_path_obj / f"{input_filename}_preprocessed.csv")
            
            # 保存数据
            data.to_csv(output_path)
            logger.info(f"💾 预处理后的数据已保存: {output_path}")
        else:
            # 如果没有指定输出目录，尝试使用输入文件所在目录
            input_dir = Path(file_path).parent
            output_path = str(input_dir / "preprocessed_data.csv")
            data.to_csv(output_path)
            logger.info(f"💾 预处理后的数据已保存: {output_path}")
        
        return {
            "status": "success",
            "preprocessed_data": data.to_dict(orient='index'),
            "output_path": output_path,
            "output_file": output_path,  # 别名，用于数据流传递
            "file_path": output_path,  # 另一个别名，确保兼容性
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

