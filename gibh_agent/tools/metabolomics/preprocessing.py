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
        
        # 🔥 CRITICAL FIX: 保留分组列（非数值列或低基数数值列）
        non_numeric_cols = df.select_dtypes(exclude=[np.number]).columns.tolist()
        numeric_cols_all = df.select_dtypes(include=[np.number]).columns.tolist()
        group_cols = []
        id_like = ("patient", "sample", "subject", "specimen", "id")
        for col in non_numeric_cols:
            unique_count = df[col].nunique()
            cn = str(col).lower().replace(" ", "").replace("_", "")
            if any(tok in cn for tok in id_like) and unique_count > max(10, len(df) // 3):
                continue
            if 2 <= unique_count <= 100:
                group_cols.append(col)
        for col in numeric_cols_all:
            unique_count = df[col].nunique()
            if 2 <= unique_count <= 5:
                group_cols.append(col)
        
        # 提取数值列（代谢物列）
        numeric_cols = df.select_dtypes(include=[np.number]).columns.tolist()
        
        # 🔥 CRITICAL FIX: 保留分组列和数值列
        columns_to_keep = group_cols + numeric_cols
        data = df[columns_to_keep].copy()
        
        logger.info(f"📊 [Preprocess] 保留分组列: {group_cols}")
        logger.info(f"📊 [Preprocess] 处理数值列: {len(numeric_cols)} 个代谢物")
        
        # 1. 缺失值填充
        if missing_imputation == "min":
            data = data.fillna(data.min())
        elif missing_imputation == "median":
            data = data.fillna(data.median())
        elif missing_imputation == "mean":
            data = data.fillna(data.mean())
        else:  # "zero"
            data = data.fillna(0)
        
        # 🔥 CRITICAL FIX: 只对数值列进行转换和标准化
        # 分离分组列和数值列
        numeric_data = data[numeric_cols].copy()
        group_data = data[group_cols].copy() if group_cols else pd.DataFrame()
        
        # 2. Log2 转换（仅对数值列）
        if log_transform:
            numeric_data = numeric_data.apply(lambda x: np.log2(x + 1))
        
        # 3. 标准化（仅对数值列）
        if standardize:
            scaler = StandardScaler()
            numeric_scaled = scaler.fit_transform(numeric_data)
            numeric_data = pd.DataFrame(
                numeric_scaled,
                index=numeric_data.index,
                columns=numeric_data.columns
            )
        
        # 🔥 CRITICAL FIX: 合并分组列和数值列
        if group_cols:
            data = pd.concat([group_data, numeric_data], axis=1)
            # 确保列顺序：分组列在前，数值列在后
            data = data[group_cols + numeric_cols]
        else:
            data = numeric_data
        
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
        
        # 🔥 TASK 3: Ensure absolute path is returned
        # Convert to absolute path to ensure data flow works correctly
        output_path_absolute = os.path.abspath(output_path)
        if output_path_absolute != output_path:
            logger.info(f"🔄 [Preprocess] 转换为绝对路径: {output_path} -> {output_path_absolute}")
        
        return {
            "status": "success",
            "preprocessed_data": data.to_dict(orient='index'),
            "output_path": output_path_absolute,  # 🔥 CRITICAL: Return absolute path
            "output_file": output_path_absolute,  # 别名，用于数据流传递
            "file_path": output_path_absolute,  # 另一个别名，确保兼容性
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

