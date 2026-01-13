"""
文件检查工具
"""
import logging
from typing import Dict, Any
from pathlib import Path
import pandas as pd

from ...core.tool_registry import registry

logger = logging.getLogger(__name__)


@registry.register(
    name="file_inspect",
    description="Inspects a file and returns basic metadata: file size, type, number of rows/columns (for tabular files), and a preview of the first few rows.",
    category="General",
    output_type="json"
)
def inspect_file(
    file_path: str,
    preview_rows: int = 5
) -> Dict[str, Any]:
    """
    检查文件基本信息
    
    Args:
        file_path: 文件路径
        preview_rows: 预览行数（默认 5）
    
    Returns:
        包含文件元数据的字典
    """
    try:
        path = Path(file_path)
        
        if not path.exists():
            return {
                "status": "error",
                "error": f"文件不存在: {file_path}"
            }
        
        metadata = {
            "status": "success",
            "file_path": str(path),
            "file_size": path.stat().st_size,
            "file_type": path.suffix,
            "exists": True
        }
        
        # 如果是 CSV 文件，读取预览
        if path.suffix.lower() == '.csv':
            df = pd.read_csv(file_path, nrows=preview_rows)
            metadata.update({
                "rows": len(df),
                "columns": list(df.columns),
                "preview": df.head(preview_rows).to_dict(orient='records')
            })
        
        return metadata
    
    except Exception as e:
        logger.error(f"❌ 文件检查失败: {e}", exc_info=True)
        return {
            "status": "error",
            "error": str(e)
        }

