"""
STED-EC 单细胞时空轨迹推断工作流（JiekaiLab 特色科研流程）。

四步骤 DAG：数据与元数据校验 -> 时间序列标准化 -> moscot 轨迹推断 -> 轨迹可视化。
"""
from typing import Dict, Any, List
from .base import BaseWorkflow


class STEDECWorkflow(BaseWorkflow):
    """
    STED-EC 工作流：单细胞时间序列最优传输轨迹推断。

    步骤顺序：
    1. sted_ec_data_validation - 数据与元数据校验
    2. sted_ec_time_series_formatting - 时间序列标准化
    3. sted_ec_moscot_trajectory - 最优传输轨迹推断
    4. sted_ec_plot_trajectory - 轨迹 UMAP/流形图可视化
    """

    def get_name(self) -> str:
        """工作流名称，用于注册表与快车道路由。"""
        return "STED_EC"

    def get_description(self) -> str:
        """工作流描述（商业化去敏，与前端技能广场一致）。"""
        return "基于最优传输理论 (Optimal Transport) 的单细胞时空轨迹推断工具。通过整合单细胞转录组与时间序列信息，精准重构细胞发育、分化或疾病演进的动态轨迹，揭示细胞命运决定的关键调控机制。"

    def get_steps_dag(self) -> Dict[str, List[str]]:
        """四步骤 DAG：顺序串联。"""
        return {
            "sted_ec_data_validation": [],
            "sted_ec_time_series_formatting": ["sted_ec_data_validation"],
            "sted_ec_moscot_trajectory": ["sted_ec_time_series_formatting"],
            "sted_ec_plot_trajectory": ["sted_ec_moscot_trajectory"],
        }

    def get_step_metadata(self, step_id: str) -> Dict[str, Any]:
        """步骤元数据：名称、描述、tool_id、默认参数。"""
        metadata_map = {
            "sted_ec_data_validation": {
                "name": "数据与元数据校验",
                "description": "加载 h5ad，校验 time_key 与 cell_type_key 是否存在并解析实际列名。",
                "tool_id": "sted_ec_data_validation",
                "default_params": {"h5ad_path": "<user_input>"},
            },
            "sted_ec_time_series_formatting": {
                "name": "时间序列标准化",
                "description": "对上一步校验后的 h5ad 执行时间序列清洗与格式转换，输出标准化 h5ad。",
                "tool_id": "sted_ec_time_series_formatting",
                "default_params": {
                    "h5ad_path": "<sted_ec_data_validation>",
                    "time_key": "day",
                    "cell_type_key": None,
                },
            },
            "sted_ec_moscot_trajectory": {
                "name": "最优传输轨迹推断",
                "description": "基于 moscot 对单细胞时间序列数据执行最优传输轨迹推断。",
                "tool_id": "sted_ec_moscot_trajectory",
                "default_params": {
                    "h5ad_path": "<sted_ec_time_series_formatting>",
                    "time_key": "day",
                    "cell_type_key": None,
                },
            },
            "sted_ec_plot_trajectory": {
                "name": "轨迹可视化",
                "description": "根据推断结果绘制轨迹 UMAP 或流形图。",
                "tool_id": "sted_ec_plot_trajectory",
                "default_params": {
                    "trajectory_data_path": "<sted_ec_moscot_trajectory>",
                    "output_plot_path": "<output_dir>/sted_ec_trajectory.png",
                    "plot_type": "umap",
                },
            },
        }
        if step_id not in metadata_map:
            raise ValueError(f"未知的步骤ID: {step_id}")
        return metadata_map[step_id]
