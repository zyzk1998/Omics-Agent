"""
STED-EC 双通道工作流（物理隔离）

- **STED_EC**（通道 A / 演示保底）：四步 DAG — 校验 → 时序标准化 → moscot → 可视化。
- **SPATIOTEMPORAL_DYNAMICS**（通道 B / 完全体）：六步工具 DAG（至通路富集）；专家解读由 Orchestrator + BaseAgent 统一生成。
"""
from typing import Dict, Any, List
from .base import BaseWorkflow


class STEDECWorkflow(BaseWorkflow):
    """
    通道 A：STED_EC 基础单细胞轨迹推断（四步，技能广场 / 领导演示绑定）。

    步骤顺序：
    1. sted_ec_data_validation - 数据与元数据校验
    2. sted_ec_time_series_formatting - 时间序列标准化
    3. sted_ec_moscot_trajectory - 最优传输轨迹推断
    4. sted_ec_plot_trajectory - 时空动力学图可视化（写出含 UMAP 的 h5ad）
    """

    def get_name(self) -> str:
        return "STED_EC"

    def get_display_title(self) -> str:
        return "STED_EC 分析流程"

    def get_description(self) -> str:
        return (
            "基于最优传输（moscot）的 **STED_EC 基础分析流程**（四步）："
            "数据校验、时间序列标准化、轨迹推断与结果可视化。"
            "适合演示与快速轨迹推断；完整驱动基因 / 通路 / 专家报告请使用自然语言「时空动力学」等触发完全体。"
        )

    def get_steps_dag(self) -> Dict[str, List[str]]:
        return {
            "sted_ec_data_validation": [],
            "sted_ec_time_series_formatting": ["sted_ec_data_validation"],
            "sted_ec_moscot_trajectory": ["sted_ec_time_series_formatting"],
            "sted_ec_plot_trajectory": ["sted_ec_moscot_trajectory"],
        }

    def get_step_metadata(self, step_id: str) -> Dict[str, Any]:
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
                "name": "时空动力学图可视化",
                "description": "根据推断结果绘制时空 UMAP、细胞类型演化等图表（与第三步「最优传输轨迹推断」相区分：本步仅作图与汇总）。",
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


class SpatiotemporalDynamicsWorkflow(BaseWorkflow):
    """
    通道 B：单细胞时空动力学分析完全体（六步；专家解读由 Orchestrator + BaseAgent 统一 Reporting 生成，非 DAG 节点）。
    """

    def get_name(self) -> str:
        return "SPATIOTEMPORAL_DYNAMICS"

    def get_display_title(self) -> str:
        return "单细胞时空动力学分析"

    def get_description(self) -> str:
        return (
            "单细胞时空动力学 **完全体** 流水线（工具六步）：校验 → 标准化 → moscot 推断 → 可视化 → "
            "驱动基因挖掘 → 通路富集（gseapy Enrichr）。"
            "执行结束后由系统与 RNA/代谢组相同的 **AnalysisSummary / diagnosis** 机制生成专家 Markdown。"
        )

    def get_steps_dag(self) -> Dict[str, List[str]]:
        return {
            "sted_ec_data_validation": [],
            "sted_ec_time_series_formatting": ["sted_ec_data_validation"],
            "sted_ec_moscot_trajectory": ["sted_ec_time_series_formatting"],
            "sted_ec_plot_trajectory": ["sted_ec_moscot_trajectory"],
            "sted_ec_driver_gene_extraction": ["sted_ec_plot_trajectory"],
            "sted_ec_pathway_enrichment": ["sted_ec_driver_gene_extraction"],
        }

    def get_step_metadata(self, step_id: str) -> Dict[str, Any]:
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
                "name": "时空动力学图可视化",
                "description": "根据推断结果绘制时空 UMAP、细胞类型演化等图表（与第三步「最优传输轨迹推断」相区分：本步仅作图与汇总）。",
                "tool_id": "sted_ec_plot_trajectory",
                "default_params": {
                    "trajectory_data_path": "<sted_ec_moscot_trajectory>",
                    "output_plot_path": "<output_dir>/sted_ec_trajectory.png",
                    "plot_type": "umap",
                },
            },
            "sted_ec_driver_gene_extraction": {
                "name": "轨迹驱动基因挖掘",
                "description": "使用 Scanpy rank_genes_groups (wilcoxon) 按时间点或细胞类型分组提取候选驱动基因，输出 driver_genes.csv 供 GSEA 等下游使用。",
                "tool_id": "sted_ec_driver_gene_extraction",
                "default_params": {
                    "h5ad_path": "<sted_ec_plot_trajectory>",
                    "top_n": 100,
                    "groupby_mode": "auto",
                },
            },
            "sted_ec_pathway_enrichment": {
                "name": "通路富集分析",
                "description": "读取 driver_genes.csv 基因名，使用 gseapy.enrichr（默认 KEGG + GO Biological Process）生成富集表与条形图/气泡图。",
                "tool_id": "sted_ec_pathway_enrichment",
                "default_params": {
                    "driver_genes_csv": "<sted_ec_driver_gene_extraction>",
                    "organism": "Human",
                },
            },
        }
        if step_id not in metadata_map:
            raise ValueError(f"未知的步骤ID: {step_id}")
        return metadata_map[step_id]
