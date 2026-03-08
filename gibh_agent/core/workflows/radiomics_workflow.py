"""
Radiomics workflow (medical imaging) — full clinical pipeline.

DAG (8 steps):
  radiomics_data_validation -> load_image -> preprocess -> preview_slice
                                            \\ -> extract_features -> radiomics_model_comparison -> calc_score -> viz_score
"""
import logging
from pathlib import Path
from typing import Dict, Any, List, Optional

from .base import BaseWorkflow

logger = logging.getLogger(__name__)


class RadiomicsWorkflow(BaseWorkflow):
    """
    Radiomics (medical imaging) workflow — full clinical pipeline.

    Steps:
    0. radiomics_data_validation - Data & ROI validation (秀肌肉)
    1. load_image                - Load NIfTI/DICOM (radiomics_load_medical_image)
    2. preprocess                - Resample + normalize (radiomics_preprocessing)
    3. preview_slice             - Export mid-slice PNG (radiomics_plot_mid_slice)
    4. extract_features         - PyRadiomics feature extraction (radiomics_extract_features)
    5. radiomics_model_comparison - Multi-algorithm ROC comparison (LR/SVM/RF)
    6. calc_score                - Rad-Score + risk probability (calculate_rad_score)
    7. viz_score                 - Bar/gauge chart of Rad-Score vs reference (plot_rad_score)
    """

    def get_name(self) -> str:
        return "Radiomics"

    def get_description(self) -> str:
        return (
            "影像组学完整临床流程：加载 -> 预处理(重采样+归一化) -> 预览 -> 特征提取 -> Rad-Score 计算 -> 评分可视化"
        )

    def get_steps_dag(self) -> Dict[str, List[str]]:
        return {
            "radiomics_data_validation": [],
            "load_image": ["radiomics_data_validation"],
            "preprocess": ["load_image"],
            "preview_slice": ["preprocess"],
            "extract_features": ["preprocess"],
            "radiomics_model_comparison": ["extract_features"],
            "calc_score": ["radiomics_model_comparison"],
            "viz_score": ["calc_score"],
        }

    def get_step_metadata(self, step_id: str) -> Dict[str, Any]:
        metadata_map = {
            "radiomics_data_validation": {
                "name": "数据与 ROI 校验",
                "description": "快速校验特征 CSV 或影像/目录：样本量、特征维度、NaN/Inf",
                "tool_id": "radiomics_data_validation",
                "default_params": {"data_path": ""},
            },
            "load_image": {
                "name": "加载医学影像",
                "description": "使用 SimpleITK 加载 NIfTI/DICOM 影像，返回尺寸与间距等元数据",
                "tool_id": "radiomics_load_medical_image",
                "default_params": {
                    "image_path": "",
                    "mask_path": "",
                },
            },
            "preprocess": {
                "name": "预处理",
                "description": "重采样至 1×1×1 mm 并强度归一化，输出预处理后的 .nii.gz 供后续特征提取",
                "tool_id": "radiomics_preprocessing",
                "default_params": {
                    "image_path": "",
                    "output_path": "",
                    "spacing_mm": [1.0, 1.0, 1.0],
                    "normalize": True,
                },
            },
            "preview_slice": {
                "name": "预览中间层",
                "description": "导出 3D 体积的中间层为 PNG，便于前端展示",
                "tool_id": "radiomics_plot_mid_slice",
                "default_params": {
                    "image_path": "",
                    "output_path": "",
                    "mask_path": "",
                },
            },
            "extract_features": {
                "name": "提取影像组学特征",
                "description": "基于 PyRadiomics 提取形状、一阶与纹理特征，输出 CSV；mask 可选",
                "tool_id": "radiomics_extract_features",
                "default_params": {
                    "image_path": "",
                    "mask_path": "",
                    "output_path": "",
                    "config_path": "",
                },
            },
            "radiomics_model_comparison": {
                "name": "多算法诊断效能对比",
                "description": "LR/SVM/RF 同划分下训练，绘制 ROC 曲线并标注 AUC。需通过参数或 AI 推荐指定 label_col（分类标签列）。",
                "tool_id": "radiomics_model_comparison",
                "default_params": {
                    "features_csv": "",
                    "output_plot_path": "<output_dir>/radiomics_roc_comparison.png",
                    "label_col": "",
                },
            },
            "calc_score": {
                "name": "计算 Rad-Score",
                "description": "基于预定义特征签名计算 Rad-Score 与风险概率（Sigmoid）",
                "tool_id": "calculate_rad_score",
                "default_params": {
                    "features_csv_path": "",
                    "output_path": "",
                    "sigmoid_scale": 1.0,
                },
            },
            "viz_score": {
                "name": "评分可视化",
                "description": "绘制患者 Rad-Score 与参考人群对比的柱状图，保存 PNG",
                "tool_id": "plot_rad_score",
                "default_params": {
                    "rad_score": None,
                    "rad_score_csv_path": "",
                    "output_path": "",
                    "reference_mean": 0.0,
                    "reference_std": 1.0,
                },
            },
        }
        return metadata_map.get(step_id, {})

    def generate_template(
        self,
        target_steps: Optional[List[str]] = None,
        file_metadata: Optional[Dict[str, Any]] = None,
    ) -> Dict[str, Any]:
        """Generate workflow template; fill image_path from file_metadata; leave mask_path and chained outputs for executor."""
        if target_steps is None:
            target_steps = list(self.steps_dag.keys())
        resolved_steps = self.resolve_dependencies(target_steps)
        if not resolved_steps:
            resolved_steps = list(self.steps_dag.keys())

        file_path = (file_metadata or {}).get("file_path")
        placeholder = "<PENDING_UPLOAD>" if not file_path else file_path
        base_dir = Path(file_path).parent if file_path else Path(".")
        output_preview = str(base_dir / "radiomics_preview.png")
        output_preprocessed = str(base_dir / "radiomics_preprocessed.nii.gz")
        output_features_csv = str(base_dir / "radiomics_features.csv")
        output_rad_score_json = str(base_dir / "rad_score.json")
        output_viz = str(base_dir / "radiomics_rad_score.png")

        steps = []
        for step_id in resolved_steps:
            step_meta = self.get_step_metadata(step_id)
            params = dict(step_meta.get("default_params", {}))
            tool_id = step_meta.get("tool_id", step_id)

            if step_id == "radiomics_data_validation":
                params["data_path"] = placeholder
            elif step_id == "load_image":
                params["image_path"] = placeholder
                params["mask_path"] = ""
            elif step_id == "preprocess":
                params["image_path"] = placeholder
                params["output_path"] = output_preprocessed if file_path else ""
                params["spacing_mm"] = [1.0, 1.0, 1.0]
                params["normalize"] = True
            elif step_id == "preview_slice":
                params["image_path"] = placeholder  # executor may substitute preprocessed_path
                params["output_path"] = output_preview if file_path else ""
                params["mask_path"] = ""
            elif step_id == "extract_features":
                params["image_path"] = placeholder
                params["mask_path"] = ""
                params["output_path"] = output_features_csv if file_path else ""
                params["config_path"] = ""
            elif step_id == "radiomics_model_comparison":
                params["features_csv"] = "<extract_features>"
                params["output_plot_path"] = str(base_dir / "radiomics_roc_comparison.png") if file_path else "<output_dir>/radiomics_roc_comparison.png"
                params["label_col"] = params.get("label_col", "")
            elif step_id == "calc_score":
                params["features_csv_path"] = "<extract_features>"
                params["output_path"] = output_rad_score_json if file_path else ""
                params["sigmoid_scale"] = 1.0
            elif step_id == "viz_score":
                params["rad_score"] = None
                params["rad_score_csv_path"] = "<calc_score>"  # executor 从 calc_score 的 output_path 填充
                params["output_path"] = output_viz if file_path else ""
                params["reference_mean"] = 0.0
                params["reference_std"] = 1.0

            steps.append({
                "id": step_id,
                "step_id": step_id,
                "tool_id": tool_id,
                "name": step_meta.get("name", step_id),
                "step_name": step_meta.get("name", step_id),
                "description": step_meta.get("description", ""),
                "desc": (step_meta.get("description", ""))[:100],
                "selected": True,
                "params": params,
            })

        workflow_name = self._generate_workflow_name(target_steps, file_metadata)
        is_template = not file_path
        return {
            "type": "workflow_config",
            "workflow_data": {
                "workflow_name": workflow_name,
                "name": workflow_name,
                "steps": steps,
            },
            "file_paths": [file_path] if file_path else [],
            "template_mode": is_template,
        }
