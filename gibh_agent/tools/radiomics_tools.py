"""
影像组学工具：校验在本地；PyRadiomics / LASSO / ROC 等重计算经 HTTP 调用 worker-pyskills（TaaS）。

路径约定：入参为已解析的绝对路径（通常位于 UPLOAD_DIR）；禁止在主 Agent 进程内加载 pyradiomics/sklearn 重依赖。
"""
from __future__ import annotations

import logging
import os
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import httpx

from ..core.executor import WorkflowExecutor
from ..core.tool_registry import registry
from ..core.file_inspector import resolve_omics_paths
from ..core.utils import safe_tool_execution

logger = logging.getLogger(__name__)

PYSKILLS_BASE_URL = (os.getenv("PYSKILLS_BASE_URL") or "http://worker-pyskills:8000").rstrip("/")
_RAD_TIMEOUT = float(os.getenv("PYSKILLS_RAD_TIMEOUT", "600"))


def _radiomics_paths_from_resolved(resolved: Dict[str, List[str]]) -> Tuple[List[str], Optional[str], Optional[str]]:
    images = resolved.get("images") or []
    masks = resolved.get("masks") or []
    tables = resolved.get("tables") or []
    unknown = resolved.get("unknown") or []
    all_paths: List[str] = images + masks + tables + unknown
    image_path: Optional[str] = images[0] if images else None
    mask_path: Optional[str] = masks[0] if masks else None
    if tables and not images and not masks:
        return all_paths, tables[0], None
    if images and not masks and tables:
        return all_paths, tables[0], None
    if not image_path and unknown:
        image_path = unknown[0]
    if not mask_path and len(unknown) >= 2:
        mask_path = unknown[1]
    return all_paths, image_path, mask_path


def _resolve_file_path(file_path: str) -> str:
    ex = WorkflowExecutor()
    return ex._resolve_file_path(file_path)


def _post_radiomics(subpath: str, body: Dict[str, Any]) -> Dict[str, Any]:
    path = subpath if subpath.startswith("/") else f"/{subpath}"
    url = f"{PYSKILLS_BASE_URL}/api/radiomics{path}"
    try:
        with httpx.Client(timeout=_RAD_TIMEOUT) as client:
            r = client.post(url, json=body)
            try:
                data = r.json() if r.content else {}
            except Exception:
                data = {"detail": r.text[:2000]}
    except Exception as e:
        logger.exception("radiomics worker POST %s", path)
        return {
            "status": "error",
            "message": f"无法连接 worker-pyskills（{url}）。请确认 PYSKILLS_BASE_URL 与容器网络。详情：{e}",
        }
    if r.status_code >= 400:
        det = data.get("detail")
        if isinstance(det, list):
            msg = "; ".join(str(x.get("msg", x)) for x in det)
        else:
            msg = str(det or data)
        return {"status": "error", "message": msg}
    if data.get("status") != "success":
        return {"status": "error", "message": data.get("message", str(data))}
    return {"status": "success", **data}


@registry.register(
    name="radiomics_data_validation",
    description="Fast validation: check feature CSV or image/mask dir for sample count, feature dim, NaN/Inf. For DAG visibility only.",
    category="Radiomics",
    output_type="json",
)
@safe_tool_execution
def radiomics_data_validation(data_path: str) -> Dict[str, Any]:
    """
    极速影像组学数据校验：读取特征 CSV 或图像/目录（本地轻量 pandas/numpy）。
    """
    try:
        import numpy as np
        import pandas as pd
    except ImportError:
        return {"status": "error", "message": "pandas 与 numpy 为必需依赖"}
    resolved = resolve_omics_paths(data_path)
    paths, image_path, mask_path = _radiomics_paths_from_resolved(resolved)
    if not paths:
        return {"status": "error", "message": f"路径为空或无效: {data_path}"}
    if len(paths) == 1 and (resolved.get("images") or resolved.get("unknown")) and not resolved.get("tables"):
        p = Path(paths[0])
        if p.is_file() and p.suffix.lower() != ".csv" and (
            p.suffix.lower() in (".nii", ".gz") or ".nii" in p.name.lower() or p.suffix.lower() == ".dcm"
        ):
            return {
                "status": "error",
                "message": "影像组学分析需要成对的原图(Image)和掩膜(Mask)文件，请同时上传或选择这两个文件。",
            }
    if image_path and not Path(image_path).exists():
        return {"status": "error", "message": f"原图路径不存在: {image_path}"}
    if mask_path and not Path(mask_path).exists():
        return {"status": "error", "message": f"掩膜路径不存在: {mask_path}"}
    first_path = image_path or mask_path or paths[0]
    p = Path(first_path)
    if not p.exists():
        return {"status": "error", "message": f"路径不存在: {first_path}"}
    try:
        if p.is_file() and p.suffix.lower() == ".csv":
            df = pd.read_csv(p, nrows=10000)
            n_samples, n_features = df.shape
            numeric = df.select_dtypes(include=[np.number])
            if numeric.size > 0:
                nan_count = numeric.isnull().sum().sum()
                inf_count = np.isinf(numeric.values).sum()
                has_issue = nan_count > 0 or inf_count > 0
            else:
                nan_count = inf_count = 0
                has_issue = False
            msg = (
                f"影像组学数据校验通过。共包含 {n_samples} 个病例，{n_features} 个高维特征。"
                + ("未发现异常值。" if not has_issue else f"发现 NaN={nan_count}, Inf={inf_count}，请预处理。")
            )
            return {
                "status": "success",
                "message": msg,
                "n_samples": int(n_samples),
                "n_features": int(n_features),
                "data_path": str(p.resolve()),
            }
        if p.is_file() and (p.suffix.lower() in (".nii", ".gz") or ".nii" in p.name.lower() or p.suffix.lower() == ".dcm"):
            import SimpleITK as sitk

            img = sitk.ReadImage(str(p))
            arr = sitk.GetArrayFromImage(img)
            n_voxels = arr.size
            has_nan = bool(np.isnan(arr).any())
            has_inf = bool(np.isinf(arr).any())
            msg = (
                f"影像组学数据校验通过。影像尺寸 {arr.shape}，体素数 {n_voxels}。"
                + ("未发现异常值。" if not (has_nan or has_inf) else "发现 NaN/Inf，请检查影像。")
            )
            return {
                "status": "success",
                "message": msg,
                "n_voxels": int(n_voxels),
                "shape": list(arr.shape),
                "data_path": str(p.resolve()),
            }
        if p.is_dir():
            files = list(p.glob("*"))
            n_files = len([f for f in files if f.is_file()])
            return {
                "status": "success",
                "message": f"影像组学数据校验通过。目录内共 {n_files} 个文件。",
                "n_files": n_files,
                "data_path": str(p.resolve()),
            }
        return {"status": "error", "message": f"不支持的路径类型或扩展名: {data_path}"}
    except Exception as e:
        logger.warning("radiomics_data_validation failed: %s", e)
        return {"status": "error", "message": str(e)}


@registry.register(
    name="extract_radiomics_features",
    description="TaaS: PyRadiomics 单病例特征提取（HTTP→worker-pyskills）。需原图+掩膜绝对路径；产出写入 /uploads/results/radiomics/。",
    category="Radiomics",
    output_type="file_path",
)
@safe_tool_execution
def extract_radiomics_features(
    image_path: str,
    mask_path: str,
    output_path: Optional[str] = None,
) -> Dict[str, Any]:
    """
    调用 worker ``POST /api/radiomics/extract`` 单对模式。

    Args:
        image_path: 影像绝对路径，或逗号/分号拼接（经 resolve_omics_paths 拆解）。
        mask_path: 掩膜绝对路径（拼接模式下可被覆盖）。
        output_path: 保留兼容，TaaS 结果路径由 worker 分配。

    Returns:
        含 ``features_csv_path``、``features_csv_url``、``n_cases`` 等。
    """
    _ = output_path
    ip, mp = str(image_path), str(mask_path)
    if ip and ("," in ip or ";" in ip):
        resolved = resolve_omics_paths(ip)
        _, ip, mp2 = _radiomics_paths_from_resolved(resolved)
        mp = mp2 or mp
    if not ip or not mp:
        return {"status": "error", "message": "需要有效的 image_path 与 mask_path。"}
    try:
        ip_abs = _resolve_file_path(ip)
        mp_abs = _resolve_file_path(mp)
    except Exception as e:
        return {"status": "error", "message": str(e)}
    data = _post_radiomics(
        "/extract",
        {"image_path": ip_abs, "mask_path": mp_abs, "case_id": Path(ip_abs).stem},
    )
    if data.get("status") != "success":
        return data
    return {
        "status": "success",
        "message": data.get("message", ""),
        "output_path": data.get("features_csv_path"),
        "features_csv_path": data.get("features_csv_path"),
        "features_csv_url": data.get("features_csv_url"),
        "n_features": data.get("n_feature_cols"),
        "n_cases": data.get("n_cases"),
        "image_path": ip_abs,
        "mask_path": mp_abs,
    }


@registry.register(
    name="radiomics_extract_manifest",
    description="TaaS: 按 manifest CSV 批量 PyRadiomics（列含 case_id,image_path,mask_path）。结果 features.csv 在 /uploads/results/radiomics/。",
    category="Radiomics",
    output_type="file_path",
)
@safe_tool_execution
def radiomics_extract_manifest(file_path: str) -> Dict[str, Any]:
    """
    Args:
        file_path: manifest CSV 的绝对路径（须位于共享上传目录下）。

    Returns:
        聚合特征表路径与公开 URL。
    """
    try:
        manifest_abs = _resolve_file_path(file_path)
    except Exception as e:
        return {"status": "error", "message": str(e)}
    data = _post_radiomics("/extract", {"manifest_path": manifest_abs})
    if data.get("status") != "success":
        return data
    return {
        "status": "success",
        "message": data.get("message", ""),
        "output_path": data.get("features_csv_path"),
        "features_csv_path": data.get("features_csv_path"),
        "features_csv_url": data.get("features_csv_url"),
        "n_cases": data.get("n_cases"),
        "n_feature_cols": data.get("n_feature_cols"),
    }


@registry.register(
    name="radiomics_lasso_select",
    description="TaaS: L1 逻辑回归 CV 特征筛选；输出 selected_features.csv 与 lasso_path.png（worker-pyskills）。",
    category="Radiomics",
    output_type="mixed",
)
@safe_tool_execution
def radiomics_lasso_select(
    file_path: str,
    label_col: str,
    cv_folds: int = 5,
) -> Dict[str, Any]:
    """
    Args:
        file_path: 多样本特征表 CSV 绝对路径。
        label_col: 二分类标签列名。
        cv_folds: 交叉验证折数（3–15）。

    Returns:
        筛选后表路径、系数图路径及 URL。
    """
    try:
        csv_abs = _resolve_file_path(file_path)
    except Exception as e:
        return {"status": "error", "message": str(e)}
    data = _post_radiomics(
        "/lasso",
        {"features_csv_path": csv_abs, "label_col": label_col, "cv_folds": int(cv_folds)},
    )
    if data.get("status") != "success":
        return data
    return {
        "status": "success",
        "message": data.get("message", ""),
        "selected_features_csv_path": data.get("selected_features_csv_path"),
        "selected_features_csv_url": data.get("selected_features_csv_url"),
        "lasso_path_png_path": data.get("lasso_path_png_path"),
        "lasso_path_png_url": data.get("lasso_path_png_url"),
        "n_selected": data.get("n_selected"),
    }


@registry.register(
    name="radiomics_model_comparison",
    description="TaaS: 逻辑回归 + 训练/测试划分 + ROC（worker-pyskills）。输出 metrics.csv 与 roc_curve.png。",
    category="Radiomics",
    output_type="mixed",
)
@safe_tool_execution
def radiomics_model_comparison(
    features_csv: str,
    output_plot_path: str,
    label_col: Optional[str] = None,
) -> Dict[str, Any]:
    """
    Args:
        features_csv: 特征 CSV 绝对路径（可为 LASSO 筛选后的表）。
        output_plot_path: 兼容保留；实际 ROC 路径由 worker 写入 results/radiomics。
        label_col: 分类标签列；若为空则尝试约定列名（与旧版行为一致）。

    Returns:
        ROC 图、指标表路径及 AUC 等。
    """
    _ = output_plot_path
    try:
        csv_abs = _resolve_file_path(features_csv)
    except Exception as e:
        return {"status": "error", "message": str(e)}

    lc = (label_col or "").strip()
    if not lc:
        try:
            import pandas as pd

            df = pd.read_csv(csv_abs, nrows=5)
            candidate_cols = [
                "label", "target", "class", "y",
                "group", "status", "diagnosis", "type", "condition", "cohort",
            ]
            cand_lower = {c.lower() for c in candidate_cols}
            candidates = [c for c in df.columns if str(c).strip().lower() in cand_lower]
            candidates = [c for c in candidates if df[c].nunique() >= 2]
            if len(candidates) == 1:
                lc = candidates[0]
            elif len(candidates) > 1:
                return {
                    "status": "error",
                    "message": f"存在多个可能标签列 {candidates}，请显式传入 label_col。",
                }
        except Exception as e:
            return {"status": "error", "message": f"无法推断 label_col: {e}"}
    if not lc:
        return {"status": "error", "message": "请指定 label_col（分类标签列名）。"}

    data = _post_radiomics("/model", {"features_csv_path": csv_abs, "label_col": lc})
    if data.get("status") != "success":
        return data
    return {
        "status": "success",
        "message": data.get("message", ""),
        "plot_path": data.get("roc_curve_png_path"),
        "roc_curve_png_path": data.get("roc_curve_png_path"),
        "roc_curve_png_url": data.get("roc_curve_png_url"),
        "metrics_csv_path": data.get("metrics_csv_path"),
        "metrics_csv_url": data.get("metrics_csv_url"),
        "features_csv": csv_abs,
        "features_csv_path": csv_abs,
        "accuracy": data.get("accuracy"),
        "roc_auc": data.get("roc_auc"),
    }
