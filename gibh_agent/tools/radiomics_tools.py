"""
Atomic Radiomics tools (pyradiomics, SimpleITK).
Extract radiomics features from image and mask. Requires pyradiomics>=3.0.1 and SimpleITK.
"""
import logging
from pathlib import Path
from typing import Dict, Any, Optional

from ..core.tool_registry import registry

logger = logging.getLogger(__name__)


@registry.register(
    name="radiomics_data_validation",
    description="Fast validation: check feature CSV or image/mask dir for sample count, feature dim, NaN/Inf. For DAG visibility only.",
    category="Radiomics",
    output_type="json",
)
def radiomics_data_validation(data_path: str) -> Dict[str, Any]:
    """
    极速影像组学数据校验：读取特征 CSV 或图像/目录，检查样本量、特征维度及 NaN/Inf。
    不修改数据，不写回文件。
    """
    try:
        import pandas as pd
        import numpy as np
    except ImportError:
        return {"status": "error", "error": "pandas and numpy are required"}
    p = Path(data_path)
    if not p.exists():
        return {"status": "error", "error": f"路径不存在: {data_path}"}
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
                + ("未发现异常值，内存预分配完成。" if not has_issue else f"发现 NaN={nan_count}, Inf={inf_count}，请预处理。")
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
                + ("未发现异常值，内存预分配完成。" if not (has_nan or has_inf) else "发现 NaN/Inf，请检查影像。")
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
                "message": f"影像组学数据校验通过。目录内共 {n_files} 个文件，内存预分配完成。",
                "n_files": n_files,
                "data_path": str(p.resolve()),
            }
        return {"status": "error", "error": f"不支持的路径类型或扩展名: {data_path}"}
    except Exception as e:
        logger.warning("radiomics_data_validation failed: %s", e)
        return {"status": "error", "error": str(e)}


@registry.register(
    name="radiomics_model_comparison",
    description="Train LR, SVM, RF on same train/test split, plot ROC curves with AUC in legend. Input: multi-sample features CSV. Requires explicit label_col (or column named label/target/class/y with 2+ classes).",
    category="Radiomics",
    output_type="mixed",
)
def radiomics_model_comparison(
    features_csv: str,
    output_plot_path: str,
    label_col: Optional[str] = None,
) -> Dict[str, Any]:
    """
    多算法诊断效能对比：同一划分下训练 LR、SVM、RF，在同一坐标系绘制 ROC 并标注 AUC。
    标签列：优先使用参数 label_col（可由大模型推荐）；未提供时仅当存在且仅存在一列约定名（label/target/class/y）且取值≥2 类时使用，否则报错并提示用户指定。
    """
    try:
        import pandas as pd
        import numpy as np
        from sklearn.model_selection import train_test_split
        from sklearn.linear_model import LogisticRegression
        from sklearn.svm import SVC
        from sklearn.ensemble import RandomForestClassifier
        from sklearn.metrics import roc_curve, auc
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError as e:
        return {"status": "error", "error": f"sklearn/matplotlib required: {e}"}
    p = Path(features_csv)
    if not p.is_file():
        return {"status": "error", "error": f"Features CSV 不存在: {features_csv}"}
    try:
        df = pd.read_csv(p)
        if df.shape[0] < 4:
            return {"status": "skipped", "error": "ROC 对比需要多样本 CSV（至少 4 行）。当前为单样本或样本过少。"}
        # 隐患 1 修复：不盲猜最后一列；显式 > 严格约定名 > 报错
        if label_col and str(label_col).strip():
            if str(label_col).strip() not in df.columns:
                return {
                    "status": "error",
                    "error": f"指定的标签列 '{label_col}' 不在 CSV 中。可用列: {list(df.columns[:20])}",
                    "message": "请在工作流参数或数据诊断后由 AI 推荐中指定正确的分类标签列 (label_col)。",
                    "user_message": "未找到指定的分类标签列，请指定用于 ROC 分析的标签列（如诊断结果、分组）。",
                }
            label_col = str(label_col).strip()
        else:
            label_col = None
            strict_names = ("label", "target", "class", "y")
            candidates = [
                c for c in df.columns
                if c in strict_names or c.lower() in strict_names
            ]
            candidates = [c for c in candidates if df[c].nunique() >= 2]
            if len(candidates) == 1:
                label_col = candidates[0]
                logger.warning("未指定 label_col，使用唯一约定名列: %s（建议在工作流中显式指定）", label_col)
            elif len(candidates) > 1:
                return {
                    "status": "error",
                    "error": f"存在多个可能的标签列 {candidates}，请通过参数 label_col 指定其一。",
                    "message": "请在工作流参数中指定分类标签列 (label_col)。",
                    "user_message": "检测到多列可能为分类标签，请明确指定用于 ROC 分析的标签列。",
                }
        if label_col is None:
            return {
                "status": "error",
                "error": "未指定分类标签列且未找到约定名列（label/target/class/y，且取值≥2 类）。",
                "message": "请在工作流参数或数据诊断后由 AI 推荐中指定 label_col。",
                "user_message": "无法自动识别分类标签列。请在参数中指定用于 ROC 分析的标签列（如诊断、分组），以保证结果科学有效。",
            }
        y = df[label_col].astype(int).values
        n_classes = len(np.unique(y))
        if n_classes < 2:
            return {
                "status": "error",
                "error": "ROC 分析需要至少两个有效分类（当前标签列唯一值不足 2）。请检查分组/诊断列是否正确。",
                "message": "请检查是否上传了正确的分组或诊断标签，且至少包含两类样本。",
                "user_message": "当前标签列仅有一个类别，无法进行 ROC 分析。请确认分组列 (label_col) 选择正确并包含至少两类样本。",
            }
        X = df.drop(columns=[label_col]).select_dtypes(include=[np.number])
        if X.shape[1] < 2:
            return {
                "status": "error",
                "error": "特征列不足 2 列，无法训练分类模型。",
                "user_message": "特征矩阵有效列数不足，请检查特征 CSV 或上游特征提取步骤。",
            }
        X = X.fillna(X.median())
        X = X.replace([np.inf, -np.inf], np.nan).fillna(X.median())
        X_arr = X.values
        random_state = 42
        X_train, X_test, y_train, y_test = train_test_split(X_arr, y, test_size=0.25, random_state=random_state)
        models = [
            ("LR", LogisticRegression(max_iter=500, random_state=random_state)),
            ("SVM", SVC(probability=True, random_state=random_state)),
            ("RF", RandomForestClassifier(n_estimators=50, random_state=random_state)),
        ]
        fig, ax = plt.subplots(1, 1, figsize=(6, 5))
        for name, clf in models:
            clf.fit(X_train, y_train)
            y_prob = clf.predict_proba(X_test)[:, 1]
            fpr, tpr, _ = roc_curve(y_test, y_prob)
            roc_auc = auc(fpr, tpr)
            ax.plot(fpr, tpr, lw=2, label=f"{name} (AUC = {roc_auc:.2f})")
        ax.plot([0, 1], [0, 1], "k--", lw=1)
        ax.set_xlabel("False Positive Rate")
        ax.set_ylabel("True Positive Rate")
        ax.set_title("Multi-Algorithm ROC Comparison")
        ax.legend(loc="lower right", fontsize=9)
        ax.grid(True, alpha=0.3)
        plt.tight_layout()
        out = Path(output_plot_path)
        out.parent.mkdir(parents=True, exist_ok=True)
        if out.suffix.lower() not in (".png", ".jpg", ".pdf"):
            out = out.parent / (out.stem + ".png")
        fig.savefig(str(out), bbox_inches="tight", dpi=150)
        plt.close()
        return {
            "status": "success",
            "message": "多算法 ROC 对比图已保存",
            "plot_path": str(out.resolve()),
            "features_csv": str(p.resolve()),
            "features_csv_path": str(p.resolve()),
        }
    except Exception as e:
        logger.exception("radiomics_model_comparison failed: %s", e)
        return {"status": "error", "error": str(e), "features_csv": features_csv}


@registry.register(
    name="extract_radiomics_features",
    description="Extract radiomics features from a medical image and its segmentation mask. Returns a CSV path with feature values (PyRadiomics).",
    category="Radiomics",
    output_type="file_path",
)
def extract_radiomics_features(
    image_path: str,
    mask_path: str,
    output_path: Optional[str] = None,
) -> Dict[str, Any]:
    """
    Extract radiomics features using PyRadiomics.

    Args:
        image_path: Path to the imaging file (e.g. NIfTI, NRRD).
        mask_path: Path to the segmentation mask (same geometry as image).
        output_path: Optional path for output CSV; if not set, writes next to image with _radiomics.csv suffix.

    Returns:
        Dict with status, output_path (CSV), n_features, and optional error.
    """
    try:
        from radiomics import featureextractor
        import pandas as pd
    except ImportError as e:
        logger.warning("pyradiomics not installed: %s", e)
        return {"status": "error", "error": "pyradiomics is required. Install with: pip install pyradiomics>=3.0.1 SimpleITK"}
    pi, pm = Path(image_path), Path(mask_path)
    if not pi.is_file():
        return {"status": "error", "error": f"Image file not found: {image_path}"}
    if not pm.is_file():
        return {"status": "error", "error": f"Mask file not found: {mask_path}"}
    if output_path is None:
        output_path = str(pi.parent / (pi.stem + "_radiomics.csv"))
    out = Path(output_path)
    out.parent.mkdir(parents=True, exist_ok=True)
    try:
        extractor = featureextractor.RadiomicsFeatureExtractor()
        result = extractor.execute(str(pi), str(pm))
        if not result:
            return {"status": "error", "error": "PyRadiomics returned no features."}
        # Drop diagnostic keys if present; keep feature names
        feature_dict = {k: v for k, v in result.items() if not k.startswith("diagnostics")}
        df = pd.DataFrame([feature_dict])
        df.to_csv(out, index=False)
        return {
            "status": "success",
            "output_path": str(out.resolve()),
            "n_features": len(feature_dict),
            "image_path": str(pi.resolve()),
            "mask_path": str(pm.resolve()),
        }
    except Exception as e:
        logger.exception("extract_radiomics_features failed: %s", e)
        return {"status": "error", "error": str(e)}
