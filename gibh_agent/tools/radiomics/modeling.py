"""
Radiomics modeling: preprocessing, Rad-Score calculation, and score visualization.

Preprocessing (resampling + normalization) before feature extraction.
Pre-defined signature: weighted sum of key features -> Rad-Score -> Sigmoid risk probability.
"""
import csv
import json
import logging
from pathlib import Path
from typing import Dict, Any, Optional, List, Union

from ...core.tool_registry import registry

logger = logging.getLogger(__name__)

# Pre-defined signature weights (feature name -> weight) for Rad-Score demo
DEFAULT_SIGNATURE = {
    "original_firstorder_Entropy": 0.25,
    "original_firstorder_Mean": -0.15,
    "original_shape_Sphericity": 0.20,
    "original_glcm_Contrast": 0.15,
    "original_glrlm_LongRunHighGrayLevelEmphasis": 0.10,
    "original_glszm_SizeZoneNonUniformity": 0.15,
}


@registry.register(
    name="radiomics_preprocessing",
    description="Resample image to isotropic spacing (e.g. 1x1x1 mm) and normalize intensity using SimpleITK. Output preprocessed .nii.gz for downstream feature extraction.",
    category="Radiomics",
    output_type="json",
)
def _parse_spacing_mm(spacing_mm: Optional[Union[str, List[float]]]) -> List[float]:
    """Parse spacing_mm from string '1,1,1' or list [1,1,1] to [float, float, float]."""
    if spacing_mm is None:
        return [1.0, 1.0, 1.0]
    if isinstance(spacing_mm, str):
        s = spacing_mm.strip()
        if not s:
            return [1.0, 1.0, 1.0]
        parts = [x.strip() for x in s.split(",") if x.strip()]
        if len(parts) != 3:
            raise ValueError("spacing_mm must be 3 values (e.g. '1,1,1' or [1,1,1]).")
        return [float(x) for x in parts]
    if isinstance(spacing_mm, (list, tuple)) and len(spacing_mm) == 3:
        return [float(x) for x in spacing_mm]
    raise ValueError("spacing_mm must be [x,y,z] or string 'x,y,z'.")


def radiomics_preprocessing(
    image_path: str,
    output_path: Optional[str] = None,
    spacing_mm: Optional[Union[str, List[float]]] = None,
    normalize: bool = True,
) -> Dict[str, Any]:
    """
    Preprocess medical image: resample to target spacing and optional intensity normalization.

    Args:
        image_path: Path to input NIfTI/DICOM.
        output_path: Path for preprocessed .nii.gz; default next to image as radiomics_preprocessed.nii.gz.
        spacing_mm: Target spacing [x,y,z] in mm; string "1,1,1" or list [1.0,1.0,1.0]; default [1,1,1].
        normalize: If True, scale intensity to [0, 1] (min-max) after resampling.

    Returns:
        Dict with status, preprocessed_path, spacing, and optional error.
    """
    try:
        import SimpleITK as sitk
        import numpy as np
    except ImportError as e:
        logger.warning("SimpleITK/numpy not installed: %s", e)
        return {"status": "error", "error": "SimpleITK and numpy are required."}

    p = Path(image_path)
    if not p.exists():
        return {"status": "error", "error": f"Image not found: {image_path}"}

    try:
        spacing_mm = _parse_spacing_mm(spacing_mm)
    except ValueError as e:
        return {"status": "error", "error": str(e)}

    try:
        img = sitk.ReadImage(str(p))
        size = list(img.GetSize())
        orig_spacing = list(img.GetSpacing())

        # Resample to target spacing (spacing_mm already parsed to [float, float, float])
        target_spacing = list(spacing_mm)
        target_size = [
            int(round(size[i] * orig_spacing[i] / target_spacing[i]))
            for i in range(3)
        ]
        target_size = [max(1, t) for t in target_size]

        resampled = sitk.Resample(
            img,
            target_size,
            sitk.Transform(),
            sitk.sitkLinear,
            img.GetOrigin(),
            target_spacing,
            img.GetDirection(),
            0.0,
            img.GetPixelID(),
        )

        if normalize:
            arr = sitk.GetArrayFromImage(resampled)
            min_v, max_v = float(arr.min()), float(arr.max())
            if max_v > min_v:
                arr = (arr - min_v) / (max_v - min_v)
            else:
                arr = arr * 0.0
            arr = arr.astype(np.float32)
            norm_img = sitk.GetImageFromArray(arr)
            norm_img.SetOrigin(resampled.GetOrigin())
            norm_img.SetSpacing(resampled.GetSpacing())
            norm_img.SetDirection(resampled.GetDirection())
            resampled = norm_img

        out_p = Path(output_path) if (output_path and str(output_path).strip()) else (p.parent / "radiomics_preprocessed.nii.gz")
        out_p = out_p.resolve()
        if out_p.suffix != ".gz" and out_p.suffix != ".nii":
            out_p = out_p.parent / (out_p.stem + ".nii.gz")
        out_p.parent.mkdir(parents=True, exist_ok=True)
        sitk.WriteImage(resampled, str(out_p))

        return {
            "status": "success",
            "preprocessed_path": str(out_p),
            "image_path": image_path,
            "spacing_mm": target_spacing,
            "size_after": list(resampled.GetSize()),
        }
    except Exception as e:
        logger.exception("radiomics_preprocessing failed: %s", e)
        return {"status": "error", "error": str(e), "image_path": image_path}


@registry.register(
    name="calculate_rad_score",
    description="Compute a Rad-Score from extracted radiomics features using a pre-defined signature (weighted sum), then risk probability via Sigmoid. Output JSON/CSV with score and risk.",
    category="Radiomics",
    output_type="json",
)
def calculate_rad_score(
    features_csv_path: str,
    output_path: Optional[str] = None,
    signature: Optional[Dict[str, float]] = None,
    sigmoid_scale: float = 1.0,
) -> Dict[str, Any]:
    """
    Calculate Rad-Score: weighted sum of signature features, then risk = sigmoid(score).

    Args:
        features_csv_path: Path to CSV from radiomics_extract_features (columns: feature, value).
        output_path: Optional path for result JSON or CSV; default next to input as rad_score.json.
        signature: Feature name -> weight; if None, use DEFAULT_SIGNATURE.
        sigmoid_scale: Scale factor inside sigmoid (score * scale); default 1.0.

    Returns:
        Dict with status, rad_score, risk_probability, output_path, and optional error.
    """
    if not features_csv_path or not str(features_csv_path).strip():
        return {
            "status": "skipped",
            "error": "未提供特征 CSV（依赖上一步「提取影像组学特征」）。请先完成特征提取；若上一步因缺少 pyradiomics 失败，请重建 Docker 镜像后重试。",
            "features_csv_path": features_csv_path or "",
        }
    csv_p = Path(features_csv_path).resolve()
    if not csv_p.exists() or csv_p.is_dir():
        return {
            "status": "error",
            "error": f"Features CSV 不存在或为目录: {features_csv_path}。请先成功执行「提取影像组学特征」（需安装 pyradiomics）。",
            "features_csv_path": features_csv_path,
        }

    sig = signature or DEFAULT_SIGNATURE

    try:
        features = {}
        with open(csv_p, "r", encoding="utf-8") as f:
            reader = csv.reader(f)
            next(reader, None)  # header
            for row in reader:
                if len(row) >= 2:
                    try:
                        features[row[0].strip()] = float(row[1])
                    except ValueError:
                        pass

        # Weighted sum (only keys present in CSV)
        score = 0.0
        used = []
        for name, w in sig.items():
            if name in features:
                score += w * features[name]
                used.append(name)

        if not used:
            return {
                "status": "error",
                "error": "No signature features found in CSV. Expected keys like original_firstorder_Entropy.",
                "features_csv_path": features_csv_path,
            }

        # Risk probability: sigmoid(score * scale)
        import math
        risk = 1.0 / (1.0 + math.exp(-score * sigmoid_scale))

        out_p = Path(output_path) if (output_path and str(output_path).strip()) else (csv_p.parent / "rad_score.json")
        out_p = out_p.resolve()
        if out_p.is_dir():
            out_p = out_p / "rad_score.json"
        out_p.parent.mkdir(parents=True, exist_ok=True)
        result = {
            "rad_score": round(score, 6),
            "risk_probability": round(risk, 6),
            "signature_used": used,
            "features_csv_path": features_csv_path,
        }
        with open(out_p, "w", encoding="utf-8") as f:
            json.dump(result, f, indent=2, ensure_ascii=False)

        # Also write a tiny CSV for compatibility
        csv_out = out_p.with_suffix(".csv")
        if csv_out == csv_p:
            csv_out = out_p.parent / "rad_score_result.csv"
        with open(csv_out, "w", newline="", encoding="utf-8") as f:
            w = csv.writer(f)
            w.writerow(["metric", "value"])
            w.writerow(["rad_score", result["rad_score"]])
            w.writerow(["risk_probability", result["risk_probability"]])

        # 追加：特征权重棒棒糖图 + 特征相关性聚类热图（不替换原有逻辑）
        lollipop_path = None
        clustermap_path = None
        out_dir = out_p.parent
        try:
            import matplotlib
            matplotlib.use("Agg")
            import matplotlib.pyplot as plt
            import numpy as np
            # 特征权重棒棒糖图：Top 10–15 按 |weight| 排序
            sorted_sig = sorted(sig.items(), key=lambda x: -abs(x[1]))[:15]
            if sorted_sig:
                names = [n[:40] for n, _ in sorted_sig]
                weights = [w for _, w in sorted_sig]
                fig, ax = plt.subplots(figsize=(7, max(4, len(names) * 0.28)))
                y_pos = range(len(names))
                ax.hlines(y_pos, 0, weights, color="lightgray", lw=1.5)
                colors = ["#e74c3c" if w > 0 else "#3498db" for w in weights]
                ax.scatter(weights, y_pos, c=colors, s=60, zorder=2)
                ax.set_yticks(y_pos)
                ax.set_yticklabels(names, fontsize=8)
                ax.set_xlabel("Feature weight (signature)")
                ax.set_title("Feature importance (Rad-Score signature)")
                ax.axvline(x=0, color="gray", linestyle="--", alpha=0.7)
                plt.tight_layout()
                lollipop_path = str(out_dir / "radiomics_feature_lollipop.png")
                fig.savefig(lollipop_path, bbox_inches="tight", dpi=150)
                plt.close()
        except Exception as ex:
            logger.warning("特征权重棒棒糖图生成失败: %s", ex)
        try:
            import pandas as pd
            import numpy as np
            import seaborn as sns
            df_wide = pd.read_csv(csv_p)
            if df_wide.shape[0] >= 2 and df_wide.shape[1] >= 10:
                numeric = df_wide.select_dtypes(include=[np.number])
                if numeric.shape[1] >= 10:
                    var = numeric.var()
                    top_cols = var.nlargest(min(50, len(var))).index.tolist()
                    sub = numeric[top_cols]
                    cor = sub.corr()
                    g = sns.clustermap(cor, cmap="coolwarm", center=0, figsize=(10, 8), yticklabels=True, xticklabels=True)
                    g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), fontsize=7)
                    g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), fontsize=7, rotation=45, ha="right")
                    clustermap_path = str(out_dir / "radiomics_feature_clustermap.png")
                    g.savefig(clustermap_path, bbox_inches="tight", dpi=150)
                    plt.close()
        except Exception as ex:
            logger.warning("特征相关性聚类热图生成失败（可能为单样本 CSV）: %s", ex)

        out_result = {
            "status": "success",
            "rad_score": result["rad_score"],
            "risk_probability": result["risk_probability"],
            "output_path": str(out_p),
            "output_csv_path": str(csv_out),
            "features_csv_path": features_csv_path,
        }
        if lollipop_path:
            out_result["lollipop_path"] = lollipop_path
        if clustermap_path:
            out_result["clustermap_path"] = clustermap_path
        return out_result
    except Exception as e:
        logger.exception("calculate_rad_score failed: %s", e)
        return {"status": "error", "error": str(e), "features_csv_path": features_csv_path}


@registry.register(
    name="plot_rad_score",
    description="Draw a bar or gauge chart of the patient Rad-Score vs reference population. Saves PNG for UI display.",
    category="Radiomics",
    output_type="json",
)
def plot_rad_score(
    rad_score: Optional[Union[float, str]] = None,
    rad_score_csv_path: Optional[str] = None,
    output_path: Optional[str] = None,
    reference_mean: Union[float, str] = 0.0,
    reference_std: Union[float, str] = 1.0,
) -> Dict[str, Any]:
    """
    Visualize Rad-Score: bar chart (patient vs reference distribution) or gauge.

    Args:
        rad_score: Direct score value; if None, read from rad_score_csv_path or JSON.
        rad_score_csv_path: Path to rad_score result CSV/JSON from calculate_rad_score.
        output_path: Path for output PNG; default next to CSV or radiomics_rad_score.png.
        reference_mean: Mock reference population mean for demo (numeric or string).
        reference_std: Mock reference population std for demo (numeric or string).

    Returns:
        Dict with status, output_path, rad_score, and optional error.
    """
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import numpy as np
    except ImportError as e:
        logger.warning("matplotlib/numpy not installed: %s", e)
        return {"status": "error", "error": "matplotlib and numpy are required."}

    try:
        ref_mean = float(reference_mean)
        ref_std = float(reference_std)
    except (TypeError, ValueError):
        ref_mean, ref_std = 0.0, 1.0

    # 空字符串视为未提供
    score_val = rad_score
    if isinstance(score_val, str) and not score_val.strip():
        score_val = None
    if score_val is None and rad_score_csv_path and str(rad_score_csv_path).strip():
        p = Path(rad_score_csv_path)
        if p.exists():
            if p.suffix.lower() == ".json":
                with open(p, "r", encoding="utf-8") as f:
                    data = json.load(f)
                    score_val = data.get("rad_score")
            else:
                with open(p, "r", encoding="utf-8") as f:
                    reader = csv.reader(f)
                    next(reader, None)
                    for row in reader:
                        if len(row) >= 2 and row[0].strip().lower() == "rad_score":
                            try:
                                score_val = float(row[1])
                                break
                            except ValueError:
                                pass
    if score_val is None:
        return {
            "status": "skipped",
            "error": "未提供 Rad-Score（依赖上一步「计算 Rad-Score」或特征 CSV）。请先完成特征提取与 Rad-Score 计算，或安装 pyradiomics 后重跑流程。",
        }
    try:
        score_val = float(score_val)
    except (TypeError, ValueError):
        return {
            "status": "skipped",
            "error": "Rad-Score 非数值（上一步可能未产出）。请先成功执行「计算 Rad-Score」或安装 pyradiomics 后重跑。",
            "rad_score_raw": str(score_val)[:50],
        }

    out_p = Path(output_path) if (output_path and str(output_path).strip()) else None
    if out_p is None and rad_score_csv_path and str(rad_score_csv_path).strip():
        p_csv = Path(rad_score_csv_path).resolve()
        if p_csv.is_file():
            out_p = p_csv.parent / "radiomics_rad_score.png"
    if out_p is None:
        out_p = Path("radiomics_rad_score.png").resolve()
    out_p = out_p.resolve()
    if out_p.is_dir():
        out_p = out_p / "radiomics_rad_score.png"
    if out_p.suffix.lower() not in (".png", ".jpg", ".jpeg"):
        out_p = out_p.parent / (out_p.name + ".png")
    out_p.parent.mkdir(parents=True, exist_ok=True)

    try:
        # Mock reference: a few samples around reference_mean
        np.random.seed(42)
        ref_scores = np.random.normal(ref_mean, ref_std, 50)
        fig, ax = plt.subplots(1, 1, figsize=(6, 4))
        ax.hist(ref_scores, bins=12, alpha=0.6, color="steelblue", edgecolor="white", label="Reference (mock)")
        ax.axvline(score_val, color="red", linewidth=2, label=f"Patient Rad-Score = {score_val:.3f}")
        ax.set_xlabel("Rad-Score")
        ax.set_ylabel("Count")
        ax.set_title("Rad-Score vs Reference Population")
        ax.legend()
        ax.grid(True, alpha=0.3)
        plt.tight_layout()
        plt.savefig(str(out_p), dpi=100, bbox_inches="tight")
        plt.close()

        return {
            "status": "success",
            "output_path": str(out_p),
            "rad_score": float(score_val),
        }
    except Exception as e:
        logger.exception("plot_rad_score failed: %s", e)
        return {"status": "error", "error": str(e), "rad_score": score_val}
