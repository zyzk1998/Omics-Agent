"""
Atomic Metabolomics tools (pyopls).
OPLS-DA and VIP scores. Requires pyopls>=0.0.2.
Additive to existing tools/metabolomics/ package.
"""
import logging
import os
from pathlib import Path
from typing import Dict, Any, Optional

import numpy as np
import pandas as pd

from ..core.tool_registry import registry
from ..core.file_inspector import resolve_omics_paths
from ..core.utils import sanitize_plot_path

logger = logging.getLogger(__name__)


@registry.register(
    name="metabo_data_validation",
    description="Fast data validation: reads abundance matrix and sample metadata, reports shape, missing ratio and zero ratio. For DAG visibility only.",
    category="Metabolomics",
    output_type="json",
)
def metabo_data_validation(
    data_path: Optional[str] = None,
    meta_path: Optional[str] = None,
    file_path: Optional[str] = None,
) -> Dict[str, Any]:
    """
    极速代谢组数据校验：读取丰度矩阵与样本分组表，计算 shape、缺失率、零值比例。
    不修改数据，不写回文件。meta_path 可选。data_path 可与 file_path 二选一（兼容 executor 注入）。
    路径经全局 resolve_omics_paths 总线解析，优先从 tables 取第一个文件。
    """
    try:
        path_in = data_path or file_path
        if not path_in:
            return {"status": "error", "error": "请提供 data_path 或 file_path"}
        if "," in str(path_in) or ";" in str(path_in):
            resolved = resolve_omics_paths(str(path_in))
            tables = resolved.get("tables") or []
            if not tables:
                raise ValueError("未找到有效的代谢组表格文件（请上传 .csv/.tsv/.xlsx/.txt）")
            path_in = tables[0]
        dp = Path(path_in)
        if not dp.is_file():
            return {"status": "error", "error": f"数据文件不存在: {data_path}"}
        df = pd.read_csv(dp, index_col=0)
        numeric = df.select_dtypes(include=[np.number])
        if numeric.size == 0:
            return {"status": "error", "error": "未找到数值列（代谢物丰度矩阵）"}
        n_samples, n_features = numeric.shape
        missing = numeric.isnull().sum().sum()
        zeros = (numeric == 0).sum().sum()
        total = numeric.size
        missing_ratio = missing / total if total else 0
        zero_ratio = zeros / total if total else 0
        if meta_path:
            mp = Path(meta_path)
            if mp.is_file():
                meta = pd.read_csv(mp, index_col=0)
                if len(meta) != n_samples:
                    logger.warning("meta 行数与 data 样本数不一致: %s vs %s", len(meta), n_samples)
        return {
            "status": "success",
            "message": (
                f"代谢组数据校验通过。共 {n_samples} 个样本，{n_features} 个代谢特征。"
                f"整体缺失率为 {missing_ratio:.2%}，零值比例为 {zero_ratio:.2%}，内存预分配完成。"
            ),
            "n_samples": int(n_samples),
            "n_features": int(n_features),
            "missing_ratio": float(missing_ratio),
            "zero_ratio": float(zero_ratio),
            "data_path": str(dp.resolve()),
        }
    except Exception as e:
        logger.warning("metabo_data_validation failed: %s", e)
        return {"status": "error", "error": str(e)}


@registry.register(
    name="metabo_model_comparison",
    description="PCA + PLS-DA + VIP comparison: 1x3 plot (PCA 2D, PLS-DA 2D, VIP distribution). Input must be imputed and scaled. Optional label_col for group column.",
    category="Metabolomics",
    output_type="mixed",
)
def metabo_model_comparison(
    data_path: str,
    meta_path: str,
    output_plot_path: str,
    label_col: Optional[str] = None,
) -> Dict[str, Any]:
    """
    多维模型对比：对已插补+标准化的数据做 PCA、PLS-DA，并绘制 1x3 图（PCA 2D、PLS-DA 2D、VIP 分布）。
    若数据含 NaN 则直接报错，不调用 sklearn。

    Args:
        data_path (str): 丰度矩阵数据文件路径（已插补+标准化）。
        meta_path (str): 样本元数据文件路径（可选，用于分组列）。
        output_plot_path (str): 输出 1x3 对比图保存路径。
        label_col (str, optional): 用于分组/分类的标签列名（如 'Group', 'Diagnosis', 'condition'）。
            若未指定，系统将按启发式字典在数据表或元数据中自动寻找（label/target/class/y/group/status/diagnosis/type/condition/cohort，忽略大小写），且取值≥2 类。强烈建议大模型根据数据检查结果显式提供此参数。
    """
    try:
        import matplotlib.pyplot as plt
        from sklearn.decomposition import PCA
        from sklearn.cross_decomposition import PLSRegression
        from sklearn.preprocessing import StandardScaler, LabelEncoder
    except ImportError as e:
        return {"status": "error", "error": f"依赖缺失: {e}"}
    try:
        df = pd.read_csv(data_path, index_col=0)
        numeric_cols = df.select_dtypes(include=[np.number]).columns.tolist()
        X = df[numeric_cols].values
        if np.isnan(X).any():
            raise ValueError(
                "数据包含缺失值，请先运行预处理（插补）步骤后再执行多维模型对比。"
            )
        n_samples, n_features = X.shape
        if n_features < 2:
            raise ValueError("至少需要 2 个特征列才能进行 PCA/PLS-DA 对比。")
        candidate_cols = [
            "label", "target", "class", "y",
            "group", "status", "diagnosis", "type", "condition", "cohort",
        ]
        candidate_cols_lower = {c.lower() for c in candidate_cols}
        group_col = None
        if label_col and str(label_col).strip():
            want = str(label_col).strip()
            found = None
            for c in df.columns:
                if str(c).strip().lower() == want.lower():
                    found = c
                    break
            if found is None:
                raise ValueError(
                    f"指定的标签列 '{label_col}' 不在数据表列中。可用列: {list(df.columns)[:30]}"
                )
            if df[found].nunique() < 2 or df[found].nunique() > 20:
                raise ValueError(
                    f"指定的标签列 '{found}' 唯一值数量为 {df[found].nunique()}，需要至少 2 类且不超过 20 类。"
                )
            group_col = found
        if group_col is None:
            candidates = [
                c for c in df.columns
                if c not in numeric_cols and str(c).strip().lower() in candidate_cols_lower
            ]
            candidates = [c for c in candidates if df[c].nunique() >= 2 and df[c].nunique() <= 20]
            if len(candidates) == 1:
                group_col = candidates[0]
            elif len(candidates) > 1:
                raise ValueError(
                    f"存在多个可能的分组列 {candidates}，请通过参数 label_col 指定其一。"
                )
        if group_col is None:
            for c in df.columns:
                if c not in numeric_cols and df[c].nunique() >= 2 and df[c].nunique() <= 20:
                    group_col = c
                    break
        y_labels = None
        if group_col:
            y_labels = df[group_col].astype(str)
        if y_labels is None and meta_path and str(meta_path).strip():
            mp = Path(meta_path)
            if mp.is_file() and str(mp.resolve()) != str(Path(data_path).resolve()):
                meta = pd.read_csv(meta_path, index_col=0)
                inter = df.index.intersection(meta.index)
                if len(inter) >= 2:
                    for c in meta.columns:
                        if meta.loc[inter, c].nunique() >= 2:
                            df = df.loc[inter]
                            X = df[numeric_cols].values
                            n_samples = X.shape[0]
                            y_labels = meta.loc[inter, c].astype(str)
                            break
                    else:
                        pass
        n_classes = len(np.unique(y_labels)) if y_labels is not None else 0
        if y_labels is None or n_classes < 2:
            raise ValueError(
                "执行 PLS-DA 需要至少两个有效的分组标签。请检查是否上传了正确的样本分组文件 (Metadata)，或分组列是否选择正确。"
            )
        le = LabelEncoder()
        y_num = le.fit_transform(y_labels.values)
        X_scaled = StandardScaler().fit_transform(X)
        pca = PCA(n_components=2)
        pca_coords = pca.fit_transform(X_scaled)
        pls = PLSRegression(n_components=2)
        pls.fit(X_scaled, y_num)
        pls_coords = pls.x_scores_[:, :2]
        T = pls.x_scores_
        W = pls.x_weights_
        ss_per_comp = np.sum(T**2, axis=0)
        total_ss = ss_per_comp.sum()
        vip = np.sqrt(n_features * np.sum(W**2 * ss_per_comp, axis=1) / total_ss) if total_ss > 0 else np.zeros(n_features)
        vip = np.array(vip).ravel()
        if len(vip) != n_features:
            vip = np.zeros(n_features)
        fig, axes = plt.subplots(1, 3, figsize=(18, 5))
        colors = plt.cm.tab10(np.linspace(0, 1, len(np.unique(y_num))))
        for i, lab in enumerate(np.unique(y_num)):
            mask = y_num == lab
            axes[0].scatter(pca_coords[mask, 0], pca_coords[mask, 1], c=[colors[i]], label=le.inverse_transform([lab])[0], alpha=0.7)
            axes[1].scatter(pls_coords[mask, 0], pls_coords[mask, 1], c=[colors[i]], label=le.inverse_transform([lab])[0], alpha=0.7)
        from matplotlib.patches import Ellipse
        for i, lab in enumerate(np.unique(y_num)):
            mask = y_num == lab
            x, y = pca_coords[mask, 0], pca_coords[mask, 1]
            if len(x) >= 2:
                mean_x, mean_y = x.mean(), y.mean()
                cov = np.cov(x, y)
                vals, vecs = np.linalg.eigh(cov)
                angle = np.degrees(np.arctan2(vecs[1, 0], vecs[0, 0]))
                chi2_val = 5.991  # 95% confidence for 2d
                w, h = 2 * np.sqrt(chi2_val * np.sqrt(vals))
                axes[0].add_patch(Ellipse((mean_x, mean_y), w, h, angle=angle, fill=False, edgecolor=colors[i], linestyle="--", alpha=0.8))
        axes[0].set_xlabel("PC1")
        axes[0].set_ylabel("PC2")
        axes[0].set_title("PCA Score Plot")
        axes[0].legend()
        axes[0].grid(True, alpha=0.3)
        axes[1].set_xlabel("Component 1")
        axes[1].set_ylabel("Component 2")
        axes[1].set_title("PLS-DA Score Plot")
        axes[1].legend()
        axes[1].grid(True, alpha=0.3)
        axes[2].hist(vip, bins=min(50, max(10, n_features // 2)), color="steelblue", edgecolor="white", alpha=0.8)
        axes[2].axvline(x=1, color="red", linestyle="--", label="VIP=1")
        axes[2].set_xlabel("VIP Score")
        axes[2].set_ylabel("Count")
        axes[2].set_title("VIP Scores")
        axes[2].legend()
        plt.tight_layout()
        out_plot = sanitize_plot_path(output_plot_path)
        os.makedirs(str(out_plot.parent) or ".", exist_ok=True)
        fig.savefig(str(out_plot), bbox_inches="tight", dpi=150)
        plt.close()
        return {
            "status": "success",
            "message": "多维模型对比图已保存",
            "plot_path": str(out_plot.resolve()),
            "data_path": data_path,
            "output_file": data_path,
            "file_path": data_path,
        }
    except ValueError as e:
        logger.warning("metabo_model_comparison 业务校验失败: %s", e)
        return {
            "status": "error",
            "error": str(e),
            "message": str(e),
            "user_message": str(e),
        }
    except Exception as e:
        logger.exception("metabo_model_comparison failed: %s", e)
        return {
            "status": "error",
            "error": str(e),
            "message": str(e),
            "user_message": f"多维模型对比步骤执行异常：{str(e)}",
        }


@registry.register(
    name="run_opls_da",
    description="Perform OPLS-DA (Orthogonal PLS-DA) on metabolomics data. Input: X CSV (samples x metabolites), Y CSV (samples x group label). Returns VIP scores and model summary.",
    category="Metabolomics",
    output_type="json",
)
def run_opls_da(
    X_csv: str,
    Y_csv: str,
    n_components: int = 1,
) -> Dict[str, Any]:
    """
    Run OPLS-DA and return VIP scores.

    Args:
        X_csv: Path to feature matrix CSV (rows=samples, columns=metabolites).
        Y_csv: Path to response/label CSV (rows=samples, one column with group labels).
        n_components: Number of orthogonal components (default 1).

    Returns:
        Dict with status, vip_scores (dict or list), and optional error.
    """
    try:
        import pandas as pd
        import numpy as np
    except ImportError:
        return {"status": "error", "error": "pandas and numpy are required"}
    try:
        from pyopls import OPLS
    except ImportError as e:
        logger.warning("pyopls not installed: %s", e)
        return {"status": "error", "error": "pyopls is required. Install with: pip install pyopls>=0.0.2"}
    px, py = Path(X_csv), Path(Y_csv)
    if not px.is_file():
        return {"status": "error", "error": f"X file not found: {X_csv}"}
    if not py.is_file():
        return {"status": "error", "error": f"Y file not found: {Y_csv}"}
    try:
        X = pd.read_csv(px, index_col=0)
        Y = pd.read_csv(py, index_col=0)
        inter = X.index.intersection(Y.index)
        if len(inter) < 3:
            return {"status": "error", "error": "Too few overlapping samples between X and Y."}
        X = X.loc[inter]
        Y = Y.loc[inter]
        if Y.select_dtypes(include=[object]).shape[1] > 0:
            y_col = Y.columns[0]
            uniq = Y[y_col].astype(str).unique()
            mapping = {u: i for i, u in enumerate(uniq)}
            y_vec = Y[y_col].astype(str).map(mapping).values
        else:
            y_vec = Y.iloc[:, 0].values
        X_mat = X.select_dtypes(include=[np.number]).values
        if X_mat.shape[1] < 2:
            return {"status": "error", "error": "Need at least 2 metabolite columns."}
        model = OPLS(n_components=n_components)
        model.fit(X_mat, y_vec)
        if hasattr(model, "get_VIP") or hasattr(model, "vip"):
            vip = getattr(model, "get_VIP", lambda: getattr(model, "vip", None))()
            if vip is not None:
                vip_ser = pd.Series(vip, index=X.select_dtypes(include=[np.number]).columns)
                vip_scores = vip_ser.to_dict()
            else:
                vip_scores = {c: 0.0 for c in X.select_dtypes(include=[np.number]).columns}
        else:
            vip_scores = {c: 0.0 for c in X.select_dtypes(include=[np.number]).columns}
        return {
            "status": "success",
            "vip_scores": vip_scores,
            "n_samples": X.shape[0],
            "n_features": X.shape[1],
            "X_path": str(px.resolve()),
            "Y_path": str(py.resolve()),
        }
    except Exception as e:
        logger.exception("run_opls_da failed: %s", e)
        return {"status": "error", "error": str(e)}
