"""
Atomic Proteomics tools (statsmodels, pandas).
Missing value imputation and batch T-tests for protein data.
"""
import logging
from pathlib import Path
from typing import Dict, Any, Optional

from ..core.tool_registry import registry

logger = logging.getLogger(__name__)


@registry.register(
    name="impute_proteomics_missing",
    description="Impute missing values in a proteomics abundance DataFrame. Methods: 'min', 'mean', 'median', 'knn'.",
    category="Proteomics",
    output_type="json",
)
def impute_proteomics_missing(
    df: Any,
    method: str = "min",
) -> Dict[str, Any]:
    """
    Impute missing values in proteomics data.
    df can be a file path (CSV) or a dict that can be converted to DataFrame (e.g. from previous step).

    Args:
        df: Path to CSV (index=proteins, columns=samples) or DataFrame-like dict.
        method: 'min' (per-column min), 'mean', 'median', or 'knn' (if available).

    Returns:
        Dict with status, imputed dataframe (as records), shape, and optional error.
    """
    try:
        import pandas as pd
        import numpy as np
    except ImportError:
        return {"status": "error", "error": "pandas and numpy are required"}
    try:
        if isinstance(df, str):
            path = Path(df)
            if path.suffix.lower() in (".csv", ".txt"):
                data = pd.read_csv(path, index_col=0)
            else:
                return {"status": "error", "error": "Unsupported file format; use CSV."}
        elif isinstance(df, dict):
            data = pd.DataFrame(df)
        else:
            return {"status": "error", "error": "df must be a file path (str) or a dict"}
        numeric = data.select_dtypes(include=[np.number])
        if numeric.empty:
            return {"status": "error", "error": "No numeric columns found."}
        if method == "min":
            fill = numeric.min().min()
            if pd.isna(fill):
                fill = 0
            imputed = numeric.fillna(fill)
        elif method == "mean":
            imputed = numeric.fillna(numeric.mean())
        elif method == "median":
            imputed = numeric.fillna(numeric.median())
        elif method == "knn":
            try:
                from sklearn.impute import KNNImputer
                imputer = KNNImputer(n_neighbors=5)
                imputed = pd.DataFrame(imputer.fit_transform(numeric), index=numeric.index, columns=numeric.columns)
            except ImportError:
                return {"status": "error", "error": "knn requires scikit-learn. Using median instead."}
                imputed = numeric.fillna(numeric.median())
        else:
            return {"status": "error", "error": f"Unknown method: {method}. Use min, mean, median, or knn."}
        n_imputed = numeric.isna().sum().sum()
        return {
            "status": "success",
            "dataframe": imputed.head(2000).to_dict(orient="index"),
            "shape": {"rows": imputed.shape[0], "cols": imputed.shape[1]},
            "method": method,
            "n_values_imputed": int(n_imputed),
        }
    except Exception as e:
        logger.exception("impute_proteomics_missing failed: %s", e)
        return {"status": "error", "error": str(e)}


@registry.register(
    name="run_proteomics_ttest",
    description="Run batch two-sample T-test for protein abundances by group. Input: DataFrame (or CSV path) and group column name.",
    category="Proteomics",
    output_type="json",
)
def run_proteomics_ttest(
    df: Any,
    group_col: str,
) -> Dict[str, Any]:
    """
    Batch T-test: for each protein, compare two groups defined by group_col.

    Args:
        df: Path to CSV or DataFrame-like dict (rows=proteins, columns include sample IDs and optionally group_col if in same table).
        group_col: Column name containing group labels (e.g. 'group' with values 'A','B').

    Returns:
        Dict with status, results (protein, pvalue, fold_change, etc.), and optional error.
    """
    try:
        import pandas as pd
        import numpy as np
        from scipy import stats
    except ImportError as e:
        logger.warning("scipy/pandas not installed: %s", e)
        return {"status": "error", "error": "pandas and scipy are required"}
    try:
        if isinstance(df, str):
            data = pd.read_csv(Path(df), index_col=0)
        elif isinstance(df, dict):
            data = pd.DataFrame(df)
        else:
            return {"status": "error", "error": "df must be a file path (str) or a dict"}
        if group_col not in data.columns:
            return {"status": "error", "error": f"Column '{group_col}' not found. Columns: {list(data.columns)}"}
        groups = data[group_col].dropna().unique()
        if len(groups) != 2:
            return {"status": "error", "error": f"Exactly two groups required; found: {groups.tolist()}"}
        g1, g2 = groups[0], groups[1]
        numeric_cols = data.select_dtypes(include=[np.number]).columns.tolist()
        if not numeric_cols:
            return {"status": "error", "error": "No numeric columns for T-test."}
        X = data[numeric_cols]
        idx1 = data[group_col] == g1
        idx2 = data[group_col] == g2
        pvals, fc = [], []
        for col in numeric_cols:
            a, b = X.loc[idx1, col], X.loc[idx2, col]
            a, b = a.dropna(), b.dropna()
            if len(a) < 2 or len(b) < 2:
                pvals.append(np.nan)
                fc.append(np.nan)
                continue
            t, p = stats.ttest_ind(a, b)
            m1, m2 = a.mean(), b.mean()
            fold = m2 / m1 if m1 != 0 else np.nan
            pvals.append(p)
            fc.append(fold)
        out = pd.DataFrame({
            "protein": numeric_cols,
            "pvalue": pvals,
            "fold_change": fc,
            "group1": g1,
            "group2": g2,
        })
        out = out.dropna(subset=["pvalue"])
        return {
            "status": "success",
            "results": out.head(2000).to_dict(orient="records"),
            "group_col": group_col,
            "n_proteins": len(out),
        }
    except Exception as e:
        logger.exception("run_proteomics_ttest failed: %s", e)
        return {"status": "error", "error": str(e)}
