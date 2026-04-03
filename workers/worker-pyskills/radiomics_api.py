"""
影像组学 TaaS：PyRadiomics 批量提取、LASSO 筛选、逻辑回归 + ROC。
仅读写 /app/uploads 下路径；产物写入 /app/uploads/results/radiomics/<run_id>/。
对外 JSON 仅含 /uploads/results/radiomics/... 相对 URL（与 nginx 挂载一致）。
"""
from __future__ import annotations

import csv
import logging
import os
import uuid
from pathlib import Path
from typing import Any, Dict, List, Optional

from fastapi import APIRouter, HTTPException
from pydantic import BaseModel, Field

logger = logging.getLogger(__name__)

router = APIRouter(prefix="/api/radiomics", tags=["radiomics"])

UPLOADS_ROOT = Path("/app/uploads")
RESULTS_RAD = Path("/app/uploads/results/radiomics")
PUBLIC_BASE = "/uploads/results/radiomics"


def _rad_out_dir(run_id: str) -> Path:
    d = RESULTS_RAD / run_id
    d.mkdir(parents=True, exist_ok=True)
    return d


def _public_url(run_id: str, filename: str) -> str:
    return f"{PUBLIC_BASE}/{run_id}/{filename}"


def secure_resolve_path(path: str) -> Path:
    """
    解析客户端传入路径，禁止目录穿越：规范化后必须落在 /app/uploads/ 之下。
    使用 os.path.abspath / realpath；越界返回 403。
    """
    if not path or not str(path).strip():
        raise HTTPException(status_code=400, detail="路径为空")
    raw = str(path).strip()
    expanded = os.path.expanduser(raw)
    if os.path.isabs(expanded):
        candidate = os.path.abspath(expanded)
    else:
        candidate = os.path.abspath(os.path.join(str(UPLOADS_ROOT), expanded))
    try:
        resolved = os.path.realpath(candidate)
    except OSError as e:
        logger.warning("secure_resolve_path realpath failed: %s", e)
        raise HTTPException(status_code=403, detail="无效路径，拒绝访问") from e
    normalized = resolved.replace("\\", "/")
    root_norm = os.path.realpath(str(UPLOADS_ROOT)).replace("\\", "/")
    if root_norm.endswith("/"):
        root_norm = root_norm.rstrip("/")
    prefix = root_norm + "/"
    if normalized != root_norm and not normalized.startswith(prefix):
        raise HTTPException(
            status_code=403,
            detail="拒绝访问：路径必须位于 /app/uploads/ 目录下",
        )
    return Path(resolved)


class RadiomicsExtractRequest(BaseModel):
    """批量提取：manifest CSV（case_id,image_path,mask_path）或单对影像。"""

    manifest_path: Optional[str] = Field(default=None, description="manifest CSV 绝对路径（/app/uploads/...）")
    image_path: Optional[str] = None
    mask_path: Optional[str] = None
    case_id: Optional[str] = Field(default=None, description="单对模式下的病例 ID")


class RadiomicsLassoRequest(BaseModel):
    features_csv_path: str
    label_col: str
    cv_folds: int = Field(default=5, ge=3, le=15)
    random_state: int = 42
    max_iter: int = Field(default=2000, ge=100, le=20000)


class RadiomicsModelRequest(BaseModel):
    features_csv_path: str
    label_col: str
    test_size: float = Field(default=0.25, gt=0.05, lt=0.5)
    random_state: int = 42


def _read_manifest_rows(manifest_path: Path) -> List[Dict[str, str]]:
    rows: List[Dict[str, str]] = []
    with open(manifest_path, newline="", encoding="utf-8", errors="replace") as f:
        reader = csv.DictReader(f)
        if not reader.fieldnames:
            raise HTTPException(status_code=400, detail="manifest 无表头")

        def g(row: Dict[str, str], *keys: str) -> str:
            for k in keys:
                v = row.get(k)
                if v is not None and str(v).strip():
                    return str(v).strip()
            return ""

        for raw in reader:
            row = {str(k).strip().lower(): (str(v).strip() if v is not None else "") for k, v in raw.items()}
            cid = g(row, "case_id", "id", "sample_id", "patient_id")
            img = g(row, "image_path", "image", "img")
            msk = g(row, "mask_path", "mask")
            if not img or not msk:
                continue
            if not cid:
                cid = f"case_{len(rows)}"
            rows.append({"case_id": cid, "image_path": img, "mask_path": msk})
    if not rows:
        raise HTTPException(status_code=400, detail="manifest 无有效行（需 image_path + mask_path）")
    return rows


@router.post("/extract")
async def radiomics_extract(req: RadiomicsExtractRequest) -> Dict[str, Any]:
    try:
        import pandas as pd
        from radiomics import featureextractor
    except ImportError as e:
        raise HTTPException(status_code=503, detail=f"PyRadiomics 未安装: {e}") from e

    run_id = uuid.uuid4().hex[:16]
    out_dir = _rad_out_dir(run_id)
    out_csv = out_dir / "features.csv"

    pairs: List[Dict[str, str]] = []
    if req.manifest_path:
        mp = secure_resolve_path(req.manifest_path)
        if not mp.is_file():
            raise HTTPException(status_code=400, detail=f"manifest 不存在: {req.manifest_path}")
        pairs = _read_manifest_rows(mp)
    elif req.image_path and req.mask_path:
        cid = (req.case_id or "case_0").strip() or "case_0"
        pairs = [{"case_id": cid, "image_path": req.image_path.strip(), "mask_path": req.mask_path.strip()}]
    else:
        raise HTTPException(status_code=400, detail="请提供 manifest_path，或同时提供 image_path 与 mask_path")

    extractor = featureextractor.RadiomicsFeatureExtractor()
    all_rows: List[Dict[str, Any]] = []

    for pr in pairs:
        img_p = secure_resolve_path(pr["image_path"])
        msk_p = secure_resolve_path(pr["mask_path"])
        if not img_p.is_file():
            raise HTTPException(status_code=400, detail=f"影像不存在: {img_p}")
        if not msk_p.is_file():
            raise HTTPException(status_code=400, detail=f"掩膜不存在: {msk_p}")
        try:
            result = extractor.execute(str(img_p), str(msk_p))
        except Exception as ex:
            logger.exception("PyRadiomics execute")
            raise HTTPException(status_code=500, detail=f"特征提取失败 ({pr['case_id']}): {ex}") from ex
        feat = {k: v for k, v in result.items() if not str(k).startswith("diagnostics")}
        row: Dict[str, Any] = {"case_id": pr["case_id"]}
        row.update(feat)
        all_rows.append(row)

    pd.DataFrame(all_rows).to_csv(out_csv, index=False)
    return {
        "status": "success",
        "message": f"已提取 {len(all_rows)} 例特征",
        "run_id": run_id,
        "features_csv_url": _public_url(run_id, "features.csv"),
        "features_csv_path": str(out_csv),
        "n_cases": len(all_rows),
        "n_feature_cols": len(all_rows[0]) - 1 if all_rows else 0,
    }


@router.post("/lasso")
async def radiomics_lasso(req: RadiomicsLassoRequest) -> Dict[str, Any]:
    try:
        import numpy as np
        import pandas as pd
        from sklearn.linear_model import LogisticRegressionCV
        from sklearn.preprocessing import StandardScaler
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError as e:
        raise HTTPException(status_code=503, detail=f"依赖缺失: {e}") from e

    fp = secure_resolve_path(req.features_csv_path)
    if not fp.is_file():
        raise HTTPException(status_code=400, detail=f"特征 CSV 不存在: {req.features_csv_path}")

    df = pd.read_csv(fp)
    want = str(req.label_col).strip().lower()
    label_col = None
    for c in df.columns:
        if str(c).strip().lower() == want:
            label_col = c
            break
    if label_col is None:
        raise HTTPException(status_code=400, detail=f"找不到标签列: {req.label_col}")

    y_raw = df[label_col]
    try:
        y = y_raw.astype(int).values
    except Exception:
        y = pd.Categorical(y_raw).codes

    X_df = df.drop(columns=[label_col]).select_dtypes(include=[np.number])
    if X_df.shape[1] < 2:
        raise HTTPException(status_code=400, detail="数值特征列不足，无法进行 LASSO")

    X = X_df.values.astype(float)
    X = np.nan_to_num(X, nan=np.nanmedian(X, axis=0), posinf=0.0, neginf=0.0)
    scaler = StandardScaler()
    Xs = scaler.fit_transform(X)

    cs = np.logspace(-3, 1, 12)
    try:
        n_splits = max(2, min(req.cv_folds, max(2, len(y) - 1)))
        model = LogisticRegressionCV(
            penalty="l1",
            solver="saga",
            cv=n_splits,
            Cs=list(cs),
            max_iter=req.max_iter,
            random_state=req.random_state,
            scoring="roc_auc",
            n_jobs=None,
        )
        model.fit(Xs, y)
    except Exception as ex:
        logger.exception("LASSO-CV")
        raise HTTPException(status_code=500, detail=f"LASSO 训练失败: {ex}") from ex

    coef = np.asarray(model.coef_).ravel()
    names = list(X_df.columns)
    order = np.argsort(np.abs(coef))[::-1]
    selected_idx = [i for i in order if abs(coef[i]) > 1e-8]
    if not selected_idx:
        selected_idx = order[: min(20, len(order))].tolist()

    selected_names = [names[i] for i in selected_idx]
    out_cols = ["case_id"] if "case_id" in df.columns else []
    out_cols = [c for c in out_cols if c in df.columns]
    keep = out_cols + [c for c in selected_names if c in df.columns] + [label_col]
    keep = list(dict.fromkeys(keep))
    sub = df[[c for c in keep if c in df.columns]]

    run_id = uuid.uuid4().hex[:16]
    out_dir = _rad_out_dir(run_id)
    sel_path = out_dir / "selected_features.csv"
    sub.to_csv(sel_path, index=False)

    png_path = out_dir / "lasso_path.png"
    try:
        top_k = min(40, len(names))
        idx = order[:top_k]
        labels = [names[i][:48] for i in idx]
        vals = coef[idx]
        fig_h = max(4.0, 0.22 * top_k)
        fig, ax = plt.subplots(figsize=(8, fig_h))
        colors = ["#2563eb" if v >= 0 else "#dc2626" for v in vals]
        ax.barh(labels[::-1], vals[::-1], color=colors[::-1])
        ax.axvline(0, color="#64748b", lw=0.8)
        ax.set_xlabel("L1-logistic coefficient")
        ax.set_title("LASSO 非零系数（按 |coef| 排序，至多 40 项）")
        ax.grid(True, axis="x", alpha=0.3)
        fig.tight_layout()
        fig.savefig(str(png_path), dpi=150, bbox_inches="tight")
        plt.close(fig)
    except Exception as ex:
        logger.warning("lasso_path plot failed: %s", ex)
        png_path.write_bytes(b"")

    return {
        "status": "success",
        "message": f"LASSO 筛选保留 {len(selected_names)} 个特征",
        "run_id": run_id,
        "selected_features_csv_url": _public_url(run_id, "selected_features.csv"),
        "lasso_path_png_url": _public_url(run_id, "lasso_path.png"),
        "selected_features_csv_path": str(sel_path),
        "lasso_path_png_path": str(out_dir / "lasso_path.png"),
        "n_selected": len(selected_names),
        "selected_feature_names": selected_names[:50],
    }


@router.post("/model")
async def radiomics_model(req: RadiomicsModelRequest) -> Dict[str, Any]:
    try:
        import numpy as np
        import pandas as pd
        from sklearn.linear_model import LogisticRegression
        from sklearn.metrics import accuracy_score, roc_auc_score, roc_curve
        from sklearn.model_selection import train_test_split
        from sklearn.preprocessing import StandardScaler
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError as e:
        raise HTTPException(status_code=503, detail=f"依赖缺失: {e}") from e

    fp = secure_resolve_path(req.features_csv_path)
    if not fp.is_file():
        raise HTTPException(status_code=400, detail=f"特征 CSV 不存在: {req.features_csv_path}")

    df = pd.read_csv(fp)
    want = str(req.label_col).strip().lower()
    label_col = None
    for c in df.columns:
        if str(c).strip().lower() == want:
            label_col = c
            break
    if label_col is None:
        raise HTTPException(status_code=400, detail=f"找不到标签列: {req.label_col}")

    y_raw = df[label_col]
    try:
        y = y_raw.astype(int).values
    except Exception:
        y = pd.Categorical(y_raw).codes

    if len(np.unique(y)) < 2:
        raise HTTPException(status_code=400, detail="标签需至少两类")

    X_df = df.drop(columns=[label_col]).select_dtypes(include=[np.number])
    if X_df.shape[1] < 1:
        raise HTTPException(status_code=400, detail="无数值特征列")

    X = np.nan_to_num(X_df.values.astype(float), nan=np.nanmedian(X_df.values.astype(float), axis=0))
    X_train, X_test, y_train, y_test = train_test_split(
        X, y, test_size=req.test_size, random_state=req.random_state, stratify=y
    )
    scaler = StandardScaler()
    X_train = scaler.fit_transform(X_train)
    X_test = scaler.transform(X_test)

    clf = LogisticRegression(max_iter=2000, random_state=req.random_state)
    clf.fit(X_train, y_train)
    y_prob = clf.predict_proba(X_test)[:, 1]
    y_pred = clf.predict(X_test)
    acc = float(accuracy_score(y_test, y_pred))
    try:
        auc_v = float(roc_auc_score(y_test, y_prob))
    except Exception:
        auc_v = float("nan")

    run_id = uuid.uuid4().hex[:16]
    out_dir = _rad_out_dir(run_id)
    metrics_path = out_dir / "metrics.csv"
    with open(metrics_path, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["metric", "value"])
        w.writerow(["accuracy", f"{acc:.6f}"])
        w.writerow(["roc_auc", f"{auc_v:.6f}" if auc_v == auc_v else ""])
        w.writerow(["n_train", len(y_train)])
        w.writerow(["n_test", len(y_test)])
        w.writerow(["n_features", X_df.shape[1]])

    fig, ax = plt.subplots(figsize=(6, 5))
    try:
        fpr, tpr, _ = roc_curve(y_test, y_prob)
        ax.plot(fpr, tpr, lw=2, label=f"LR (AUC = {auc_v:.3f})" if auc_v == auc_v else "LR")
    except Exception as ex:
        logger.warning("roc_curve: %s", ex)
        ax.text(0.1, 0.5, "ROC 不可用（标签或预测异常）")
    ax.plot([0, 1], [0, 1], "k--", lw=1)
    ax.set_xlabel("FPR")
    ax.set_ylabel("TPR")
    ax.set_title("ROC — Radiomics model")
    ax.legend(loc="lower right")
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    roc_path = out_dir / "roc_curve.png"
    fig.savefig(str(roc_path), dpi=150, bbox_inches="tight")
    plt.close(fig)

    return {
        "status": "success",
        "message": "模型评估完成",
        "run_id": run_id,
        "metrics_csv_url": _public_url(run_id, "metrics.csv"),
        "roc_curve_png_url": _public_url(run_id, "roc_curve.png"),
        "metrics_csv_path": str(metrics_path),
        "roc_curve_png_path": str(roc_path),
        "accuracy": acc,
        "roc_auc": auc_v,
    }
