#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
EIS → DRT（弛豫时间分布）离散 Tikhonov + 非负最小二乘 — 供 ``drt_tool.eis_drt_analysis`` 子进程调用。

随仓库分发；依赖 numpy / scipy / pandas / matplotlib（主 API 镜像已包含）。

退出码：0 成功；1 失败（错误信息在 stderr）。
"""
from __future__ import annotations

import argparse
import json
import re
import sys
import traceback
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
import pandas as pd


def _norm_header(name: str) -> str:
    s = str(name).strip().replace("\ufeff", "")
    s = s.replace("′", "'").replace("’", "'").replace("‘", "'")
    s = s.replace("″", '"').replace('"', '"')
    s = s.lower()
    s = re.sub(r"\s+", "", s)
    return s


def _read_csv_flexible(path: Path) -> pd.DataFrame:
    """尝试编码与分隔符；兼容技能包内中文表头、分号分隔等。"""
    last_err: Optional[Exception] = None
    for enc in ("utf-8-sig", "utf-8", "gbk", "gb2312"):
        for sep in (",", ";", "\t"):
            try:
                kwargs = dict(encoding=enc, sep=sep, engine="python")
                try:
                    kwargs["on_bad_lines"] = "skip"
                    df = pd.read_csv(path, **kwargs)
                except TypeError:
                    df = pd.read_csv(path, encoding=enc, sep=sep, engine="python")
                if df.shape[1] >= 3 and df.shape[0] >= 3:
                    df.columns = [str(c).strip() for c in df.columns]
                    return df
            except Exception as e:
                last_err = e
                continue
    try:
        df = pd.read_csv(path)
        df.columns = [str(c).strip() for c in df.columns]
        return df
    except Exception as e:
        raise ValueError(f"无法解析 CSV（{path}）：{last_err or e}") from e


def _column_classify(name: str) -> Optional[str]:
    """返回 'freq' | 'zreal' | 'zimag' | None（保守匹配，避免误伤）。"""
    raw = str(name).strip()
    rl = raw.lower()
    n = _norm_header(name)
    # 虚部：双撇号 / Zimag / 虚部（必须先于单撇号 Z'，否则 Z'' 会误匹配实部）
    if (
        "z''" in rl
        or 'z"' in rl
        or "″" in raw
        or "zimag" in n
        or "z_imag" in n
        or "-z_imag" in n
        or "虚部" in raw
    ):
        return "zimag"
    # 实部：仅单撇 Z'（排除已判为 Z''）、Zreal、实部
    rp = raw.replace("′", "'")
    if (
        ("z'" in rp.lower() and "z''" not in rp.lower())
        or "zreal" in n
        or "z_real" in n
        or "实部" in raw
    ):
        return "zreal"
    # 频率
    if (
        "freq" in n
        or "频率" in raw
        or n.endswith("/hz")
        or n.endswith("hz")
        or n == "f"
    ):
        return "freq"
    return None


def _guess_numeric_triplet(df: pd.DataFrame) -> Tuple[str, str, str]:
    """三列数值表：用数量级跨度识别频率列，其余两列为阻抗。"""
    numeric_cols: List[Tuple[str, np.ndarray]] = []
    for c in df.columns:
        v = pd.to_numeric(df[c], errors="coerce").to_numpy(dtype=float)
        if np.isfinite(v).sum() >= max(8, len(df) // 3):
            numeric_cols.append((c, v))
    if len(numeric_cols) < 3:
        raise ValueError("数值列不足 3 列，无法启发式推断。")

    best_f: Optional[str] = None
    best_span = -1.0
    for c, v in numeric_cols:
        vv = v[np.isfinite(v)]
        pos = vv[vv > 0]
        if pos.size < 8:
            continue
        span = float(np.nanmax(np.log10(pos)) - np.nanmin(np.log10(pos)))
        # 频率跨多个十倍频；阻抗原侧常见 |Z| 量级相对稳定
        if span > best_span and np.median(pos) < 1e12:
            best_span = span
            best_f = c
    if best_f is None:
        raise ValueError("启发式推断失败（未找到合理的频率列）。")

    others = [c for c, _ in numeric_cols if c != best_f]
    if len(others) < 2:
        raise ValueError("启发式推断失败（阻抗列不足）。")
    # 默认保持 CSV 原始列顺序（多数仪器：Freq, Z', Z''）
    ordered = [c for c in df.columns if c in others][:2]
    if len(ordered) == 2:
        return best_f, ordered[0], ordered[1]
    return best_f, others[0], others[1]


def _find_columns(df: pd.DataFrame) -> Tuple[str, str, str]:
    buckets: Dict[str, List[str]] = {"freq": [], "zreal": [], "zimag": []}
    for col in df.columns:
        cat = _column_classify(col)
        if cat:
            buckets[cat].append(col)
    if len(buckets["freq"]) == 1 and len(buckets["zreal"]) == 1 and len(buckets["zimag"]) == 1:
        return buckets["freq"][0], buckets["zreal"][0], buckets["zimag"][0]

    try:
        return _guess_numeric_triplet(df)
    except ValueError:
        pass

    raise ValueError(
        "无法在 CSV 中识别频率与 Z'、Z'' 列。请使用常见表头（Freq/Hz、Z'/ohm、Z''/ohm 等）。"
        f" 当前列: {list(df.columns)}"
    )


def _estimate_r_inf(freq_hz: np.ndarray, z_real: np.ndarray) -> float:
    idx = np.argsort(freq_hz)
    f = freq_hz[idx]
    zr = z_real[idx]
    k = max(3, len(f) // 10)
    return float(np.mean(zr[-k:]))


def _second_diff_matrix(n: int) -> np.ndarray:
    if n < 3:
        return np.zeros((0, n))
    D = np.zeros((n - 2, n))
    for i in range(n - 2):
        D[i, i] = 1.0
        D[i, i + 1] = -2.0
        D[i, i + 2] = 1.0
    return D


def run_drt(
    freq_hz: np.ndarray,
    z_complex: np.ndarray,
    *,
    n_tau: int = 80,
    regularization_lambda: float = 0.1,
) -> Tuple[np.ndarray, np.ndarray, float, np.ndarray]:
    omega = 2 * np.pi * np.asarray(freq_hz, dtype=float)
    z = np.asarray(z_complex, dtype=complex)

    om_min = max(np.min(omega), 1e-12)
    om_max = max(np.max(omega), om_min * 1.001)
    tau_min = max(1.0 / om_max / 10.0, 1e-15)
    tau_max = 1.0 / om_min * 10.0
    tau = np.logspace(np.log10(tau_min), np.log10(tau_max), n_tau)

    r_inf = _estimate_r_inf(freq_hz, z.real)
    y = z - r_inf

    log_tau = np.log(tau)
    dln = np.gradient(log_tau)

    omega_col = omega.reshape(-1, 1)
    tau_row = tau.reshape(1, -1)
    K = dln.reshape(1, -1) / (1.0 + 1j * omega_col * tau_row)

    Ar = np.vstack([K.real, K.imag])
    yr = np.concatenate([y.real, y.imag])

    D2 = _second_diff_matrix(n_tau)
    lam = float(max(regularization_lambda, 1e-12))
    if D2.shape[0] > 0:
        X = np.vstack([Ar, np.sqrt(lam) * D2])
        yvec = np.concatenate([yr, np.zeros(D2.shape[0])])
    else:
        X = Ar
        yvec = yr

    n_col = X.shape[1]
    lb = np.zeros(n_col)
    ub = np.full(n_col, np.inf)

    try:
        from scipy.optimize import lsq_linear

        sol = lsq_linear(X, yvec, bounds=(lb, ub), method="trf", verbose=0)
        g = sol.x
    except Exception:
        try:
            from scipy.optimize import nnls

            g, _ = nnls(X, yvec)
        except Exception:
            g, _, _, _ = np.linalg.lstsq(X, yvec, rcond=None)
    g = np.maximum(np.asarray(g, dtype=float), 0.0)

    z_fit = r_inf + K @ g
    return tau, g, r_inf, z_fit


def _rmse(a: np.ndarray, b: np.ndarray) -> float:
    return float(np.sqrt(np.mean(np.abs(a - b) ** 2)))


def plot_drt(tau: np.ndarray, gamma: np.ndarray, out_png: Path) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(8, 4.5))
    ax.plot(np.log10(tau), gamma, "b-", lw=1.5)
    ax.set_xlabel(r"$\log_{10}(\tau\,/\,\mathrm{s})$")
    ax.set_ylabel(r"$G(\ln\tau)$ (arb. u.)")
    ax.set_title("DRT (Tikhonov / nonnegative LS)")
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(out_png, dpi=150)
    plt.close(fig)


def plot_nyquist(z_meas: np.ndarray, z_fit: np.ndarray, out_png: Path) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(5.5, 5.5))
    ax.plot(z_meas.real, -z_meas.imag, "ko", ms=3, label="measured (-Z'')")
    ax.plot(z_fit.real, -z_fit.imag, "r-", lw=1.2, label="fit (-Z'')")
    ax.set_xlabel(r"$Z'$ ($\Omega$)")
    ax.set_ylabel(r"$-Z''$ ($\Omega$)")
    ax.legend()
    ax.grid(True, alpha=0.3)
    try:
        ax.set_aspect("equal", adjustable="datalim")
    except Exception:
        pass
    fig.tight_layout()
    fig.savefig(out_png, dpi=150)
    plt.close(fig)


def main() -> int:
    ap = argparse.ArgumentParser(description="EIS DRT analysis")
    ap.add_argument("--input", required=True)
    ap.add_argument("--output", required=True)
    ap.add_argument("--lambda", dest="regularization_lambda", type=float, default=0.1)
    ap.add_argument("--method", default="tikhonov")
    args = ap.parse_args()

    try:
        inp = Path(args.input).expanduser().resolve()
        out = Path(args.output).expanduser().resolve()
        out.mkdir(parents=True, exist_ok=True)

        df = _read_csv_flexible(inp)
        fc, zr_c, zi_c = _find_columns(df)
        freq = pd.to_numeric(df[fc], errors="coerce").to_numpy()
        z_real = pd.to_numeric(df[zr_c], errors="coerce").to_numpy()
        z_imag = pd.to_numeric(df[zi_c], errors="coerce").to_numpy()
        mask = np.isfinite(freq) & np.isfinite(z_real) & np.isfinite(z_imag) & (freq > 0)
        freq = freq[mask]
        z_real = z_real[mask]
        z_imag = z_imag[mask]
        if len(freq) < 5:
            print("ERROR: 有效频点过少（<5）。", file=sys.stderr)
            return 1

        z_complex = z_real + 1j * z_imag

        tau, gamma, r_inf, z_fit = run_drt(
            freq,
            z_complex,
            n_tau=80,
            regularization_lambda=args.regularization_lambda,
        )

        err = _rmse(z_complex, z_fit)

        summary: Dict[str, Any] = {
            "input": str(inp),
            "output_dir": str(out),
            "n_points": int(freq.size),
            "r_inf_ohm": r_inf,
            "rmse_complex_fit": err,
            "regularization_lambda": args.regularization_lambda,
            "columns_detected": {"frequency": fc, "z_real": zr_c, "z_imag": zi_c},
            "tau_s": tau.tolist(),
            "gamma": gamma.tolist(),
            "freq_hz": freq.tolist(),
            "z_meas": {"real": z_real.tolist(), "imag": z_imag.tolist()},
            "z_fit": {"real": z_fit.real.tolist(), "imag": z_fit.imag.tolist()},
        }

        json_path = out / "drt_summary.json"
        with open(json_path, "w", encoding="utf-8") as f:
            json.dump(summary, f, ensure_ascii=False, indent=2)

        plot_drt(tau, gamma, out / "drt_distribution.png")
        plot_nyquist(z_complex, z_fit, out / "nyquist_fit.png")

        print(f"OK: wrote {json_path}", file=sys.stderr)
        return 0
    except Exception:
        traceback.print_exc(file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
