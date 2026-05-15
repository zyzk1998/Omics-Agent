#!/usr/bin/env python3
"""
protease_kinetics.py

Protease kinetics analysis based on fluorescence readings.
Fits Michaelis-Menten model to substrate concentration vs initial reaction rate data.

Dependencies: numpy, scipy, matplotlib

Usage:
    python protease_kinetics.py --json-input '{"time_points":[0,10,20,...],"substrate_concentrations":[10,20,...],"fluorescence_data":[[...],[...]],"enzyme_concentration":0.1}' --output-dir /tmp/results
"""

import argparse
import base64
import json
import os
import sys
import warnings
from typing import Dict, List, Optional, Any, Tuple

import numpy as np
from scipy.optimize import curve_fit

# matplotlib backend for headless
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


# ---------------------------------------------------------------------------
# Michaelis-Menten Model
# ---------------------------------------------------------------------------
def michaelis_menten(S, Vmax, Km):
    """Michaelis-Menten equation: v = Vmax * [S] / (Km + [S])"""
    return (Vmax * S) / (Km + S)


def lineweaver_burk(S, v):
    """Lineweaver-Burk double reciprocal: 1/v vs 1/[S]"""
    # Handle zero values
    S_lb = np.array(S, dtype=float)
    v_lb = np.array(v, dtype=float)
    # Filter out zero or near-zero values
    mask = (S_lb > 1e-12) & (v_lb > 1e-12)
    inv_S = 1.0 / S_lb[mask]
    inv_v = 1.0 / v_lb[mask]
    return inv_S, inv_v


# ---------------------------------------------------------------------------
# Initial Rate Calculation from Fluorescence Data
# ---------------------------------------------------------------------------
def calculate_initial_rates(
    time_points: List[float],
    fluorescence_data: List[List[float]],
    substrate_concentrations: List[float],
    enzyme_concentration: float,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, Dict]:
    """
    Calculate initial reaction rates from fluorescence vs time data.
    
    For each substrate concentration, fit linear regression to the early time points
    to obtain initial rate (slope of fluorescence vs time).
    
    Returns:
        S: substrate concentrations (array)
        v: initial rates (array)
        v_err: standard errors of slopes (array)
        details: dict with per-substrate fitting details
    """
    time = np.array(time_points, dtype=float)
    n_substrates = len(substrate_concentrations)
    
    rates = []
    rate_errors = []
    details = {
        "per_substrate": [],
        "method": "linear regression on all time points",
        "time_points_used": list(time),
    }
    
    for i, (S_val, fluo_row) in enumerate(zip(substrate_concentrations, fluorescence_data)):
        fluo = np.array(fluo_row, dtype=float)
        
        if len(fluo) != len(time):
            raise ValueError(
                f"Fluorescence row {i} has {len(fluo)} values but {len(time)} time points"
            )
        
        # Linear fit: fluorescence = slope * time + intercept
        # slope = initial rate (in fluorescence units per second)
        # Use all points for robustness, or early half if curvature suspected
        slope, intercept, r_value, p_value, std_err = _linregress(time, fluo)
        
        rate = slope
        rate_err = std_err
        
        rates.append(rate)
        rate_errors.append(rate_err)
        
        details["per_substrate"].append({
            "substrate_concentration": S_val,
            "slope": float(rate),
            "intercept": float(intercept),
            "r_squared": float(r_value ** 2),
            "p_value": float(p_value),
            "std_error": float(rate_err),
            "fluorescence_range": [float(min(fluo)), float(max(fluo))],
        })
    
    S = np.array(substrate_concentrations, dtype=float)
    v = np.array(rates, dtype=float)
    v_err = np.array(rate_errors, dtype=float)
    
    return S, v, v_err, details


def _linregress(x, y):
    """Simple linear regression. Returns slope, intercept, r, p, stderr."""
    n = len(x)
    x_mean = np.mean(x)
    y_mean = np.mean(y)
    
    ss_xx = np.sum((x - x_mean) ** 2)
    ss_xy = np.sum((x - x_mean) * (y - y_mean))
    
    if ss_xx == 0:
        return 0.0, y_mean, 0.0, 1.0, 0.0
    
    slope = ss_xy / ss_xx
    intercept = y_mean - slope * x_mean
    
    # Predicted values
    y_pred = slope * x + intercept
    ss_res = np.sum((y - y_pred) ** 2)
    ss_tot = np.sum((y - y_mean) ** 2)
    
    r_value = np.sqrt(1 - ss_res / ss_tot) if ss_tot > 0 else 0.0
    if slope < 0:
        r_value = -abs(r_value)
    
    # Standard error of slope
    df = n - 2
    if df > 0 and ss_xx > 0:
        mse = ss_res / df
        stderr = np.sqrt(mse / ss_xx)
    else:
        stderr = 0.0
    
    # p-value (two-tailed t-test for slope != 0)
    if stderr > 1e-12:
        t_stat = slope / stderr
        # Approximate p-value using normal distribution
        p_value = 2 * (1 - _normal_cdf(abs(t_stat)))
    else:
        p_value = 0.0 if abs(slope) > 1e-12 else 1.0
    
    return slope, intercept, r_value, p_value, stderr


def _normal_cdf(x):
    """Approximate normal CDF."""
    import math
    return 0.5 * (1 + math.erf(x / math.sqrt(2)))


# ---------------------------------------------------------------------------
# Michaelis-Menten Fitting
# ---------------------------------------------------------------------------
def fit_michaelis_menten(
    S: np.ndarray,
    v: np.ndarray,
    v_err: Optional[np.ndarray] = None,
) -> Dict[str, Any]:
    """
    Fit Michaelis-Menten model to substrate-velocity data.
    
    Uses scipy.optimize.curve_fit for non-linear least squares.
    Also performs Lineweaver-Burk linearization as validation.
    
    Returns dict with fitted parameters and statistics.
    """
    # Filter out non-positive rates
    mask = v > 1e-12
    if np.sum(mask) < 2:
        raise ValueError("Need at least 2 positive initial rates to fit MM model")
    
    S_fit = S[mask]
    v_fit = v[mask]
    
    # Initial parameter estimates
    Vmax_est = np.max(v_fit)
    Km_est = S_fit[np.argmin(np.abs(v_fit - Vmax_est / 2))] if len(S_fit) > 0 else np.median(S_fit)
    
    # Bounds: Vmax > 0, Km > 0
    bounds = ([1e-12, 1e-12], [np.inf, np.inf])
    p0 = [Vmax_est, Km_est]
    
    try:
        if v_err is not None and np.any(v_err[mask] > 1e-12):
            sigma = np.where(v_err[mask] > 1e-12, v_err[mask], 1.0)
            popt, pcov = curve_fit(
                michaelis_menten,
                S_fit,
                v_fit,
                p0=p0,
                bounds=bounds,
                sigma=sigma,
                absolute_sigma=True,
                maxfev=10000,
            )
        else:
            popt, pcov = curve_fit(
                michaelis_menten,
                S_fit,
                v_fit,
                p0=p0,
                bounds=bounds,
                maxfev=10000,
            )
        
        Vmax, Km = popt
        
        # Standard errors from covariance matrix
        perr = np.sqrt(np.diag(pcov))
        Vmax_err, Km_err = perr
        
        # Calculate R^2
        v_pred = michaelis_menten(S_fit, Vmax, Km)
        ss_res = np.sum((v_fit - v_pred) ** 2)
        ss_tot = np.sum((v_fit - np.mean(v_fit)) ** 2)
        r_squared = 1 - ss_res / ss_tot if ss_tot > 0 else 0.0
        
        # Lineweaver-Burk validation
        inv_S, inv_v = lineweaver_burk(S_fit, v_fit)
        lb_details = None
        if len(inv_S) >= 2:
            lb_slope, lb_intercept, lb_r, lb_p, lb_stderr = _linregress(inv_S, inv_v)
            lb_Vmax = 1.0 / lb_intercept if abs(lb_intercept) > 1e-12 else float('inf')
            lb_Km = lb_slope / lb_intercept if abs(lb_intercept) > 1e-12 else float('inf')
            lb_details = {
                "Vmax_lb": float(lb_Vmax),
                "Km_lb": float(lb_Km),
                "r_squared_lb": float(lb_r ** 2),
            }
        
        return {
            "Vmax": float(Vmax),
            "Vmax_std_error": float(Vmax_err),
            "Km": float(Km),
            "Km_std_error": float(Km_err),
            "r_squared": float(r_squared),
            "num_points_fitted": int(len(S_fit)),
            "lineweaver_burk": lb_details,
        }
    
    except Exception as e:
        raise RuntimeError(f"Michaelis-Menten fitting failed: {str(e)}")


def calculate_kcat(Vmax: float, enzyme_concentration: float) -> float:
    """kcat = Vmax / [E]"""
    if enzyme_concentration <= 0:
        return float('inf')
    return Vmax / enzyme_concentration


def calculate_kcat_Km(kcat: float, Km: float) -> float:
    """Catalytic efficiency = kcat / Km"""
    if Km <= 0:
        return float('inf')
    return kcat / Km


# ---------------------------------------------------------------------------
# Visualization
# ---------------------------------------------------------------------------
def plot_michaelis_menten(
    S: np.ndarray,
    v: np.ndarray,
    v_err: Optional[np.ndarray],
    fit_params: Dict,
    output_path: str,
    title: str = "Michaelis-Menten Kinetics",
) -> str:
    """Generate MM fit curve plot."""
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    
    Vmax = fit_params["Vmax"]
    Km = fit_params["Km"]
    
    # --- Plot 1: MM Curve ---
    ax1 = axes[0]
    
    # Smooth curve for fit
    S_smooth = np.linspace(0, max(S) * 1.2, 200)
    v_smooth = michaelis_menten(S_smooth, Vmax, Km)
    
    ax1.plot(S_smooth, v_smooth, "b-", linewidth=2, label=f"Fit: Vmax={Vmax:.3f}, Km={Km:.3f}")
    
    # Data points with error bars
    if v_err is not None and np.any(v_err > 1e-12):
        ax1.errorbar(S, v, yerr=v_err, fmt="ro", capsize=4, label="Data")
    else:
        ax1.scatter(S, v, c="red", s=60, zorder=5, label="Data")
    
    # Km dashed line
    ax1.axhline(Vmax / 2, color="gray", linestyle="--", alpha=0.5)
    ax1.axvline(Km, color="gray", linestyle="--", alpha=0.5)
    ax1.annotate(f"Km = {Km:.2f}", xy=(Km, Vmax / 2), xytext=(Km * 1.3, Vmax * 0.4),
                 fontsize=9, arrowprops=dict(arrowstyle="->", color="gray"))
    
    ax1.set_xlabel("Substrate Concentration [S] (μM)", fontsize=11)
    ax1.set_ylabel("Initial Rate v (a.u./s)", fontsize=11)
    ax1.set_title("Michaelis-Menten Fit", fontsize=12)
    ax1.legend(loc="lower right")
    ax1.set_xlim(left=0)
    ax1.set_ylim(bottom=0)
    ax1.grid(True, alpha=0.3)
    
    # --- Plot 2: Lineweaver-Burk ---
    ax2 = axes[1]
    
    inv_S, inv_v = lineweaver_burk(S, v)
    if len(inv_S) >= 2:
        # Fit line
        inv_S_range = np.linspace(min(inv_S) * 0.8, max(inv_S) * 1.2, 100)
        lb_slope, lb_intercept, _, _, _ = _linregress(inv_S, inv_v)
        inv_v_line = lb_slope * inv_S_range + lb_intercept
        ax2.plot(inv_S_range, inv_v_line, "g-", linewidth=2, label="LB linear fit")
        ax2.scatter(inv_S, inv_v, c="red", s=60, zorder=5, label="Data")
        
        ax2.set_xlabel("1 / [S] (μM⁻¹)", fontsize=11)
        ax2.set_ylabel("1 / v (s/a.u.)", fontsize=11)
        ax2.set_title("Lineweaver-Burk Plot", fontsize=12)
        ax2.legend()
        ax2.grid(True, alpha=0.3)
    else:
        ax2.text(0.5, 0.5, "Insufficient data\nfor Lineweaver-Burk", 
                 ha="center", va="center", transform=ax2.transAxes)
    
    fig.suptitle(title, fontsize=13, fontweight="bold")
    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close()
    
    return output_path


def plot_fluorescence_timecourses(
    time_points: List[float],
    fluorescence_data: List[List[float]],
    substrate_concentrations: List[float],
    output_path: str,
) -> str:
    """Plot fluorescence vs time for each substrate concentration."""
    fig, ax = plt.subplots(figsize=(9, 6))
    
    time = np.array(time_points)
    colors = plt.cm.viridis(np.linspace(0, 1, len(substrate_concentrations)))
    
    for i, (S_val, fluo_row) in enumerate(zip(substrate_concentrations, fluorescence_data)):
        fluo = np.array(fluo_row)
        ax.plot(time, fluo, "o-", color=colors[i], linewidth=1.5, 
                markersize=4, label=f"[{S_val} μM]")
    
    ax.set_xlabel("Time (s)", fontsize=11)
    ax.set_ylabel("Fluorescence (a.u.)", fontsize=11)
    ax.set_title("Fluorescence Time Courses", fontsize=12)
    ax.legend(title="Substrate", loc="upper left")
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close()
    
    return output_path


# ---------------------------------------------------------------------------
# Report Generation
# ---------------------------------------------------------------------------
def generate_report(
    input_data: Dict,
    rate_details: Dict,
    fit_params: Dict,
    enzyme_concentration: float,
    output_path: str,
) -> str:
    """Generate detailed analysis report."""
    lines = []
    lines.append("=" * 70)
    lines.append("PROTEASE KINETICS ANALYSIS REPORT")
    lines.append("=" * 70)
    lines.append("")
    
    # Input summary
    lines.append("INPUT PARAMETERS")
    lines.append("-" * 70)
    lines.append(f"Enzyme concentration: {enzyme_concentration} μM")
    lines.append(f"Number of substrate concentrations: {len(input_data.get('substrate_concentrations', []))}")
    lines.append(f"Time points: {input_data.get('time_points', [])}")
    lines.append("")
    
    # Initial rates
    lines.append("INITIAL RATE CALCULATION")
    lines.append("-" * 70)
    lines.append(f"Method: {rate_details.get('method', 'linear regression')}")
    lines.append("")
    for detail in rate_details.get("per_substrate", []):
        lines.append(f"  [S] = {detail['substrate_concentration']} μM:")
        lines.append(f"    Slope (v) = {detail['slope']:.4f} ± {detail['std_error']:.4f} a.u./s")
        lines.append(f"    R² = {detail['r_squared']:.4f}, p = {detail['p_value']:.4e}")
    lines.append("")
    
    # Fitted parameters
    lines.append("MICHAELIS-MENTEN FIT RESULTS")
    lines.append("-" * 70)
    lines.append(f"Vmax = {fit_params['Vmax']:.4f} ± {fit_params['Vmax_std_error']:.4f} a.u./s")
    lines.append(f"Km   = {fit_params['Km']:.4f} ± {fit_params['Km_std_error']:.4f} μM")
    lines.append(f"R² (fit) = {fit_params['r_squared']:.4f}")
    lines.append(f"Points used: {fit_params['num_points_fitted']}")
    lines.append("")
    
    # Catalytic parameters
    kcat = calculate_kcat(fit_params['Vmax'], enzyme_concentration)
    kcat_Km = calculate_kcat_Km(kcat, fit_params['Km'])
    
    lines.append("CATALYTIC PARAMETERS")
    lines.append("-" * 70)
    lines.append(f"kcat = Vmax / [E] = {kcat:.4f} s⁻¹")
    lines.append(f"kcat/Km = {kcat_Km:.4f} μM⁻¹·s⁻¹")
    lines.append("")
    
    # Lineweaver-Burk validation
    lb = fit_params.get("lineweaver_burk")
    if lb:
        lines.append("LINEWEAVER-BURK VALIDATION")
        lines.append("-" * 70)
        lines.append(f"Vmax (LB) = {lb['Vmax_lb']:.4f} a.u./s")
        lines.append(f"Km (LB)   = {lb['Km_lb']:.4f} μM")
        lines.append(f"R² (LB)   = {lb['r_squared_lb']:.4f}")
        lines.append("")
    
    lines.append("=" * 70)
    lines.append("END OF REPORT")
    lines.append("=" * 70)
    
    with open(output_path, "w") as f:
        f.write("\n".join(lines))
    
    return output_path


def generate_json_output(
    input_data: Dict,
    rate_details: Dict,
    fit_params: Dict,
    enzyme_concentration: float,
    output_path: str,
) -> str:
    """Generate machine-readable JSON output."""
    kcat = calculate_kcat(fit_params['Vmax'], enzyme_concentration)
    kcat_Km = calculate_kcat_Km(kcat, fit_params['Km'])
    
    result = {
        "input": {
            "enzyme_concentration_uM": float(enzyme_concentration),
            "substrate_concentrations_uM": input_data.get("substrate_concentrations", []),
            "time_points_s": input_data.get("time_points", []),
        },
        "initial_rates": rate_details,
        "michaelis_menten_fit": fit_params,
        "catalytic_parameters": {
            "kcat_s-1": float(kcat),
            "kcat_over_Km_uM-1_s-1": float(kcat_Km),
        },
    }
    
    with open(output_path, "w") as f:
        json.dump(result, f, indent=2)
    
    return output_path


# ---------------------------------------------------------------------------
# Main Pipeline
# ---------------------------------------------------------------------------
def run_analysis(input_data: Dict[str, Any], output_dir: str) -> Dict[str, Any]:
    """
    Main analysis pipeline.
    
    Returns dict with:
        - research_log: analysis log string
        - mm_plot_path: Michaelis-Menten plot path
        - timecourse_plot_path: fluorescence timecourse plot path
        - report_path: detailed report path
        - json_path: machine-readable JSON path
        - error: error message if any
    """
    log_lines = []
    
    try:
        # 1. Parse input
        log_lines.append("Step 1: Parsing input data...")
        time_points = input_data.get("time_points", [])
        substrate_concentrations = input_data.get("substrate_concentrations", [])
        fluorescence_data = input_data.get("fluorescence_data", [])
        enzyme_concentration = float(input_data.get("enzyme_concentration", 0.1))
        
        if not time_points:
            raise ValueError("time_points is required")
        if not substrate_concentrations:
            raise ValueError("substrate_concentrations is required")
        if not fluorescence_data:
            raise ValueError("fluorescence_data is required")
        if len(fluorescence_data) != len(substrate_concentrations):
            raise ValueError(
                f"fluorescence_data has {len(fluorescence_data)} rows but "
                f"{len(substrate_concentrations)} substrate concentrations"
            )
        
        log_lines.append(f"  - {len(substrate_concentrations)} substrate concentrations")
        log_lines.append(f"  - {len(time_points)} time points")
        log_lines.append(f"  - Enzyme concentration: {enzyme_concentration} μM")
        
        os.makedirs(output_dir, exist_ok=True)
        
        # 2. Calculate initial rates
        log_lines.append("")
        log_lines.append("Step 2: Calculating initial reaction rates from fluorescence data...")
        S, v, v_err, rate_details = calculate_initial_rates(
            time_points, fluorescence_data, substrate_concentrations, enzyme_concentration
        )
        
        for detail in rate_details["per_substrate"]:
            log_lines.append(
                f"  [S]={detail['substrate_concentration']} μM: "
                f"v={detail['slope']:.4f}±{detail['std_error']:.4f} a.u./s, R²={detail['r_squared']:.4f}"
            )
        
        # 3. Fit Michaelis-Menten
        log_lines.append("")
        log_lines.append("Step 3: Fitting Michaelis-Menten model...")
        fit_params = fit_michaelis_menten(S, v, v_err)
        
        log_lines.append(f"  Vmax = {fit_params['Vmax']:.4f} ± {fit_params['Vmax_std_error']:.4f} a.u./s")
        log_lines.append(f"  Km   = {fit_params['Km']:.4f} ± {fit_params['Km_std_error']:.4f} μM")
        log_lines.append(f"  R²   = {fit_params['r_squared']:.4f}")
        
        # 4. Catalytic parameters
        kcat = calculate_kcat(fit_params['Vmax'], enzyme_concentration)
        kcat_Km = calculate_kcat_Km(kcat, fit_params['Km'])
        log_lines.append("")
        log_lines.append("Step 4: Calculating catalytic parameters...")
        log_lines.append(f"  kcat     = {kcat:.4f} s⁻¹")
        log_lines.append(f"  kcat/Km  = {kcat_Km:.4f} μM⁻¹·s⁻¹")
        
        # 5. Generate plots
        log_lines.append("")
        log_lines.append("Step 5: Generating plots...")
        
        mm_plot_path = os.path.join(output_dir, "michaelis_menten_plot.png")
        plot_michaelis_menten(S, v, v_err, fit_params, mm_plot_path)
        log_lines.append(f"  - MM plot: {mm_plot_path}")
        
        tc_plot_path = os.path.join(output_dir, "fluorescence_timecourses.png")
        plot_fluorescence_timecourses(time_points, fluorescence_data, substrate_concentrations, tc_plot_path)
        log_lines.append(f"  - Timecourse plot: {tc_plot_path}")
        
        # 6. Generate reports
        log_lines.append("")
        log_lines.append("Step 6: Generating reports...")
        
        report_path = os.path.join(output_dir, "kinetics_report.txt")
        generate_report(input_data, rate_details, fit_params, enzyme_concentration, report_path)
        log_lines.append(f"  - Report: {report_path}")
        
        json_path = os.path.join(output_dir, "kinetics_data.json")
        generate_json_output(input_data, rate_details, fit_params, enzyme_concentration, json_path)
        log_lines.append(f"  - JSON: {json_path}")
        
        log_lines.append("")
        log_lines.append("Analysis complete!")
        
        return {
            "research_log": "\n".join(log_lines),
            "mm_plot_path": mm_plot_path,
            "timecourse_plot_path": tc_plot_path,
            "report_path": report_path,
            "json_path": json_path,
            "Vmax": fit_params["Vmax"],
            "Km": fit_params["Km"],
            "kcat": kcat,
            "kcat_over_Km": kcat_Km,
            "r_squared": fit_params["r_squared"],
            "error": None,
        }
    
    except Exception as e:
        log_lines.append(f"ERROR: {str(e)}")
        return {
            "research_log": "\n".join(log_lines),
            "mm_plot_path": None,
            "timecourse_plot_path": None,
            "report_path": None,
            "json_path": None,
            "error": str(e),
        }


def _png_file_to_base64(path: Optional[str]) -> Optional[str]:
    if not path or not os.path.isfile(path):
        return None
    with open(path, "rb") as f:
        return base64.b64encode(f.read()).decode("ascii")


def embed_mm_fit_plot_base64(result: Dict[str, Any]) -> Dict[str, Any]:
    """将 Michaelis-Menten 拟合曲线图以 Base64 写入结果。"""
    if not isinstance(result, dict):
        return result
    b64 = _png_file_to_base64(result.get("mm_plot_path"))
    if b64:
        result["michaelis_menten_fit_plot_png_base64"] = b64
    return result


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(description="Protease Kinetics Analysis")
    parser.add_argument("--json-input", type=str, help="JSON string with input data")
    parser.add_argument("--output-dir", type=str, default="/tmp/protease_kinetics", help="Output directory")
    parser.add_argument("--output-json", type=str, help="Write full result as JSON to this path")
    
    args = parser.parse_args()
    
    if args.json_input:
        input_data = json.loads(args.json_input)
    else:
        print("Error: Provide --json-input", file=sys.stderr)
        sys.exit(1)
    
    result = run_analysis(input_data, args.output_dir)
    result = embed_mm_fit_plot_base64(result)
    
    # Print result as JSON
    print(json.dumps(result, indent=2, default=str))
    
    if args.output_json:
        with open(args.output_json, "w") as f:
            json.dump(result, f, indent=2, default=str)


if __name__ == "__main__":
    main()
