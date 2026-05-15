#!/usr/bin/env python3
"""
cell-cycle-phase: Estimate cell-cycle phase durations from EdU/BrdU dual-pulse flow-cytometry data.

Usage:
    python cell_cycle_analysis.py --data '{"time_points":[...],"edu_positive":[...],...}' \
        --initial '{"g1_duration":6,"s_duration":8,"g2m_duration":4,"death_rate":0.02}'
"""

import argparse
import json
import sys
import traceback
from typing import Dict, List, Optional

try:
    import numpy as np
    from scipy.optimize import minimize, Bounds
    SCIPY_AVAILABLE = True
except ImportError:
    SCIPY_AVAILABLE = False


def parse_json_input(data_str: str) -> Dict:
    """Parse flow-cytometry JSON data."""
    try:
        data = json.loads(data_str)
    except json.JSONDecodeError as e:
        raise ValueError(f"Invalid JSON in --data: {e}")
    
    required_keys = ['time_points', 'edu_positive', 'brdu_positive', 'double_positive']
    for key in required_keys:
        if key not in data:
            raise ValueError(f"Missing required key in data: '{key}'")
    
    n = len(data['time_points'])
    for key in required_keys:
        if len(data[key]) != n:
            raise ValueError(f"Length mismatch: '{key}' has {len(data[key])} elements, expected {n}")
    
    return data


def parse_initial_estimates(est_str: str) -> Dict:
    """Parse initial estimate JSON."""
    try:
        est = json.loads(est_str)
    except json.JSONDecodeError as e:
        raise ValueError(f"Invalid JSON in --initial: {e}")
    
    defaults = {'g1_duration': 6.0, 's_duration': 8.0, 'g2m_duration': 4.0, 'death_rate': 0.02}
    for k, v in defaults.items():
        est.setdefault(k, v)
    
    return est


def forward_model(time_points: List[float], g1: float, s: float, g2m: float,
                  death_rate: float, eta_edu: float, eta_brd: float, eta_double: float,
                  t_offset: float = 0.0) -> Dict[str, List[float]]:
    """
    Continuous dual-labeling forward model with time offset.
    
    Cells are exposed to both EdU and BrdU starting at t = -t_offset
    (relative to the first sampling time point). A cell becomes positive
    once it enters S phase during the labeling period.
    """
    tc = g1 + s + g2m
    if tc <= 0:
        return {'edu': [0.0]*len(time_points), 'brdu': [0.0]*len(time_points), 'double': [0.0]*len(time_points)}
    
    edu_pred = []
    brdu_pred = []
    double_pred = []
    
    for t in time_points:
        tau = t + t_offset  # effective labeling time
        
        # Core progression: fraction of cells that have entered S by time tau
        if tau <= 0:
            progress = s / tc  # cells already in S at labeling start
        elif tau < g1:
            progress = (tau + s) / tc
        elif tau < g1 + g2m:
            progress = (g1 + s + (tau - g1)) / tc
        elif tau < tc:
            progress = (s + tau) / tc
        else:
            progress = 1.0
        
        # Mortality attenuation
        mortality_factor = max(0.0, 1.0 - death_rate * tau)
        
        f_t = progress * mortality_factor
        
        edu_pred.append(max(0.0, min(1.0, eta_edu * f_t)))
        brdu_pred.append(max(0.0, min(1.0, eta_brd * f_t)))
        double_pred.append(max(0.0, min(1.0, eta_double * f_t)))
    
    return {'edu': edu_pred, 'brdu': brdu_pred, 'double': double_pred}


def loss_function(params: np.ndarray, time_points: List[float],
                  edu_obs: List[float], brdu_obs: List[float], double_obs: List[float]) -> float:
    """Compute weighted sum of squared residuals."""
    g1, s, g2m, death_rate, eta_edu, eta_brd, eta_double, t_offset = params
    
    preds = forward_model(time_points, g1, s, g2m, death_rate, eta_edu, eta_brd, eta_double, t_offset)
    
    # Scale predictions to 0-100 (observations are already in 0-100)
    edu_pred = [v * 100.0 for v in preds['edu']]
    brdu_pred = [v * 100.0 for v in preds['brdu']]
    double_pred = [v * 100.0 for v in preds['double']]
    
    # Observations are already percentages (0-100), no scaling needed
    edu_obs_s = edu_obs
    brdu_obs_s = brdu_obs
    double_obs_s = double_obs
    
    # Weighted residuals: double-positive is usually most informative
    weights = {'edu': 1.0, 'brdu': 1.0, 'double': 1.5}
    
    sse = 0.0
    for i in range(len(time_points)):
        sse += weights['edu'] * (edu_pred[i] - edu_obs_s[i])**2
        sse += weights['brdu'] * (brdu_pred[i] - brdu_obs_s[i])**2
        sse += weights['double'] * (double_pred[i] - double_obs_s[i])**2
    
    # Penalize extreme t_offset to avoid overfitting
    sse += 0.1 * (t_offset ** 2)
    
    return sse


def estimate_phases(data: Dict, initial: Dict) -> Dict:
    """Run optimization to estimate cell-cycle parameters."""
    if not SCIPY_AVAILABLE:
        raise RuntimeError("scipy and numpy are required. Install with: pip install scipy numpy")
    
    time_points = data['time_points']
    edu_obs = data['edu_positive']
    brdu_obs = data['brdu_positive']
    double_obs = data['double_positive']
    
    # Initial parameter vector
    g1_0 = float(initial['g1_duration'])
    s_0 = float(initial['s_duration'])
    g2m_0 = float(initial['g2m_duration'])
    d_0 = float(initial['death_rate'])
    
    # Initial labeling efficiencies from data at final time point
    final_edu = edu_obs[-1] / 100.0
    final_brdu = brdu_obs[-1] / 100.0
    final_double = double_obs[-1] / 100.0
    
    # Rough estimate of saturation f(t) at final time point
    f_est = 0.95  # assume near saturation
    
    eta_edu_0 = min(1.5, max(0.1, final_edu / f_est))
    eta_brd_0 = min(1.5, max(0.1, final_brdu / f_est))
    eta_double_0 = min(1.5, max(0.05, final_double / f_est))
    
    t_offset_0 = 0.0
    
    x0 = np.array([g1_0, s_0, g2m_0, d_0, eta_edu_0, eta_brd_0, eta_double_0, t_offset_0])
    
    # Bounds: G1, S, G2M > 0.1; death_rate 0..0.5; efficiencies 0.01..1.5; t_offset 0..20
    lower = np.array([0.1, 0.1, 0.1, 0.0, 0.01, 0.01, 0.01, 0.0])
    upper = np.array([50.0, 50.0, 30.0, 0.5, 1.5, 1.5, 1.5, 20.0])
    bounds = Bounds(lower, upper)
    
    # Optimization
    result = minimize(
        loss_function,
        x0,
        args=(time_points, edu_obs, brdu_obs, double_obs),
        method='L-BFGS-B',
        bounds=bounds,
        options={'maxiter': 1000, 'ftol': 1e-12, 'gtol': 1e-12}
    )
    
    g1, s, g2m, d, eta_edu, eta_brd, eta_double, t_offset = result.x
    tc = g1 + s + g2m
    
    # Predictions with fitted parameters
    preds = forward_model(time_points, g1, s, g2m, d, eta_edu, eta_brd, eta_double, t_offset)
    
    # Compute fit statistics
    def r_squared(obs, pred):
        obs_arr = np.array(obs)
        pred_arr = np.array([p * 100.0 for p in pred])
        ss_res = np.sum((obs_arr - pred_arr) ** 2)
        ss_tot = np.sum((obs_arr - np.mean(obs_arr)) ** 2)
        if ss_tot == 0:
            return 1.0
        return 1.0 - ss_res / ss_tot
    
    r2_edu = r_squared(edu_obs, preds['edu'])
    r2_brdu = r_squared(brdu_obs, preds['brdu'])
    r2_double = r_squared(double_obs, preds['double'])
    
    return {
        'success': result.success,
        'message': result.message,
        'final_loss': result.fun,
        'iterations': result.nit,
        'g1_duration': float(g1),
        's_duration': float(s),
        'g2m_duration': float(g2m),
        'total_cycle_time': float(tc),
        'death_rate': float(d),
        'death_rate_percent': float(d * 100.0),
        'eta_edu': float(eta_edu),
        'eta_brd': float(eta_brd),
        'eta_double': float(eta_double),
        't_offset': float(t_offset),
        'phase_fractions': {
            'g1': float(g1 / tc),
            's': float(s / tc),
            'g2m': float(g2m / tc),
        },
        'predictions': {
            'time_points': time_points,
            'edu_positive': [float(v * 100.0) for v in preds['edu']],
            'brdu_positive': [float(v * 100.0) for v in preds['brdu']],
            'double_positive': [float(v * 100.0) for v in preds['double']],
        },
        'fit_statistics': {
            'r2_edu': float(r2_edu),
            'r2_brdu': float(r2_brdu),
            'r2_double': float(r2_double),
            'r2_average': float((r2_edu + r2_brdu + r2_double) / 3.0),
        },
        'initial_estimates': initial,
        'observed_data': data,
        'warnings': [],
    }


def format_log_report(result: Dict) -> str:
    """Format a detailed text analysis log."""
    lines = []
    
    lines.append("=" * 70)
    lines.append("CELL CYCLE PHASE DURATION ANALYSIS LOG")
    lines.append("=" * 70)
    lines.append("")
    
    # Section 1: Input Data
    lines.append("-" * 70)
    lines.append("SECTION 1: INPUT DATA")
    lines.append("-" * 70)
    lines.append("")
    
    data = result['observed_data']
    lines.append("Flow Cytometry Measurements (percent positive cells):")
    lines.append("")
    lines.append(f"{'Time (h)':>10s} {'EdU+':>12s} {'BrdU+':>12s} {'Double+':>12s}")
    lines.append("-" * 50)
    for i in range(len(data['time_points'])):
        lines.append(f"{data['time_points'][i]:>10.1f} {data['edu_positive'][i]:>11.1f}% {data['brdu_positive'][i]:>11.1f}% {data['double_positive'][i]:>11.1f}%")
    lines.append("")
    
    init = result['initial_estimates']
    lines.append("Initial Estimates Provided by User:")
    lines.append(f"  G1 duration:    {init['g1_duration']:.2f} h")
    lines.append(f"  S duration:     {init['s_duration']:.2f} h")
    lines.append(f"  G2/M duration:  {init['g2m_duration']:.2f} h")
    lines.append(f"  Death rate:     {init['death_rate']:.4f} ({init['death_rate']*100:.2f}% per hour)")
    lines.append("")
    
    # Section 2: Optimization Process
    lines.append("-" * 70)
    lines.append("SECTION 2: OPTIMIZATION PROCESS")
    lines.append("-" * 70)
    lines.append("")
    lines.append(f"Optimization Algorithm: L-BFGS-B (bounded nonlinear optimization)")
    lines.append(f"Number of Parameters:   8 (G1, S, G2M, death_rate, eta_edu, eta_brd, eta_double, t_offset)")
    lines.append(f"Convergence Status:     {'SUCCESS' if result['success'] else 'WARNING - may not have converged'}")
    lines.append(f"Optimizer Message:      {result['message']}")
    lines.append(f"Iterations Performed:   {result['iterations']}")
    lines.append(f"Final Loss (SSE):       {result['final_loss']:.4f}")
    lines.append("")
    
    # Section 3: Estimated Results
    lines.append("-" * 70)
    lines.append("SECTION 3: ESTIMATED CELL CYCLE PARAMETERS")
    lines.append("-" * 70)
    lines.append("")
    
    tc = result['total_cycle_time']
    lines.append(f"Total Cell Cycle Time (Tc): {tc:.2f} hours")
    lines.append("")
    lines.append("Phase Durations:")
    lines.append(f"  G1 phase:   {result['g1_duration']:>8.2f} h  ({result['phase_fractions']['g1']*100:>6.2f}% of cycle)")
    lines.append(f"  S phase:    {result['s_duration']:>8.2f} h  ({result['phase_fractions']['s']*100:>6.2f}% of cycle)")
    lines.append(f"  G2/M phase: {result['g2m_duration']:>8.2f} h  ({result['phase_fractions']['g2m']*100:>6.2f}% of cycle)")
    lines.append("")
    
    lines.append("Estimated Mortality:")
    d = result['death_rate']
    lines.append(f"  Death rate: {d:.4f} per hour")
    lines.append(f"  Death rate: {result['death_rate_percent']:.2f}% per hour")
    if d > 1e-6:
        lines.append(f"  Half-life:  {0.693 / d:.2f} hours (if exponential)")
    else:
        lines.append("  Half-life:  N/A (death rate is effectively zero)")
    lines.append("")
    
    lines.append("Estimated Labeling Parameters:")
    lines.append(f"  EdU efficiency:    {result['eta_edu']*100:.1f}%")
    lines.append(f"  BrdU efficiency:   {result['eta_brd']*100:.1f}%")
    lines.append(f"  Double-positive eff: {result['eta_double']*100:.1f}%")
    if result.get('t_offset', 0) > 0.5:
        lines.append(f"  Labeling started {result['t_offset']:.1f} h before first measurement")
    lines.append("")
    
    # Check for warnings
    warnings = []
    if result['eta_edu'] > 1.0 or result['eta_brd'] > 1.0 or result['eta_double'] > 1.0:
        warnings.append("Effective labeling efficiency > 100%. This may indicate quiescent/non-cycling cells or data that exceeds the model's theoretical maximum.")
    if result['fit_statistics']['r2_average'] < 0.5:
        warnings.append("Model fit quality is poor (R² < 0.5). Check data quality or experimental design.")
    if result['s_duration'] / result['total_cycle_time'] > 0.8:
        warnings.append("S phase occupies >80% of cycle. This is biologically unusual; verify input data.")
    
    if warnings:
        lines.append("⚠️  WARNINGS:")
        for w in warnings:
            lines.append(f"  • {w}")
        lines.append("")
    
    # Section 4: Fit Quality
    lines.append("-" * 70)
    lines.append("SECTION 4: MODEL FIT QUALITY")
    lines.append("-" * 70)
    lines.append("")
    
    stats = result['fit_statistics']
    lines.append(f"R-squared (EdU):      {stats['r2_edu']:.4f}")
    lines.append(f"R-squared (BrdU):     {stats['r2_brdu']:.4f}")
    lines.append(f"R-squared (Double):   {stats['r2_double']:.4f}")
    lines.append(f"R-squared (Average):  {stats['r2_average']:.4f}")
    lines.append("")
    
    # Predicted vs Observed table
    lines.append("Predicted vs. Observed Values:")
    lines.append("")
    lines.append(f"{'Time':>6s} {'EdU(obs)':>10s} {'EdU(pred)':>10s} {'BrdU(obs)':>10s} {'BrdU(pred)':>10s} {'Dbl(obs)':>10s} {'Dbl(pred)':>10s}")
    lines.append("-" * 72)
    preds = result['predictions']
    for i in range(len(preds['time_points'])):
        lines.append(f"{preds['time_points'][i]:>6.1f} "
                     f"{data['edu_positive'][i]:>9.1f}% {preds['edu_positive'][i]:>9.1f}% "
                     f"{data['brdu_positive'][i]:>9.1f}% {preds['brdu_positive'][i]:>9.1f}% "
                     f"{data['double_positive'][i]:>9.1f}% {preds['double_positive'][i]:>9.1f}%")
    lines.append("")
    
    # Section 5: Interpretation
    lines.append("-" * 70)
    lines.append("SECTION 5: BIOLOGICAL INTERPRETATION")
    lines.append("-" * 70)
    lines.append("")
    
    g1_pct = result['phase_fractions']['g1'] * 100
    s_pct = result['phase_fractions']['s'] * 100
    g2m_pct = result['phase_fractions']['g2m'] * 100
    
    lines.append("Phase Distribution Analysis:")
    if g1_pct > 50:
        lines.append("  • G1 phase occupies >50% of the cycle, suggesting a long quiescent/G1 period.")
    elif g1_pct < 20:
        lines.append("  • Short G1 phase; cells rapidly commit to DNA replication.")
    else:
        lines.append("  • G1 duration is within typical mammalian range (30-60% of cycle).")
    
    if s_pct > 40:
        lines.append("  • Extended S phase may indicate slow DNA replication or replication stress.")
    elif s_pct < 15:
        lines.append("  • Very short S phase suggests rapid DNA synthesis (e.g., embryonic cells).")
    else:
        lines.append("  • S phase duration is within normal range for most mammalian cells.")
    
    if g2m_pct > 25:
        lines.append("  • Long G2/M may indicate checkpoint arrest or slow mitotic progression.")
    else:
        lines.append("  • G2/M duration appears normal.")
    
    lines.append("")
    
    lines.append("Mortality Assessment:")
    if d < 0.005:
        lines.append("  • Very low death rate (<0.5%/h). Cell population is highly viable.")
    elif d < 0.02:
        lines.append("  • Moderate death rate (0.5-2%/h). Typical for standard culture conditions.")
    elif d < 0.05:
        lines.append("  • Elevated death rate (2-5%/h). Consider cytotoxic stress or apoptosis.")
    else:
        lines.append("  • High death rate (>5%/h). Significant cell loss; check culture conditions.")
    
    lines.append("")
    
    lines.append("Model Caveats:")
    lines.append("  1. Assumes steady-state population and constant phase durations.")
    lines.append("  2. Labeling efficiency may vary across cell subpopulations.")
    lines.append("  3. Mortality is modeled as a constant first-order process.")
    lines.append("  4. For best accuracy, time points should span at least one full cell cycle.")
    lines.append("")
    
    lines.append("=" * 70)
    lines.append("END OF ANALYSIS LOG")
    lines.append("=" * 70)
    
    return "\n".join(lines)


def main():
    parser = argparse.ArgumentParser(description="Cell cycle phase duration estimation from EdU/BrdU flow data")
    parser.add_argument("--data", required=True, help="JSON string with time_points, edu_positive, brdu_positive, double_positive")
    parser.add_argument("--initial", required=True, help="JSON string with g1_duration, s_duration, g2m_duration, death_rate")
    parser.add_argument("--output-format", default="text", choices=["text", "json"], help="Output format")
    
    args = parser.parse_args()
    
    try:
        data = parse_json_input(args.data)
        initial = parse_initial_estimates(args.initial)
        result = estimate_phases(data, initial)
        
        if args.output_format == "json":
            print(json.dumps(result, indent=2))
        else:
            print(format_log_report(result))
    
    except Exception as e:
        print(f"ERROR: {e}", file=sys.stderr)
        traceback.print_exc(file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()