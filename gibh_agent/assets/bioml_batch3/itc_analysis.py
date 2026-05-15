#!/usr/bin/env python3
"""
ITC (Isothermal Titration Calorimetry) Analysis Tool
One-site binding model fitting to extract thermodynamic parameters:
Kd, dH, dS, dG, and n (stoichiometry)

Usage:
    python itc_analysis.py --data "[[inj1, vol1, heat1], [inj2, vol2, heat2], ...]" \
                          --protein-conc 1e-6 \
                          --ligand-conc 10e-6 \
                          --cell-volume 1.4e-3 \
                          --temperature 25.0

    python itc_analysis.py --file data.csv \
                          --protein-conc 1e-6 \
                          --ligand-conc 10e-6 \
                          --cell-volume 1.4e-3 \
                          --temperature 25.0
"""

import sys
import argparse
import json
import csv
import io
import base64
from typing import Optional

import numpy as np
from scipy.optimize import curve_fit, minimize
from scipy import stats

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

# Constants
R = 8.314  # J/(mol·K)
T_KELVIN_25C = 298.15


def parse_data_array(data_str):
    """Parse JSON-like array string into numpy array."""
    try:
        data = json.loads(data_str)
        return np.array(data, dtype=float)
    except Exception as e:
        raise ValueError(f"Failed to parse data array: {e}")


def load_csv_data(filepath):
    """Load ITC data from CSV/TSV file. Expects 3 columns: injection_number, volume, heat."""
    data = []
    with open(filepath, 'r') as f:
        # Try to detect delimiter
        sample = f.read(1024)
        f.seek(0)
        delimiter = '\t' if '\t' in sample else ','
        reader = csv.reader(f, delimiter=delimiter)
        
        for row in reader:
            if len(row) == 0:
                continue
            # Skip header rows
            try:
                vals = [float(x.strip()) for x in row[:3]]
                data.append(vals)
            except ValueError:
                continue  # Skip non-numeric rows
                
    return np.array(data, dtype=float)


def calculate_bound_fraction(Pt, Lt, Kd, n):
    """
    Calculate the fraction of bound sites for one-site binding model.
    
    Pt: total protein concentration in cell
    Lt: total ligand concentration in cell
    Kd: dissociation constant
    n: stoichiometry (binding sites per protein)
    
    Returns: fraction of bound sites theta = [PL] / (n*Pt)
    """
    P_eff = n * Pt
    
    if P_eff <= 0 or Lt < 0 or Kd <= 0:
        return 0.0
    
    a = 1.0
    b = -(P_eff + Lt + Kd)
    c = P_eff * Lt
    
    discriminant = b**2 - 4*a*c
    
    if discriminant < 0:
        return min(1.0, Lt / P_eff) if P_eff > 0 else 0.0
    
    sqrt_disc = np.sqrt(discriminant)
    PL = (-b - sqrt_disc) / (2*a)
    PL = np.clip(PL, 0, min(P_eff, Lt))
    
    theta = PL / P_eff if P_eff > 0 else 0.0
    theta = np.clip(theta, 0, 1)
    
    return theta


def injection_heat(injections, volumes, P0, L0, V_cell, Kd, dH, n, 
                   dilution_heat=0.0):
    """
    Calculate heat per injection (differential heat).
    
    Returns: array of heat per injection (same units as dH * concentration * volume)
    """
    N = len(injections)
    Q_inj = np.zeros(N)
    
    Pt = P0
    Lt = 0.0
    V_total = V_cell
    Q_prev = 0.0
    
    for i in range(N):
        vol_inj = volumes[i]
        
        # Update ligand concentration (accounting for dilution + addition)
        Lt = Lt * (1 - vol_inj / V_total) + L0 * vol_inj / V_total
        
        # Calculate bound fraction
        theta = calculate_bound_fraction(Pt, Lt, Kd, n)
        
        # Cumulative heat = fraction_bound * total_moles_protein * dH
        # Total moles of protein = P0 * V_cell (approximately constant)
        Q_cum = theta * n * Pt * V_cell * dH
        
        # Injection heat = cumulative difference + dilution correction
        dilution_corr = dilution_heat * vol_inj / V_total
        
        Q_inj[i] = Q_cum - Q_prev + dilution_corr
        Q_prev = Q_cum
        
    return Q_inj


def generate_initial_guesses(injections, volumes, heat_observed, P0, L0, V_cell):
    """Generate robust initial guesses for Kd, dH, and n."""
    N = len(heat_observed)
    
    # Estimate total heat and check sign
    total_heat = np.sum(heat_observed)
    abs_heat = np.abs(heat_observed)
    max_heat = np.max(abs_heat) if len(abs_heat) > 0 else 1e-15
    
    # Estimate n: default to 1.0 (most common for protein-ligand binding)
    # Only deviate if there's clear evidence of non-1:1 stoichiometry
    n_guess = 1.0
    
    # Estimate dH from total heat (assuming n=1 initially)
    moles_P = P0 * V_cell
    if moles_P > 0 and abs(total_heat) > 1e-15:
        dH_guess = total_heat / moles_P  # n=1 assumption
    else:
        dH_guess = -1e5 if total_heat < 0 else 1e5
    
    # Sanity check: typical ITC dH is -200 to +200 kJ/mol
    if abs(dH_guess) > 5e5:
        dH_guess = -1e5 if total_heat < 0 else 1e5
    elif abs(dH_guess) < 1e3:
        dH_guess = -1e5 if total_heat < 0 else 1e5
    
    # Estimate Kd from titration curvature
    cum_heat = np.cumsum(heat_observed)
    
    if abs(total_heat) > 1e-15 and len(cum_heat) > 3:
        # Find the injection where cumulative heat is ~50% of total
        half_fraction = 0.5
        target_heat = total_heat * half_fraction
        
        diff = np.abs(cum_heat - target_heat)
        half_idx = np.argmin(diff)
        
        if half_idx >= 0:
            # Ligand concentration at this point
            Lt_half = np.sum(L0 * volumes[:half_idx+1] / V_cell)
            
            # Check if we see saturation
            late_heat_avg = np.mean(abs_heat[-3:]) if N >= 3 else max_heat
            has_saturation = late_heat_avg < 0.25 * max_heat
            
            if has_saturation:
                # With saturation: Kd is roughly comparable to Lt at half-point
                # For 1:1 binding: Kd ≈ Lt_half when P0 << Kd
                # More generally: use a fraction of Lt_half
                Kd_guess = max(Lt_half / 3, 1e-12)
            else:
                # No saturation: Kd is larger than total ligand added
                Kd_guess = max(Lt_half * 2, 1e-12)
        else:
            Kd_guess = 1e-6
    else:
        Kd_guess = 1e-6
    
    # Final sanity checks
    Kd_guess = np.clip(Kd_guess, 1e-15, 1e-1)
    dH_guess = np.clip(dH_guess, -1e6, 1e6)
    
    return Kd_guess, dH_guess, n_guess


def fit_itc_model(injections, volumes, heat_observed, P0, L0, V_cell,
                  temperature_celsius=25.0,
                  Kd_guess=None, dH_guess=None, n_guess=None,
                  fit_dilution=False):
    """
    Fit one-site binding model to ITC data.
    
    Returns: fitted parameters and statistics
    """
    # Generate initial guesses
    Kd_est, dH_est, n_est = generate_initial_guesses(
        injections, volumes, heat_observed, P0, L0, V_cell
    )
    
    Kd_guess = Kd_guess if Kd_guess is not None else Kd_est
    dH_guess = dH_guess if dH_guess is not None else dH_est
    n_guess = n_guess if n_guess is not None else n_est
    
    N = len(heat_observed)
    
    # Parameter bounds
    # Kd: 1e-15 to 1e-1 M
    # dH: -1e6 to 1e6 J/mol
    # n: 0.1 to 10
    
    if fit_dilution and N >= 5:
        # Fit with dilution heat as additional parameter
        p0 = [Kd_guess, dH_guess, n_guess, 0.0]
        bounds = ([1e-15, -1e6, 0.1, -1e6], [1e-1, 1e6, 10.0, 1e6])
        
        def fit_func(x, Kd, dH, n, dH_dil):
            return injection_heat(injections, volumes, P0, L0, V_cell, Kd, dH, n, dH_dil)
        
        n_params = 4
    else:
        p0 = [Kd_guess, dH_guess, n_guess]
        bounds = ([1e-15, -1e6, 0.1], [1e-1, 1e6, 10.0])
        
        def fit_func(x, Kd, dH, n):
            return injection_heat(injections, volumes, P0, L0, V_cell, Kd, dH, n)
        
        n_params = 3
    
    popt = None
    pcov = None
    curve_fit_ok = False
    
    try:
        # Try curve_fit first with improved initial guesses
        popt, pcov = curve_fit(fit_func, injections, heat_observed, 
                               p0=p0, bounds=bounds, maxfev=50000,
                               method='trf', ftol=1e-15, xtol=1e-15, gtol=1e-15)
        
        # Check if curve_fit result is good
        heat_test = fit_func(injections, *popt)
        ss_res_test = np.sum((heat_observed - heat_test)**2)
        ss_tot_test = np.sum((heat_observed - np.mean(heat_observed))**2)
        if ss_tot_test > 0 and ss_res_test / ss_tot_test < 0.5:
            curve_fit_ok = True
            
    except Exception as e:
        pass
    
    if not curve_fit_ok or popt is None:
        # Fallback: use differential evolution for global optimization
        try:
            from scipy.optimize import differential_evolution
            
            def objective(params):
                predicted = fit_func(injections, *params)
                residuals = heat_observed - predicted
                return np.sum(residuals**2)
            
            # Set up bounds for differential evolution
            de_bounds = [(bounds[0][i], bounds[1][i]) for i in range(n_params)]
            
            # Run differential evolution
            result_de = differential_evolution(objective, de_bounds, 
                                              maxiter=2000, seed=42,
                                              polish=True,
                                              tol=1e-12,
                                              atol=1e-12)
            
            if result_de.success:
                popt = result_de.x
                # Approximate covariance
                pcov = np.eye(n_params) * np.nan
            else:
                return {'error': f"Optimization failed: {result_de.message}", 'converged': False}
                
        except ImportError:
            if popt is None:
                return {'error': "Curve fit failed and differential evolution not available", 'converged': False}
    
    # Extract parameters
    Kd_fit, dH_fit, n_fit = popt[0], popt[1], popt[2]
    dH_dil_fit = popt[3] if fit_dilution and N >= 5 else 0.0
    
    # Calculate parameter standard errors
    try:
        perr = np.sqrt(np.diag(pcov))
    except:
        perr = [np.nan] * n_params
    
    # Calculate fitted values
    if fit_dilution and N >= 5:
        heat_fitted = injection_heat(injections, volumes, P0, L0, V_cell, 
                                     Kd_fit, dH_fit, n_fit, dH_dil_fit)
    else:
        heat_fitted = injection_heat(injections, volumes, P0, L0, V_cell, 
                                     Kd_fit, dH_fit, n_fit)
    
    # Calculate residuals
    residuals = heat_observed - heat_fitted
    ss_res = np.sum(residuals**2)
    ss_tot = np.sum((heat_observed - np.mean(heat_observed))**2)
    
    # R-squared
    r_squared = 1 - ss_res / ss_tot if ss_tot != 0 else 0
    
    # Degrees of freedom
    n_points = len(heat_observed)
    dof = n_points - n_params
    
    # Reduced chi-square
    if dof > 0:
        chi2_reduced = ss_res / dof
    else:
        chi2_reduced = np.inf
    
    # Calculate thermodynamic parameters
    T_kelvin = temperature_celsius + 273.15
    Ka = 1.0 / Kd_fit if Kd_fit > 0 else np.inf
    
    if Ka > 0 and Ka != np.inf:
        dG = -R * T_kelvin * np.log(Ka)
        # Handle very large/small Ka values
        if not np.isfinite(dG):
            dG = np.nan
    else:
        dG = np.nan
    
    dS = (dH_fit - dG) / T_kelvin if np.isfinite(dG) else np.nan
    
    # RMSE
    rmse = np.sqrt(np.mean(residuals**2))
    
    # Standard errors for derived parameters
    Kd_err = perr[0] if len(perr) > 0 else np.nan
    dH_err = perr[1] if len(perr) > 1 else np.nan
    n_err = perr[2] if len(perr) > 2 else np.nan
    
    return {
        'Kd': float(Kd_fit),
        'Kd_err': float(Kd_err),
        'dH': float(dH_fit),
        'dH_err': float(dH_err),
        'n': float(n_fit),
        'n_err': float(n_err),
        'Ka': float(Ka),
        'dG': float(dG),
        'dS': float(dS),
        'T_kelvin': float(T_kelvin),
        'r_squared': float(r_squared),
        'chi2_reduced': float(chi2_reduced),
        'rmse': float(rmse),
        'residuals': residuals.tolist(),
        'heat_fitted': heat_fitted.tolist(),
        'heat_observed': heat_observed.tolist(),
        'injections': injections.tolist(),
        'volumes': volumes.tolist(),
        'converged': True,
        'dof': int(dof),
        'dilution_heat': float(dH_dil_fit) if fit_dilution else None,
        'method': 'curve_fit' if curve_fit_ok else 'differential_evolution'
    }


def format_result(result):
    """Format analysis result as human-readable text."""
    if not result.get('converged', False):
        return f"ITC Analysis Failed: {result.get('error', 'Unknown error')}"
    
    lines = []
    lines.append("=" * 60)
    lines.append("ITC One-Site Binding Model Analysis Results")
    lines.append("=" * 60)
    lines.append("")
    
    # Fitted parameters
    lines.append("Fitted Parameters:")
    lines.append(f"  Kd (Dissociation constant)  = {result['Kd']:.4e} M")
    lines.append(f"  Kd standard error           = {result['Kd_err']:.4e} M")
    lines.append(f"  ΔH (Enthalpy change)        = {result['dH']:.2f} J/mol = {result['dH']/1000:.2f} kJ/mol")
    lines.append(f"  ΔH standard error           = {result['dH_err']:.2f} J/mol")
    lines.append(f"  n (Stoichiometry)           = {result['n']:.3f}")
    lines.append(f"  n standard error            = {result['n_err']:.3f}")
    lines.append("")
    
    # Derived thermodynamic parameters
    lines.append("Derived Thermodynamic Parameters:")
    lines.append(f"  Ka (Association constant)   = {result['Ka']:.4e} M⁻¹")
    lines.append(f"  ΔG (Gibbs free energy)      = {result['dG']:.2f} J/mol = {result['dG']/1000:.2f} kJ/mol")
    lines.append(f"  ΔS (Entropy change)         = {result['dS']:.2f} J/(mol·K) = {result['dS']/1000:.3f} kJ/(mol·K)")
    lines.append(f"  -TΔS                        = {-result['T_kelvin'] * result['dS']:.2f} J/mol")
    lines.append("")
    
    # Fit quality
    lines.append("Fit Quality:")
    lines.append(f"  R² (coefficient of determination) = {result['r_squared']:.4f}")
    lines.append(f"  Reduced χ²                       = {result['chi2_reduced']:.4f}")
    lines.append(f"  RMSE (root mean square error)     = {result['rmse']:.4f}")
    lines.append(f"  Degrees of freedom               = {result['dof']}")
    lines.append("")
    
    # Data summary
    lines.append("Data Summary:")
    lines.append(f"  Number of injections: {len(result['injections'])}")
    lines.append(f"  Temperature: {result['T_kelvin']:.2f} K ({result['T_kelvin']-273.15:.1f} °C)")
    lines.append("")
    
    lines.append("=" * 60)
    
    return "\n".join(lines)


def itc_fit_plot_png_base64(result: dict) -> Optional[str]:
    """将一结合位点拟合的观测/模型差分热曲线绘为 PNG，返回 Base64（ASCII），失败则 None。"""
    if not result.get("converged"):
        return None
    try:
        inj = np.asarray(result["injections"], dtype=float)
        y_obs = np.asarray(result["heat_observed"], dtype=float)
        y_fit = np.asarray(result["heat_fitted"], dtype=float)
        fig, ax = plt.subplots(figsize=(7.5, 4.0))
        ax.plot(inj, y_obs, "o", ms=5, label="Observed ΔQ")
        ax.plot(inj, y_fit, "-", lw=2, label="Fitted (one-site)")
        ax.set_xlabel("Injection #")
        ax.set_ylabel("Differential heat")
        ax.legend(loc="best")
        ax.grid(True, alpha=0.3)
        ax.set_title("ITC: observed vs one-site model")
        buf = io.BytesIO()
        fig.savefig(buf, format="png", dpi=120, bbox_inches="tight")
        plt.close(fig)
        return base64.b64encode(buf.getvalue()).decode("ascii")
    except Exception:
        return None


def main():
    parser = argparse.ArgumentParser(description='ITC Data Analysis - One-Site Binding Model')
    parser.add_argument('--data', type=str, help='JSON array of [injection, volume, heat] triples')
    parser.add_argument('--file', type=str, help='Path to CSV/TSV file with injection, volume, heat columns')
    parser.add_argument('--protein-conc', type=float, required=True, help='Protein concentration in cell (M)')
    parser.add_argument('--ligand-conc', type=float, required=True, help='Ligand concentration in syringe (M)')
    parser.add_argument('--cell-volume', type=float, default=1.4e-3, help='Cell volume in L (default: 1.4 mL)')
    parser.add_argument('--temperature', type=float, default=25.0, help='Temperature in Celsius (default: 25)')
    parser.add_argument('--volume-unit', type=str, default='L', choices=['L', 'mL', 'uL', 'µL'],
                        help='Unit for injection volumes (default: L)')
    parser.add_argument('--output-format', type=str, default='text', choices=['text', 'json'], 
                        help='Output format')
    parser.add_argument('--output-file', type=str, help='Save output to file')
    
    args = parser.parse_args()
    
    # Load data
    if args.data:
        raw_data = parse_data_array(args.data)
    elif args.file:
        raw_data = load_csv_data(args.file)
    else:
        print("Error: Must provide either --data or --file", file=sys.stderr)
        sys.exit(1)
    
    if raw_data.shape[1] < 3:
        print("Error: Data must have at least 3 columns (injection, volume, heat)", file=sys.stderr)
        sys.exit(1)
    
    injections = raw_data[:, 0]
    volumes = raw_data[:, 1]
    heat = raw_data[:, 2]
    
    # Convert volume units to liters
    volume_conversion = {
        'L': 1.0,
        'mL': 1e-3,
        'uL': 1e-6,
        'µL': 1e-6
    }
    volumes = volumes * volume_conversion.get(args.volume_unit, 1.0)
    
    # Check minimum data requirements
    if len(injections) < 4:
        print("Warning: At least 4 injections recommended for reliable fitting.", file=sys.stderr)
        print(f"Current data has {len(injections)} injection(s). Results may be unreliable.", file=sys.stderr)
    
    # Run analysis
    result = fit_itc_model(
        injections=injections,
        volumes=volumes,
        heat_observed=heat,
        P0=args.protein_conc,
        L0=args.ligand_conc,
        V_cell=args.cell_volume,
        temperature_celsius=args.temperature
    )
    
    # Output
    if args.output_format == 'json':
        payload = dict(result)
        b64p = itc_fit_plot_png_base64(result)
        if b64p:
            payload["fit_plot_png_base64"] = b64p
        output = json.dumps(payload, indent=2, default=str)
    else:
        output = format_result(result)
    
    if args.output_file:
        with open(args.output_file, 'w') as f:
            f.write(output)
        print(f"Results saved to {args.output_file}")
    else:
        print(output)
    
    return 0 if result.get('converged', False) else 1


if __name__ == '__main__':
    sys.exit(main())
