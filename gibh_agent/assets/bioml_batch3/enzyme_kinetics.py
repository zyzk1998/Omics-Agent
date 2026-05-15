#!/usr/bin/env python3
"""
enzyme_kinetics.py

Enzyme kinetics experiment analysis tool.
Simulates enzyme reactions under various substrate concentrations and modulators,
generates time-course data, substrate kinetics curves, and modulator dose-response curves.
Fits Michaelis-Menten parameters and outputs CSV data files.

Dependencies: numpy, scipy, pandas, matplotlib

Usage:
    python enzyme_kinetics.py --json-input '{"enzyme_name":"Trypsin",...}' --output-dir /tmp/results
"""

import argparse
import base64
import csv
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
# Constants
# ---------------------------------------------------------------------------
DEFAULT_TIME_POINTS = [0, 5, 10, 15, 20, 30, 45, 60]
NOISE_LEVEL = 0.03  # 3% relative noise


# ---------------------------------------------------------------------------
# Kinetic Models
# ---------------------------------------------------------------------------
def michaelis_menten(S, Vmax, Km):
    """v = Vmax * [S] / (Km + [S])"""
    return (Vmax * S) / (Km + S)


def competitive_inhibition(S, Vmax, Km, I, Ki):
    """Competitive inhibition: v = Vmax*[S] / (Km*(1+[I]/Ki) + [S])"""
    return (Vmax * S) / (Km * (1 + I / Ki) + S)


def noncompetitive_inhibition(S, Vmax, Km, I, Ki):
    """Non-competitive inhibition: v = Vmax*[S] / ((Km+[S])*(1+[I]/Ki))"""
    return (Vmax * S) / ((Km + S) * (1 + I / Ki))


def uncompetitive_inhibition(S, Vmax, Km, I, Ki):
    """Uncompetitive inhibition: v = Vmax*[S] / (Km + [S]*(1+[I]/Ki))"""
    return (Vmax * S) / (Km + S * (1 + I / Ki))


def mixed_inhibition(S, Vmax, Km, I, Ki, alpha=1.0):
    """Mixed inhibition with alpha factor"""
    return (Vmax * S) / (Km * (1 + I / Ki) + S * (1 + I / (alpha * Ki)))


def activation_vmax(S, Vmax, Km, A, Ka):
    """Activator increases Vmax: Vmax_eff = Vmax * (1 + [A]/(Ka + [A]))"""
    Vmax_eff = Vmax * (1 + A / (Ka + A)) if Ka > 0 else Vmax
    return (Vmax_eff * S) / (Km + S)


def activation_km(S, Vmax, Km, A, Ka):
    """Activator decreases Km: Km_eff = Km / (1 + [A]/Ka)"""
    Km_eff = Km / (1 + A / Ka) if Ka > 0 else Km
    return (Vmax * S) / (Km_eff + S)


def allosteric_michaelis_menten(S, Vmax, Km, h):
    """Hill equation for allosteric enzymes: v = Vmax*[S]^h / (Km^h + [S]^h)"""
    return (Vmax * S**h) / (Km**h + S**h)


def allosteric_activation(S, Vmax, Km, h, A, Ka):
    """Allosteric activator: decreases Km and increases cooperativity"""
    if Ka <= 0 or A <= 0:
        return allosteric_michaelis_menten(S, Vmax, Km, h)
    Km_eff = Km / (1 + A / Ka)
    h_eff = h * (1 + A / (2 * Ka))
    return (Vmax * S**h_eff) / (Km_eff**h_eff + S**h_eff)


def allosteric_inhibition(S, Vmax, Km, h, I, Ki):
    """Allosteric inhibitor: increases Km and decreases cooperativity"""
    if Ki <= 0 or I <= 0:
        return allosteric_michaelis_menten(S, Vmax, Km, h)
    Km_eff = Km * (1 + I / Ki)
    h_eff = h / (1 + I / (2 * Ki))
    return (Vmax * S**h_eff) / (Km_eff**h_eff + S**h_eff)




def competitive_product_inhibition(S, P, Vmax, Km, KiP):
    """Product inhibition (competitive): v = Vmax*[S] / (Km*(1+[P]/KiP) + [S])"""
    return (Vmax * S) / (Km * (1 + P / KiP) + S)


def noncompetitive_product_inhibition(S, P, Vmax, Km, KiP):
    """Product inhibition (non-competitive): v = Vmax*[S] / ((Km+[S])*(1+[P]/KiP))"""
    return (Vmax * S) / ((Km + S) * (1 + P / KiP))


def uncompetitive_product_inhibition(S, P, Vmax, Km, KiP):
    """Product inhibition (uncompetitive): v = Vmax*[S] / (Km + [S]*(1+[P]/KiP))"""
    return (Vmax * S) / (Km + S * (1 + P / KiP))


def arrhenius_temperature_factor(T, T_ref, Ea):
    """
    Arrhenius temperature dependence of kcat.
    factor = exp(-Ea/R * (1/T - 1/T_ref))
    T in Kelvin, Ea in J/mol, R = 8.314 J/(mol*K)
    """
    R = 8.314
    return np.exp(-Ea / R * (1.0 / T - 1.0 / T_ref))


def bell_shaped_pH_activity(pH, pH_opt, width):
    """
    Bell-shaped pH activity profile.
    activity = exp(-((pH - pH_opt)/width)^2)
    """
    return np.exp(-((pH - pH_opt) / width) ** 2)


def mwc_model(S, Vmax, K_R, K_T, L, c, n):
    """
    Monod-Wyman-Changeux (MWC) model for allosteric enzymes.
    
    Parameters:
        S: substrate concentration
        Vmax: maximum rate
        K_R: dissociation constant for R state
        K_T: dissociation constant for T state
        L: allosteric constant (T0/R0)
        c: K_R / K_T (cooperativity parameter)
        n: number of subunits
    
    The fraction of R state:
        R_frac = (1 + S/K_R)^n / ((1 + S/K_R)^n + L * (1 + S/K_T)^n)
    
    Rate:
        v = Vmax * R_frac * S / (K_R + S)
    
    Note: Simplified to use R_frac as effective fraction of active enzyme.
    """
    term_R = (1 + S / K_R) ** n
    term_T = (1 + S / K_T) ** n
    R_frac = term_R / (term_R + L * term_T)
    return Vmax * R_frac * S / (K_R + S)


def mwc_with_activator(S, Vmax, K_R, K_T, L, c, n, A, Ka):
    """MWC model with activator: activator binds R state, reducing L"""
    if A <= 0 or Ka <= 0:
        return mwc_model(S, Vmax, K_R, K_T, L, c, n)
    L_eff = L / (1 + A / Ka)
    return mwc_model(S, Vmax, K_R, K_T, L_eff, c, n)


def mwc_with_inhibitor(S, Vmax, K_R, K_T, L, c, n, I, Ki):
    """MWC model with inhibitor: inhibitor binds T state, increasing L"""
    if I <= 0 or Ki <= 0:
        return mwc_model(S, Vmax, K_R, K_T, L, c, n)
    L_eff = L * (1 + I / Ki)
    return mwc_model(S, Vmax, K_R, K_T, L_eff, c, n)

def _modulated_rate(S, Vmax, Km, modulator_type, modulator_conc, Ki_or_Ka,
                     hill_n=1.0, product_conc=0.0, KiP=None, T_kelvin=None,
                     T_ref=None, Ea=None, pH=None, pH_opt=None, pH_width=None,
                     mwc_params=None, **kwargs):
    """Calculate reaction rate with modulator. Unified dispatcher."""
    # Temperature correction
    if T_kelvin is not None and T_ref is not None and Ea is not None:
        Vmax = Vmax * arrhenius_temperature_factor(T_kelvin, T_ref, Ea)
    
    # pH correction
    if pH is not None and pH_opt is not None and pH_width is not None:
        Vmax = Vmax * bell_shaped_pH_activity(pH, pH_opt, pH_width)
    
    # Product inhibition
    if product_conc > 0 and KiP is not None and KiP > 0:
        if modulator_type == "competitive_product_inhibition":
            return competitive_product_inhibition(S, product_conc, Vmax, Km, KiP)
        elif modulator_type == "noncompetitive_product_inhibition":
            return noncompetitive_product_inhibition(S, product_conc, Vmax, Km, KiP)
        elif modulator_type == "uncompetitive_product_inhibition":
            return uncompetitive_product_inhibition(S, product_conc, Vmax, Km, KiP)
        else:
            # Default: competitive product inhibition
            return competitive_product_inhibition(S, product_conc, Vmax, Km, KiP)
    
    # MWC allosteric model
    if mwc_params is not None:
        K_R = mwc_params.get("K_R", Km)
        K_T = mwc_params.get("K_T", Km * 100)
        L = mwc_params.get("L", 100)
        c = mwc_params.get("c", 0.01)
        n = mwc_params.get("n", 4)
        if modulator_type == "mwc_activator":
            return mwc_with_activator(S, Vmax, K_R, K_T, L, c, n, modulator_conc, Ki_or_Ka)
        elif modulator_type == "mwc_inhibitor":
            return mwc_with_inhibitor(S, Vmax, K_R, K_T, L, c, n, modulator_conc, Ki_or_Ka)
        else:
            return mwc_model(S, Vmax, K_R, K_T, L, c, n)
    
    if modulator_conc <= 0 or modulator_type == "none":
        if hill_n != 1.0:
            return allosteric_michaelis_menten(S, Vmax, Km, hill_n)
        return michaelis_menten(S, Vmax, Km)
    if modulator_type in ("competitive_inhibition", "inhibitor"):
        return competitive_inhibition(S, Vmax, Km, modulator_conc, Ki_or_Ka)
    elif modulator_type == "noncompetitive_inhibition":
        return noncompetitive_inhibition(S, Vmax, Km, modulator_conc, Ki_or_Ka)
    elif modulator_type == "uncompetitive_inhibition":
        return uncompetitive_inhibition(S, Vmax, Km, modulator_conc, Ki_or_Ka)
    elif modulator_type == "mixed_inhibition":
        return mixed_inhibition(S, Vmax, Km, modulator_conc, Ki_or_Ka)
    elif modulator_type in ("activation_vmax", "activator"):
        return activation_vmax(S, Vmax, Km, modulator_conc, Ki_or_Ka)
    elif modulator_type == "activation_km":
        return activation_km(S, Vmax, Km, modulator_conc, Ki_or_Ka)
    elif modulator_type == "allosteric_activation":
        return allosteric_activation(S, Vmax, Km, hill_n, modulator_conc, Ki_or_Ka)
    elif modulator_type == "allosteric_inhibition":
        return allosteric_inhibition(S, Vmax, Km, hill_n, modulator_conc, Ki_or_Ka)
    else:
        if hill_n != 1.0:
            return allosteric_michaelis_menten(S, Vmax, Km, hill_n)
        return michaelis_menten(S, Vmax, Km)


# ---------------------------------------------------------------------------
# Simulation Engine
# ---------------------------------------------------------------------------
def simulate_time_course(
    time_points: np.ndarray,
    S0: float,
    E_conc: float,
    Vmax: float,
    Km: float,
    modulator_name: Optional[str] = None,
    modulator_conc: float = 0.0,
    modulator_type: str = "none",
    Ki_or_Ka: float = 1.0,
    noise: bool = True,
    product_inhibition: bool = False,
    KiP: float = None,
    **env_kwargs,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Simulate substrate and product concentrations over time via numerical integration.
    
    Uses Euler method with adaptive small steps.
    Returns (substrate_concentrations, product_concentrations) at each time point.
    """
    dt = 0.01  # small time step in minutes
    S = S0
    P = 0.0
    S_history = [S]
    P_history = [P]
    t_current = 0.0
    
    for t_target in time_points[1:]:
        while t_current < t_target:
            step = min(dt, t_target - t_current)
            
            # Calculate current rate
            v = _modulated_rate(S, Vmax, Km, modulator_type, modulator_conc, Ki_or_Ka,
                                product_conc=P if product_inhibition else 0.0,
                                KiP=KiP, **env_kwargs)
            
            # Euler step: dS/dt = -v, dP/dt = +v
            S = max(0.0, S - v * step)
            P = P + v * step
            t_current += step
        
        S_history.append(S)
        P_history.append(P)
    
    S_result = np.array(S_history)
    P_result = np.array(P_history)
    
    # Add noise
    if noise:
        S_result = S_result * (1 + np.random.normal(0, NOISE_LEVEL, len(S_result)))
        S_result = np.maximum(S_result, 0.0)
        P_result = P_result * (1 + np.random.normal(0, NOISE_LEVEL, len(P_result)))
        P_result = np.maximum(P_result, 0.0)
    
    return S_result, P_result


def simulate_initial_rates(
    substrate_concentrations: List[float],
    E_conc: float,
    Vmax: float,
    Km: float,
    modulator_name: Optional[str] = None,
    modulator_conc: float = 0.0,
    modulator_type: str = "none",
    Ki_or_Ka: float = 1.0,
    noise: bool = True,
    **env_kwargs,
) -> np.ndarray:
    """Simulate initial reaction rates at various substrate concentrations."""
    rates = []
    for S in substrate_concentrations:
        v = _modulated_rate(S, Vmax, Km, modulator_type, modulator_conc, Ki_or_Ka, **env_kwargs)
        
        if noise:
            v = v * (1 + np.random.normal(0, NOISE_LEVEL))
            v = max(v, 1e-6)
        
        rates.append(v)
    
    return np.array(rates)


def simulate_dose_response(
    modulator_concentrations: List[float],
    S_ref: float,
    E_conc: float,
    Vmax: float,
    Km: float,
    modulator_type: str = "none",
    Ki_or_Ka: float = 1.0,
    noise: bool = True,
    **env_kwargs,
) -> np.ndarray:
    """Simulate activity at various modulator concentrations (fixed substrate)."""
    activities = []
    for conc in modulator_concentrations:
        v = _modulated_rate(S_ref, Vmax, Km, modulator_type, conc, Ki_or_Ka, **env_kwargs)
        
        if noise:
            v = v * (1 + np.random.normal(0, NOISE_LEVEL))
            v = max(v, 1e-6)
        
        activities.append(v)
    
    return np.array(activities)


# ---------------------------------------------------------------------------
# Parameter Estimation
# ---------------------------------------------------------------------------
def fit_michaelis_menten(S_data: np.ndarray, v_data: np.ndarray) -> Dict[str, Any]:
    """Fit MM model to substrate-velocity data."""
    # Filter positive
    mask = v_data > 1e-12
    if np.sum(mask) < 3:
        return {"Vmax": None, "Km": None, "r_squared": 0.0, "error": "Insufficient data"}
    
    S_fit = S_data[mask]
    v_fit = v_data[mask]
    
    Vmax_est = np.max(v_fit)
    Km_est = S_fit[np.argmin(np.abs(v_fit - Vmax_est / 2))] if len(S_fit) > 0 else np.median(S_fit)
    
    bounds = ([1e-12, 1e-12], [np.inf, np.inf])
    p0 = [Vmax_est, Km_est]
    
    try:
        popt, pcov = curve_fit(
            michaelis_menten, S_fit, v_fit, p0=p0, bounds=bounds, maxfev=10000
        )
        Vmax, Km = popt
        perr = np.sqrt(np.diag(pcov))
        
        v_pred = michaelis_menten(S_fit, Vmax, Km)
        ss_res = np.sum((v_fit - v_pred) ** 2)
        ss_tot = np.sum((v_fit - np.mean(v_fit)) ** 2)
        r_squared = 1 - ss_res / ss_tot if ss_tot > 0 else 0.0
        
        return {
            "Vmax": float(Vmax),
            "Vmax_std_error": float(perr[0]),
            "Km": float(Km),
            "Km_std_error": float(perr[1]),
            "r_squared": float(r_squared),
            "num_points": int(len(S_fit)),
        }
    except Exception as e:
        return {"Vmax": None, "Km": None, "r_squared": 0.0, "error": str(e)}


# ---------------------------------------------------------------------------
# Visualization
# ---------------------------------------------------------------------------
def plot_substrate_kinetics(
    substrate_concentrations: np.ndarray,
    rates: np.ndarray,
    fit_params: Optional[Dict],
    output_path: str,
    title: str = "Substrate Kinetics",
) -> str:
    """Plot substrate kinetics with MM fit."""
    fig, ax = plt.subplots(figsize=(8, 6))
    
    ax.scatter(substrate_concentrations, rates, c="red", s=80, zorder=5, label="Simulated data")
    
    if fit_params and fit_params.get("Vmax") is not None:
        Vmax = fit_params["Vmax"]
        Km = fit_params["Km"]
        S_smooth = np.linspace(0, max(substrate_concentrations) * 1.2, 200)
        v_smooth = michaelis_menten(S_smooth, Vmax, Km)
        ax.plot(S_smooth, v_smooth, "b-", linewidth=2, 
                label=f"Fit: Vmax={Vmax:.2f}, Km={Km:.2f} μM")
        
        ax.axhline(Vmax / 2, color="gray", linestyle="--", alpha=0.5)
        ax.axvline(Km, color="gray", linestyle="--", alpha=0.5)
        ax.annotate(f"Km = {Km:.1f}", xy=(Km, Vmax / 2), 
                     xytext=(Km * 1.3, Vmax * 0.4), fontsize=9)
    
    ax.set_xlabel("Substrate Concentration [S] (μM)", fontsize=11)
    ax.set_ylabel("Initial Rate v (a.u./min)", fontsize=11)
    ax.set_title(title, fontsize=12)
    ax.legend()
    ax.set_xlim(left=0)
    ax.set_ylim(bottom=0)
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close()
    return output_path


def plot_time_courses(
    time_points: np.ndarray,
    substrate_concentrations: List[float],
    time_course_data: Dict[str, np.ndarray],
    output_path: str,
    title: str = "Time Courses",
    product_data: Dict[str, np.ndarray] = None,
) -> str:
    """Plot substrate consumption and product formation over time."""
    if product_data:
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))
    else:
        fig, ax1 = plt.subplots(figsize=(10, 6))
        ax2 = None
    
    colors = plt.cm.viridis(np.linspace(0, 1, len(substrate_concentrations)))
    
    for i, (S0, label) in enumerate(zip(substrate_concentrations, time_course_data.keys())):
        data = time_course_data[label]
        pct = (data / S0) * 100
        ax1.plot(time_points, pct, "o-", color=colors[i], linewidth=1.5,
                markersize=4, label=f"[{S0}] μM")
    
    ax1.set_xlabel("Time (min)", fontsize=11)
    ax1.set_ylabel("Substrate Remaining (%)", fontsize=11)
    ax1.set_title(f"{title} - Substrate", fontsize=12)
    ax1.legend(title="Substrate", loc="upper right")
    ax1.set_ylim(0, 105)
    ax1.grid(True, alpha=0.3)
    
    if ax2 and product_data:
        for i, (S0, label) in enumerate(zip(substrate_concentrations, time_course_data.keys())):
            p_data = product_data[label]
            ax2.plot(time_points, p_data, "s--", color=colors[i], linewidth=1.5,
                    markersize=4, label=f"[{S0}] μM")
        ax2.set_xlabel("Time (min)", fontsize=11)
        ax2.set_ylabel("Product Formed (μM)", fontsize=11)
        ax2.set_title(f"{title} - Product", fontsize=12)
        ax2.legend(title="Substrate", loc="lower right")
        ax2.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close()
    return output_path


def plot_dose_response(
    modulator_concentrations: np.ndarray,
    activities: np.ndarray,
    modulator_name: str,
    modulator_type: str,
    output_path: str,
) -> str:
    """Plot modulator dose-response curve."""
    fig, ax = plt.subplots(figsize=(8, 6))
    
    # For log scale, filter out zero
    conc_plot = np.where(modulator_concentrations <= 0, 1e-3, modulator_concentrations)
    
    ax.semilogx(conc_plot, activities, "bo-", linewidth=2, markersize=6)
    ax.scatter(conc_plot, activities, c="red", s=60, zorder=5)
    
    ax.set_xlabel(f"{modulator_name} Concentration (μM)", fontsize=11)
    ax.set_ylabel("Activity (a.u./min)", fontsize=11)
    ax.set_title(f"{modulator_type.title()} Dose-Response: {modulator_name}", fontsize=12)
    ax.grid(True, alpha=0.3, which="both")
    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close()
    return output_path


# ---------------------------------------------------------------------------
# CSV Data Generation
# ---------------------------------------------------------------------------
def write_time_course_csv(
    time_points: np.ndarray,
    substrate_concentrations: List[float],
    time_course_data: Dict[str, np.ndarray],
    output_path: str,
    product_data: Dict[str, np.ndarray] = None,
) -> str:
    """Write time-course data to CSV (substrate + optional product)."""
    with open(output_path, "w", newline="") as f:
        writer = csv.writer(f)
        # Header
        header = ["Time_min"]
        for S0 in substrate_concentrations:
            header.append(f"S_{S0}uM_remaining")
        if product_data:
            for S0 in substrate_concentrations:
                header.append(f"P_{S0}uM_produced")
        writer.writerow(header)
        # Data
        for i, t in enumerate(time_points):
            row = [t]
            for S0 in substrate_concentrations:
                label = f"S_{S0}"
                row.append(f"{time_course_data[label][i]:.4f}")
            if product_data:
                for S0 in substrate_concentrations:
                    label = f"S_{S0}"
                    row.append(f"{product_data[label][i]:.4f}")
            writer.writerow(row)
    return output_path


def write_substrate_kinetics_csv(
    substrate_concentrations: List[float],
    rates: np.ndarray,
    output_path: str,
) -> str:
    """Write substrate kinetics data to CSV."""
    with open(output_path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["Substrate_uM", "Initial_Rate_a.u._per_min"])
        for S, v in zip(substrate_concentrations, rates):
            writer.writerow([S, f"{v:.6f}"])
    return output_path


def write_dose_response_csv(
    modulator_concentrations: List[float],
    activities: np.ndarray,
    modulator_name: str,
    output_path: str,
) -> str:
    """Write dose-response data to CSV."""
    with open(output_path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow([f"{modulator_name}_uM", "Activity_a.u._per_min"])
        for conc, act in zip(modulator_concentrations, activities):
            writer.writerow([conc, f"{act:.6f}"])
    return output_path


# ---------------------------------------------------------------------------
# Report Generation
# ---------------------------------------------------------------------------
def generate_report(
    enzyme_name: str,
    enzyme_concentration: float,
    substrate_concentrations: List[float],
    fit_params: Dict,
    modulators: Dict,
    output_path: str,
    temperature: float = None,
    pH: float = None,
    product_inhibition: bool = False,
    KiP: float = None,
    mwc_params: Dict = None,
) -> str:
    """Generate detailed analysis report."""
    lines = []
    lines.append("=" * 70)
    lines.append("ENZYME KINETICS EXPERIMENT ANALYSIS REPORT")
    lines.append("=" * 70)
    lines.append("")
    
    lines.append(f"Enzyme: {enzyme_name}")
    lines.append(f"Enzyme Concentration: {enzyme_concentration} nM")
    lines.append(f"Substrate Concentrations: {substrate_concentrations} μM")
    if temperature is not None:
        lines.append(f"Temperature: {temperature}°C")
    if pH is not None:
        lines.append(f"pH: {pH}")
    lines.append("")
    
    if fit_params and fit_params.get("Vmax") is not None:
        lines.append("MICHAELIS-MENTEN FIT RESULTS")
        lines.append("-" * 70)
        lines.append(f"Vmax = {fit_params['Vmax']:.4f} ± {fit_params.get('Vmax_std_error', 0):.4f} a.u./min")
        lines.append(f"Km   = {fit_params['Km']:.4f} ± {fit_params.get('Km_std_error', 0):.4f} μM")
        lines.append(f"R²   = {fit_params['r_squared']:.4f}")
        lines.append(f"Points fitted: {fit_params.get('num_points', 0)}")
        
        # kcat
        E_uM = enzyme_concentration / 1000.0  # nM -> μM
        kcat = fit_params['Vmax'] / E_uM if E_uM > 0 else 0
        lines.append(f"kcat = {kcat:.2f} min⁻¹")
        lines.append(f"kcat/Km = {kcat / fit_params['Km']:.4f} μM⁻¹·min⁻¹")
        lines.append("")
    
    if mwc_params:
        lines.append("MWC ALLOSTERIC MODEL")
        lines.append("-" * 70)
        for k, v in mwc_params.items():
            lines.append(f"  {k}: {v}")
        lines.append("")
    
    if product_inhibition and KiP is not None:
        lines.append("PRODUCT INHIBITION")
        lines.append("-" * 70)
        lines.append(f"  KiP = {KiP:.2f} μM")
        lines.append("")
    
    if modulators:
        lines.append("MODULATORS")
        lines.append("-" * 70)
        for name, concs in modulators.items():
            lines.append(f"  {name}: concentrations {concs} μM")
        lines.append("")
    
    lines.append("OUTPUT FILES")
    lines.append("-" * 70)
    lines.append("  - time_course_data.csv: Substrate consumption + product formation over time")
    lines.append("  - substrate_kinetics.csv: Initial rate vs [S]")
    lines.append("  - dose_response_*.csv: Modulator dose-response curves")
    lines.append("  - Plots: PNG files for each analysis")
    lines.append("")
    
    lines.append("=" * 70)
    lines.append("END OF REPORT")
    lines.append("=" * 70)
    
    with open(output_path, "w") as f:
        f.write("\n".join(lines))
    return output_path


def generate_json_metadata(
    enzyme_name: str,
    enzyme_concentration: float,
    substrate_concentrations: List[float],
    time_points: List[float],
    fit_params: Dict,
    modulators: Dict,
    modulator_config: Dict,
    file_paths: Dict,
    output_path: str,
    temperature: float = None,
    pH: float = None,
    product_inhibition: bool = False,
    KiP: float = None,
    mwc_params: Dict = None,
) -> str:
    """Generate JSON metadata file."""
    E_uM = enzyme_concentration / 1000.0
    kcat = fit_params.get("Vmax", 0) / E_uM if E_uM > 0 and fit_params.get("Vmax") else 0
    
    result = {
        "enzyme_name": enzyme_name,
        "enzyme_concentration_nM": enzyme_concentration,
        "substrate_concentrations_uM": substrate_concentrations,
        "time_points_min": time_points,
        "environmental_conditions": {
            "temperature_C": temperature,
            "pH": pH,
            "product_inhibition": product_inhibition,
            "KiP_uM": KiP,
            "mwc_params": mwc_params,
        },
        "michaelis_menten_fit": fit_params,
        "catalytic_parameters": {
            "kcat_min-1": float(kcat) if fit_params.get("Vmax") else None,
            "kcat_over_Km_uM-1_min-1": float(kcat / fit_params["Km"]) if fit_params.get("Km") else None,
        },
        "modulators": modulators,
        "modulator_config": modulator_config,
        "output_files": file_paths,
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
    
    Simulates enzyme kinetics data and generates:
      - Time-course CSV (substrate consumption)
      - Substrate kinetics CSV (v vs [S])
      - Dose-response CSVs (one per modulator)
      - Plots for each dataset
      - Analysis report and JSON metadata
    """
    log_lines = []
    np.random.seed(42)  # Reproducible simulations
    
    try:
        # 1. Parse input
        log_lines.append("Step 1: Parsing input data...")
        enzyme_name = input_data.get("enzyme_name", "Unknown")
        substrate_concentrations = input_data.get("substrate_concentrations", [])
        enzyme_concentration_nM = float(input_data.get("enzyme_concentration", 50))
        modulators = input_data.get("modulators", {})
        modulator_config = input_data.get("modulator_config", {})
        time_points = input_data.get("time_points", DEFAULT_TIME_POINTS)
        
        if not substrate_concentrations:
            raise ValueError("substrate_concentrations is required")
        
        log_lines.append(f"  Enzyme: {enzyme_name} ({enzyme_concentration_nM} nM)")
        log_lines.append(f"  Substrates: {len(substrate_concentrations)} concentrations")
        log_lines.append(f"  Time points: {len(time_points)} points")
        if modulators:
            log_lines.append(f"  Modulators: {len(modulators)}")
        
        os.makedirs(output_dir, exist_ok=True)
        
        # Default kinetic parameters (could be estimated from data in future)
        # Using reasonable defaults for a typical protease
        Vmax_base = 10.0  # a.u./min
        Km_base = max(substrate_concentrations) / 5.0  # μM, estimated
        
        time_arr = np.array(time_points, dtype=float)
        substrate_arr = np.array(substrate_concentrations, dtype=float)
        
        file_paths = {}
        
        # 2. Parse environmental parameters
        env_params = {}
        temperature = input_data.get("temperature")
        pH = input_data.get("pH")
        product_inhibition = input_data.get("product_inhibition", False)
        KiP = input_data.get("KiP")
        mwc_params = input_data.get("mwc_params")
        
        if temperature is not None:
            T_k = float(temperature) + 273.15
            env_params["T_kelvin"] = T_k
            env_params["T_ref"] = 298.15
            env_params["Ea"] = 50000.0  # J/mol, typical for enzymes
            log_lines.append(f"  Temperature: {temperature}°C (T={T_k:.1f}K)")
        if pH is not None:
            env_params["pH"] = float(pH)
            env_params["pH_opt"] = 7.5
            env_params["pH_width"] = 2.0
            log_lines.append(f"  pH: {pH} (opt={env_params['pH_opt']}, width={env_params['pH_width']})")
        if product_inhibition:
            env_params["product_inhibition"] = True
            env_params["KiP"] = KiP if KiP is not None else Km_base * 5.0
            log_lines.append(f"  Product inhibition: KiP={env_params['KiP']:.2f} μM")
        if mwc_params:
            env_params["mwc_params"] = mwc_params
            log_lines.append(f"  MWC model: {mwc_params}")
        
        # 3. Simulate time courses
        log_lines.append("")
        log_lines.append("Step 2: Simulating time-course data...")
        time_course_data = {}
        product_data = {}
        for S0 in substrate_concentrations:
            label = f"S_{S0}"
            S_remaining, P_produced = simulate_time_course(
                time_arr, S0, enzyme_concentration_nM, Vmax_base, Km_base,
                **env_params
            )
            time_course_data[label] = S_remaining
            product_data[label] = P_produced
        
        tc_csv = os.path.join(output_dir, "time_course_data.csv")
        write_time_course_csv(time_arr, substrate_concentrations, time_course_data, tc_csv, product_data)
        file_paths["time_course_csv"] = tc_csv
        log_lines.append(f"  - Time-course CSV: {tc_csv}")
        
        tc_plot = os.path.join(output_dir, "time_courses.png")
        plot_time_courses(time_arr, substrate_concentrations, time_course_data, tc_plot, product_data=product_data)
        file_paths["time_course_plot"] = tc_plot
        log_lines.append(f"  - Time-course plot: {tc_plot}")
        
        # 3. Simulate substrate kinetics (no modulator)
        log_lines.append("")
        log_lines.append("Step 3: Simulating substrate kinetics...")
        rates = simulate_initial_rates(substrate_concentrations, enzyme_concentration_nM, Vmax_base, Km_base, **env_params)
        
        sk_csv = os.path.join(output_dir, "substrate_kinetics.csv")
        write_substrate_kinetics_csv(substrate_concentrations, rates, sk_csv)
        file_paths["substrate_kinetics_csv"] = sk_csv
        log_lines.append(f"  - Substrate kinetics CSV: {sk_csv}")
        
        # 4. Fit Michaelis-Menten
        log_lines.append("")
        log_lines.append("Step 4: Fitting Michaelis-Menten model...")
        fit_params = fit_michaelis_menten(substrate_arr, rates)
        
        if fit_params.get("Vmax") is not None:
            log_lines.append(f"  Vmax = {fit_params['Vmax']:.4f} a.u./min")
            log_lines.append(f"  Km   = {fit_params['Km']:.4f} μM")
            log_lines.append(f"  R²   = {fit_params['r_squared']:.4f}")
        else:
            log_lines.append(f"  Fit failed: {fit_params.get('error', 'unknown error')}")
        
        sk_plot = os.path.join(output_dir, "substrate_kinetics.png")
        plot_substrate_kinetics(substrate_arr, rates, fit_params if fit_params.get("Vmax") else None, sk_plot,
                                title=f"{enzyme_name} Substrate Kinetics")
        file_paths["substrate_kinetics_plot"] = sk_plot
        log_lines.append(f"  - Substrate kinetics plot: {sk_plot}")
        
        # 5. Simulate modulator dose-responses
        if modulators:
            log_lines.append("")
            log_lines.append("Step 5: Simulating modulator dose-responses...")
            
            # Use middle substrate concentration as reference for dose-response
            S_ref = np.median(substrate_concentrations)
            
            for mod_name, mod_concs in modulators.items():
                # Resolve modulator type and parameters
                cfg = modulator_config.get(mod_name, {})
                mod_type = cfg.get("type", None)
                if mod_type is None:
                    # Fallback: infer from name
                    mod_type = "competitive_inhibition" if "inhibitor" in mod_name.lower() or "inhib" in mod_name.lower() else "activation_vmax"
                
                Ki_or_Ka = cfg.get("Ki", cfg.get("Ka", Km_base * 2.0 if "inhib" in mod_type else Km_base * 0.5))
                hill_n = cfg.get("hill_coefficient", 2.0)
                
                mod_concs_arr = np.array(mod_concs, dtype=float)
                activities = simulate_dose_response(
                    mod_concs, S_ref, enzyme_concentration_nM, Vmax_base, Km_base,
                    modulator_type=mod_type, Ki_or_Ka=Ki_or_Ka, **env_params
                )
                
                dr_csv = os.path.join(output_dir, f"dose_response_{mod_name.replace(' ', '_')}.csv")
                write_dose_response_csv(mod_concs, activities, mod_name, dr_csv)
                file_paths[f"dose_response_{mod_name}_csv"] = dr_csv
                log_lines.append(f"  - {mod_name} ({mod_type}) dose-response CSV: {dr_csv}")
                
                dr_plot = os.path.join(output_dir, f"dose_response_{mod_name.replace(' ', '_')}.png")
                plot_dose_response(mod_concs_arr, activities, mod_name, mod_type, dr_plot)
                file_paths[f"dose_response_{mod_name}_plot"] = dr_plot
                log_lines.append(f"  - {mod_name} ({mod_type}) dose-response plot: {dr_plot}")
        
        # 6. Generate report
        log_lines.append("")
        log_lines.append("Step 6: Generating report and metadata...")
        
        report_path = os.path.join(output_dir, "kinetics_report.txt")
        generate_report(enzyme_name, enzyme_concentration_nM, substrate_concentrations,
                        fit_params, modulators, report_path,
                        temperature=temperature, pH=pH,
                        product_inhibition=product_inhibition, KiP=KiP,
                        mwc_params=mwc_params)
        file_paths["report"] = report_path
        log_lines.append(f"  - Report: {report_path}")
        
        metadata_path = os.path.join(output_dir, "experiment_metadata.json")
        generate_json_metadata(enzyme_name, enzyme_concentration_nM, substrate_concentrations,
                               time_points, fit_params, modulators, modulator_config, file_paths,
                               metadata_path, temperature=temperature, pH=pH,
                               product_inhibition=product_inhibition, KiP=KiP,
                               mwc_params=mwc_params)
        file_paths["metadata"] = metadata_path
        log_lines.append(f"  - Metadata JSON: {metadata_path}")
        
        log_lines.append("")
        log_lines.append("Analysis complete!")
        
        return {
            "research_log": "\n".join(log_lines),
            "enzyme_name": enzyme_name,
            "fit_params": fit_params,
            "output_files": file_paths,
            "error": None,
        }
    
    except Exception as e:
        log_lines.append(f"ERROR: {str(e)}")
        return {
            "research_log": "\n".join(log_lines),
            "enzyme_name": input_data.get("enzyme_name", "Unknown"),
            "fit_params": {},
            "output_files": {},
            "error": str(e),
        }


def _png_file_to_base64(path: Optional[str]) -> Optional[str]:
    if not path or not os.path.isfile(path):
        return None
    with open(path, "rb") as f:
        return base64.b64encode(f.read()).decode("ascii")


def embed_kinetics_fit_plot_base64(result: Dict[str, Any]) -> Dict[str, Any]:
    """将底物动力学 MM 拟合图以 Base64 写入结果，供上层 Markdown 嵌入。"""
    if not isinstance(result, dict):
        return result
    ofs = result.get("output_files") or {}
    b64 = _png_file_to_base64(ofs.get("substrate_kinetics_plot"))
    if b64:
        result["substrate_kinetics_fit_plot_png_base64"] = b64
    return result


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(description="Enzyme Kinetics Experiment Analysis")
    parser.add_argument("--json-input", type=str, help="JSON string with input data")
    parser.add_argument("--output-dir", type=str, default="/tmp/enzyme_kinetics", help="Output directory")
    parser.add_argument("--output-json", type=str, help="Write full result as JSON to this path")
    
    args = parser.parse_args()
    
    if args.json_input:
        input_data = json.loads(args.json_input)
    else:
        print("Error: Provide --json-input", file=sys.stderr)
        sys.exit(1)
    
    result = run_analysis(input_data, args.output_dir)
    result = embed_kinetics_fit_plot_base64(result)
    
    # Print result as JSON
    print(json.dumps(result, indent=2, default=str))
    
    if args.output_json:
        with open(args.output_json, "w") as f:
            json.dump(result, f, indent=2, default=str)


if __name__ == "__main__":
    main()
