#!/usr/bin/env python3
"""
circular-dichroism: CD spectroscopy analysis for secondary structure and thermal stability.

Usage:
    python cd_analysis.py --data '{"sample_name":"MyProtein","sample_type":"protein","wavelength_data":[190,200,210,220],"cd_signal_data":[0.1,-0.2,0.3,-0.4]}' --output-dir /tmp
"""

import argparse
import json
import sys
import os
import csv
import math
from typing import Dict, List, Optional, Tuple, Any

# =============================================================================
# REFERENCE SPECTRA FOR PROTEIN SECONDARY STRUCTURE
# =============================================================================

# Reference CD spectra (normalized) for common secondary structure elements
# Values are approximate mean residue ellipticity [θ] in deg·cm²·dmol⁻¹
# Wavelengths in nm

PROTEIN_REFERENCES = {
    "alpha_helix": {
        190: 72000, 195: 55000, 200: 30000, 205: 10000,
        208: -30000, 210: -33000, 215: -36000, 220: -38000,
        222: -40000, 225: -35000, 230: -25000, 235: -15000, 240: -8000,
    },
    "beta_sheet": {
        190: 35000, 195: 45000, 200: 40000, 205: 25000,
        210: 10000, 215: -5000, 218: -20000, 220: -22000,
        225: -18000, 230: -12000, 235: -8000, 240: -5000,
    },
    "beta_turn": {
        190: 20000, 195: 25000, 200: 15000, 205: 0,
        210: -10000, 215: -18000, 220: -15000, 225: -10000,
        230: -5000, 235: -2000, 240: 0,
    },
    "random_coil": {
        190: -20000, 195: -30000, 200: -40000, 205: -35000,
        210: -25000, 215: -15000, 220: -8000, 225: -5000,
        230: -3000, 235: -2000, 240: -1000,
    },
    "poly_proline_ii": {
        190: 50000, 195: 40000, 200: 20000, 205: 0,
        210: -15000, 215: -25000, 220: -20000, 225: -15000,
        230: -10000, 235: -5000, 240: -2000,
    },
}

# =============================================================================
# NUCLEIC ACID REFERENCE SIGNATURES
# =============================================================================

NUCLEIC_SIGNATURES = {
    "B_DNA": {
        "description": "Right-handed double helix (Watson-Crick)",
        "positive_peak": (275, 280),
        "negative_peak": (245, 250),
        "crossover": (258, 262),
    },
    "A_DNA": {
        "description": "Right-handed double helix (dehydrated)",
        "positive_peak": (260, 265),
        "negative_peak": (210, 215),
        "crossover": (240, 245),
    },
    "Z_DNA": {
        "description": "Left-handed double helix",
        "positive_peak": (260, 265),
        "negative_peak": (290, 295),
        "crossover": (270, 275),
    },
    "G_quadruplex_parallel": {
        "description": "Parallel G-quadruplex",
        "positive_peak": (260, 265),
        "negative_peak": (240, 245),
        "crossover": (250, 255),
    },
    "G_quadruplex_antiparallel": {
        "description": "Antiparallel G-quadruplex",
        "positive_peak": (295, 300),
        "negative_peak": (265, 270),
        "crossover": (280, 285),
    },
    "G_quadruplex_hybrid": {
        "description": "Hybrid (3+1) G-quadruplex",
        "positive_peak": (290, 295),
        "negative_peak": (240, 245),
        "crossover": (265, 270),
    },
    "RNA_A_form": {
        "description": "A-form RNA duplex",
        "positive_peak": (265, 270),
        "negative_peak": (210, 215),
        "crossover": (240, 245),
    },
    "RNA_single_strand": {
        "description": "Single-stranded RNA",
        "positive_peak": (270, 275),
        "negative_peak": (245, 250),
        "crossover": (260, 265),
    },
    "triplex_DNA": {
        "description": "DNA triplex",
        "positive_peak": (280, 285),
        "negative_peak": (215, 220),
        "crossover": (250, 255),
    },
}

# =============================================================================
# PARSING
# =============================================================================

def parse_input(data_str: str) -> Dict:
    try:
        data = json.loads(data_str)
    except json.JSONDecodeError as e:
        raise ValueError(f"Invalid JSON: {e}")
    
    required = ["sample_name", "sample_type", "wavelength_data", "cd_signal_data"]
    for key in required:
        if key not in data:
            raise ValueError(f"Missing required key: '{key}'")
    
    sample_type = data["sample_type"].lower()
    if sample_type not in ("protein", "nucleic_acid"):
        raise ValueError(f"sample_type must be 'protein' or 'nucleic_acid', got '{sample_type}'")
    data["sample_type"] = sample_type
    
    wl = data["wavelength_data"]
    cd = data["cd_signal_data"]
    if not isinstance(wl, list) or not isinstance(cd, list):
        raise ValueError("wavelength_data and cd_signal_data must be lists")
    if len(wl) != len(cd):
        raise ValueError(f"Wavelength data ({len(wl)} points) and CD signal data ({len(cd)} points) must have same length")
    if len(wl) < 3:
        raise ValueError("At least 3 data points required for meaningful analysis")
    
    try:
        data["wavelength_data"] = [float(x) for x in wl]
        data["cd_signal_data"] = [float(x) for x in cd]
    except (ValueError, TypeError):
        raise ValueError("All wavelength and CD signal values must be numeric")
    
    if "temperature_data" in data and data["temperature_data"] is not None:
        if "thermal_cd_data" not in data or data["thermal_cd_data"] is None:
            raise ValueError("thermal_cd_data is required when temperature_data is provided")
        temp = data["temperature_data"]
        t_cd = data["thermal_cd_data"]
        if len(temp) != len(t_cd):
            raise ValueError(f"Temperature data ({len(temp)}) and thermal CD data ({len(t_cd)}) must match")
        if len(temp) < 4:
            raise ValueError("At least 4 thermal data points required for Tm fitting")
        data["temperature_data"] = [float(x) for x in temp]
        data["thermal_cd_data"] = [float(x) for x in t_cd]
    
    return data

# =============================================================================
# UTILITY FUNCTIONS
# =============================================================================

def interpolate_ref(ref_dict: Dict[float, float], wavelength: float) -> float:
    wavelengths = sorted(ref_dict.keys())
    if wavelength <= wavelengths[0]:
        return ref_dict[wavelengths[0]]
    if wavelength >= wavelengths[-1]:
        return ref_dict[wavelengths[-1]]
    for i in range(len(wavelengths) - 1):
        if wavelengths[i] <= wavelength <= wavelengths[i+1]:
            w1, w2 = wavelengths[i], wavelengths[i+1]
            v1, v2 = ref_dict[w1], ref_dict[w2]
            return v1 + (v2 - v1) * (wavelength - w1) / (w2 - w1)
    return ref_dict[wavelengths[-1]]

def normalize_spectrum(wavelengths: List[float], signals: List[float]) -> List[float]:
    max_abs = max(abs(x) for x in signals) if any(signals) else 1.0
    if max_abs == 0:
        return [0.0] * len(signals)
    return [x / max_abs for x in signals]

def compute_similarity(spectrum1: List[float], spectrum2: List[float]) -> float:
    if len(spectrum1) != len(spectrum2):
        return 0.0
    dot = sum(a * b for a, b in zip(spectrum1, spectrum2))
    norm1 = math.sqrt(sum(a * a for a in spectrum1))
    norm2 = math.sqrt(sum(b * b for b in spectrum2))
    if norm1 == 0 or norm2 == 0:
        return 0.0
    return dot / (norm1 * norm2)

def compute_rmse(observed: List[float], predicted: List[float]) -> float:
    n = len(observed)
    if n == 0:
        return float('inf')
    return math.sqrt(sum((o - p) ** 2 for o, p in zip(observed, predicted)) / n)

# =============================================================================
# PROTEIN SECONDARY STRUCTURE ANALYSIS
# =============================================================================

def analyze_protein_cd(wavelengths: List[float], signals: List[float]) -> Dict:
    norm_signals = normalize_spectrum(wavelengths, signals)
    ref_spectra = {}
    for struct_name, ref_dict in PROTEIN_REFERENCES.items():
        ref_at_wl = [interpolate_ref(ref_dict, wl) for wl in wavelengths]
        ref_spectra[struct_name] = normalize_spectrum(wavelengths, ref_at_wl)
    
    similarities = {}
    for struct_name, ref_norm in ref_spectra.items():
        sim = compute_similarity(norm_signals, ref_norm)
        similarities[struct_name] = max(0.0, sim)
    
    total_sim = sum(similarities.values())
    if total_sim > 0:
        fractions = {k: v / total_sim for k, v in similarities.items()}
    else:
        fractions = {k: 0.2 for k in similarities}
    
    total_frac = sum(fractions.values())
    if total_frac > 0:
        fractions = {k: v / total_frac * 100 for k, v in fractions.items()}
    
    dominant = max(fractions, key=fractions.get)
    dominant_pct = fractions[dominant]
    
    features = analyze_spectral_features(wavelengths, signals)
    refined = refine_protein_estimates(fractions, features, wavelengths, signals)
    
    return {
        "fractions": refined,
        "dominant_structure": dominant.replace("_", " ").title(),
        "dominant_percentage": round(refined.get(dominant, 0), 1),
        "spectral_features": features,
        "confidence": "medium",
    }

def analyze_spectral_features(wavelengths: List[float], signals: List[float]) -> Dict:
    features = {}
    max_idx = max(range(len(signals)), key=lambda i: signals[i])
    min_idx = min(range(len(signals)), key=lambda i: signals[i])
    features["max_wavelength"] = wavelengths[max_idx]
    features["max_signal"] = signals[max_idx]
    features["min_wavelength"] = wavelengths[min_idx]
    features["min_signal"] = signals[min_idx]
    features["amplitude"] = signals[max_idx] - signals[min_idx]
    
    has_208_222_doublet = False
    has_218_minimum = False
    has_200_minimum = False
    
    for wl, sig in zip(wavelengths, signals):
        if 206 <= wl <= 210:
            features["signal_at_208"] = sig
        if 220 <= wl <= 224:
            features["signal_at_222"] = sig
        if 216 <= wl <= 220:
            features["signal_at_218"] = sig
        if 198 <= wl <= 202:
            features["signal_at_200"] = sig
    
    if "signal_at_208" in features and "signal_at_222" in features:
        if features["signal_at_208"] < -0.05 and features["signal_at_222"] < -0.05:
            has_208_222_doublet = True
    if "signal_at_218" in features and features["signal_at_218"] < -0.05:
        has_218_minimum = True
    if "signal_at_200" in features and features["signal_at_200"] < -0.05:
        has_200_minimum = True
    
    features["has_208_222_doublet"] = has_208_222_doublet
    features["has_218_minimum"] = has_218_minimum
    features["has_200_minimum"] = has_200_minimum
    
    return features

def refine_protein_estimates(fractions: Dict[str, float], features: Dict, 
                             wavelengths: List[float], signals: List[float]) -> Dict[str, float]:
    refined = dict(fractions)
    
    if features.get("has_208_222_doublet"):
        refined["alpha_helix"] = min(95, refined.get("alpha_helix", 0) * 1.3 + 10)
        others = [k for k in refined if k != "alpha_helix"]
        total_others = sum(refined.get(k, 0) for k in others)
        if total_others > 0:
            scale = max(0, 100 - refined["alpha_helix"]) / total_others
            for k in others:
                refined[k] *= scale
    
    if features.get("has_218_minimum") and not features.get("has_208_222_doublet"):
        refined["beta_sheet"] = min(80, refined.get("beta_sheet", 0) * 1.2 + 5)
        others = [k for k in refined if k != "beta_sheet"]
        total_others = sum(refined.get(k, 0) for k in others)
        if total_others > 0:
            scale = max(0, 100 - refined["beta_sheet"]) / total_others
            for k in others:
                refined[k] *= scale
    
    if features.get("has_200_minimum") and not features.get("has_208_222_doublet"):
        refined["random_coil"] = min(70, refined.get("random_coil", 0) * 1.2 + 5)
        others = [k for k in refined if k != "random_coil"]
        total_others = sum(refined.get(k, 0) for k in others)
        if total_others > 0:
            scale = max(0, 100 - refined["random_coil"]) / total_others
            for k in others:
                refined[k] *= scale
    
    total = sum(refined.values())
    if total > 0:
        refined = {k: round(v / total * 100, 1) for k, v in refined.items()}
    
    return refined

# =============================================================================
# NUCLEIC ACID STRUCTURE ANALYSIS
# =============================================================================

def analyze_nucleic_cd(wavelengths: List[float], signals: List[float]) -> Dict:
    features = analyze_spectral_features(wavelengths, signals)
    scores = {}
    
    for struct_name, signature in NUCLEIC_SIGNATURES.items():
        score = 0.0
        pos_wl_range = signature["positive_peak"]
        for wl, sig in zip(wavelengths, signals):
            if pos_wl_range[0] - 5 <= wl <= pos_wl_range[1] + 5:
                if sig > 0:
                    score += sig * 2
        neg_wl_range = signature["negative_peak"]
        for wl, sig in zip(wavelengths, signals):
            if neg_wl_range[0] - 5 <= wl <= neg_wl_range[1] + 5:
                if sig < 0:
                    score += abs(sig) * 2
        scores[struct_name] = score
    
    if scores:
        best_match = max(scores, key=scores.get)
        best_score = scores[best_match]
    else:
        best_match = "unknown"
        best_score = 0
    
    has_290_negative = any(288 <= wl <= 292 and sig < -0.05 for wl, sig in zip(wavelengths, signals))
    has_260_positive = any(258 <= wl <= 262 and sig > 0.05 for wl, sig in zip(wavelengths, signals))
    has_240_negative = any(238 <= wl <= 242 and sig < -0.05 for wl, sig in zip(wavelengths, signals))
    has_295_positive = any(293 <= wl <= 297 and sig > 0.05 for wl, sig in zip(wavelengths, signals))
    
    structure_type = best_match.replace("_", " ").title()
    confidence = "low"
    
    if has_295_positive and has_260_positive:
        structure_type = "Hybrid G-quadruplex (3+1)"
        confidence = "high"
    elif has_295_positive:
        structure_type = "Antiparallel G-quadruplex"
        confidence = "high"
    elif has_260_positive and has_240_negative and not has_290_negative:
        structure_type = "Parallel G-quadruplex"
        confidence = "high"
    elif has_290_negative and has_260_positive:
        structure_type = "Z-DNA (Left-handed)"
        confidence = "high"
    elif has_260_positive and has_240_negative:
        if any(275 <= wl <= 280 and sig > 0 for wl, sig in zip(wavelengths, signals)):
            structure_type = "B-DNA (Right-handed)"
            confidence = "high"
        else:
            structure_type = "A-DNA or A-RNA"
            confidence = "medium"
    elif any(265 <= wl <= 270 and sig > 0 for wl, sig in zip(wavelengths, signals)):
        if any(210 <= wl <= 215 and sig < 0 for wl, sig in zip(wavelengths, signals)):
            structure_type = "A-DNA / RNA A-form"
            confidence = "high"
        else:
            structure_type = "Single-stranded RNA"
            confidence = "medium"
    
    content = estimate_nucleic_content(structure_type, signals)
    
    return {
        "structure_type": structure_type,
        "confidence": confidence,
        "spectral_features": features,
        "scores": {k.replace("_", " ").title(): round(v, 2) for k, v in scores.items()},
        "structure_content": content,
    }

def estimate_nucleic_content(structure_type: str, signals: List[float]) -> Dict:
    if "B-DNA" in structure_type:
        return {"double_stranded": 85, "single_stranded": 10, "other": 5}
    elif "A-DNA" in structure_type or "A-form" in structure_type:
        return {"double_stranded": 80, "single_stranded": 15, "other": 5}
    elif "Z-DNA" in structure_type:
        return {"z_form": 70, "b_form": 25, "other": 5}
    elif "G-quadruplex" in structure_type:
        return {"quadruplex": 75, "unfolded": 20, "other": 5}
    elif "Single-stranded" in structure_type:
        return {"single_stranded": 90, "double_stranded": 5, "other": 5}
    else:
        return {"unknown": 100}

# =============================================================================
# THERMAL STABILITY ANALYSIS
# =============================================================================

def analyze_thermal_stability(temperatures: List[float], thermal_signals: List[float]) -> Dict:
    max_sig = max(thermal_signals)
    min_sig = min(thermal_signals)
    amplitude = max_sig - min_sig
    
    if amplitude < 0.01:
        return {
            "tm": None,
            "tm_error": None,
            "transition_amplitude": amplitude,
            "baseline_start": thermal_signals[0],
            "baseline_end": thermal_signals[-1],
            "analysis_method": "insufficient_signal_change",
            "folded_fraction": [],
            "reliability": "unreliable",
            "warning": "Thermal signal change too small (<0.01) to determine Tm reliably",
        }
    
    n = len(temperatures)
    if n >= 2:
        slope = (thermal_signals[-1] - thermal_signals[0]) / (temperatures[-1] - temperatures[0])
    else:
        slope = 0
    
    if slope > 0:
        signal_folded = min_sig
        signal_unfolded = max_sig
    else:
        signal_folded = max_sig
        signal_unfolded = min_sig
    
    folded_fractions = []
    for sig in thermal_signals:
        if signal_folded != signal_unfolded:
            f = (sig - signal_unfolded) / (signal_folded - signal_unfolded)
        else:
            f = 0.5
        f = max(0.0, min(1.0, f))
        folded_fractions.append(f)
    
    tm = None
    tm_idx = None
    for i in range(len(folded_fractions) - 1):
        if (folded_fractions[i] >= 0.5 and folded_fractions[i+1] < 0.5) or \
           (folded_fractions[i] < 0.5 and folded_fractions[i+1] >= 0.5):
            f1, f2 = folded_fractions[i], folded_fractions[i+1]
            t1, t2 = temperatures[i], temperatures[i+1]
            tm = t1 + (0.5 - f1) * (t2 - t1) / (f2 - f1)
            tm_idx = i
            break
    
    if tm is None:
        closest_idx = min(range(len(folded_fractions)), key=lambda i: abs(folded_fractions[i] - 0.5))
        tm = temperatures[closest_idx]
        tm_idx = closest_idx
    
    if tm_idx is not None and 1 <= tm_idx < len(folded_fractions) - 1:
        dt = temperatures[tm_idx + 1] - temperatures[tm_idx - 1]
        df = folded_fractions[tm_idx + 1] - folded_fractions[tm_idx - 1]
        if dt > 0 and abs(df) > 0.001:
            slope_at_tm = abs(df / dt)
            R = 8.314
            Tm_k = tm + 273.15
            delta_h = 4 * R * Tm_k * Tm_k * slope_at_tm
        else:
            delta_h = None
    else:
        delta_h = None
    
    t_90 = None
    t_10 = None
    for i in range(len(folded_fractions) - 1):
        if (folded_fractions[i] >= 0.9 and folded_fractions[i+1] < 0.9) or \
           (folded_fractions[i] < 0.9 and folded_fractions[i+1] >= 0.9):
            t_90 = temperatures[i]
        if (folded_fractions[i] >= 0.1 and folded_fractions[i+1] < 0.1) or \
           (folded_fractions[i] < 0.1 and folded_fractions[i+1] >= 0.1):
            t_10 = temperatures[i]
    
    if t_90 is not None and t_10 is not None:
        transition_width = abs(t_10 - t_90)
    else:
        transition_width = None
    
    if amplitude < 0.05:
        reliability = "low"
    elif transition_width is not None and transition_width < 15:
        reliability = "high"
    elif transition_width is not None and transition_width < 30:
        reliability = "medium"
    else:
        reliability = "low"
    
    return {
        "tm": round(tm, 1) if tm is not None else None,
        "tm_error": 2.0 if reliability == "high" else (5.0 if reliability == "medium" else 10.0),
        "transition_amplitude": round(amplitude, 4),
        "baseline_start": round(thermal_signals[0], 4),
        "baseline_end": round(thermal_signals[-1], 4),
        "analysis_method": "two_state_linear_baseline",
        "folded_fraction": [round(f, 3) for f in folded_fractions],
        "reliability": reliability,
        "delta_h_vant_hoff": round(delta_h / 1000, 1) if delta_h else None,
        "transition_width": round(transition_width, 1) if transition_width else None,
        "warning": None if reliability != "low" else "Low thermal signal amplitude or broad transition. Tm estimate may be unreliable.",
    }

# =============================================================================
# CSV GENERATION
# =============================================================================

def generate_cd_csv(wavelengths: List[float], signals: List[float], 
                      sample_name: str, sample_type: str, output_path: str) -> str:
    with open(output_path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["Sample", sample_name])
        writer.writerow(["Type", sample_type])
        writer.writerow([])
        writer.writerow(["Wavelength_nm", "CD_Signal"])
        for wl, sig in zip(wavelengths, signals):
            writer.writerow([wl, sig])
    return output_path

def generate_thermal_csv(temperatures: List[float], thermal_signals: List[float],
                           folded_fractions: List[float], sample_name: str, output_path: str) -> str:
    with open(output_path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["Sample", sample_name])
        writer.writerow(["Analysis", "Thermal_Denaturation"])
        writer.writerow([])
        writer.writerow(["Temperature_C", "CD_Signal", "Folded_Fraction"])
        for t, sig, ff in zip(temperatures, thermal_signals, folded_fractions):
            writer.writerow([t, sig, ff])
    return output_path

# =============================================================================
# REPORT GENERATION
# =============================================================================

def generate_report(data: Dict, cd_analysis: Dict, thermal_analysis: Optional[Dict],
                    cd_csv_path: str, thermal_csv_path: Optional[str]) -> str:
    lines = []
    lines.append("=" * 75)
    lines.append("CIRCULAR DICHROISM (CD) SPECTROSCOPY ANALYSIS REPORT")
    lines.append("=" * 75)
    lines.append("")
    
    lines.append("-" * 75)
    lines.append("SECTION 1: SAMPLE INFORMATION")
    lines.append("-" * 75)
    lines.append("")
    lines.append(f"  Sample Name:        {data['sample_name']}")
    lines.append(f"  Sample Type:        {data['sample_type'].replace('_', ' ').title()}")
    lines.append(f"  Data Points:        {len(data['wavelength_data'])} (wavelength range: {min(data['wavelength_data']):.1f} - {max(data['wavelength_data']):.1f} nm)")
    if data.get('temperature_data'):
        lines.append(f"  Thermal Data:       Yes ({len(data['temperature_data'])} points, {min(data['temperature_data']):.1f} - {max(data['temperature_data']):.1f} °C)")
    else:
        lines.append(f"  Thermal Data:       No")
    lines.append("")
    
    lines.append("-" * 75)
    lines.append("SECTION 2: CD SPECTRUM ANALYSIS")
    lines.append("-" * 75)
    lines.append("")
    
    lines.append("Raw CD Data:")
    lines.append(f"  {'Wavelength (nm)':<18} {'CD Signal':>12}")
    lines.append(f"  {'-'*18} {'-'*12}")
    for wl, sig in zip(data['wavelength_data'], data['cd_signal_data']):
        lines.append(f"  {wl:<18.1f} {sig:>12.4f}")
    lines.append("")
    
    features = cd_analysis.get("spectral_features", {})
    lines.append("Spectral Features:")
    lines.append(f"  Maximum signal:     {features.get('max_signal', 0):.4f} at {features.get('max_wavelength', 0):.1f} nm")
    lines.append(f"  Minimum signal:     {features.get('min_signal', 0):.4f} at {features.get('min_wavelength', 0):.1f} nm")
    lines.append(f"  Spectral amplitude: {features.get('amplitude', 0):.4f}")
    lines.append("")
    
    if data['sample_type'] == 'protein':
        lines.append("Secondary Structure Estimation:")
        lines.append("  (Based on reference spectrum comparison and wavelength markers)")
        lines.append("")
        lines.append(f"  Dominant structure: {cd_analysis['dominant_structure']} ({cd_analysis['dominant_percentage']}%)")
        lines.append("")
        
        fractions = cd_analysis.get("fractions", {})
        lines.append("  Estimated composition:")
        for struct, pct in sorted(fractions.items(), key=lambda x: -x[1]):
            struct_name = struct.replace("_", " ").title()
            bar = "█" * int(pct / 5)
            lines.append(f"    {struct_name:<20} {pct:>6.1f}%  {bar}")
        lines.append("")
        
        lines.append("Interpretation:")
        dom = cd_analysis['dominant_structure'].lower()
        if 'alpha' in dom:
            lines.append("  • Strong α-helical signature detected. Typical of coiled-coil proteins,")
            lines.append("    membrane proteins, or proteins with extensive helix bundles.")
        elif 'beta' in dom and 'turn' not in dom:
            lines.append("  • Predominant β-sheet structure. Common in immunoglobulins,")
            lines.append("    amyloid fibrils, and many globular proteins.")
        elif 'random' in dom or 'coil' in dom:
            lines.append("  • High random coil content suggests natively unfolded or intrinsically")
            lines.append("    disordered protein regions. May indicate partial denaturation.")
        elif 'poly_proline' in dom:
            lines.append("  • PPII helix signature detected. Characteristic of extended conformations")
            lines.append("    and often seen in collagen-like or proline-rich regions.")
        else:
            lines.append("  • Mixed secondary structure content.")
        
        if features.get("has_208_222_doublet"):
            lines.append("  • Double minimum at ~208 and ~222 nm confirms α-helical content.")
        if features.get("has_218_minimum"):
            lines.append("  • Minimum near 218 nm suggests β-sheet contribution.")
        if features.get("has_200_minimum"):
            lines.append("  • Strong negative signal near 200 nm indicates disordered/random coil regions.")
        lines.append("")
        
    else:
        lines.append("Nucleic Acid Structure Analysis:")
        lines.append("")
        lines.append(f"  Predicted structure:  {cd_analysis['structure_type']}")
        lines.append(f"  Confidence:           {cd_analysis['confidence'].title()}")
        lines.append("")
        
        content = cd_analysis.get("structure_content", {})
        if content:
            lines.append("  Estimated content:")
            for key, pct in content.items():
                struct_name = key.replace("_", " ").title()
                lines.append(f"    {struct_name:<20} {pct}%")
            lines.append("")
        
        scores = cd_analysis.get("scores", {})
        if scores:
            lines.append("  Structure match scores:")
            for struct, score in sorted(scores.items(), key=lambda x: -x[1]):
                lines.append(f"    {struct:<25} {score:.3f}")
            lines.append("")
        
        lines.append("Interpretation:")
        st = cd_analysis['structure_type']
        if 'B-DNA' in st:
            lines.append("  • B-DNA signature: Positive band ~275 nm, negative ~245 nm.")
            lines.append("    Standard right-handed Watson-Crick duplex conformation.")
        elif 'A-DNA' in st or 'A-form' in st:
            lines.append("  • A-form signature: Positive band ~260 nm, negative ~210 nm.")
            lines.append("    Common in dehydrated conditions or RNA duplexes.")
        elif 'Z-DNA' in st:
            lines.append("  • Z-DNA signature: Negative band ~290 nm, positive ~260 nm.")
            lines.append("    Left-handed conformation, often stabilized by high salt or negative supercoiling.")
        elif 'G-quadruplex' in st:
            lines.append("  • G-quadruplex signature detected. Characteristic of guanine-rich sequences.")
            lines.append("    Structure depends on strand orientation and loop connectivity.")
        elif 'Single-stranded' in st:
            lines.append("  • Single-stranded signature: Reduced ellipticity compared to duplex forms.")
        lines.append("")
    
    if thermal_analysis:
        lines.append("-" * 75)
        lines.append("SECTION 3: THERMAL STABILITY ANALYSIS")
        lines.append("-" * 75)
        lines.append("")
        
        tm = thermal_analysis.get("tm")
        if tm is not None:
            lines.append(f"  Melting Temperature (Tm):     {tm:.1f} ± {thermal_analysis.get('tm_error', 'N/A')} °C")
            lines.append(f"  Transition Amplitude:         {thermal_analysis['transition_amplitude']:.4f}")
            lines.append(f"  Baseline (start):             {thermal_analysis['baseline_start']:.4f}")
            lines.append(f"  Baseline (end):               {thermal_analysis['baseline_end']:.4f}")
            lines.append(f"  Analysis Method:              {thermal_analysis['analysis_method']}")
            
            dh = thermal_analysis.get("delta_h_vant_hoff")
            if dh:
                lines.append(f"  ΔH (van't Hoff):              {dh:.1f} kJ/mol")
            
            tw = thermal_analysis.get("transition_width")
            if tw:
                lines.append(f"  Transition Width (10%-90%):   {tw:.1f} °C")
            
            lines.append(f"  Reliability:                  {thermal_analysis['reliability'].title()}")
            lines.append("")
            
            lines.append("Thermal Denaturation Data:")
            lines.append(f"  {'Temperature (°C)':<18} {'CD Signal':>12} {'Folded Fraction':>16}")
            lines.append(f"  {'-'*18} {'-'*12} {'-'*16}")
            for t, sig, ff in zip(data['temperature_data'], data['thermal_cd_data'], 
                                  thermal_analysis['folded_fraction']):
                lines.append(f"  {t:<18.1f} {sig:>12.4f} {ff:>16.3f}")
            lines.append("")
            
            lines.append("Thermal Stability Interpretation:")
            if tm < 40:
                lines.append(f"  • Tm = {tm:.1f}°C indicates low thermal stability. The sample may be")
                lines.append("    partially unfolded at physiological temperature or intrinsically disordered.")
            elif tm < 60:
                lines.append(f"  • Tm = {tm:.1f}°C indicates moderate stability. Typical of many globular")
                lines.append("    proteins or nucleic acid duplexes under standard conditions.")
            elif tm < 80:
                lines.append(f"  • Tm = {tm:.1f}°C indicates good thermal stability. The structure is")
                lines.append("    well-maintained at physiological and elevated temperatures.")
            else:
                lines.append(f"  • Tm = {tm:.1f}°C indicates very high thermal stability. Typical of")
                lines.append("    hyperstable proteins, highly structured RNAs, or proteins with extensive")
                lines.append("    disulfide crosslinking.")
            
            if thermal_analysis['reliability'] == 'low':
                lines.append("  ⚠️ WARNING: Tm estimate has low reliability. Consider collecting more")
                lines.append("     temperature points near the transition region for improved accuracy.")
            
            if tw and tw < 5:
                lines.append("  • Sharp transition suggests highly cooperative two-state unfolding.")
            elif tw and tw > 20:
                lines.append("  • Broad transition suggests gradual unfolding or multiple unfolding domains.")
            
            lines.append("")
        else:
            lines.append("  Could not determine Tm. Insufficient signal change during thermal scan.")
            lines.append("")
    
    lines.append("-" * 75)
    lines.append("SECTION 4: RECOMMENDATIONS")
    lines.append("-" * 75)
    lines.append("")
    
    if data['sample_type'] == 'protein':
        lines.append("  • For improved secondary structure accuracy, collect data in the far-UV")
        lines.append("    region (190-260 nm) with at least 1 nm resolution.")
        if features.get("amplitude", 0) < 0.1:
            lines.append("  ⚠️ Low spectral amplitude detected. Verify protein concentration")
            lines.append("     (recommended: 0.1-1.0 mg/mL) and path length (0.1-1.0 mm).")
        lines.append("  • Consider using CONTIN, SELCON3, or CDSSTR algorithms for more")
        lines.append("    rigorous secondary structure deconvolution with full reference datasets.")
        if not data.get('temperature_data'):
            lines.append("  • Add thermal denaturation data to assess conformational stability.")
    else:
        lines.append("  • For nucleic acids, verify buffer conditions (salt, pH) as these strongly")
        lines.append("    influence conformational equilibria.")
        lines.append("  • Consider collecting near-UV CD (250-320 nm) and IR spectroscopy data")
        lines.append("    for complementary structural information.")
        if 'G-quadruplex' in cd_analysis.get('structure_type', ''):
            lines.append("  • G-quadruplex topology can be confirmed by NMR or X-ray crystallography.")
            lines.append("    CD alone cannot distinguish all topological variants unambiguously.")
    
    lines.append("")
    
    lines.append("-" * 75)
    lines.append("SECTION 5: OUTPUT FILES")
    lines.append("-" * 75)
    lines.append("")
    lines.append(f"  CD Spectrum CSV:      {cd_csv_path}")
    if thermal_csv_path:
        lines.append(f"  Thermal Data CSV:     {thermal_csv_path}")
    lines.append("")
    
    lines.append("=" * 75)
    lines.append("END OF REPORT")
    lines.append("=" * 75)
    
    return "\n".join(lines)

# =============================================================================
# MAIN
# =============================================================================

def main():
    parser = argparse.ArgumentParser(description="Circular Dichroism (CD) spectrum analysis")
    parser.add_argument("--data", required=True, help="JSON string with CD data parameters")
    parser.add_argument("--output-dir", default="/tmp", help="Output directory for CSV files")
    parser.add_argument("--output-format", default="text", choices=["text", "json"], help="Output format")
    parser.add_argument("--cd-csv", default="cd_spectrum.csv", help="CD spectrum CSV filename")
    parser.add_argument("--thermal-csv", default="thermal_denaturation.csv", help="Thermal data CSV filename")
    
    args = parser.parse_args()
    
    try:
        data = parse_input(args.data)
        
        if data['sample_type'] == 'protein':
            cd_analysis = analyze_protein_cd(data['wavelength_data'], data['cd_signal_data'])
        else:
            cd_analysis = analyze_nucleic_cd(data['wavelength_data'], data['cd_signal_data'])
        
        thermal_analysis = None
        if data.get('temperature_data') and data.get('thermal_cd_data'):
            thermal_analysis = analyze_thermal_stability(
                data['temperature_data'], data['thermal_cd_data']
            )
        
        os.makedirs(args.output_dir, exist_ok=True)
        cd_csv_path = os.path.join(args.output_dir, args.cd_csv)
        generate_cd_csv(data['wavelength_data'], data['cd_signal_data'],
                        data['sample_name'], data['sample_type'], cd_csv_path)
        
        thermal_csv_path = None
        if thermal_analysis:
            thermal_csv_path = os.path.join(args.output_dir, args.thermal_csv)
            generate_thermal_csv(
                data['temperature_data'], data['thermal_cd_data'],
                thermal_analysis['folded_fraction'], data['sample_name'], thermal_csv_path
            )
        
        report = generate_report(data, cd_analysis, thermal_analysis, cd_csv_path, thermal_csv_path)
        
        if args.output_format == 'json':
            output = {
                "report_log": report,
                "cd_spectrum_analysis_file_url": f"file://{cd_csv_path}",
            }
            if thermal_csv_path:
                output["thermal_analysis_file_url"] = f"file://{thermal_csv_path}"
            
            warnings = []
            if data['sample_type'] == 'protein':
                features = cd_analysis.get("spectral_features", {})
                if features.get("amplitude", 0) < 0.1:
                    warnings.append("Low CD spectral amplitude. Verify concentration and path length.")
            if thermal_analysis and thermal_analysis.get("warning"):
                warnings.append(thermal_analysis["warning"])
            
            if warnings:
                output["warning"] = "; ".join(warnings)
            
            print(json.dumps(output, indent=2, ensure_ascii=False))
        else:
            print(report)
            print(f"\nCD Spectrum CSV: {cd_csv_path}")
            if thermal_csv_path:
                print(f"Thermal Data CSV: {thermal_csv_path}")
    
    except Exception as e:
        print(f"ERROR: {e}", file=sys.stderr)
        import traceback
        traceback.print_exc(file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    main()
