#!/usr/bin/env python3
"""
immune-cell: Simulate immune cell isolation and purification from tissue samples.

Usage:
    python immune_cell_isolation.py --data '{"tissue_type":"spleen","target_cell_type":"macrophages","enzyme_type":"collagenase","digestion_time_min":60}'
"""

import argparse
import json
import sys
import os
import csv
import random
from datetime import datetime
from typing import Dict, List, Optional, Tuple

# =============================================================================
# KNOWLEDGE BASE
# =============================================================================

TISSUE_PROFILES = {
    "spleen": {
        "total_cells_per_mg": 5_000_000,
        "recommended_enzymes": ["collagenase", "liberase"],
        "typical_yield_mg": 100,
        "digestion_efficiency_base": 0.85,
        "cell_fragility": 0.15,  # higher = more fragile
    },
    "adipose": {
        "total_cells_per_mg": 500_000,
        "recommended_enzymes": ["collagenase", "collagenase_pronase"],
        "typical_yield_mg": 500,
        "digestion_efficiency_base": 0.70,
        "cell_fragility": 0.25,
    },
    "kidney": {
        "total_cells_per_mg": 3_000_000,
        "recommended_enzymes": ["collagenase", "dispase", "trypsin"],
        "typical_yield_mg": 200,
        "digestion_efficiency_base": 0.75,
        "cell_fragility": 0.20,
    },
    "liver": {
        "total_cells_per_mg": 4_000_000,
        "recommended_enzymes": ["collagenase", "liberase"],
        "typical_yield_mg": 300,
        "digestion_efficiency_base": 0.80,
        "cell_fragility": 0.18,
    },
    "lung": {
        "total_cells_per_mg": 2_000_000,
        "recommended_enzymes": ["collagenase", "dispase", "elastase"],
        "typical_yield_mg": 150,
        "digestion_efficiency_base": 0.72,
        "cell_fragility": 0.22,
    },
    "tumor": {
        "total_cells_per_mg": 1_500_000,
        "recommended_enzymes": ["collagenase", "liberase", "accutase"],
        "typical_yield_mg": 250,
        "digestion_efficiency_base": 0.65,
        "cell_fragility": 0.30,
    },
    "blood": {
        "total_cells_per_mg": 0,  # special case: per mL
        "recommended_enzymes": ["none", "ack_lysis"],
        "typical_yield_mg": 0,
        "typical_volume_ml": 10,
        "total_cells_per_ml": 5_000_000,
        "digestion_efficiency_base": 0.95,
        "cell_fragility": 0.05,
    },
    "lymph_node": {
        "total_cells_per_mg": 8_000_000,
        "recommended_enzymes": ["collagenase", "liberase", "dnase"],
        "typical_yield_mg": 50,
        "digestion_efficiency_base": 0.88,
        "cell_fragility": 0.12,
    },
    "thymus": {
        "total_cells_per_mg": 10_000_000,
        "recommended_enzymes": ["collagenase", "liberase"],
        "typical_yield_mg": 30,
        "digestion_efficiency_base": 0.90,
        "cell_fragility": 0.10,
    },
    "bone_marrow": {
        "total_cells_per_mg": 0,
        "recommended_enzymes": ["none", "collagenase"],
        "typical_yield_mg": 0,
        "typical_volume_ml": 1,
        "total_cells_per_ml": 15_000_000,
        "digestion_efficiency_base": 0.92,
        "cell_fragility": 0.08,
    },
}

ENZYME_PROFILES = {
    "collagenase": {
        "digestion_rate": 1.0,
        "cytotoxicity": 0.10,
        "optimal_time_min": 45,
        "max_time_min": 90,
    },
    "liberase": {
        "digestion_rate": 1.2,
        "cytotoxicity": 0.08,
        "optimal_time_min": 30,
        "max_time_min": 60,
    },
    "dispase": {
        "digestion_rate": 0.9,
        "cytotoxicity": 0.12,
        "optimal_time_min": 50,
        "max_time_min": 100,
    },
    "trypsin": {
        "digestion_rate": 1.5,
        "cytotoxicity": 0.20,
        "optimal_time_min": 10,
        "max_time_min": 20,
    },
    "pronase": {
        "digestion_rate": 1.3,
        "cytotoxicity": 0.25,
        "optimal_time_min": 20,
        "max_time_min": 40,
    },
    "collagenase_pronase": {
        "digestion_rate": 1.4,
        "cytotoxicity": 0.22,
        "optimal_time_min": 25,
        "max_time_min": 50,
    },
    "elastase": {
        "digestion_rate": 1.1,
        "cytotoxicity": 0.15,
        "optimal_time_min": 35,
        "max_time_min": 70,
    },
    "accutase": {
        "digestion_rate": 0.8,
        "cytotoxicity": 0.05,
        "optimal_time_min": 10,
        "max_time_min": 20,
    },
    "dnase": {
        "digestion_rate": 0.5,
        "cytotoxicity": 0.02,
        "optimal_time_min": 15,
        "max_time_min": 30,
    },
    "ack_lysis": {
        "digestion_rate": 0.6,
        "cytotoxicity": 0.15,
        "optimal_time_min": 5,
        "max_time_min": 10,
    },
    "none": {
        "digestion_rate": 0.3,
        "cytotoxicity": 0.0,
        "optimal_time_min": 0,
        "max_time_min": 10,
    },
}

CELL_TYPE_PROFILES = {
    "macrophages": {
        "frequency_in_tissue": 0.05,
        "macs_antibody": "CD11b",
        "macs_purity_base": 0.90,
        "macs_recovery_base": 0.75,
        "alternative_antibodies": ["F4/80", "CD68", "CD115"],
        "density_g_ml": 1.077,
    },
    "t_cells": {
        "frequency_in_tissue": 0.25,
        "macs_antibody": "CD3",
        "macs_purity_base": 0.92,
        "macs_recovery_base": 0.80,
        "alternative_antibodies": ["CD4", "CD8", "TCR"],
        "density_g_ml": 1.075,
    },
    "b_cells": {
        "frequency_in_tissue": 0.15,
        "macs_antibody": "CD19",
        "macs_purity_base": 0.91,
        "macs_recovery_base": 0.78,
        "alternative_antibodies": ["CD20", "B220", "CD79a"],
        "density_g_ml": 1.075,
    },
    "nk_cells": {
        "frequency_in_tissue": 0.05,
        "macs_antibody": "CD56",
        "macs_purity_base": 0.88,
        "macs_recovery_base": 0.70,
        "alternative_antibodies": ["NKp46", "NKG2D", "CD16"],
        "density_g_ml": 1.076,
    },
    "dendritic_cells": {
        "frequency_in_tissue": 0.02,
        "macs_antibody": "CD11c",
        "macs_purity_base": 0.85,
        "macs_recovery_base": 0.65,
        "alternative_antibodies": ["MHCII", "CD103", "XCR1"],
        "density_g_ml": 1.078,
    },
    "neutrophils": {
        "frequency_in_tissue": 0.30,
        "macs_antibody": "Ly6G",
        "macs_purity_base": 0.93,
        "macs_recovery_base": 0.82,
        "alternative_antibodies": ["CD11b", "Gr1", "CXCR2"],
        "density_g_ml": 1.080,
    },
    "eosinophils": {
        "frequency_in_tissue": 0.03,
        "macs_antibody": "SiglecF",
        "macs_purity_base": 0.87,
        "macs_recovery_base": 0.68,
        "alternative_antibodies": ["CCR3", "CD11b", "IL5Ra"],
        "density_g_ml": 1.082,
    },
    "mast_cells": {
        "frequency_in_tissue": 0.01,
        "macs_antibody": "cKit",
        "macs_purity_base": 0.80,
        "macs_recovery_base": 0.60,
        "alternative_antibodies": ["FcεRI", "CD117", "Tryptase"],
        "density_g_ml": 1.074,
    },
    "monocytes": {
        "frequency_in_tissue": 0.08,
        "macs_antibody": "CD14",
        "macs_purity_base": 0.89,
        "macs_recovery_base": 0.72,
        "alternative_antibodies": ["CD16", "Ly6C", "CCR2"],
        "density_g_ml": 1.077,
    },
    "platelets": {
        "frequency_in_tissue": 0.0,
        "macs_antibody": "CD41",
        "macs_purity_base": 0.95,
        "macs_recovery_base": 0.85,
        "alternative_antibodies": ["CD61", "CD42b", "GPVI"],
        "density_g_ml": 1.060,
    },
}


# =============================================================================
# PARSING
# =============================================================================

def parse_input(data_str: str) -> Dict:
    """Parse JSON input and validate."""
    try:
        data = json.loads(data_str)
    except json.JSONDecodeError as e:
        raise ValueError(f"Invalid JSON: {e}")
    
    required = ["tissue_type", "target_cell_type"]
    for key in required:
        if key not in data:
            raise ValueError(f"Missing required key: '{key}'")
    
    # Normalize cell type
    cell_type = data["target_cell_type"].lower().replace(" ", "_").replace("-", "_")
    data["target_cell_type"] = cell_type
    
    # Set defaults
    data.setdefault("enzyme_type", "collagenase")
    data.setdefault("digestion_time_min", 45)
    data.setdefault("macs_antibody", None)
    
    # Validate tissue
    tissue = data["tissue_type"].lower()
    if tissue not in TISSUE_PROFILES:
        valid = ", ".join(TISSUE_PROFILES.keys())
        raise ValueError(f"Unknown tissue type '{tissue}'. Valid: {valid}")
    data["tissue_type"] = tissue
    
    # Validate enzyme
    enzyme = data["enzyme_type"].lower()
    if enzyme not in ENZYME_PROFILES:
        valid = ", ".join(ENZYME_PROFILES.keys())
        raise ValueError(f"Unknown enzyme '{enzyme}'. Valid: {valid}")
    data["enzyme_type"] = enzyme
    
    # Validate cell type
    if cell_type not in CELL_TYPE_PROFILES:
        valid = ", ".join(CELL_TYPE_PROFILES.keys())
        raise ValueError(f"Unknown cell type '{cell_type}'. Valid: {valid}")
    
    # Auto-select MACS antibody if not provided
    if data["macs_antibody"] is None:
        data["macs_antibody"] = CELL_TYPE_PROFILES[cell_type]["macs_antibody"]
    
    return data


# =============================================================================
# SIMULATION ENGINE
# =============================================================================

def simulate_isolation(params: Dict) -> Dict:
    """Run the full isolation and purification simulation."""
    tissue = TISSUE_PROFILES[params["tissue_type"]]
    enzyme = ENZYME_PROFILES[params["enzyme_type"]]
    cell_type = CELL_TYPE_PROFILES[params["target_cell_type"]]
    
    digestion_time = params["digestion_time_min"]
    macs_antibody = params["macs_antibody"]
    
    # Determine sample amount
    if params["tissue_type"] in ("blood", "bone_marrow"):
        sample_amount = tissue.get("typical_volume_ml", 10)
        total_cells_start = tissue["total_cells_per_ml"] * sample_amount
        unit = "mL"
    else:
        sample_amount = tissue.get("typical_yield_mg", 100)
        total_cells_start = tissue["total_cells_per_mg"] * sample_amount
        unit = "mg"
    
    # Step 1: Tissue harvesting and mincing
    harvest_loss = random.uniform(0.90, 0.98)
    cells_after_harvest = int(total_cells_start * harvest_loss)
    
    # Step 2: Enzymatic digestion
    # Digestion efficiency depends on time vs optimal
    time_factor = 1.0 - abs(digestion_time - enzyme["optimal_time_min"]) / enzyme["max_time_min"]
    time_factor = max(0.3, min(1.0, time_factor))
    
    digestion_eff = tissue["digestion_efficiency_base"] * enzyme["digestion_rate"] * time_factor
    digestion_eff = min(0.95, digestion_eff)
    
    # Cytotoxicity increases with time
    overtime_ratio = max(0, (digestion_time - enzyme["optimal_time_min"]) / enzyme["max_time_min"])
    cytotoxicity = enzyme["cytotoxicity"] * (1 + overtime_ratio)
    
    cells_after_digestion = int(cells_after_harvest * digestion_eff)
    viability_after_digestion = max(0.55, 1.0 - cytotoxicity - tissue["cell_fragility"] * 0.5)
    
    # Step 3: Filtration and washing
    filter_recovery = random.uniform(0.85, 0.95)
    cells_after_filter = int(cells_after_digestion * filter_recovery)
    
    # Step 4: Density gradient (optional, for blood/bone marrow)
    if params["tissue_type"] in ("blood", "bone_marrow"):
        gradient_recovery = random.uniform(0.60, 0.80)
        cells_after_gradient = int(cells_after_filter * gradient_recovery)
        gradient_viability = random.uniform(0.90, 0.98)
        viability_after_gradient = viability_after_digestion * gradient_viability
        use_gradient = True
    else:
        cells_after_gradient = cells_after_filter
        viability_after_gradient = viability_after_digestion
        use_gradient = False
    
    # Step 5: RBC lysis (for spleen, blood, bone marrow)
    if params["tissue_type"] in ("spleen", "blood", "bone_marrow"):
        rbc_lysis_eff = random.uniform(0.85, 0.95)
        cells_after_rbc = int(cells_after_gradient * 0.60)  # RBCs removed, only nucleated remain
        rbc_ratio = random.uniform(0.40, 0.60)  # RBC fraction removed
        use_rbc_lysis = True
    else:
        cells_after_rbc = cells_after_gradient
        rbc_ratio = 0.0
        use_rbc_lysis = False
    
    # Step 6: MACS enrichment
    macs_purity = cell_type["macs_purity_base"] * random.uniform(0.95, 1.05)
    macs_purity = min(0.98, macs_purity)
    
    macs_recovery = cell_type["macs_recovery_base"] * random.uniform(0.90, 1.10)
    macs_recovery = min(0.95, macs_recovery)
    
    target_cells_pre_macs = int(cells_after_rbc * cell_type["frequency_in_tissue"])
    target_cells_after_macs = int(target_cells_pre_macs * macs_recovery)
    
    non_target_contaminants = int(cells_after_rbc * (1 - cell_type["frequency_in_tissue"]) * (1 - macs_purity) * 0.1)
    total_cells_after_macs = target_cells_after_macs + non_target_contaminants
    
    final_purity = target_cells_after_macs / total_cells_after_macs if total_cells_after_macs > 0 else 0
    final_viability = viability_after_gradient * random.uniform(0.92, 0.98)
    
    # Calculate overall metrics
    overall_recovery = target_cells_after_macs / total_cells_start if total_cells_start > 0 else 0
    enrichment_factor = final_purity / cell_type["frequency_in_tissue"] if cell_type["frequency_in_tissue"] > 0 else 0
    
    return {
        "params": params,
        "sample_amount": sample_amount,
        "unit": unit,
        "total_cells_start": total_cells_start,
        "cells_after_harvest": cells_after_harvest,
        "cells_after_digestion": cells_after_digestion,
        "viability_after_digestion": viability_after_digestion,
        "cells_after_filter": cells_after_filter,
        "cells_after_gradient": cells_after_gradient,
        "viability_after_gradient": viability_after_gradient,
        "use_gradient": use_gradient,
        "cells_after_rbc": cells_after_rbc,
        "rbc_ratio": rbc_ratio,
        "use_rbc_lysis": use_rbc_lysis,
        "target_cells_pre_macs": target_cells_pre_macs,
        "target_cells_after_macs": target_cells_after_macs,
        "non_target_contaminants": non_target_contaminants,
        "total_cells_after_macs": total_cells_after_macs,
        "final_purity": final_purity,
        "final_viability": final_viability,
        "overall_recovery": overall_recovery,
        "enrichment_factor": enrichment_factor,
        "macs_purity": macs_purity,
        "macs_recovery": macs_recovery,
        "digestion_efficiency": digestion_eff,
        "cytotoxicity": cytotoxicity,
        "time_factor": time_factor,
        "harvest_loss": harvest_loss,
        "filter_recovery": filter_recovery,
    }


# =============================================================================
# CSV GENERATION
# =============================================================================

def generate_csv(result: Dict, output_path: str) -> str:
    """Generate cell count and viability CSV file."""
    headers = [
        "Step", "Total_Cells", "Target_Cells", "Viability_%",
        "Purity_%", "Recovery_%", "Notes"
    ]
    
    rows = [
        ["Initial Tissue", result["total_cells_start"], "—", "—", "—", "—", f"{result['sample_amount']} {result['unit']}"],
        ["Harvest/Mince", result["cells_after_harvest"], "—", f"{result['harvest_loss']*100:.1f}%", "—", "—", "Mechanical disruption"],
        ["Digestion", result["cells_after_digestion"], "—", f"{result['viability_after_digestion']*100:.1f}%", "—", f"{result['digestion_efficiency']*100:.1f}%", f"{result['params']['enzyme_type']}, {result['params']['digestion_time_min']} min"],
        ["Filtration", result["cells_after_filter"], "—", f"{result['viability_after_digestion']*100:.1f}%", "—", f"{result['filter_recovery']*100:.1f}%", "70 μm cell strainer"],
    ]
    
    if result["use_gradient"]:
        rows.append(["Density Gradient", result["cells_after_gradient"], "—", f"{result['viability_after_gradient']*100:.1f}%", "—", f"{result['cells_after_gradient']/result['cells_after_filter']*100:.1f}%", "Ficoll-Paque separation"])
    
    if result["use_rbc_lysis"]:
        rows.append(["RBC Lysis", result["cells_after_rbc"], "—", f"{result['viability_after_gradient']*100:.1f}%", "—", f"{(1-result['rbc_ratio'])*100:.1f}%", f"Removed {result['rbc_ratio']*100:.1f}% RBCs"])
    
    rows.append([
        "Pre-MACS", result["cells_after_rbc"], result["target_cells_pre_macs"],
        f"{result['viability_after_gradient']*100:.1f}%",
        f"{result['params']['target_cell_type'].replace('_', ' ').title()} {CELL_TYPE_PROFILES[result['params']['target_cell_type']]['frequency_in_tissue']*100:.1f}%",
        "—",
        "Before enrichment"
    ])
    
    rows.append([
        "MACS (+)", result["total_cells_after_macs"], result["target_cells_after_macs"],
        f"{result['final_viability']*100:.1f}%",
        f"{result['final_purity']*100:.1f}%",
        f"{result['macs_recovery']*100:.1f}%",
        f"Anti-{result['params']['macs_antibody']} microbeads"
    ])
    
    rows.append([
        "MACS (flow-through)",
        result["cells_after_rbc"] - result["total_cells_after_macs"],
        result["target_cells_pre_macs"] - result["target_cells_after_macs"],
        f"{result['viability_after_gradient']*100:.1f}%",
        f"{(1-result['final_purity'])*100:.1f}%",
        f"{(1-result['macs_recovery'])*100:.1f}%",
        "Depleted fraction"
    ])
    
    with open(output_path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(headers)
        writer.writerows(rows)
    
    return output_path


# =============================================================================
# LOG FORMATTING
# =============================================================================

def format_log_report(result: Dict, csv_path: str) -> str:
    """Format a detailed text analysis log."""
    p = result["params"]
    tissue_name = p["tissue_type"].replace("_", " ").title()
    cell_name = p["target_cell_type"].replace("_", " ").title()
    enzyme_name = p["enzyme_type"].title()
    
    lines = []
    lines.append("=" * 75)
    lines.append("IMMUNE CELL ISOLATION & PURIFICATION SIMULATION LOG")
    lines.append("=" * 75)
    lines.append(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    lines.append("")
    
    # Section 1: Input Parameters
    lines.append("-" * 75)
    lines.append("SECTION 1: INPUT PARAMETERS")
    lines.append("-" * 75)
    lines.append("")
    lines.append(f"  Tissue Type:          {tissue_name}")
    lines.append(f"  Target Cell Type:     {cell_name}")
    lines.append(f"  Digestion Enzyme:     {enzyme_name}")
    lines.append(f"  Digestion Time:       {p['digestion_time_min']} minutes")
    lines.append(f"  MACS Antibody:        Anti-{p['macs_antibody']} microbeads")
    lines.append(f"  Sample Amount:        {result['sample_amount']} {result['unit']} (typical)")
    lines.append("")
    
    # Tissue profile
    tissue_prof = TISSUE_PROFILES[p["tissue_type"]]
    lines.append("  Tissue Characteristics:")
    lines.append(f"    • Recommended enzymes: {', '.join(tissue_prof['recommended_enzymes'])}")
    lines.append(f"    • Cell fragility score:  {tissue_prof['cell_fragility']*100:.0f}%")
    lines.append(f"    • Digestion efficiency:  {tissue_prof['digestion_efficiency_base']*100:.0f}% (base)")
    lines.append("")
    
    # Cell type profile
    cell_prof = CELL_TYPE_PROFILES[p["target_cell_type"]]
    lines.append("  Target Cell Characteristics:")
    lines.append(f"    • Natural frequency:     {cell_prof['frequency_in_tissue']*100:.1f}% in tissue")
    lines.append(f"    • Default MACS antibody: {cell_prof['macs_antibody']}")
    lines.append(f"    • Alternative markers:   {', '.join(cell_prof['alternative_antibodies'])}")
    lines.append(f"    • Expected MACS purity:  {cell_prof['macs_purity_base']*100:.0f}%")
    lines.append(f"    • Expected recovery:     {cell_prof['macs_recovery_base']*100:.0f}%")
    lines.append("")
    
    # Section 2: Step-by-Step Simulation
    lines.append("-" * 75)
    lines.append("SECTION 2: STEP-BY-STEP SIMULATION")
    lines.append("-" * 75)
    lines.append("")
    
    # Step 1
    lines.append("STEP 1: TISSUE HARVESTING & MINCING")
    lines.append(f"  • Sample:              {result['sample_amount']} {result['unit']} of {tissue_name}")
    lines.append(f"  • Estimated cells:     {result['total_cells_start']:,}")
    lines.append(f"  • Harvest recovery:    {result['harvest_loss']*100:.1f}%")
    lines.append(f"  • Cells after harvest: {result['cells_after_harvest']:,}")
    lines.append("")
    
    # Step 2
    lines.append("STEP 2: ENZYMATIC DIGESTION")
    lines.append(f"  • Enzyme:              {enzyme_name}")
    lines.append(f"  • Digestion time:      {p['digestion_time_min']} min (optimal: {ENZYME_PROFILES[p['enzyme_type']]['optimal_time_min']} min)")
    lines.append(f"  • Time efficiency:     {result['time_factor']*100:.1f}%")
    lines.append(f"  • Cytotoxicity:        {result['cytotoxicity']*100:.1f}%")
    lines.append(f"  • Digestion efficiency: {result['digestion_efficiency']*100:.1f}%")
    lines.append(f"  • Cells released:      {result['cells_after_digestion']:,}")
    lines.append(f"  • Viability:           {result['viability_after_digestion']*100:.1f}%")
    if result['cytotoxicity'] > 0.20:
        lines.append("  ⚠️ WARNING: High cytotoxicity. Consider reducing digestion time or switching to a gentler enzyme.")
    lines.append("")
    
    # Step 3
    lines.append("STEP 3: FILTRATION & WASHING")
    lines.append(f"  • Filter:              70 μm cell strainer")
    lines.append(f"  • Recovery:            {result['filter_recovery']*100:.1f}%")
    lines.append(f"  • Cells after filter:  {result['cells_after_filter']:,}")
    lines.append("")
    
    # Step 4 (conditional)
    if result["use_gradient"]:
        lines.append("STEP 4: DENSITY GRADIENT CENTRIFUGATION")
        lines.append(f"  • Medium:              Ficoll-Paque (1.077 g/mL)")
        grad_recov = result['cells_after_gradient'] / result['cells_after_filter'] * 100
        lines.append(f"  • Recovery:            {grad_recov:.1f}%")
        lines.append(f"  • Cells recovered:     {result['cells_after_gradient']:,}")
        lines.append(f"  • Viability:           {result['viability_after_gradient']*100:.1f}%")
        lines.append("")
    
    # Step 5 (conditional)
    if result["use_rbc_lysis"]:
        lines.append("STEP 5: RBC LYSIS")
        lines.append(f"  • Method:              ACK lysis buffer")
        lines.append(f"  • RBC fraction:        {result['rbc_ratio']*100:.1f}% (removed)")
        lines.append(f"  • Nucleated cells:     {result['cells_after_rbc']:,}")
        lines.append("")
    
    # Step 6
    lines.append("STEP 6: MACS (MAGNETIC-ACTIVATED CELL SORTING)")
    lines.append(f"  • Strategy:            Positive selection")
    lines.append(f"  • Antibody:            Anti-{p['macs_antibody']} microbeads")
    lines.append(f"  • Expected purity:     {result['macs_purity']*100:.1f}%")
    lines.append(f"  • Expected recovery:   {result['macs_recovery']*100:.1f}%")
    lines.append("")
    
    # Section 3: Final Results
    lines.append("-" * 75)
    lines.append("SECTION 3: FINAL RESULTS")
    lines.append("-" * 75)
    lines.append("")
    
    lines.append("Pre-MACS Cell Population:")
    lines.append(f"  • Total cells:         {result['cells_after_rbc']:,}")
    lines.append(f"  • Target cells:        {result['target_cells_pre_macs']:,} ({cell_prof['frequency_in_tissue']*100:.1f}% of total)")
    lines.append(f"  • Viability:           {result['viability_after_gradient']*100:.1f}%")
    lines.append("")
    
    lines.append("Post-MACS (+) Fraction:")
    lines.append(f"  • Total cells:         {result['total_cells_after_macs']:,}")
    lines.append(f"  • Target cells:        {result['target_cells_after_macs']:,}")
    lines.append(f"  • Purity:              {result['final_purity']*100:.1f}%")
    lines.append(f"  • Viability:           {result['final_viability']*100:.1f}%")
    lines.append(f"  • MACS recovery:       {result['macs_recovery']*100:.1f}%")
    lines.append("")
    
    lines.append("Post-MACS (flow-through) Fraction:")
    lines.append(f"  • Total cells:         {result['cells_after_rbc'] - result['total_cells_after_macs']:,}")
    lines.append(f"  • Target cells lost:   {result['target_cells_pre_macs'] - result['target_cells_after_macs']:,}")
    lines.append("")
    
    lines.append("Overall Performance Metrics:")
    lines.append(f"  • Overall recovery:    {result['overall_recovery']*100:.2f}% (from starting tissue)")
    lines.append(f"  • Enrichment factor:   {result['enrichment_factor']:.1f}x")
    lines.append(f"  • Final yield:         {result['target_cells_after_macs']:,} {cell_name}")
    lines.append("")
    
    # Section 4: Recommendations
    lines.append("-" * 75)
    lines.append("SECTION 4: RECOMMENDATIONS & QUALITY ASSESSMENT")
    lines.append("-" * 75)
    lines.append("")
    
    # Quality flags
    quality_flags = []
    if result['final_purity'] < 0.85:
        quality_flags.append(f"  ⚠️ Low purity ({result['final_purity']*100:.1f}%). Consider a two-step MACS or FACS sorting.")
    if result['final_viability'] < 0.80:
        quality_flags.append(f"  ⚠️ Low viability ({result['final_viability']*100:.1f}%). Reduce digestion time or use gentler conditions.")
    if result['overall_recovery'] < 0.01:
        quality_flags.append(f"  ⚠️ Very low overall recovery ({result['overall_recovery']*100:.2f}%). Verify tissue quality and cell frequency.")
    if result['cytotoxicity'] > 0.20:
        quality_flags.append(f"  ⚠️ High cytotoxicity from digestion. Consider switching from {enzyme_name} to Liberase or Accutase.")
    if p['enzyme_type'] not in tissue_prof['recommended_enzymes']:
        quality_flags.append(f"  ⚠️ {enzyme_name} is not typically recommended for {tissue_name}. Consider: {', '.join(tissue_prof['recommended_enzymes'])}.")
    if p['digestion_time_min'] > ENZYME_PROFILES[p['enzyme_type']]['max_time_min']:
        quality_flags.append(f"  ⚠️ Digestion time ({p['digestion_time_min']} min) exceeds recommended maximum ({ENZYME_PROFILES[p['enzyme_type']]['max_time_min']} min).")
    
    if quality_flags:
        lines.append("Quality Flags:")
        for flag in quality_flags:
            lines.append(flag)
        lines.append("")
    else:
        lines.append("✅ No major quality concerns detected.")
        lines.append("")
    
    # Recommendations
    lines.append("Optimization Suggestions:")
    if result['final_purity'] < 0.90:
        lines.append(f"  • For higher purity: Use Anti-{cell_prof['macs_antibody']} MACS followed by FACS with {cell_prof['alternative_antibodies'][0]}.")
    if result['final_viability'] < 0.85:
        lines.append(f"  • For higher viability: Reduce digestion to {max(15, ENZYME_PROFILES[p['enzyme_type']]['optimal_time_min'] - 10)} min and keep on ice.")
    if result['overall_recovery'] < 0.05:
        lines.append(f"  • For higher yield: Increase starting tissue to {result['sample_amount'] * 2:.0f} {result['unit']}.")
    lines.append(f"  • Recommended medium: RPMI 1640 + 10% FBS + 1% Pen/Strep for {cell_name} culture.")
    lines.append(f"  • Expected doubling time: ~{random.randint(12, 48)} hours (varies by activation state).")
    lines.append("")
    
    # Section 5: File Output
    lines.append("-" * 75)
    lines.append("SECTION 5: OUTPUT FILES")
    lines.append("-" * 75)
    lines.append("")
    lines.append(f"  CSV Report:            {csv_path}")
    lines.append(f"  • Contains:            Cell counts, viability, purity, and recovery per step")
    lines.append(f"  • Rows:                All isolation steps + MACS fractions")
    lines.append("")
    
    lines.append("=" * 75)
    lines.append("END OF SIMULATION LOG")
    lines.append("=" * 75)
    
    return "\n".join(lines)


# =============================================================================
# MAIN
# =============================================================================

def main():
    parser = argparse.ArgumentParser(description="Immune cell isolation & purification simulator")
    parser.add_argument("--data", required=True, help="JSON string with simulation parameters")
    parser.add_argument("--output-csv", default="/tmp/immune_cell_isolation.csv", help="Output CSV path")
    parser.add_argument("--output-format", default="text", choices=["text", "json"], help="Output format")
    parser.add_argument("--seed", type=int, default=None, help="Random seed for reproducibility")
    
    args = parser.parse_args()
    
    if args.seed is not None:
        random.seed(args.seed)
    
    try:
        params = parse_input(args.data)
        result = simulate_isolation(params)
        csv_path = generate_csv(result, args.output_csv)
        
        if args.output_format == "json":
            output = {
                "params": params,
                "result": {k: v for k, v in result.items() if k != "params"},
                "csv_path": csv_path,
            }
            print(json.dumps(output, indent=2))
        else:
            print(format_log_report(result, csv_path))
            print(f"\nCSV file saved to: {csv_path}")
    
    except Exception as e:
        print(f"ERROR: {e}", file=sys.stderr)
        import traceback
        traceback.print_exc(file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
