# -*- coding: utf-8 -*-
"""
生物医药大类 · 输出效果预览（demo_visualization）精修覆盖。

合并顺序：Phase1 → Phase2 → Bulk → **本模块**（仅覆盖 demo_visualization 等显式字段）。
"""
from __future__ import annotations

from typing import Any, Dict

from gibh_agent.db.skill_detail_demo_visualizations import (
    DEMO_VIZ_ADMET_SCREENING,
    DEMO_VIZ_ANTIBODY_HUMANNESS,
    DEMO_VIZ_BEPIPRED_ENHANCED,
    DEMO_VIZ_CD_SPECTRUM,
    DEMO_VIZ_CRISPR_EDIT,
    DEMO_VIZ_ENZYME_KINETICS,
    DEMO_VIZ_IMMUNE_CELL_ISOLATION,
    DEMO_VIZ_MHC_EPITOPE,
    DEMO_VIZ_MOLECULAR_DOCKING,
    DEMO_VIZ_PROTEIN_HOMOLOGY,
    DEMO_VIZ_RNAFOLD,
)
from gibh_agent.db.skill_detail_demo_visualizations_biopharma_remainder import (
    BIOPHARMA_REMAINDER_VIZ,
)

# tool_id → 需覆盖的 detailed_spec 字段（当前以 demo_visualization 为主）
SKILL_DETAILED_SPECS_BIOPHARMA_VIZ: Dict[str, Dict[str, Any]] = {
    "bepipred3_prediction": {
        "demo_visualization": DEMO_VIZ_BEPIPRED_ENHANCED,
    },
    "rnafold_analysis": {
        "demo_visualization": DEMO_VIZ_RNAFOLD,
    },
    "protein_homology_structure_assessment": {
        "demo_visualization": DEMO_VIZ_PROTEIN_HOMOLOGY,
    },
    "crispr_cas9_simulation": {
        "demo_visualization": DEMO_VIZ_CRISPR_EDIT,
    },
    "antibody_humanness_oasis_evaluation": {
        "demo_visualization": DEMO_VIZ_ANTIBODY_HUMANNESS,
    },
    "circular_dichroism_analysis": {
        "demo_visualization": DEMO_VIZ_CD_SPECTRUM,
    },
    "enzyme_kinetics_analysis": {
        "demo_visualization": DEMO_VIZ_ENZYME_KINETICS,
    },
    "mhc_epitope_search": {
        "demo_visualization": DEMO_VIZ_MHC_EPITOPE,
    },
    "immune_cell_isolation_simulation": {
        "demo_visualization": DEMO_VIZ_IMMUNE_CELL_ISOLATION,
    },
    "admet_property_screening": {
        "demo_visualization": DEMO_VIZ_ADMET_SCREENING,
    },
    "diffdock_pose_preview": {
        "demo_visualization": DEMO_VIZ_MOLECULAR_DOCKING,
    },
}

# Phase 1.3：剩余带 [Skill_Route] 的生物医药技能
for _tid, _html in BIOPHARMA_REMAINDER_VIZ.items():
    SKILL_DETAILED_SPECS_BIOPHARMA_VIZ[_tid] = {"demo_visualization": _html}
