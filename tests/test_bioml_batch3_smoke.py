# -*- coding: utf-8 -*-
"""第三批生物医药 batch3 工具冒烟。"""
from __future__ import annotations

import asyncio
import importlib.util
import json
from pathlib import Path

from gibh_agent.tools import load_all_tools
from gibh_agent.core.tool_registry import registry
import pytest

from gibh_agent.tools.bioml_batch3_tools import (
    antibody_humanness_oasis_evaluation,
    circular_dichroism_analysis,
    enzyme_kinetics_analysis,
    immune_cell_isolation_simulation,
    itc_binding_thermodynamic_analysis,
    protease_kinetics_analysis,
    protein_homology_structure_assessment,
    protein_sequence_conservation_analysis,
    rna_secondary_structure_analysis,
    cell_cycle_phase_duration_estimation,
    _asset_script,
)


def test_registry_has_bioml_tools():
    load_all_tools()
    for name in (
        "crispr_cas9_simulation",
        "immune_cell_isolation_simulation",
        "rna_secondary_structure_analysis",
        "circular_dichroism_analysis",
        "itc_binding_thermodynamic_analysis",
        "cell_cycle_phase_duration_estimation",
        "protease_kinetics_analysis",
        "enzyme_kinetics_analysis",
        "protein_sequence_conservation_analysis",
        "protein_homology_structure_assessment",
        "antibody_humanness_oasis_evaluation",
    ):
        assert registry.get_tool(name), name


def test_asset_scripts_exist():
    for fn in (
        "immune_cell_isolation.py",
        "rna_structure.py",
        "cd_analysis.py",
        "crispr_cas9.py",
        "itc_analysis.py",
        "cell_cycle_analysis.py",
        "protease_kinetics.py",
        "enzyme_kinetics.py",
        "sequence_conservation.py",
        "protein_assess.py",
        "oasis_analysis.py",
    ):
        assert _asset_script(fn).is_file(), fn


def test_immune_cell_smoke():
    sim = json.dumps(
        {
            "tissue_type": "spleen",
            "target_cell_type": "macrophages",
            "enzyme_type": "collagenase",
            "digestion_time_min": 45,
        },
        ensure_ascii=False,
    )
    r = asyncio.run(immune_cell_isolation_simulation(simulation_json=sim, output_format="text", random_seed=1))
    assert r.get("status") == "success", r
    assert "markdown" in r and len((r.get("markdown") or "")) > 20


def test_rna_structure_json_smoke():
    sj = json.dumps({"structure": "..((..))", "sequence": "AAGGCCUU"}, ensure_ascii=False)
    r = asyncio.run(rna_secondary_structure_analysis(structure_json=sj, file_path="", result_format="json"))
    assert r.get("status") == "success", r
    md = r.get("markdown") or ""
    assert "RNA" in md or "json" in md.lower(), md[:400]
    od = r.get("output_dir")
    assert od and Path(od).is_dir()


def test_cd_smoke():
    spec = {
        "sample_name": "Demo",
        "sample_type": "protein",
        "wavelength_data": [190, 200, 210, 220],
        "cd_signal_data": [0.1, -0.2, 0.3, -0.4],
    }
    r = asyncio.run(circular_dichroism_analysis(spectrum_json=json.dumps(spec), output_format="text"))
    assert r.get("status") == "success", r
    assert (r.get("markdown") or "").strip(), r


def test_itc_smoke():
    inj = json.dumps([[1, 2e-6, -1200], [2, 2e-6, -2100], [3, 2e-6, -2800], [4, 2e-6, -3100], [5, 2e-6, -2900]])
    r = asyncio.run(
        itc_binding_thermodynamic_analysis(
            protein_conc=1e-6,
            ligand_conc=10e-6,
            itc_injections_json=inj,
            file_path="",
        )
    )
    assert r.get("status") == "success", r
    md = r.get("markdown") or ""
    assert "ITC" in md and ("base64" in md.lower() or "data:image" in md), md[:500]


def test_cell_cycle_smoke():
    fc = json.dumps(
        {
            "time_points": [0, 2, 4, 6, 8],
            "edu_positive": [5, 15, 28, 35, 38],
            "brdu_positive": [3, 12, 30, 40, 42],
            "double_positive": [1, 8, 18, 25, 28],
        }
    )
    r = asyncio.run(cell_cycle_phase_duration_estimation(flow_cytometry_json=fc, initial_estimates_json=""))
    assert r.get("status") == "success", r
    assert "G1" in (r.get("markdown") or ""), r


def test_protease_smoke():
    kin = json.dumps(
        {
            "time_points": [0, 30, 60, 90],
            "substrate_concentrations": [10, 40],
            "fluorescence_data": [[100, 150, 180, 200], [100, 130, 150, 165]],
            "enzyme_concentration": 0.1,
        }
    )
    r = asyncio.run(protease_kinetics_analysis(kinetics_input_json=kin))
    assert r.get("status") == "success", r
    md = r.get("markdown") or ""
    assert "蛋白酶" in md and ("data:image" in md or "base64" in md.lower()), md[:600]


def test_enzyme_smoke():
    kin = json.dumps(
        {
            "enzyme_name": "DemoEnzyme",
            "enzyme_concentration": 50,
            "substrate_concentrations": [5, 10, 25, 50, 100],
        }
    )
    r = asyncio.run(enzyme_kinetics_analysis(kinetics_input_json=kin))
    assert r.get("status") == "success", r
    md = r.get("markdown") or ""
    assert "酶动力学" in md and ("data:image" in md or "base64" in md.lower()), md[:600]


def test_sequence_conservation_smoke():
    aj = json.dumps({"protein_sequences": ["MKTLLL", "MKTLLL", "MKALLL"]})
    r = asyncio.run(protein_sequence_conservation_analysis(alignment_json=aj))
    assert r.get("status") == "success", r
    assert "保守" in (r.get("markdown") or ""), r


def test_protein_homology_test_mode_smoke():
  r = asyncio.run(
      protein_homology_structure_assessment(uniprot_id="O95292", test_mode=True, top_n=3)
  )
  assert r.get("status") == "success", r
  md = r.get("markdown") or ""
  assert "同源" in md or "UniProt" in md or "O95292" in md, md[:800]


def test_oasis_json_param_smoke():
    aj = json.dumps(
        {
            "heavy_chain": "EVQLVESGGGLVQPGGSLRLSCAASGFTFSSYAMSWVRQAPGKGLEWVSAISGSGGSTYYADSVKGRFTISRDNSKNTLYLQMNSLRAEDTAVYYCARDYGDYWGQGTLVTVSS",
            "light_chain": "DIQMTQSPSSLSASVGDRVTITCRASQDVNTAVAWYQQKPGKAPKLLIYSASFLYSGVPSRFSGSRSGTDFTLTISSLQPEDFATYYCQQHYTTPPTFGQGTKVEIK",
        }
    )
    r = asyncio.run(
        antibody_humanness_oasis_evaluation(
            heavy_chain="",
            light_chain="",
            antibody_sequences_json=aj,
        )
    )
    assert r.get("status") == "success", r


@pytest.mark.skipif(importlib.util.find_spec("promb") is None, reason="promb not installed")
def test_oasis_humanness_smoke():
    hc = "EVQLVESGGGLVQPGGSLRLSCAASGFTFSSYAMSWVRQAPGKGLEWVSAISGSGGSTYYADSVKGRFTISRDNSKNTLYLQMNSLRAEDTAVYYCARDYGDYWGQGTLVTVSS"
    lc = "DIQMTQSPSSLSASVGDRVTITCRASQDVNTAVAWYQQKPGKAPKLLIYSASFLYSGVPSRFSGSRSGTDFTLTISSLQPEDFATYYCQQHYTTPPTFGQGTKVEIK"
    r = asyncio.run(antibody_humanness_oasis_evaluation(heavy_chain=hc, light_chain=lc, antibody_name="DemoAb"))
    assert r.get("status") == "success", r
    assert "OASis" in (r.get("markdown") or "") or "人源性" in (r.get("markdown") or ""), r
