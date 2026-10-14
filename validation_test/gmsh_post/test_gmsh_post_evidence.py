from __future__ import annotations

import json
import math
from pathlib import Path


HERE = Path(__file__).resolve().parent
SCHEMA = "radia.validation.gmsh-post.v1"


def _load(name: str, case: str) -> dict:
    result = json.loads((HERE / name).read_text(encoding="utf-8"))
    assert result["schema"] == SCHEMA
    assert result["case"] == case
    assert result["source_notebook"] == f"docs/gmsh_post/{case}.ipynb"
    return result


def test_fieldline_evidence_preserves_analytic_and_flux_checks():
    result = _load("em_fieldlines_results.json", "em_fieldlines")["results"]
    assert result["case_a"]["B0_rel_err"] < 1.0e-3
    assert result["case_a"]["psi_rel_drift"] < 1.0e-2
    assert result["case_a"]["n_stream_lines"] > 0
    assert result["case_b"]["n_stream_lines"] > 0
    assert len(result["case_b"]["iso_levels_T"]) >= 3


def test_gallery_evidence_preserves_solver_and_postprocessing_checks():
    result = _load("em_post_gallery_results.json", "em_post_gallery")
    assert result["iron_enhancement"] > 1.0
    assert result["solve_outer_iterations"] > 0
    assert 0.0 < result["selection_fraction"] < 1.0
    assert 0.0 < result["raycast_cad_covered_fraction"] < 1.0
    assert result["lic_outline_segments"] > 0
    assert all(math.isfinite(value) and value > 0.0 for value in result["sweep_domain_max_T"])


def test_particle_orbit_evidence_preserves_tracking_checks():
    result = _load("em_particle_orbits_results.json", "em_particle_orbits")
    rays = result["dispersion"]["rays"]
    assert len(rays) >= 3
    assert max(ray["ke_drift_rel"] for ray in rays) < 1.0e-8
    quadrupole = result["quadrupole"]
    assert abs(quadrupole["focus_mean_m"] / quadrupole["f_thick_m"] - 1.0) < 0.02
    assert quadrupole["focus_spread_m"] < 1.0e-3
    edge_rows = result["edge_focusing"]["rows"]
    assert min(abs(row["ratio"] - 1.0) for row in edge_rows) < 1.0e-3
