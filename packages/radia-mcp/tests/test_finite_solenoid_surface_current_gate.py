import copy
import json
import math

from radia_mcp.radia_ngsolve.finite_solenoid_surface_current_gate import finite_solenoid_surface_current_gate
from radia_mcp.radia_ngsolve.server import finite_solenoid_surface_current_gate as mcp_gate


def good_summary():
    z = [float(value) for value in range(-40, 41)]
    mu0 = 4.0 * math.pi * 1.0e-7
    analytic = [mu0 / 2.0 * ((value + 5.0) / math.sqrt(100 + (value + 5.0) ** 2) - (value - 5.0) / math.sqrt(100 + (value - 5.0) ** 2)) for value in z]
    auto = [value * 1.018 for value in analytic]
    best = [value * 1.004 for value in analytic]
    return {
        "geometry": {"radius": 10.0, "length": 10.0, "center_z": 0.0},
        "units": {"length": "mm", "field": "T", "surface_current_density": "A/m"},
        "current_density_a_per_m": 1.0,
        "profile_z": z,
        "levels": [
            {"label": "source_auto", "tet_count": 10000, "minimum_mesh_quality": 0.2, "Bx_T": [0.0] * len(z), "By_T": [0.0] * len(z), "Bz_T": auto},
            {"label": "candidate_fine", "tet_count": 140000, "minimum_mesh_quality": 0.18, "Bx_T": [0.0] * len(z), "By_T": [0.0] * len(z), "Bz_T": [value * 1.01 for value in analytic]},
            {"label": "candidate_best", "tet_count": 42000, "minimum_mesh_quality": 0.19, "Bx_T": [0.0] * len(z), "By_T": [0.0] * len(z), "Bz_T": best},
        ],
        "linearity": {"current_density_a_per_m": [1.0, 2.0, -1.0], "Bz_T": [best, [2.0 * value for value in best], [-value for value in best]]},
        "gate_tolerances": {
            "minimum_samples": 81,
            "maximum_center_relative_error": 0.02,
            "maximum_profile_error_normalized_by_center": 0.02,
            "maximum_rms_error_normalized_by_center": 0.005,
            "maximum_mirror_symmetry_relative": 0.015,
            "maximum_transverse_leakage_relative": 0.025,
            "maximum_linearity_relative_error": 1.0e-10,
            "minimum_baseline_to_best_improvement": 1.5,
        },
    }


def test_accepts_analytic_profile_and_signed_linearity():
    result = finite_solenoid_surface_current_gate(good_summary())
    assert result["status"] == "ok"
    assert result["metrics"]["best_level"]["label"] == "candidate_best"
    assert result["metrics"]["best_is_largest_tet_mesh"] is False
    assert json.loads(mcp_gate(json.dumps(good_summary())))["status"] == "ok"


def test_rejects_missing_negative_current_reversal():
    bad = copy.deepcopy(good_summary())
    bad["linearity"]["Bz_T"][2] = bad["linearity"]["Bz_T"][0]
    result = finite_solenoid_surface_current_gate(bad)
    assert result["status"] == "needs_attention"
    assert result["checks"]["positive_and_negative_scaling_is_linear"] is False
    assert result["checks"]["center_field_sign_tracks_current"] is False
