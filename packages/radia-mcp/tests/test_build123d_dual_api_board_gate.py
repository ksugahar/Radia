import copy
import json

from radia_mcp.build123d.dual_api_board_gate import dual_api_perforated_board_gate
from radia_mcp.build123d.server import build123d_dual_api_perforated_board_gate as mcp_gate


def _api_row(name: str) -> dict:
    return {
        "implementation": name,
        "source_example": "circuit_board.py" if name == "builder" else "circuit_board_algebra.py",
        "source_sha256": ("a" if name == "builder" else "b") * 64,
        "source_execution_mode": "exact_source_with_display_stub",
        "shape_valid_access": "property",
        "valid": True,
        "solid_count": 1,
        "vertex_count": 190,
        "edge_count": 285,
        "face_count": 97,
        "volume": 5767.5000452165295,
        "surface_area": 5173.203492338984,
        "bbox_size": [70.0, 30.0, 3.0],
        "analytic_relative_error": 0.0,
        "self_roundtrip_relative_error": 6.4e-15,
        "step_sha256": ("c" if name == "builder" else "d") * 64,
    }


def _summary() -> dict:
    bias = 1.510238421649028e-6
    external = []
    for mode in ("noheal", "default_heal"):
        for api in ("builder", "algebra"):
            external.append(
                {
                    "implementation": api,
                    "import_mode": mode,
                    "volume_count": 1,
                    "surface_count": 97,
                    "volume_relative_error": bias,
                    "signed_volume_relative_bias": bias,
                    "surface_area_relative_error": 8.8e-16,
                    "bbox_max_absolute_error": 1.8e-14,
                }
            )
    return {
        "version": "0.10.0",
        "upstream_commit": "72cecc0f8950507e6cf4d150b0413ca4c7561d6d",
        "board_dimensions": {"length": 70.0, "width": 30.0, "height": 3.0},
        "hole_contract": {
            "radius1_boundary_half_holes": 31,
            "radius1_full_holes_after_overlap_and_large_hole_replacement": 25,
            "radius2_full_holes": 4,
        },
        "analytic_volume": 5767.5000452165295,
        "api_rows": [_api_row("builder"), _api_row("algebra")],
        "external_rows": external,
        "external_volume_bias_classification": "systematic_kernel_mass_property_bias",
        "external_volume_bias_tolerance_basis": "dual_api_dual_import_mode_consistency",
        "external_execution_mode": "python_api_headless",
        "persistent_gui_started": False,
        "external_process_exit_code": 2,
        "startup_diagnostics": [
            "Could not open file: <install>/plugins",
            "Could not open file: -commandplugindir",
        ],
        "script_error_lines": [],
        "external_result_artifact_fresh": True,
    }


def test_accepts_dual_api_dual_import_mode_crosscheck():
    result = dual_api_perforated_board_gate(_summary())
    assert result["status"] == "ok"
    assert result["metrics"]["topology"] == [190, 285, 97]
    assert json.loads(mcp_gate(json.dumps(_summary())))["status"] == "ok"


def test_rejects_topology_drift_and_unclassified_volume_bias():
    bad = copy.deepcopy(_summary())
    bad["api_rows"][1]["face_count"] = 96
    bad["external_volume_bias_classification"] = "ignored_tolerance"
    result = dual_api_perforated_board_gate(bad)
    assert result["status"] == "needs_attention"
    assert set(result["issues"]) >= {
        "api_topology_matches",
        "external_bias_classified_not_hidden",
    }


def test_rejects_source_identity_and_external_mode_asymmetry():
    bad = copy.deepcopy(_summary())
    bad["upstream_commit"] = ""
    bad["api_rows"][0]["source_execution_mode"] = "rewritten_harness"
    bad["external_rows"].pop()
    result = dual_api_perforated_board_gate(bad)
    assert result["status"] == "needs_attention"
    assert set(result["issues"]) >= {
        "upstream_commit_recorded",
        "exact_source_display_stub_mode",
        "two_external_import_modes_per_api",
    }
