from __future__ import annotations

import copy
import json
import math

from radia_mcp.cubit.server import cubit_mapped_boundary_layer_shell_gate


def _summary() -> dict:
    first = 0.03
    growth = 1.2
    cumulative = [
        sum(first * growth**index for index in range(count))
        for count in range(1, 4)
    ]
    return {
        "source_journal": "mapped_boundary_layer.jou",
        "source_sha256": "a" * 64,
        "execution_mode": "python_api_headless",
        "headless_flags": ["-nographics", "-batch"],
        "persistent_gui_started": False,
        "batch_wrapper_mode": "single_line_compile_wrapper",
        "direct_multiline_batch_rejected": True,
        "result_artifact_fresh": True,
        "process_exit_code": 3,
        "startup_diagnostics": [
            "ERROR: Could not open file: install/bin/plugins",
            "ERROR: Could not open file: -commandplugindir",
            "ERROR: Could not open file: -nojournal",
        ],
        "script_error_lines": [],
        "element_counts": {"hex": 92, "pyramid": 0, "wedge": 0, "tet": 0},
        "boundary_layer": {
            "outer_radius": 1.0,
            "first_height": first,
            "growth": growth,
            "layers": 3,
            "radial_node_levels": [0.2, *(1.0 - value for value in cumulative), 1.0],
            "radial_shell_element_counts": {
                "wall_layer_1": 20,
                "wall_layer_2": 20,
                "wall_layer_3": 20,
                "core": 32,
            },
        },
        "quality": {
            "scaled_jacobian": {"count": 92, "min": 0.80},
            "shape": {"count": 92, "min": 0.20},
        },
        "cad_volume_before_scale": math.pi,
        "analytic_volume_before_scale": math.pi,
        "cad_volume_after_scale": math.pi * 1.0e-9,
        "unit_scale": 1.0e-3,
        "coordinate_scale_max_abs_error": 0.0,
    }


def _gate(summary: dict) -> dict:
    return json.loads(cubit_mapped_boundary_layer_shell_gate(summary))


def test_accepts_nodal_shells_and_allowlisted_headless_diagnostics():
    result = _gate(_summary())
    assert result["status"] == "ok"
    assert result["metrics"]["boundary_layer_shell_counts"] == [20, 20, 20]
    assert result["launcher_classification"] == "allowlisted_startup_diagnostic_with_clean_script"
    assert max(result["metrics"]["radial_level_absolute_errors"]) < 1.0e-12


def test_rejects_missing_shell_even_when_hex_quality_is_good():
    summary = copy.deepcopy(_summary())
    levels = summary["boundary_layer"]["radial_node_levels"]
    del levels[min(range(len(levels)), key=lambda index: abs(levels[index] - 0.934))]
    summary["boundary_layer"]["radial_shell_element_counts"]["wall_layer_2"] = 0
    result = _gate(summary)
    assert result["status"] == "needs_attention"
    assert set(result["issues"]) >= {
        "every_requested_radial_interface_present",
        "every_boundary_layer_shell_occupied",
        "all_hexes_classified_by_shell_or_core",
    }


def test_rejects_unrelated_process_error_and_direct_multiline_batch():
    summary = _summary()
    summary["startup_diagnostics"] = ["ERROR: mesh generation failed"]
    summary["script_error_lines"] = ["SyntaxError: list was never closed"]
    summary["batch_wrapper_mode"] = "direct_multiline_batch"
    summary["direct_multiline_batch_rejected"] = False
    result = _gate(summary)
    assert result["status"] == "needs_attention"
    assert set(result["issues"]) >= {
        "single_line_compile_wrapper_recorded",
        "direct_multiline_batch_failure_recorded",
        "script_error_lines_empty",
        "process_exit_semantics_acceptable",
    }
