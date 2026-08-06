"""Replay the saved TEAM 28 CoilBuilder/eddy-bubble validation evidence."""

from __future__ import annotations

import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parent.parent
RESULT = HERE / "team28_coilbuilder_eddy_bubble_results.json"
EXCHANGE = HERE / "team28_coilbuilder_hcurl_eddy_cln.json"
MATLAB_RESULT = HERE / "team28_coilbuilder_matlab_results.json"


def test_team28_coilbuilder_eddy_bubble_artifact_is_complete():
    payload = json.loads(RESULT.read_text(encoding="utf-8"))

    assert next(iter(payload)) == "radia_version"
    assert payload["schema"] == "cae-ai-lab.solver-run.v1"
    assert payload["pass"] is True
    assert payload["checks"]["validation_passed"] is True
    assert all(payload["checks"].values())
    assert payload["details"]["reduced_model"]["evrs_rank"] == 3
    assert payload["details"]["observables"]["upward_force_relative_error"] < 0.01
    assert (
        payload["details"]["coil_builder"]["field_cross_check"]["flux_density_relative_l2"] < 1.0e-3
    )
    for relative_path in payload["result_files"]:
        assert (REPO_ROOT / relative_path).is_file()


def test_team28_matlab_exchange_identifies_coilbuilder_source():
    payload = json.loads(EXCHANGE.read_text(encoding="utf-8"))

    assert payload["schema"] == "radia.hcurl.eddy_cln.exchange.v1"
    assert payload["state_order"] == 3
    assert payload["port_count"] == 1
    assert payload["metadata"]["coil_source"] == "radia.coil_builder.CoilBuilder"
    assert payload["metadata"]["winding_directions"] == "counter-wound"


def test_team28_matlab_replay_preserves_cross_language_parity():
    payload = json.loads(MATLAB_RESULT.read_text(encoding="utf-8"))

    assert next(iter(payload)) == "radia_version"
    assert payload["schema"] == "cae-ai-lab.solver-run.v1"
    assert payload["pass"] is True
    assert payload["checks"]["validation_passed"] is True
    assert all(payload["checks"]["details"].values())
    assert payload["comparisons"]["mex_vs_python_relative_error"] < 1.0e-10
    assert payload["comparisons"]["upward_force_reference_relative_error"] < 0.01
    assert payload["comparisons"]["force_over_current_squared_relative_spread"] < 1.0e-12
    assert [row["coil_current_A"] for row in payload["current_sweep"]] == [
        10,
        20,
        30,
    ]
