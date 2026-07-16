import copy
import json

import pytest

from radia_mcp.radia_ngsolve.nonlinear_inductance_sweep_gate import (
    nonlinear_inductance_sweep_gate,
)
from radia_mcp.radia_ngsolve.server import nonlinear_inductance_sweep_gate as mcp_gate


def _summary():
    runs = []
    levels = (
        (0.5, 2.0, 2.5),
        (2.0, 3.0, 3.6),
        (10.0, 2.0, 0.8),
        (50.0, 0.8, 0.2),
    )
    for current, apparent, incremental in levels:
        for replay in (1, 2):
            matrix = [[apparent, -0.5 * apparent], [-0.5 * apparent, apparent]]
            tangent = [
                [incremental, -0.5 * incremental],
                [-0.5 * incremental, incremental],
            ]
            flux = [apparent * current, -0.5 * apparent * current]
            energy = 0.4 * current * flux[0]
            runs.append(
                {
                    "current_A_requested": current,
                    "replay": replay,
                    "apparent_inductance_H": matrix,
                    "incremental_inductance_H": tangent,
                    "current_A": [current, 0.0],
                    "flux_linkage_Vs": flux,
                    "energy_J": energy,
                    "coenergy_J": current * flux[0] - energy,
                    "final_nonlinear_residual_log10": -7.0,
                    "result_metadata": {
                        "energy": {"run_id": 0},
                        "coenergy": {"run_id": 0},
                        "residual": {"run_id": 0},
                    },
                }
            )
    return {"runs": runs}


def test_nonlinear_inductance_sweep_accepts_crossover_duality_and_replay():
    result = nonlinear_inductance_sweep_gate(_summary())
    assert result["status"] == "ok"
    assert all(result["checks"].values())
    assert result["differential_to_apparent_primary_ratios"][0] > 1.0
    assert result["differential_to_apparent_primary_ratios"][-1] < 1.0


def test_nonlinear_inductance_sweep_rejects_wrong_global_order_and_duality():
    bad = copy.deepcopy(_summary())
    for row in bad["runs"]:
        row["incremental_inductance_H"] = row["apparent_inductance_H"]
    bad["runs"][0]["coenergy_J"] *= 0.5
    result = nonlinear_inductance_sweep_gate(bad)
    assert result["status"] == "needs_attention"
    assert result["checks"]["initial_magnetization_rise_is_observed"] is False
    assert result["checks"]["differential_to_apparent_crossover_is_observed"] is False
    assert result["checks"]["all_run_identities_and_matrices_close"] is False


def test_nonlinear_inductance_sweep_rejects_asymmetric_incremental_matrix():
    bad = copy.deepcopy(_summary())
    bad["runs"][3]["incremental_inductance_H"][0][1] *= 0.5
    result = nonlinear_inductance_sweep_gate(bad)
    assert result["status"] == "needs_attention"
    assert result["checks"]["all_run_identities_and_matrices_close"] is False


def test_nonlinear_inductance_sweep_mcp_dispatches_and_rejects_bad_shape():
    result = json.loads(mcp_gate(json.dumps(_summary())))
    assert result["status"] == "ok"
    invalid = json.loads(mcp_gate('{"runs": []}'))
    assert invalid["status"] == "invalid_input"


def test_nonlinear_inductance_sweep_rejects_indefinite_tangent_matrix():
    bad = copy.deepcopy(_summary())
    bad["runs"][0]["incremental_inductance_H"][0][0] *= -1.0
    result = nonlinear_inductance_sweep_gate(bad)
    assert result["status"] == "needs_attention"
    assert result["checks"]["all_run_identities_and_matrices_close"] is False


@pytest.mark.parametrize(
    "case_id",
    ["matrix_symmetry", "requested_current", "legendre_duality", "replay_matrix", "nonlinear_residual"],
)
def test_counterfactual_curriculum90_public(case_id):
    bad = copy.deepcopy(_summary())
    if case_id == "matrix_symmetry":
        bad["runs"][0]["incremental_inductance_H"][0][1] *= 0.5
    elif case_id == "requested_current":
        bad["runs"][0]["current_A"][0] *= 0.5
    elif case_id == "legendre_duality":
        bad["runs"][0]["coenergy_J"] *= 0.5
    elif case_id == "replay_matrix":
        bad["runs"][1]["apparent_inductance_H"][0][0] *= 1.1
    else:
        bad["runs"][0]["final_nonlinear_residual_log10"] = -2.0
    assert nonlinear_inductance_sweep_gate(bad)["status"] == "needs_attention"


def test_generalization_v3s_rejects_nonzero_open_secondary_current():
    bad = copy.deepcopy(_summary())
    bad["runs"][0]["current_A"][1] = 1.0
    assert nonlinear_inductance_sweep_gate(bad)["status"] == "needs_attention"


@pytest.mark.parametrize(
    "case_id",
    ["v4_apparent_symmetry", "v4_flux_identity", "v4_negative_energy", "v4_duplicate_current_level", "v4_saturation_reversal"],
)
def test_counterfactual_curriculum90_v4_public(case_id):
    bad = copy.deepcopy(_summary())
    if case_id == "v4_apparent_symmetry":
        bad["runs"][0]["apparent_inductance_H"][0][1] *= 0.25
    elif case_id == "v4_flux_identity":
        bad["runs"][0]["flux_linkage_Vs"][0] *= 1.2
    elif case_id == "v4_negative_energy":
        bad["runs"][0]["energy_J"] = -1.0
    elif case_id == "v4_duplicate_current_level":
        bad["runs"][2]["current_A_requested"] = bad["runs"][0]["current_A_requested"]
    else:
        bad["runs"][6]["apparent_inductance_H"][0][0] *= 10.0
    assert nonlinear_inductance_sweep_gate(bad)["status"] == "needs_attention"


def test_generalization_v5_rejects_noncanonical_replay_index():
    bad = copy.deepcopy(_summary())
    bad["runs"][0]["replay"] = 3
    assert nonlinear_inductance_sweep_gate(bad)["status"] == "needs_attention"


@pytest.mark.parametrize(
    "case_id",
    ["v6_public_incremental_matrix_asymmetry", "v6_public_result_metadata_run_mismatch"],
)
def test_generalization_v6_public(case_id):
    bad = copy.deepcopy(_summary())
    if case_id == "v6_public_incremental_matrix_asymmetry":
        bad["runs"][0]["incremental_inductance_H"][0][1] *= 0.50
    else:
        bad["runs"][0]["result_metadata"]["energy"]["run_id"] = "wrong-run"
    assert nonlinear_inductance_sweep_gate(bad)["status"] == "needs_attention"


@pytest.mark.parametrize(
    "case_id",
    [
        "v7_public_operating_point_matrix_mix",
        "v7_public_coenergy_unit_shadowing",
    ],
)
def test_generalization_v7_public(case_id):
    bad = copy.deepcopy(_summary())
    for index, row in enumerate(bad["runs"]):
        operating_point_id = f"op-{index // 2}"
        row.update(
            {
                "operating_point_id": operating_point_id,
                "apparent_matrix_operating_point_id": operating_point_id,
                "incremental_matrix_operating_point_id": operating_point_id,
                "apparent_matrix_current_A": list(row["current_A"]),
                "incremental_matrix_current_A": list(row["current_A"]),
                "reported_units": {
                    "current": "A",
                    "flux_linkage": "Vs",
                    "inductance": "H",
                    "energy": "J",
                    "coenergy": "J",
                },
                "artifact_units": {
                    "current": "A",
                    "flux_linkage": "Vs",
                    "inductance": "H",
                    "energy": "J",
                    "coenergy": "J",
                },
            }
        )
    if case_id == "v7_public_operating_point_matrix_mix":
        for row in bad["runs"][:2]:
            row["incremental_matrix_operating_point_id"] = "op-1"
            row["incremental_matrix_current_A"] = [2.0, 0.0]
    else:
        bad["runs"][0]["artifact_units"]["coenergy"] = "mJ"
    assert nonlinear_inductance_sweep_gate(bad)["status"] == "needs_attention"
