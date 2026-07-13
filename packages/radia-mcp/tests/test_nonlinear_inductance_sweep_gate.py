import copy
import json

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
