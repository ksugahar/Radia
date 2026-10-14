import json

from radia_mcp.radia_ngsolve.inductance_matrix_gate import inductance_matrix_family_gate
from radia_mcp.radia_ngsolve.server import inductance_matrix_family_gate as mcp_gate


def _cases():
    return [
        {"case_id": "open", "topology_class": "open_path", "matrix_H": [[4.0, -1.0], [-1.001, 1.0]]},
        {"case_id": "closed", "topology_class": "closed_path", "matrix_H": [[9.0, -2.7], [-2.699, 1.0]]},
    ]


def test_inductance_matrix_family_accepts_reciprocal_psd_family():
    result = inductance_matrix_family_gate(
        _cases(), expected_strongest_coupling_case="closed"
    )
    assert result["status"] == "ok"
    assert result["strongest_coupling_case"] == "closed"
    assert result["maximum_reciprocity_relative_error"] < 0.002


def test_inductance_matrix_family_rejects_nonphysical_mutual_and_wrong_ranking():
    cases = _cases()
    cases[1]["matrix_H"] = [[1.0, 2.0], [2.0, 1.0]]
    result = inductance_matrix_family_gate(
        cases, expected_strongest_coupling_case="open"
    )
    assert result["status"] == "needs_attention"
    assert result["checks"]["all_case_matrices_valid"] is False
    assert result["cases"][1]["checks"]["symmetrized_matrix_positive_semidefinite"] is False


def test_inductance_matrix_family_mcp_dispatches_and_handles_bad_shape():
    result = json.loads(mcp_gate(_cases(), "closed"))
    assert result["status"] == "ok"
    bad = json.loads(mcp_gate([{"case_id": "bad", "matrix_H": [[1.0]]}], None))
    assert bad["status"] == "invalid_input"


def _turn_cases():
    return [
        {
            "case_id": "n50",
            "matrix_H": [[1.0, 2.0], [2.0, 16.0]],
            "turns": [10, 50],
            "current_A": [1.0, 1.0],
            "flux_linkage_Vs": [3.0, 18.0],
            "energy_J": 10.5,
            "coenergy_J": 10.5,
            "replay_count": 2,
            "replay_max_relative_error": 1.0e-12,
        },
        {
            "case_id": "n100",
            "matrix_H": [[1.0, 4.0], [4.0, 64.0]],
            "turns": [10, 100],
            "current_A": [1.0, 1.0],
            "flux_linkage_Vs": [5.0, 68.0],
            "energy_J": 36.5,
            "coenergy_J": 36.5,
            "replay_count": 2,
            "replay_max_relative_error": 1.0e-12,
        },
    ]


def test_inductance_matrix_family_closes_identities_replay_and_turn_scaling():
    result = inductance_matrix_family_gate(_turn_cases())
    assert result["policy"] == "inductance_matrix_family_gate_v2"
    assert result["status"] == "ok"
    assert result["checks"]["turn_scaling_within_tolerance"] is True


def test_inductance_matrix_family_rejects_wrong_flux_energy_scaling_and_replay():
    cases = _turn_cases()
    cases[0]["flux_linkage_Vs"][0] = 4.0
    cases[0]["energy_J"] = 9.0
    cases[0]["replay_count"] = 1
    cases[1]["matrix_H"][1][1] = 50.0
    result = inductance_matrix_family_gate(cases)
    assert result["status"] == "needs_attention"
    assert result["cases"][0]["checks"]["flux_linkage_matches_matrix_current"] is False
    assert result["cases"][0]["checks"]["energy_matches_quadratic_form"] is False
    assert result["cases"][0]["checks"]["independent_replay_is_stable"] is False
    assert result["checks"]["turn_scaling_within_tolerance"] is False
