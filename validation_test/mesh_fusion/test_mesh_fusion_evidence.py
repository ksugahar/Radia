import json
from pathlib import Path


HERE = Path(__file__).resolve().parent


def _load(phase):
    return json.loads((HERE / f"results_phase{phase}.json").read_text(encoding="utf-8"))


def _strictly_decreases(values):
    return all(right < left for left, right in zip(values, values[1:]))


def test_poisson_and_em_phase_convergence():
    phase1 = _load(1)
    phase2 = _load(2)
    phase3 = _load(3)

    phase1_errors = [row["l2_err"] for row in phase1["results"]]
    phase2_errors = [row["l2_err_total"] for row in phase2["results"]]
    phase2_jumps = [row["interface_jump"] for row in phase2["results"]]
    phase3_errors = [row["l2_err_total"] for row in phase3["results"]]
    phase3_jumps = [row["interface_jump"] for row in phase3["results"]]

    derived = {
        1: {
            "l2_error_decreases": _strictly_decreases(phase1_errors),
            "finest_l2_error_lt_1e_6": phase1_errors[-1] < 1e-6,
        },
        2: {
            "l2_error_decreases": _strictly_decreases(phase2_errors),
            "interface_jump_decreases": _strictly_decreases(phase2_jumps),
            "finest_l2_error_lt_2e_7": phase2_errors[-1] < 2e-7,
        },
        3: {
            "l2_error_decreases": _strictly_decreases(phase3_errors),
            "interface_jump_decreases": _strictly_decreases(phase3_jumps),
            "finest_l2_error_lt_2e_6": phase3_errors[-1] < 2e-6,
        },
    }

    for phase, data in ((1, phase1), (2, phase2), (3, phase3)):
        assert data["schema"] == f"radia.validation.mesh_fusion.phase{phase}.v1"
        assert data["checks"] == derived[phase]
        assert all(derived[phase].values())


def test_harmonic_mortar_fourier_orthogonality():
    data = _load(4)
    ratios = [row["off_diag_ratio"] for row in data["results"]]
    derived = {"all_off_diagonal_ratios_lt_1e_6": max(ratios) < 1e-6}

    assert data["schema"] == "radia.validation.mesh_fusion.phase4.v1"
    assert data["checks"] == derived
    assert all(derived.values())


def test_shim_yoke_field_converges():
    data = _load(5)
    rows = data["results"]
    fields = [row["field_at_center"] for row in rows]
    meshes = [row["maxh_yoke"] for row in rows]
    field_changes = [abs(right - left) for left, right in zip(fields, fields[1:])]
    relative_change = abs(fields[-1] - fields[-2]) / abs(fields[-1])
    derived = {
        "field_converges_with_refinement": _strictly_decreases(field_changes),
        "last_two_relative_change_lt_1e_3": relative_change < 1e-3,
    }

    assert data["schema"] == "radia.validation.mesh_fusion.phase5.v1"
    assert data["checks"] == derived
    assert all(derived.values())
    assert _strictly_decreases(meshes)
