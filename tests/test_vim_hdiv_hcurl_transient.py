from __future__ import annotations

import json

import numpy as np
import pytest

from radia.vim import (
    AssembleHybridVIM,
    CoupleHybridVIMWithHDivMMM,
    MagnetizationBasis,
    SampledCurrentBasis,
    save_transient_artifact,
    solve_hdiv_hcurl_transient,
)


def _system():
    current = SampledCurrentBasis(
        points=np.asarray([[0.10, 0.0, 0.0], [0.20, 0.0, 0.0]]),
        weights=np.asarray([1.0e-3, 1.0e-3]),
        modes=np.asarray([[[1.0, 0.0, 0.0], [1.0, 0.0, 0.0]]]),
        kind="volume",
        names=("volume",),
    )
    magnetization = MagnetizationBasis(
        points=np.asarray([[0.0, 0.0, 0.0]]),
        weights=np.asarray([1.0e-3]),
        magnetization_modes=np.asarray([[[0.0, 0.0, 1.0]]]),
        names=("magnet",),
    )
    eddy = AssembleHybridVIM(current, sigma=1.0e6, kernel_epsilon=1.0e-3)
    return CoupleHybridVIMWithHDivMMM(
        magnetization,
        eddy,
        (current,),
        magnetic_operator=np.asarray([[2.0]]),
        mu=1.0,
        kernel_epsilon=1.0e-3,
    )


def _two_magnet_mode_system():
    current = SampledCurrentBasis(
        points=np.asarray([[0.10, 0.0, 0.0], [0.20, 0.0, 0.0]]),
        weights=np.asarray([1.0e-3, 1.0e-3]),
        modes=np.asarray([[[1.0, 0.0, 0.0], [1.0, 0.0, 0.0]]]),
        kind="volume",
        names=("volume",),
    )
    magnetization = MagnetizationBasis(
        points=np.asarray([[0.0, 0.0, 0.0]]),
        weights=np.asarray([1.0e-3]),
        magnetization_modes=np.asarray(
            [[[0.0, 0.0, 1.0]], [[1.0, 0.0, 0.0]]]
        ),
        names=("magnet_z", "magnet_x"),
    )
    eddy = AssembleHybridVIM(current, sigma=1.0e6, kernel_epsilon=1.0e-3)
    return CoupleHybridVIMWithHDivMMM(
        magnetization,
        eddy,
        (current,),
        magnetic_operator=np.eye(2),
        mu=1.0,
        kernel_epsilon=1.0e-3,
    )


def test_backward_euler_records_all_steps_and_energy_observables():
    system = _system()
    result = solve_hdiv_hcurl_transient(
        system,
        [0.0, 1.0e-3, 2.0e-3, 3.0e-3],
        magnetic_operator=np.asarray([[2.0]]),
        magnetic_rhs=lambda step, time_s, previous: np.asarray([1.0 + time_s]),
        mu=1.0,
    )

    assert result["n_steps"] == 3
    assert result["all_steps_converged"]
    assert len(result["states"]) == 3
    assert result["n_snapshots"] == 4
    assert [snapshot["step"] for snapshot in result["snapshots"]] == [0, 1, 2, 3]
    assert result["max_residual_relative_norm"] < 1.0e-10
    assert result["max_abs_energy_balance_residual_w"] < 1.0e-10
    assert result["max_energy_balance_relative_norm"] < 1.0e-10
    assert result["max_energy_balance_mixed_norm"] < 1.0
    assert result["all_energy_steps_balanced"]
    assert result["contract"]["enforce_energy_balance"]
    for state in result["states"]:
        assert state["joule_loss_w"] >= -1.0e-12
        assert state["backward_euler_dissipation_w"] >= -1.0e-12
        assert state["operator_motion_work_w"] == pytest.approx(0.0)
        assert state["energy_balance_relative_norm"] < 1.0e-10
        assert np.isfinite(state["magnetic_energy_j"])
        assert np.isfinite(state["eddy_energy_j"])
        assert np.isfinite(state["energy_balance_residual_w"])


def test_transient_accepts_a_moved_system_provider():
    calls = []

    def provider(step, time_s, previous):
        calls.append((step, time_s, tuple(previous["magnetization"])))
        return _system()

    result = solve_hdiv_hcurl_transient(
        provider,
        [0.0, 1.0e-3, 2.0e-3],
        magnetic_rhs=np.asarray([1.0]),
        mu=1.0,
    )

    assert result["n_steps"] == 2
    assert [call[0] for call in calls] == [0, 1, 2]


def test_nonzero_initial_state_counts_first_operator_motion_work():
    result = solve_hdiv_hcurl_transient(
        _system(),
        [0.0, 1.0e-3],
        magnetic_operator=lambda step, time_s, previous: np.asarray([[2.0 + step]]),
        magnetic_rhs=np.asarray([0.0]),
        initial_magnetization=np.asarray([1.0]),
        initial_eddy=np.asarray([0.0]),
        mu=1.0,
    )

    state = result["states"][0]
    assert state["operator_motion_work_w"] != pytest.approx(0.0)
    assert state["magnetic_operator_motion_work_w"] == pytest.approx(500.0)
    assert state["energy_balance_residual_w"] == pytest.approx(0.0, abs=1.0e-10)


def test_complex_surface_impedance_is_rejected_for_local_time_driver():
    with pytest.raises(ValueError, match="convolution-quadrature"):
        solve_hdiv_hcurl_transient(
            _system(),
            [0.0, 1.0e-3],
            magnetic_operator=np.asarray([[2.0]]),
            magnetic_rhs=np.asarray([1.0]),
            surface_impedance=1.0 + 2.0j,
            mu=1.0,
        )


def test_energy_inconsistent_operator_is_rejected_by_default():
    with pytest.raises(RuntimeError, match="energy balance mixed norm"):
        solve_hdiv_hcurl_transient(
            _two_magnet_mode_system(),
            [0.0, 1.0e-3],
            magnetic_operator=np.asarray([[2.0, 10.0], [-10.0, 2.0]]),
            magnetic_rhs=np.asarray([1.0, 0.0]),
            mu=1.0,
            energy_balance_absolute_tolerance=1.0e-20,
            energy_balance_relative_tolerance=1.0e-20,
        )


def test_transient_artifact_is_json_safe(tmp_path):
    result = solve_hdiv_hcurl_transient(
        _system(),
        [0.0, 1.0e-3],
        magnetic_operator=np.asarray([[2.0]]),
        magnetic_rhs=np.asarray([1.0]),
        mu=1.0,
    )
    path = tmp_path / "transient.json"
    save_transient_artifact(
        result,
        path,
        metadata={"case_id": "manufactured", "accepted": np.bool_(True)},
    )
    saved = json.loads(path.read_text(encoding="utf-8"))
    assert saved["schema"] == "cae-ai-lab.radia-vim.hdiv-hcurl-transient.v1"
    assert saved["metadata"]["case_id"] == "manufactured"
    assert saved["metadata"]["accepted"] is True
    assert saved["states"][0]["step"] == 1
