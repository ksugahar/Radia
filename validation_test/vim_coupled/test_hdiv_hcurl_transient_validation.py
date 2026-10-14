from __future__ import annotations

import json

import numpy as np

from radia.vim import (
    AssembleHybridVIM,
    CoupleHybridVIMWithHDivMMM,
    MagnetizationBasis,
    SampledCurrentBasis,
    save_transient_artifact,
    solve_hdiv_hcurl_transient,
)


def _validation_system(step: int):
    current = SampledCurrentBasis(
        points=np.asarray(
            [[0.10, 0.0, 0.0], [0.20, 0.0, 0.0], [0.10, 0.02, 0.0]],
        ),
        weights=np.asarray([1.0e-3, 1.0e-3, 1.0e-3]),
        modes=np.asarray(
            [
                [[1.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 0.0, 0.0]],
                [[0.0, 1.0, 0.0], [0.0, 1.0, 0.0], [0.0, 1.0, 0.0]],
            ]
        ),
        kind="volume",
        names=("volume_x", "volume_y"),
    )
    magnetization = MagnetizationBasis(
        points=np.asarray([[0.0, 0.0, 0.0], [0.0, 0.02, 0.0]]),
        weights=np.asarray([1.0e-3, 1.0e-3]),
        magnetization_modes=np.asarray(
            [
                [[0.0, 0.0, 1.0], [0.0, 0.0, 1.0]],
                [[1.0, 0.0, 0.0], [1.0, 0.0, 0.0]],
            ]
        ),
        names=("magnet_z", "magnet_x"),
    )
    eddy = AssembleHybridVIM(current, sigma=1.0e6, kernel_epsilon=1.0e-3)
    stiffness = np.asarray(
        [[2.0 + 0.05 * step, 0.1], [0.1, 1.0 + 0.03 * step]],
        dtype=float,
    )
    return CoupleHybridVIMWithHDivMMM(
        magnetization,
        eddy,
        (current,),
        magnetic_operator=stiffness,
        mu=1.0,
        kernel_epsilon=1.0e-3,
    )


def test_eleven_step_moving_transient_production_contract(tmp_path):
    times = np.linspace(0.0, 1.0e-2, 11)

    def provider(step, _time_s, _previous):
        return _validation_system(int(step))

    result = solve_hdiv_hcurl_transient(
        provider,
        times,
        magnetic_rhs=lambda step, time_s, _previous: np.asarray(
            [1.0 + 0.1 * time_s, 0.2],
            dtype=float,
        ),
        mu=1.0,
    )

    assert result["n_steps"] == 10
    assert result["n_snapshots"] == 11
    assert result["all_steps_converged"]
    assert result["all_energy_steps_balanced"]
    assert result["contract"]["enforce_energy_balance"]
    assert result["max_energy_balance_mixed_norm"] < 1.0
    assert all(
        state["joule_loss_w"] >= -1.0e-12
        and state["backward_euler_dissipation_w"] >= -1.0e-12
        and np.isfinite(state["energy_balance_residual_w"])
        for state in result["states"]
    )

    artifact_path = tmp_path / "moving_transient_validation.json"
    save_transient_artifact(
        result,
        artifact_path,
        metadata={
            "case_id": "manufactured_moving_transient_validation",
            "verification_class": "validation_test",
        },
    )
    artifact = json.loads(artifact_path.read_text(encoding="utf-8"))
    assert artifact["metadata"]["verification_class"] == "validation_test"
    assert artifact["states"][-1]["step"] == 10
