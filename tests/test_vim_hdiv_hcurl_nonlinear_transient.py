from __future__ import annotations

import numpy as np
import pytest

from radia.vim import (
    AssembleHybridVIM,
    CoupleHybridVIMWithHDivMMM,
    MagnetizationBasis,
    MU0,
    SampledCurrentBasis,
    solve_hdiv_hcurl_nonlinear_transient,
    solve_hdiv_hcurl_transient,
)


def _system():
    current = SampledCurrentBasis(
        points=np.asarray([[0.1, 0.0, 0.0], [0.2, 0.0, 0.0]]),
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


def _linear_bh(mu_r: float = 2.0):
    return np.asarray(
        [
            [0.0, 0.0],
            [1.0, MU0 * mu_r],
            [10.0, 10.0 * MU0 * mu_r],
        ]
    )


def _rectangular_coupling_system():
    current = SampledCurrentBasis(
        points=np.asarray([[0.1, 0.0, 0.0], [0.2, 0.0, 0.0]]),
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


def test_nonlinear_driver_accepts_rectangular_magnetic_eddy_coupling():
    result = solve_hdiv_hcurl_nonlinear_transient(
        _rectangular_coupling_system(),
        [0.0, 1.0e-3],
        bh_curve=_linear_bh(),
        demag_operator=np.eye(2),
        magnetic_rhs=np.asarray([1.0, 0.0]),
        mu=1.0,
        residual_tolerance=1.0e-10,
    )

    assert result["all_steps_converged"]
    assert result["final_magnetization"].shape == (2,)
    assert result["final_eddy"].shape == (1,)


def test_linear_bh_matches_existing_transient_driver():
    system = _system()
    times = [0.0, 1.0e-3, 2.0e-3]
    demag = np.asarray([[1.999]])
    linear = solve_hdiv_hcurl_transient(
        system,
        times,
        magnetic_operator=np.asarray([[2.0]]),
        magnetic_rhs=np.asarray([1.0]),
        mu=1.0,
    )
    nonlinear = solve_hdiv_hcurl_nonlinear_transient(
        system,
        times,
        bh_curve=_linear_bh(),
        demag_operator=demag,
        magnetic_rhs=np.asarray([1.0]),
        mu=1.0,
        residual_tolerance=1.0e-10,
    )

    assert nonlinear["all_steps_converged"]
    assert nonlinear["all_energy_steps_balanced"]
    assert nonlinear["final_magnetization"] == pytest.approx(
        linear["final_magnetization"], rel=1.0e-8, abs=1.0e-12
    )
    assert nonlinear["final_eddy"] == pytest.approx(
        linear["final_eddy"], rel=1.0e-8, abs=1.0e-12
    )


def test_saturating_bh_runs_coupled_newton_and_changes_response():
    system = _system()
    times = [0.0, 1.0e-3, 2.0e-3, 3.0e-3]
    linear = solve_hdiv_hcurl_nonlinear_transient(
        system,
        times,
        bh_curve=_linear_bh(mu_r=20.0),
        demag_operator=np.asarray([[1.0]]),
        magnetic_rhs=lambda _step, time_s, _previous: np.asarray([200.0 * time_s]),
        mu=1.0,
        residual_tolerance=1.0e-9,
    )
    saturating = solve_hdiv_hcurl_nonlinear_transient(
        system,
        times,
        bh_curve=np.asarray(
            [
                [0.0, 0.0],
                [0.1, MU0 * 2.0],
                [1.0, MU0 * 2.5],
                [10.0, MU0 * 3.0],
            ]
        ),
        demag_operator=np.asarray([[1.0]]),
        magnetic_rhs=lambda _step, time_s, _previous: np.asarray([200.0 * time_s]),
        mu=1.0,
        residual_tolerance=1.0e-9,
    )

    assert saturating["all_steps_converged"]
    assert saturating["all_energy_steps_balanced"]
    assert saturating["max_nonlinear_iterations"] >= 2
    assert not np.allclose(
        saturating["final_magnetization"], linear["final_magnetization"]
    )
    assert all(state["joule_loss_w"] >= -1.0e-12 for state in saturating["states"])


def test_nonlinear_driver_requires_demag_contract():
    with pytest.raises(ValueError, match="demag_operator"):
        solve_hdiv_hcurl_nonlinear_transient(
            _system(),
            [0.0, 1.0e-3],
            bh_curve=_linear_bh(),
            magnetic_rhs=np.asarray([1.0]),
            mu=1.0,
        )


def test_nonlinear_driver_rejects_nonmonotone_bh():
    with pytest.raises(ValueError, match="strictly increasing"):
        solve_hdiv_hcurl_nonlinear_transient(
            _system(),
            [0.0, 1.0e-3],
            bh_curve=np.asarray([[0.0, 0.0], [1.0, 2.0], [2.0, 1.0]]),
            demag_operator=np.asarray([[1.0]]),
            magnetic_rhs=np.asarray([1.0]),
            mu=1.0,
        )
