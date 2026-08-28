from concurrent.futures import ThreadPoolExecutor

import numpy as np
import pytest

from radia import _radia_pybind as _rb
from radia.vim import _solve as _solve_module


def _point_gram(seed, count):
    rng = np.random.default_rng(seed)
    points = rng.uniform(-1.0, 1.0, size=(count, 3))
    weights = rng.uniform(0.05, 0.2, size=count)
    self_energy = rng.uniform(0.5, 1.0, size=count)
    return _rb._ChargeGramHMatrix(
        points.ravel(), weights, self_energy, 1e-8, 8, 2.0
    )


def test_concurrent_chargegram_builds_keep_callback_state_isolated():
    # The pybind constructor releases the GIL.  Different sizes make a shared
    # callback-manager race fail loudly instead of accidentally reading a
    # same-shaped peer object.
    with ThreadPoolExecutor(max_workers=2) as executor:
        small_future = executor.submit(_point_gram, 1, 23)
        large_future = executor.submit(_point_gram, 2, 41)
        small = small_future.result()
        large = large_future.result()

    assert small.ndof() == 23
    assert large.ndof() == 41
    assert np.isfinite(small.matvec_sym(np.ones(23))).all()
    assert np.isfinite(large.matvec_sym(np.ones(41))).all()


def test_image_folded_negative_diagonal_is_rejected():
    # Two repeated antisymmetric images are deliberately not a valid group
    # projection.  They make the folded self-energy negative and exercise the
    # low-level safety gate independently of high-level image validation.
    x = 0.01
    tet = np.array(
        [[x, 0.0, 0.0], [1.0 + x, 0.0, 0.0],
         [x, 1.0, 0.0], [x, 0.0, 1.0]],
        dtype=np.float64,
    )
    # Assert the CONTRACT (a genuinely negative folded diagonal is refused),
    # not the wording: the message now quotes the roundoff band it fell outside.
    with pytest.raises(RuntimeError, match="broken entry oracle"):
        _rb._ChargeGramHMatrix(
            tet.ravel(), np.empty(0, dtype=np.float64), 1,
            1e-8, 8, 2.0, 1e30,
            np.array([1, 1], dtype=np.int32),
            np.array([-1.0, -1.0], dtype=np.float64),
            True, 0,
        )


def test_antisymmetric_fixed_plane_roundoff_remains_buildable():
    # A face fixed by x reflection has an exactly annihilated charge in exact
    # arithmetic.  The analytic quadrature leaves a tiny positive roundoff,
    # which remains a valid (very small) normalized diagonal.
    face = np.array(
        [[0.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
        dtype=np.float64,
    )
    gram = _rb._ChargeGramHMatrix(
        np.empty(0, dtype=np.float64), face.ravel(), 0,
        1e-8, 8, 2.0, 1e30,
        np.array([1], dtype=np.int32), np.array([-1.0], dtype=np.float64),
        True, 0,
    )
    assert gram.entry(0, 0) >= 0.0


def test_antisymmetric_fixed_plane_negative_roundoff_remains_buildable():
    # Same annihilated-charge situation as the test above, but the analytic
    # quadrature happens to land on the NEGATIVE side of zero.  The sign of the
    # residue is geometry dependent, so an exact-zero acceptance test rejects
    # every mesh whose on-plane cancellation rounds down -- measured on the
    # ESRF example-5 quarter models, that was all five (TET and HEX alike, at
    # d ~ -5e-7 against an O(1) diagonal scale).  The accept band is relative
    # to the largest diagonal actually present, so an O(1) companion charge has
    # to be in the operator for the band to mean anything.
    on_plane = np.array(
        [[0.0, -0.7053, 0.9186], [0.0, -0.4457, -0.1660],
         [0.0, 0.9758, -0.6317]],
        dtype=np.float64,
    )
    off_plane = np.array(
        [[2.1834, -0.6511, 0.1998], [1.0668, 0.8241, -0.9299],
         [1.7891, 0.3672, 0.5540]],
        dtype=np.float64,
    )
    faces = np.vstack([on_plane, off_plane]).ravel()
    gram = _rb._ChargeGramHMatrix(
        np.empty(0, dtype=np.float64), faces, 0,
        1e-8, 8, 2.0, 1e30,
        np.array([1], dtype=np.int32), np.array([-1.0], dtype=np.float64),
        True, 0,
    )
    annihilated = gram.entry(0, 0)
    companion = gram.entry(1, 1)
    assert companion > 1.0e-3
    assert abs(annihilated) <= 1.0e-12 * companion


def test_fill_exception_restores_chargegram_and_global_hacapk_state(monkeypatch):
    points = np.array(
        [[0.0, 0.0, 0.0], [0.5, 0.1, 0.0], [1.0, 0.0, 0.2],
         [0.2, 0.7, 0.1], [0.8, 0.8, 0.0]],
        dtype=np.float64,
    )
    weights = np.linspace(0.1, 0.3, len(points), dtype=np.float64)
    gram = _rb._ChargeGramHMatrix.from_sampled_laplace(
        points.ravel(), weights, 1.0e-3, 1.0e-12, 4, 2.0, False
    )
    expected_entry = weights[0] ** 2 / (4.0 * np.pi * 1.0e-3)

    monkeypatch.setenv("RADIA_HDIV_TEST_FAIL_FILL_AFTER", "0")
    with pytest.raises(RuntimeError, match="injected ChargeGram fill failure"):
        gram.build_hmatrix(eps=1.0e-12, leaf=4, eta=2.0)
    monkeypatch.delenv("RADIA_HDIV_TEST_FAIL_FILL_AFTER")

    # The failed normalized fill must not leak Ghat through the physical entry
    # oracle, and the symmetric-fill global must not poison the next PEEC build.
    assert gram.entry(0, 0) == pytest.approx(expected_entry, rel=2.0e-15)
    assert _rb._TestPEECHACApKSanity(16) < 1.0e-6


def test_nonlinear_timing_collector_keeps_latest_solver_outcome():
    _solve_module._clear_cpp_solve_timings()
    _solve_module._capture_cpp_solve_timings({
        "timings": {
            "solve_total_s": 1.25,
            "last_solve_converged": 0.0,
            "last_solve_final_relative_residual": 2e-3,
        }
    })
    _solve_module._capture_cpp_solve_timings({
        "timings": {
            "solve_total_s": 0.75,
            "last_solve_converged": 1.0,
            "last_solve_final_relative_residual": 4e-9,
        }
    })

    assert _solve_module._LAST_CPP_SOLVE_TIMINGS == {
        "solve_total_s": 2.0,
        "last_solve_converged": 1.0,
        "last_solve_final_relative_residual": 4e-9,
    }
