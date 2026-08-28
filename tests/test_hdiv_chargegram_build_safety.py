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
    for iteration in range(6):
        with ThreadPoolExecutor(max_workers=2) as executor:
            small_future = executor.submit(_point_gram, 2 * iteration + 1, 23)
            large_future = executor.submit(_point_gram, 2 * iteration + 2, 41)
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
    with pytest.raises(RuntimeError, match="beyond the image-cancellation"):
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

    # The same object must recover through a clean build. A failed rebuild has
    # no leaves, so retaining the previous sigma-active state is invalid.
    gram.build_hmatrix(eps=1.0e-12, leaf=4, eta=2.0)
    result = gram.matvec_sym(np.ones(len(points)))
    assert np.isfinite(result).all()


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
