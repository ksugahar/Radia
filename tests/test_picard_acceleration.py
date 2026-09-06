"""Constrained Anderson mixing for the real-valued material Picard loops."""
import numpy as np
import pytest

from radia.picard_acceleration import (
    ConstrainedAndersonAccelerator, estimate_contraction_rate)


def _iterate(accelerator, mapping, x0, tol=1.0e-10, maxit=500):
    x = np.asarray(x0, dtype=float)
    for iteration in range(1, maxit + 1):
        g = mapping(x)
        x_next = accelerator.step(x, g)
        if np.max(np.abs(x_next - x)) < tol:
            return x_next, iteration
        x = x_next
    return x, maxit


def _slow_linear_map(n=30):
    # A contraction with three slow modes (0.985, 0.95, 0.9) among many fast ones,
    # the structure of a saturating yoke's material map: the damped Picard
    # iteration crawls on the slow modes, the Anderson window removes them.
    rng = np.random.default_rng(7)
    q, _ = np.linalg.qr(rng.standard_normal((n, n)))
    eigenvalues = np.concatenate([[0.985, 0.95, 0.9], rng.uniform(-0.3, 0.3, n - 3)])
    A = q @ np.diag(eigenvalues) @ q.T
    b = rng.standard_normal(n)
    fixed_point = np.linalg.solve(np.eye(n) - A, b)
    return (lambda x: A @ x + b), fixed_point


def test_depth_zero_is_the_damped_picard_update():
    accelerator = ConstrainedAndersonAccelerator(depth=0, relaxation=0.3)
    x = np.array([1.0, 2.0, 3.0])
    g = np.array([2.0, 0.0, 3.0])
    expected = 0.3 * g + (1.0 - 0.3) * x
    assert np.array_equal(accelerator.step(x, g), expected)
    assert accelerator.accelerated_steps == 0
    assert accelerator.residual_history == [pytest.approx(2.0)]


def test_anderson_converges_far_faster_than_damped_picard_on_a_slow_mode():
    mapping, fixed_point = _slow_linear_map()
    x0 = np.zeros(fixed_point.size)
    plain, plain_iterations = _iterate(
        ConstrainedAndersonAccelerator(depth=0, relaxation=0.5), mapping, x0, maxit=2000)
    mixed, mixed_iterations = _iterate(
        ConstrainedAndersonAccelerator(depth=3, relaxation=0.5), mapping, x0, maxit=2000)
    assert np.allclose(mixed, fixed_point, atol=1.0e-7)
    # The lab's complex ESIM accelerator takes exactly the same count on this map
    # (checked when this module was written); the sliding window is not GMRES.
    assert mixed_iterations * 4 < plain_iterations
    assert np.max(np.abs(plain - fixed_point)) > np.max(np.abs(mixed - fixed_point))


def test_iterates_are_projected_onto_the_admissible_interval():
    accelerator = ConstrainedAndersonAccelerator(depth=2, relaxation=1.0, lower=1.0, upper=50.0)
    x = np.array([10.0, 20.0, 30.0])
    for target in ([0.5, 100.0, 30.0], [0.2, 200.0, 31.0], [0.1, 300.0, 32.0]):
        x = accelerator.step(x, np.asarray(target, dtype=float))
        assert np.all(x >= 1.0)
        assert np.all(x <= 50.0)
        assert np.all(np.isreal(x))
    assert accelerator.projections >= 1


def test_residual_growth_restarts_the_history():
    accelerator = ConstrainedAndersonAccelerator(depth=3, relaxation=0.5, restart_growth=2.0)
    x = np.zeros(2)
    accelerator.step(x, np.array([1.0, 1.0]))
    accelerator.step(x, np.array([1.0, 1.0]))
    assert accelerator.restarts == 0
    accelerator.step(x, np.array([10.0, 10.0]))
    assert accelerator.restarts == 1
    assert accelerator.stats()["restarts"] == 1


def test_non_finite_input_fails_loud():
    accelerator = ConstrainedAndersonAccelerator(depth=1)
    with pytest.raises(ValueError, match="non-finite"):
        accelerator.step(np.array([1.0, np.nan]), np.array([1.0, 1.0]))
    with pytest.raises(ValueError, match="shape"):
        accelerator.step(np.array([1.0, 1.0]), np.array([1.0]))


def test_constructor_rejects_bad_parameters():
    with pytest.raises(ValueError):
        ConstrainedAndersonAccelerator(depth=-1)
    with pytest.raises(ValueError):
        ConstrainedAndersonAccelerator(relaxation=0.0)
    with pytest.raises(ValueError):
        ConstrainedAndersonAccelerator(restart_growth=1.0)
    with pytest.raises(ValueError):
        ConstrainedAndersonAccelerator(lower=2.0, upper=1.0)


def test_contraction_rate_estimate_recovers_a_geometric_sequence():
    history = [0.5 ** k for k in range(12)]
    assert estimate_contraction_rate(history) == pytest.approx(0.5)
    assert estimate_contraction_rate([1.0, 0.0]) is None
    assert estimate_contraction_rate([]) is None
