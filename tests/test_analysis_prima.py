"""radia.analysis PEEC series model on the PRIMA congruence projection."""

import numpy as np
import pytest

from radia.analysis import PEECAnalysisSolver, UnifiedAnalysis


def _peec(n=12, seed=0):
    rng = np.random.default_rng(seed)
    a = rng.standard_normal((n, n))
    inductance = 1.0e-7 * (a @ a.T / n + np.eye(n))
    resistance = np.diag(1.0e-3 * (1.0 + rng.random(n)))
    return inductance, resistance


@pytest.mark.parametrize("order", [1, 3, 5])
def test_prima_series_model_is_exact_for_any_order(order):
    inductance, resistance = _peec()
    solver = PEECAnalysisSolver()
    solver.set_peec_matrices(inductance, resistance)
    solver.set_reduction_order(order)
    solver.build()
    reduced = solver._reduced
    q = reduced["Q"]
    np.testing.assert_allclose(q.T @ resistance @ q, np.eye(reduced["order"]), atol=1e-10)
    assert reduced["R_series"] == pytest.approx(np.trace(resistance), rel=1e-13)
    assert reduced["L_series"] == pytest.approx(inductance.sum(), rel=1e-13)

    frequencies = np.array([0.0, 1.0e3, 1.0e5])
    result = solver.solve_frequency(frequencies)
    expected = np.trace(resistance) + 2j * np.pi * frequencies * inductance.sum()
    np.testing.assert_allclose(result.impedance, expected, rtol=1e-13)


def test_transient_converges_to_the_exact_series_rl_step_response():
    inductance, resistance = _peec()
    r_total, l_total = np.trace(resistance), inductance.sum()
    tau = l_total / r_total
    errors = []
    for steps in (500, 1000, 2000):
        analysis = UnifiedAnalysis()
        analysis.set_peec_model(inductance, resistance, reduction_order=4)
        time = np.linspace(0.0, 5.0 * tau, steps + 1)
        result = analysis.transient(time, lambda _t: 1.0, "voltage")
        exact = (1.0 - np.exp(-time / tau)) / r_total
        errors.append(np.max(np.abs(result.current - exact)) * r_total)
    assert errors[-1] < 1.0e-3
    # backward Euler is first order
    assert errors[0] / errors[1] == pytest.approx(2.0, rel=0.05)
    assert errors[1] / errors[2] == pytest.approx(2.0, rel=0.05)


def test_build_fails_loudly_on_invalid_matrices():
    inductance, resistance = _peec()
    solver = PEECAnalysisSolver()
    skew = inductance.copy()
    skew[0, 1] += 1.0e-6
    solver.set_peec_matrices(skew, resistance)
    with pytest.raises(ValueError, match="symmetric"):
        solver.build()
    negative = resistance.copy()
    negative[0, 0] = -1.0
    solver.set_peec_matrices(inductance, negative)
    with pytest.raises(ValueError, match="positive definite"):
        solver.build()
    with pytest.raises(ValueError, match="positive integer"):
        solver.set_reduction_order(0)


def test_loop_star_sweep_requires_its_matrices():
    inductance, resistance = _peec()
    solver = PEECAnalysisSolver()
    solver.set_peec_matrices(inductance, resistance)
    with pytest.raises(ValueError, match="Loop-Star matrices not set"):
        solver.solve_frequency_loop_star(np.array([1.0e3]))
