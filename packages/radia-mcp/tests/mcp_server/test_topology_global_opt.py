"""Differential Evolution -- global derivative-free optimizer. Gated on MULTIMODAL / non-convex
functions with KNOWN global optima (Rastrigin/Ackley -> 0 at the origin, Rosenbrock -> (1,..,1)),
where local methods (Nelder-Mead/LM) stall in a local minimum. Deterministic via a seeded RNG.
"""
import math
import os
import sys

import pytest

np = pytest.importorskip("numpy")

_SRC = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", "src"))
if _SRC not in sys.path:
    sys.path.insert(0, _SRC)

from radia_mcp.topology_optimization.global_optimizers import (
    best_feasible_record,
    constraint_violation,
    differential_evolution,
)


def rastrigin(x):
    x = np.asarray(x)
    return 10 * len(x) + float(np.sum(x * x - 10 * np.cos(2 * np.pi * x)))


def ackley(x):
    x = np.asarray(x); d = len(x)
    return float(-20 * np.exp(-0.2 * np.sqrt(np.sum(x * x) / d))
                 - np.exp(np.sum(np.cos(2 * np.pi * x)) / d) + 20 + np.e)


def rosen(x):
    x = np.asarray(x)
    return float(np.sum(100 * (x[1:] - x[:-1] ** 2) ** 2 + (1 - x[:-1]) ** 2))


def test_rastrigin_2d_global():
    r = differential_evolution(rastrigin, [(-5.12, 5.12)] * 2, seed=12345, maxiter=600)
    assert r["fun"] < 1e-6
    assert np.max(np.abs(r["x"])) < 1e-3


def test_ackley_global():
    r = differential_evolution(ackley, [(-32.0, 32.0)] * 3, seed=12345, maxiter=600)
    assert r["fun"] < 1e-6
    assert np.max(np.abs(r["x"])) < 1e-3


def test_rosenbrock_global():
    r = differential_evolution(rosen, [(-5.0, 10.0)] * 4, seed=12345, maxiter=800)
    assert r["fun"] < 1e-6
    assert np.max(np.abs(r["x"] - 1.0)) < 1e-3


def test_rastrigin_5d_needs_bigger_budget():
    # 5-D Rastrigin reliably reached only with popsize~25 & ~2000 gens (the scaling caveat)
    r = differential_evolution(rastrigin, [(-5.12, 5.12)] * 5, popsize=25, F=0.6,
                               maxiter=2000, seed=1)
    assert r["fun"] < 1e-6
    assert np.max(np.abs(r["x"])) < 1e-3


def test_determinism_same_seed():
    a = differential_evolution(rastrigin, [(-5.12, 5.12)] * 2, seed=42, maxiter=300)
    b = differential_evolution(rastrigin, [(-5.12, 5.12)] * 2, seed=42, maxiter=300)
    assert np.allclose(a["x"], b["x"]) and a["fun"] == b["fun"]


def test_agreement_with_scipy_de():
    sp = pytest.importorskip("scipy.optimize")
    mine = differential_evolution(rastrigin, [(-5.12, 5.12)] * 2, seed=7, maxiter=600)
    theirs = sp.differential_evolution(rastrigin, [(-5.12, 5.12)] * 2, seed=1, tol=1e-10)
    assert mine["fun"] < 1e-6 and theirs.fun < 1e-6
    assert np.max(np.abs(mine["x"])) < 1e-3 and np.max(np.abs(theirs.x)) < 1e-3


def test_constraint_violation_norms():
    vals = [-1.0, 0.0, 2.0, 3.0]
    assert constraint_violation(vals) == pytest.approx(5.0)
    assert constraint_violation(vals, norm="l2") == pytest.approx(math.sqrt(13.0))
    assert constraint_violation(vals, norm="linf") == pytest.approx(3.0)
    assert constraint_violation([]) == pytest.approx(0.0)
    with pytest.raises(ValueError):
        constraint_violation(vals, norm="bad")


def test_best_feasible_record_prefers_feasible_objective():
    rows = [
        {"value": 0.1, "constraints": [0.2], "params": {"x": 0}},
        {"value": 3.0, "constraints": [-1.0, 0.0], "params": {"x": 1}},
        {"value": 4.0, "constraints": [-0.5], "params": {"x": 2}},
    ]
    best = best_feasible_record(rows)
    assert best["params"] == {"x": 1}
    assert best["feasible"] is True
    assert best["constraint_violation"] == pytest.approx(0.0)
    assert "_objective_sort" not in best

    best_max = best_feasible_record(rows, minimize=False)
    assert best_max["params"] == {"x": 2}


def test_best_feasible_record_falls_back_to_least_violation():
    rows = [
        {"value": 1.0, "constraints": [0.4], "params": {"x": 0}},
        {"value": 2.0, "constraints": [0.1, 0.2], "params": {"x": 1}},
        {"value": 0.5, "constraints": [0.3], "params": {"x": 2}},
    ]
    best = best_feasible_record(rows)
    assert best["params"] == {"x": 2}
    assert best["feasible"] is False
    assert best["constraint_violation"] == pytest.approx(0.3)

    with pytest.raises(ValueError):
        best_feasible_record([])
