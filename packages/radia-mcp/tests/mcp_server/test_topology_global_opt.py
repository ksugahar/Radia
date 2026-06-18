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

from radia_mcp.topology_optimization.global_optimizers import differential_evolution


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
