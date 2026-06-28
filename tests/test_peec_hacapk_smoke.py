"""PEEC x HACApK end-to-end smoke / golden lock.

The PEEC HACApK path (RadHACApKPEECManager -> HACApKPEECManager pybind ->
peec_hacapk_solver.PEECHACApKSolver) is a LIVE but otherwise UNTESTED route
(2026-06-26 inventory: no tests/ file exercised it; only
validation_test/solver_benchmarks/bench_peec_hacapk.py does).  This file gives the
maintenance-thin path a foothold so a regression in the H-matrix PEEC
adapter cannot land silently.

No Cubit, no NGSolve -- pure filament geometry.  Fast (<1 s).

Checks:
  (1) _TestPEECHACApKSanity(n): HACApK MatVec vs dense Ruehli L, machine
      precision agreement (the H-matrix assembly + matvec correctness).
  (2) PEECHACApKSolver public path: build H-matrix L, physical self/mutual
      inductance signs + decay, complex BiCGSTAB convergence at 100 kHz.
"""
import os
import sys

import numpy as np
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "../src/radia"))

# Skip cleanly when the C++ extension is not built (CI minimal-dep collect).
_pb = pytest.importorskip("radia._radia_pybind")
if not hasattr(_pb, "HACApKPEECManager"):
    pytest.skip("radia built without HACApKPEECManager", allow_module_level=True)

from peec_hacapk_solver import PEECHACApKSolver  # noqa: E402

SIGMA_CU = 5.8e7


@pytest.mark.parametrize("n", [4, 16, 64])
def test_hacapk_matches_dense_L(n):
    """HACApK L MatVec must match the dense Ruehli L to ~machine precision."""
    err = _pb._TestPEECHACApKSanity(n)
    assert err >= 0.0, f"sanity check failed with code {err}"
    assert err < 1e-6, f"HACApK vs dense L max_rel_err={err:.3e} too large"


def _parallel_filaments(n, length=0.05, wh=0.002, spacing=0.01):
    centers = np.zeros((n, 3))
    centers[:, 0] = np.arange(n) * spacing
    directions = np.tile([0.0, 0.0, 1.0], (n, 1))
    lengths = np.full(n, length)
    widths = np.full(n, wh)
    heights = np.full(n, wh)
    sigmas = np.full(n, SIGMA_CU)
    return centers, directions, lengths, widths, heights, sigmas


def test_solver_inductance_signs_and_decay():
    """Self-inductance > 0, parallel mutual > 0 and decaying with distance."""
    n = 12
    solver = PEECHACApKSolver(*_parallel_filaments(n))
    assert solver.N == n

    e0 = np.zeros(n)
    e0[0] = 1.0
    col0 = solver.apply_L(e0)  # column 0 of L
    assert col0[0] > 0.0, f"self-inductance L00={col0[0]:.3e} must be > 0"
    assert col0[1] > 0.0, "parallel same-sense mutual must be > 0"
    assert abs(col0[1]) > abs(col0[n - 1]), "mutual L must decay with distance"
    # 5 cm x 2 mm square filament self-inductance is tens of nH.
    assert 1e-9 < col0[0] < 1e-6, f"L00={col0[0]:.3e} H outside physical band"


def test_solver_frequency_domain_solve_converges():
    """Complex BiCGSTAB converges; driven filament carries the most current."""
    n = 12
    solver = PEECHACApKSolver(*_parallel_filaments(n))
    V = np.zeros(n, dtype=np.complex128)
    V[0] = 1.0
    res = solver.solve(100.0e3, V)
    assert res.converged, f"BiCGSTAB info={res.info} res={res.final_residual:.2e}"
    assert np.all(np.isfinite(res.x))
    assert int(np.argmax(np.abs(res.x))) == 0, "driven filament must carry max |I|"


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "-v"]))
