"""Golden lock for the HACApK H-LU preconditioner (hacapk_hlu_precond, 2026-06-26).

method-2 (SolveMomentHACApK) defaults to an element block-Jacobi preconditioner whose BiCGSTAB iteration
count GROWS with mu_r / size (the non-normal demag complement; under block-Jacobi the field-null loop modes
collapse to ~1/chi).  Opt-in rad.SolverConfig(hacapk_hlu_precond=True) replaces block-Jacobi with the full
A(chi) moment H-matrix factored ONCE (cHACApK_hlu, accum_cap=0 so the non-symmetric Schur-update factors are
not over-recompressed -- the default cap=64 deferred path gives ~5.38 round-trip error on the moment system,
cap=0 ~6e-2, near-exact as a preconditioner).  The H-LU bounds the linear iters to single/low-double digits
with the IDENTICAL converged field.  See memory mmmm-preconditioner-loop-vs-factorization +
docs/multipole_moment_mmm/ACA_MOMENT_DESIGN.

These lock (a) same field as block-Jacobi (correctness), (b) the H-LU BOUNDS the iters far below block-Jacobi
(the preconditioner win, proving it is LIVE not a no-op), (c) it works through the nonlinear Picard loop
(per-element chi), (d) the flag round-trips + defaults off.  Self-contained (mesh-less ObjHexahedron), no
NGSolve, fast.
"""
import numpy as np
import pytest
import radia as rad

MU0 = 4e-7 * np.pi
H0 = 1000.0
BH = [[0.0, 0.0], [50.0, 0.30], [150.0, 0.80], [500.0, 1.30], [2000.0, 1.70],
      [10000.0, 1.95], [80000.0, 2.10], [500000.0, 2.25]]


@pytest.fixture(autouse=True)
def _clean():
    rad.UtiDelAll(); rad.set_demag_backend("auto")
    yield
    rad.SolverConfig(hacapk_hlu_precond=False)   # never leak the global flag to other goldens
    rad.set_demag_backend("auto"); rad.UtiDelAll()


def _solve_block(hlu, mat, n=3, L=0.01, H0z=H0):
    """Solve an n x n x n pure-hex iron block in a uniform Hz with method 2; return (linear_iters, M_probe)."""
    rad.UtiDelAll(); rad.set_demag_backend("yano")
    rad.SolverConfig(hacapk_hlu_precond=bool(hlu), bicgstab_tol=1e-10)
    objs = []
    for ix in range(n):
        for iy in range(n):
            for iz in range(n):
                x0, y0, z0 = ix * L, iy * L, iz * L
                v = [[x0, y0, z0], [x0 + L, y0, z0], [x0 + L, y0 + L, z0], [x0, y0 + L, z0],
                     [x0, y0, z0 + L], [x0 + L, y0, z0 + L], [x0 + L, y0 + L, z0 + L], [x0, y0 + L, z0 + L]]
                h = rad.ObjHexahedron(v, [0, 0, 0]); rad.MatApl(h, mat()); objs.append(h)
    cont = rad.ObjCnt(objs + [rad.ObjBckg(lambda p: [0.0, 0.0, MU0 * H0z])])
    rad.Solve(cont, 1e-8, 5000, 2)
    iters = rad.GetSolveStats().get("linear_iterations")
    M = np.asarray(rad.ObjM(objs[len(objs) // 2])["magnetization"], float)
    rad.UtiDelAll()
    return iters, M


def test_hlu_matches_blockjacobi_linear():
    """method 2 linear (MatLin mu_r=1000): H-LU gives the IDENTICAL field as block-Jacobi (correctness)."""
    _, Mb = _solve_block(False, lambda: rad.MatLin(1000.0))
    _, Mh = _solve_block(True, lambda: rad.MatLin(1000.0))
    rel = np.linalg.norm(Mh - Mb) / max(np.linalg.norm(Mb), 1e-30)
    assert rel < 1e-4, f"H-LU field != block-Jacobi (rel {rel:.2e})"


def test_hlu_bounds_iters():
    """The H-LU BOUNDS the linear iters far below block-Jacobi (the win; proves the precond is LIVE)."""
    itb, _ = _solve_block(False, lambda: rad.MatLin(5000.0))
    ith, _ = _solve_block(True, lambda: rad.MatLin(5000.0))
    assert ith < itb, f"H-LU iters {ith} not < block-Jacobi {itb} -- preconditioner not LIVE / no-op"
    assert ith <= 30, f"H-LU iters {ith} not bounded (expected single/low-double digits)"


def test_hlu_nonlinear_picard():
    """Nonlinear (MatSatIsoTab, per-element chi through Picard): same field + bounded iters."""
    itb, Mb = _solve_block(False, lambda: rad.MatSatIsoTab(BH), n=2, H0z=2e4)
    ith, Mh = _solve_block(True, lambda: rad.MatSatIsoTab(BH), n=2, H0z=2e4)
    rel = np.linalg.norm(Mh - Mb) / max(np.linalg.norm(Mb), 1e-30)
    assert rel < 1e-4, f"H-LU nonlinear field != block-Jacobi (rel {rel:.2e})"
    assert ith < itb, f"H-LU nonlinear iters {ith} not < block-Jacobi {itb}"


def test_hlu_config_roundtrip_and_default_off():
    """The hacapk_hlu_precond flag round-trips through SolverConfig/GetSolverConfig and defaults OFF."""
    assert rad.GetSolverConfig()["hacapk_hlu_precond"] is False
    rad.SolverConfig(hacapk_hlu_precond=True)
    assert rad.GetSolverConfig()["hacapk_hlu_precond"] is True
    rad.SolverConfig(hacapk_hlu_precond=False)
    assert rad.GetSolverConfig()["hacapk_hlu_precond"] is False
