"""Golden lock for the multipole-moment (MMMM) method-2 H-LU safety system (2026-06-26).

method-2 (SolveMomentHACApK) H-LU preconditions a Krylov solve.  On a loop-heavy thin geometry the no-pivot
block H-LU on the NON-SYMMETRIC A(chi) can produce an UNSTABLE factor (broad-probe round-trip ~1e6-1e7, a
catastrophic amplification of the loop-charge near-null space).  A residual-based Krylov then "converges" the
H-matrix residual onto a loop-polluted solution the residual CANNOT see: BiCGSTAB diverges (true residual
~1e99) while GMRES reports tol-"converged" at a |M| 3x-47x too large (memory
mmmm-preconditioner-loop-vs-factorization #13).  Either way the old code returned a positive iteration count
the caller read as success -> a silently-wrong field (the worst failure class per "No Fallbacks: a wrong field
claiming convergence is worse than no field").

SolveMomentHACApK now has a 3-layer safety system:
 (1) H-LU FACTOR SELF-TEST -- a deterministic broad-support (pseudo-random) probe round-trip
     ||A(chi) M_H^-1 p - p||/||p|| BEFORE the Krylov loop; reject (fail loud, return -14) when it exceeds 1e3.
     This excites the loop-charge null space (which the RHS b is ~orthogonal to), so it catches the instability
     for BOTH BiCGSTAB and GMRES, instantly, before it pollutes the solve.  Calibrated: usable factors
     round-trip <=~5 (even rough high-mu_r ones), unstable ones ~1e6+ (>5 orders of magnitude apart).
 (2) EARLY-DIVERGENCE bailout -- break when the recursive/restart residual blows past 1e10*||b|| or goes
     non-finite, so any remaining divergence fails FAST instead of grinding 10000 iterations.
 (3) TRUE-RESIDUAL guard -- recompute ||b-Ax||/||b|| at the solver exit; return -13 if it exceeds tol*10, so a
     drifted recursive (BiCGSTAB) / maxiter-exhausted residual can never be reported as convergence.

These lock:
 1. healthy method-2 H-LU still converges + matches method-0 dense LU + does NOT raise (no false-fire);
 2. on a loop-heavy thin C-yoke (mu_r=1000, the realistic regime), method-2 RAISES (fail loud) -- or, if a
    future stable-H-LU fix makes it converge, matches LU -- it NEVER returns a wrong field as "converged".
    Both the BiCGSTAB and the GMRES paths are covered (the factor self-test is solver-independent).

Self-contained (mesh-less ObjHexahedron), no NGSolve.  The divergent case is a few seconds (the factor
self-test fires up front).  Lives in validation_test (the heavy tier split from CI), not tests/.
"""
import numpy as np
import pytest
import radia as rad

MU0 = 4e-7 * np.pi


@pytest.fixture(autouse=True)
def _clean():
    rad.UtiDelAll(); rad.set_demag_backend("auto")
    yield
    rad.SolverConfig(hacapk_hlu_precond=False, moment_krylov="bicgstab")  # never leak globals to other goldens
    rad.set_demag_backend("auto"); rad.UtiDelAll()


def _cube(n, L, mu_r):
    objs = []
    for ix in range(n):
        for iy in range(n):
            for iz in range(n):
                x0, y0, z0 = ix * L, iy * L, iz * L
                v = [[x0, y0, z0], [x0 + L, y0, z0], [x0 + L, y0 + L, z0], [x0, y0 + L, z0],
                     [x0, y0, z0 + L], [x0 + L, y0, z0 + L], [x0 + L, y0 + L, z0 + L], [x0, y0 + L, z0 + L]]
                h = rad.ObjHexahedron(v, [0, 0, 0]); rad.MatApl(h, rad.MatLin(mu_r)); objs.append(h)
    return objs


def _inside_cyoke(cx, cy):
    return (-0.06 <= cx <= 0.06) and (-0.06 <= cy <= 0.06) and not (-0.035 < cx < 0.035 and -0.035 < cy < 0.035) and not (cx > 0.018)


def _cyoke(nxy, nz, mu_r):
    xs = np.linspace(-0.06, 0.06, nxy + 1); zs = np.linspace(-0.02, 0.02, nz + 1)
    objs = []
    for k in range(nz):
        for j in range(nxy):
            for i in range(nxy):
                if not _inside_cyoke(0.5 * (xs[i] + xs[i + 1]), 0.5 * (xs[j] + xs[j + 1])):
                    continue
                x0, x1, y0, y1, z0, z1 = xs[i], xs[i + 1], xs[j], xs[j + 1], zs[k], zs[k + 1]
                v = [[x0, y0, z0], [x1, y0, z0], [x1, y1, z0], [x0, y1, z0],
                     [x0, y0, z1], [x1, y0, z1], [x1, y1, z1], [x0, y1, z1]]
                h = rad.ObjHexahedron(v, [0, 0, 0]); rad.MatApl(h, rad.MatLin(mu_r)); objs.append(h)
    return objs


def _solve_Mnorm(objs, bckg, method):
    cont = rad.ObjCnt(objs + [rad.ObjBckg(bckg)])
    rad.Solve(cont, 1e-8, 5000, method)
    M = np.array([rad.ObjM(h)["magnetization"] for h in objs], float)
    return float(np.linalg.norm(M))


def test_healthy_method2_matches_lu_no_raise():
    """(1) small cube, mu_r=1000: method-2 H-LU converges, matches method-0 dense LU, guard does NOT fire."""
    bckg = lambda p: [0.0, 0.0, MU0 * 1e3]
    rad.UtiDelAll(); rad.set_demag_backend("yano"); rad.SolverConfig(hacapk_hlu_precond=False, bicgstab_tol=1e-10)
    Mlu = _solve_Mnorm(_cube(3, 0.01, 1000.0), bckg, 0)
    rad.UtiDelAll(); rad.set_demag_backend("yano"); rad.SolverConfig(hacapk_hlu_precond=True, bicgstab_tol=1e-10)
    Mh = _solve_Mnorm(_cube(3, 0.01, 1000.0), bckg, 2)   # must NOT raise
    rel = abs(Mh - Mlu) / max(Mlu, 1e-30)
    assert rel < 1e-4, f"healthy method-2 H-LU field != LU (rel {rel:.2e}); guard may be false-firing"


@pytest.mark.parametrize("krylov", ["bicgstab", "gmres"])
def test_no_silent_wrong_on_loopheavy_divergence(krylov):
    """(2) loop-heavy thin C-yoke (nz=2, mu_r=1000) where the no-pivot block H-LU goes unstable.
    method-2 must RAISE (fail loud) or match LU; it must NEVER return a wrong field claiming convergence.
    The factor self-test is solver-independent, so both BiCGSTAB and GMRES fail loud here."""
    bckg = lambda p: [0.0, MU0 * 1e3, 0.0]
    NXY, NZ, MU_R = 40, 2, 1000.0   # 7296 DoF; factor round-trip ~2.7e6 >> 1e3 -> self-test fires up front
    rad.UtiDelAll(); rad.set_demag_backend("yano"); rad.SolverConfig(hacapk_hlu_precond=False, bicgstab_tol=1e-8)
    Mlu = _solve_Mnorm(_cyoke(NXY, NZ, MU_R), bckg, 0)
    rad.UtiDelAll(); rad.set_demag_backend("yano")
    rad.SolverConfig(hacapk_hlu_precond=True, bicgstab_tol=1e-8, moment_krylov=krylov)
    raised = False
    Mh = None
    try:
        Mh = _solve_Mnorm(_cyoke(NXY, NZ, MU_R), bckg, 2)
    except RuntimeError:
        raised = True   # the expected fail-loud outcome
    if not raised:
        rel = abs(Mh - Mlu) / max(Mlu, 1e-30)
        assert rel < 1e-2, (f"method-2 {krylov} returned success but |M| disagrees with LU (rel {rel:.2f}) "
                            f"-- SILENT-WRONG regression: the H-LU safety system is not firing")
