"""Golden lock for the multipole-moment (MMMM) method-2 TWO-SIDED deflated block-Jacobi preconditioner
(rad.SolverConfig(moment_deflation=True), 2026-06-26).

MMMM's system matrix A is NON-symmetric, so its RIGHT near-null space (the surface-charge loops,
GetLoopBasis) and its LEFT near-null space (the "co-loops" = orthonorm(A@loops)) are DISTINCT (the small
singular directions are the loops at sigma~1/chi; cond(A) ~ chi).  ONE-sided loop deflation does NOT condition
the complement (the left co-loops remain) -- this is why every prior one-sided loop-deflation attempt barely
helped.  TWO-sided deflation (coarse correction zc = W (Yco^T A W)^-1 Yco^T r, W=loops, Yco=co-loops) removes
BOTH near-null spaces, conditioning the complement so the block-Jacobi smoother's iteration count collapses.

Measured (loop-heavy C-yoke): plain block-Jacobi 167 / 397 / 2917 iters (mu_r 1e3@1080DoF / 1e3@1824 /
1e4@1824) -> two-sided deflated 41 / 84 / 114 (4x-25x fewer, the win grows with mu_r), same field as dense LU.

NOTE: this is a MEDIUM-N (<=~15k) constant-factor speedup over dense LU, NOT asymptotically scalable -- the
loop space is O(N)-dimensional so the coarse factor is O(N^3) / O(N^2)-storage (a true scalable loop-heavy
solver needs a multilevel/non-symmetric-AMS coarse solve, or the symmetric HDiv-VIM).  DEFAULT ON (Sugahara
2026-06-27): the target workload is loop-heavy accelerator electromagnets (C-yokes), where deflation is a
large win; set moment_deflation=False for a non-loop-heavy solve like a benchmark cube (the O(N^3) coarse
setup then costs more than it saves).  Exclusive with the H-LU precond.

Locks: (a) the flag round-trips + defaults ON (verified in a fresh process); (b) deflated method-2 matches
method-0 dense LU AND uses fewer iters than plain block-Jacobi on a loop-heavy C-yoke.
"""
import numpy as np
import pytest
import radia as rad

MU0 = 4e-7 * np.pi


@pytest.fixture(autouse=True)
def _clean():
    rad.UtiDelAll(); rad.set_demag_backend("auto")
    yield
    rad.SolverConfig(moment_deflation=True, hacapk_hlu_precond=False)   # restore the default (ON)
    rad.set_demag_backend("auto"); rad.UtiDelAll()


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


def _solve(nxy, nz, mu_r, method, deflate):
    rad.UtiDelAll(); rad.set_demag_backend("collocation_mmmm")
    rad.SolverConfig(hacapk_hlu_precond=False, moment_deflation=bool(deflate), bicgstab_tol=1e-8)
    objs = _cyoke(nxy, nz, mu_r)
    cont = rad.ObjCnt(objs + [rad.ObjBckg(lambda p: [0.0, MU0 * 1e3, 0.0])])
    rad.Solve(cont, 1e-8, 5000, method)
    M = float(np.linalg.norm(np.array([rad.ObjM(h)["magnetization"] for h in objs], float)))
    it = rad.GetSolveStats().get("linear_iterations")
    rad.UtiDelAll()
    return M, it


def test_deflation_default_on_in_fresh_process():
    """The C++ default is moment_deflation=ON -- target workload = loop-heavy accelerator electromagnets.
    Checked in a FRESH process: the autouse fixture mutates the flag in-process, so the pristine default must
    be read before any SolverConfig call."""
    import subprocess, sys
    out = subprocess.run([sys.executable, "-c",
                          "import radia as rad; print(rad.GetSolverConfig()['moment_deflation'])"],
                         capture_output=True, text=True)
    assert out.stdout.strip() == "True", f"default not ON: stdout={out.stdout!r} stderr={out.stderr[-200:]!r}"


def test_deflation_config_roundtrip():
    rad.SolverConfig(moment_deflation=True)
    assert rad.GetSolverConfig()["moment_deflation"] is True
    rad.SolverConfig(moment_deflation=False)
    assert rad.GetSolverConfig()["moment_deflation"] is False


def test_deflation_matches_lu_and_cuts_iters():
    """Loop-heavy C-yoke: two-sided deflated method-2 matches dense LU AND uses fewer iters than block-Jacobi."""
    NXY, NZ, MU_R = 16, 2, 1000.0
    Mlu, _ = _solve(NXY, NZ, MU_R, 0, False)
    Mbj, it_bj = _solve(NXY, NZ, MU_R, 2, False)     # plain block-Jacobi
    Mdf, it_df = _solve(NXY, NZ, MU_R, 2, True)       # two-sided deflated block-Jacobi
    rel = abs(Mdf - Mlu) / max(Mlu, 1e-30)
    assert rel < 1e-2, f"deflated field != dense LU (rel {rel:.2e})"
    assert it_df < it_bj, f"deflation did not cut iters ({it_df} not < block-Jacobi {it_bj}) -- precond not LIVE"
