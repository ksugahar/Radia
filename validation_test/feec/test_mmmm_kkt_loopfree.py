"""Golden lock: method-0 (dense LU) KKT/Schur loop-free collocation-MMMM solve (2026-06-30).

The collocation surface-charge (MMMM) operator has a field-null loop near-null space (cell-graph
cycles; A q = (1/chi) q).  The EXACT dense-LU solve returns A^-1 b, whose loop content is amplified
by chi -> the per-element |M| is inflated (loop-POLLUTED; field-null, so the EXTERNAL field is still
right).  rad.SolverConfig(loop_deflate=True) enforces the loop-free constraint Q^T x = 0 via a Schur
complement that reuses the EXACT A^-1 (one dense factorization of the augmented RHS [b | Q], Q =
BuildLoopBasis): S = Q^T(A^-1 Q), l = S^-1(Q^T A^-1 b), x = A^-1 b - (A^-1 Q) l.  A is NEVER modified
(penalty / deflation that modified A made the iterative method-1 solver diverge; an iterative A^-1 also
diverges on the loop near-null RHS -- the DIRECT LU is what makes this work).

These lock:
  (1) default loop_deflate is OFF (regression / backward compat),
  (2) plain method-0 LU is loop-POLLUTED (|M|max inflated),
  (3) method-0 + KKT/Schur is LOOP-FREE (|M|max drops to the physical level), for linear AND nonlinear,
  (4) the EXTERNAL field is preserved (loops ~field-null) -> the physical observable is unchanged,
  (5) method-1 (iterative) + loop_deflate FAILS LOUD (cannot resolve the loop near-null space of A).
Self-contained (mesh-less ObjHexahedron + MatLin / MatSatIsoTab), no NGSolve, fast.
"""
import numpy as np
import pytest
import radia as rad

MU0 = 4e-7 * np.pi
HEX = 0.004  # element edge length (m)
EXT = np.array([[0.03, 0, 0], [0, 0.03, 0], [0, 0, 0.03],
                [0.02, 0.02, 0.01], [-0.025, 0.01, 0.005]])  # far external probes
# soft-iron-ish saturating B-H curve (mu_r ~ 1e4 at low B, saturates ~2.2 T)
BH = [[0., 0.], [50., 0.6], [100., 1.0], [200., 1.3], [500., 1.6],
      [1000., 1.8], [3000., 1.95], [10000., 2.05], [50000., 2.2]]


@pytest.fixture(autouse=True)
def _clean():
    rad.UtiDelAll(); rad.set_demag_backend("collocation_mmmm")
    yield
    rad.SolverConfig(loop_deflate=False)
    rad.set_demag_backend("auto"); rad.UtiDelAll()


def _frame(material, N=10, bw=2):
    """A loop-heavy closed square frame of hexes (an annulus of single-layer cubes) -> many cell-graph
    loops; the field-null circulation is what loop_deflate removes."""
    es = []
    for i in range(N):
        for j in range(N):
            if (bw <= i < N - bw) and (bw <= j < N - bw):
                continue
            x0, y0, z0 = (i - N / 2) * HEX, (j - N / 2) * HEX, 0.0
            V = [[x0, y0, z0], [x0 + HEX, y0, z0], [x0 + HEX, y0 + HEX, z0], [x0, y0 + HEX, z0],
                 [x0, y0, z0 + HEX], [x0 + HEX, y0, z0 + HEX], [x0 + HEX, y0 + HEX, z0 + HEX], [x0, y0 + HEX, z0 + HEX]]
            e = rad.ObjHexahedron(V, [0., 0., 0.])
            rad.MatApl(e, rad.MatLin(1e4) if material == "lin" else rad.MatSatIsoTab(BH))
            es.append(e)
    cents = np.array([[(i - N / 2 + 0.5) * HEX, (j - N / 2 + 0.5) * HEX, 0.0]
                      for i in range(N) for j in range(N)
                      if not ((bw <= i < N - bw) and (bw <= j < N - bw))])
    return es, cents


def _solve(material, method, deflate, Hdrive=(2500., 0., 0.)):
    rad.UtiDelAll(); rad.set_demag_backend("collocation_mmmm")
    rad.SolverConfig(loop_deflate=bool(deflate))
    es, cents = _frame(material)
    cont = rad.ObjCnt(es + [rad.ObjBckg(lambda p: (MU0 * np.array(Hdrive)).tolist())])
    rad.Solve(cont, 1e-7, 4000, method)
    M = np.array(rad.Fld(cont, 'm', cents)).reshape(-1, 3)
    Bext = np.array(rad.Fld(cont, 'b', EXT)).reshape(-1, 3)
    rad.SolverConfig(loop_deflate=False); rad.UtiDelAll()
    return float(np.linalg.norm(M, axis=1).max()), Bext


def test_default_loop_deflate_off():
    assert rad.GetSolverConfig().get("loop_deflate") is False


@pytest.mark.parametrize("material", ["lin", "hinput"])
def test_method0_kkt_is_loop_free(material):
    mmax_plain, B_plain = _solve(material, 0, False)
    mmax_kkt, B_kkt = _solve(material, 0, True)
    # (2) plain LU returns the exact (loop-polluted) A^-1 b -> |M|max inflated above the physical level
    assert mmax_plain > 1.15e5, f"{material}: plain |M|max={mmax_plain:.3e} not polluted"
    # (3) KKT/Schur removes the loop circulation -> |M|max drops to the physical band (~1.05e5)
    assert 0.98e5 < mmax_kkt < 1.10e5, f"{material}: KKT |M|max={mmax_kkt:.3e} out of loop-free band"
    assert mmax_kkt < 0.92 * mmax_plain, f"{material}: KKT did not reduce |M|max ({mmax_kkt:.3e} vs {mmax_plain:.3e})"
    # (4) loops are ~field-null -> the EXTERNAL field is preserved (KKT vs plain within ~2%)
    dB = np.linalg.norm(B_kkt - B_plain) / np.linalg.norm(B_plain)
    assert dB < 2.0e-2, f"{material}: external B shifted {dB:.2e} (>2%)"


def test_method1_loop_deflate_fails_loud():
    # The iterative method-1 solver cannot resolve the loop near-null space (A^-1 Q diverges):
    # loop_deflate must fail loud rather than silently return a polluted solve.
    rad.UtiDelAll(); rad.set_demag_backend("collocation_mmmm")
    rad.SolverConfig(loop_deflate=True)
    es, _ = _frame("lin")
    cont = rad.ObjCnt(es + [rad.ObjBckg(lambda p: [0., 0., MU0 * 5e4])])
    with pytest.raises(RuntimeError):
        rad.Solve(cont, 1e-7, 4000, 1)
    rad.SolverConfig(loop_deflate=False); rad.UtiDelAll()
