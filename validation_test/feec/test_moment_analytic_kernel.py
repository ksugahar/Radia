"""Golden lock for the analytic closed-form moment kernel (moment_analytic_kernel, 2026-06-26).

The multipole-moment surface-charge kernel computes, per source face, the demag field H and field-gradient
gH at a target centroid.  Default since 2026-07-02: each face is fan-triangulated and integrated with the
CLOSED FORM (FieldGradFromChargedTriangleLocal) -- H = van Oosterom-Strackee, gH = its Mathematica-verified
symbolic gradient (the quadrupole field-gradient).  ~64x fewer kernel evals/face; EXACT for planar faces.
The 64pt (8x8) Gauss bilinear-quad path remains selectable with
rad.SolverConfig(moment_analytic_kernel=False) for cross-checks.

WIRING (2026-06-26): the analytic path is wired into the matrix-BUILD paths -- method 0 (dense
BuildCentroidFieldGrad) and method 2 (HACApK MomentSystemEntry), both delegating to the single analytic-
capable CentroidFieldGradFromFace.  The method-1 matrix-free matvec + method-2 block-Jacobi precond still
use the precomputed 64 Gauss samples (they would need the face corners cached) -- a documented follow-up.

These tests lock that (a) the analytic path is actually LIVE (changes the result vs Gauss, not a no-op),
(b) it reproduces the Gauss demag physics, (c) it produces the correct cube demag, and (d) the flag
round-trips + defaults on + does not leak.  Derivation + verification (gH rel ~ Gauss accuracy, symmetric
+ traceless to machine eps; planar quad exact): docs/multipole_moment_mmm/quadrupole_hessian_derivation.wls
+ quad_split_validation.wls.  Self-contained (mesh-less ObjHexahedron + MatLin), no NGSolve, fast.
"""
import numpy as np
import pytest
import radia as rad

MU0 = 4e-7 * np.pi
H0 = 1000.0


@pytest.fixture(autouse=True)
def _clean():
    rad.UtiDelAll(); rad.set_demag_backend("auto")
    yield
    rad.SolverConfig(moment_analytic_kernel=True)   # restore the DEFAULT (analytic, 2026-07-02 flip); never leak Gauss to other goldens
    rad.set_demag_backend("auto"); rad.UtiDelAll()


def _iron_block_extB(analytic, method, n=2, mu_r=200.0, L=0.01):
    """Solve an n x n x n iron-hex block in a uniform Hz with the chosen face kernel + method; return ext B."""
    rad.UtiDelAll(); rad.set_demag_backend("collocation_mmmm")
    rad.SolverConfig(bicgstab_tol=1e-10, moment_analytic_kernel=bool(analytic))
    objs = []
    for ix in range(n):
        for iy in range(n):
            for iz in range(n):
                x0, y0, z0 = ix * L, iy * L, iz * L
                v = [[x0, y0, z0], [x0 + L, y0, z0], [x0 + L, y0 + L, z0], [x0, y0 + L, z0],
                     [x0, y0, z0 + L], [x0 + L, y0, z0 + L], [x0 + L, y0 + L, z0 + L], [x0, y0 + L, z0 + L]]
                h = rad.ObjHexahedron(v, [0, 0, 0]); rad.MatApl(h, rad.MatLin(mu_r)); objs.append(h)
    cont = rad.ObjCnt(objs + [rad.ObjBckg(lambda p: [0.0, 0.0, MU0 * H0])])
    rad.Solve(cont, 1e-8, 3000, method)
    pts = ([0.05, 0.01, 0.01], [0.0, 0.05, 0.02], [0.012, 0.012, 0.06])
    B = np.asarray([rad.Fld(cont, "b", p) for p in pts], float)
    rad.UtiDelAll()
    return B


def test_analytic_kernel_matches_gauss_method0_dense():
    """method 0 (dense BuildCentroidFieldGrad): the analytic kernel is LIVE (changes the result vs Gauss)
    AND reproduces the Gauss demag physics.  rel>0 proves the analytic path is actually taken (guards
    against a silent no-op / false pass); rel<5e-3 proves same physics (the small residual is Gauss error,
    the analytic form being exact for these planar faces)."""
    Bg = _iron_block_extB(analytic=False, method=0)
    Ba = _iron_block_extB(analytic=True, method=0)
    rel = np.linalg.norm(Ba - Bg) / max(np.linalg.norm(Bg), 1e-30)
    assert rel > 1e-9, f"analytic flag had NO effect on method 0 (rel {rel:.2e}) -- path not wired / no-op"
    assert rel < 5e-3, f"analytic method0 external B != gauss (rel {rel:.2e})"


def test_analytic_kernel_matches_gauss_method2_hacapk():
    """method 2 (HACApK MomentSystemEntry): the analytic kernel is LIVE and reproduces the Gauss physics."""
    Bg = _iron_block_extB(analytic=False, method=2)
    Ba = _iron_block_extB(analytic=True, method=2)
    rel = np.linalg.norm(Ba - Bg) / max(np.linalg.norm(Bg), 1e-30)
    assert rel > 1e-9, f"analytic flag had NO effect on method 2 (rel {rel:.2e}) -- path not wired / no-op"
    assert rel < 5e-3, f"analytic method2 external B != gauss (rel {rel:.2e})"


def test_analytic_kernel_matches_gauss_method1_matvec():
    """method 1 (matrix-free MomentKernelMatVec6x6): the analytic kernel is LIVE and reproduces the Gauss
    physics.  This is the matvec hot path the handover flagged -- the analytic branch replaces the per-matvec
    64-sample sum with the closed form."""
    Bg = _iron_block_extB(analytic=False, method=1)
    Ba = _iron_block_extB(analytic=True, method=1)
    rel = np.linalg.norm(Ba - Bg) / max(np.linalg.norm(Bg), 1e-30)
    assert rel > 1e-9, f"analytic flag had NO effect on method 1 (rel {rel:.2e}) -- path not wired / no-op"
    assert rel < 5e-3, f"analytic method1 external B != gauss (rel {rel:.2e})"


def test_analytic_kernel_cube_demag_physical():
    """A single iron cube in uniform Hz with the analytic kernel (method 0) magnetizes with demag ~1/3 ->
    M_z ~ 3*H0, transverse ~ 0.  Confirms the analytic path produces CORRECT physics (not merely 'differs
    from Gauss')."""
    rad.UtiDelAll(); rad.set_demag_backend("collocation_mmmm")
    rad.SolverConfig(moment_analytic_kernel=True)
    L = 0.01
    v = [[-L, -L, -L], [L, -L, -L], [L, L, -L], [-L, L, -L],
         [-L, -L, L], [L, -L, L], [L, L, L], [-L, L, L]]
    h = rad.ObjHexahedron(v, [0, 0, 0]); rad.MatApl(h, rad.MatLin(1000.0))
    cont = rad.ObjCnt([h, rad.ObjBckg(lambda p: [0.0, 0.0, MU0 * H0])])
    rad.Solve(cont, 1e-6, 1000, 0)
    M = np.asarray(rad.ObjM(h)["magnetization"], float)
    rad.UtiDelAll()
    assert 2.0 * H0 < M[2] < 4.0 * H0, f"analytic cube M_z={M[2]:.1f} not ~3*H0"
    assert abs(M[0]) < 0.05 * H0 and abs(M[1]) < 0.05 * H0, f"analytic cube transverse M not ~0: {M[:2]}"


def test_analytic_kernel_config_roundtrip_and_default_on():
    """The moment_analytic_kernel flag round-trips through SolverConfig/GetSolverConfig and defaults ON
    (the 2026-07-02 deliberate flip: exact closed form + 1.5x faster method-2 H-matrix build; Gauss stays
    selectable via moment_analytic_kernel=False for cross-checks)."""
    assert rad.GetSolverConfig()["moment_analytic_kernel"] is True      # default on (2026-07-02)
    rad.SolverConfig(moment_analytic_kernel=False)
    assert rad.GetSolverConfig()["moment_analytic_kernel"] is False
    rad.SolverConfig(moment_analytic_kernel=True)
    assert rad.GetSolverConfig()["moment_analytic_kernel"] is True


def test_cross_solve_analytic_flip_respected_no_utildelall():
    """DETERMINISTIC guard for the 2026-07-02 cross-solve moment-K cache fix (commit 9120bb9d).

    method 2 caches the chi-free geometry K on radTApplication ACROSS Solve calls (optimization
    inner loops re-Solve the SAME geometry many times).  Here we Solve the SAME container TWICE
    with method 2 WITHOUT UtiDelAll in between, flipping moment_analytic_kernel between the two
    solves.  The 2nd solve MUST reflect the new kernel -- the cached K must be rebuilt.

    If the cross-solve cache validity key were missing the kernel flag (regression on part (a) of
    the fix) AND the lifecycle invalidation hooks were absent (part (b)), the 2nd solve would
    reuse the 1st solve's stale K and the two field sets would be bit-identical (rel==0).  The
    rel>1e-9 assert fails in that regression; rel<5e-3 confirms both are the same physics.

    This is the DETERMINISTIC complement to the ABA lifecycle bug (which is nondeterministic on
    LAB via heap-address reuse -- it only bit when UtiDelAll freed then an identical geometry
    rebuilt at the same address).  The 'flag baked into K must be in the key' requirement is not
    address-dependent, so this test catches a key/hook regression every run.  General rule this
    locks (bug_patterns cross-solve-cache-config-flag-key-and-lifecycle): a cross-call cache must
    key on EVERY config flag baked into the cached artifact AND invalidate at EVERY lifecycle site
    that can free the pointed-to object."""
    rad.UtiDelAll(); rad.set_demag_backend("collocation_mmmm")
    rad.SolverConfig(bicgstab_tol=1e-10)
    L = 0.01; n = 2
    objs = []
    for ix in range(n):
        for iy in range(n):
            for iz in range(n):
                x0, y0, z0 = ix * L, iy * L, iz * L
                v = [[x0, y0, z0], [x0 + L, y0, z0], [x0 + L, y0 + L, z0], [x0, y0 + L, z0],
                     [x0, y0, z0 + L], [x0 + L, y0, z0 + L], [x0 + L, y0 + L, z0 + L], [x0, y0 + L, z0 + L]]
                h = rad.ObjHexahedron(v, [0, 0, 0]); rad.MatApl(h, rad.MatLin(200.0)); objs.append(h)
    cont = rad.ObjCnt(objs + [rad.ObjBckg(lambda p: [0.0, 0.0, MU0 * H0])])
    pts = ([0.05, 0.01, 0.01], [0.0, 0.05, 0.02], [0.012, 0.012, 0.06])
    # Solve #1: GAUSS kernel, method 2 -> builds + caches the moment K (gauss) on radTApplication
    rad.SolverConfig(moment_analytic_kernel=False)
    rad.Solve(cont, 1e-8, 3000, 2)
    Bg = np.asarray([rad.Fld(cont, "b", p) for p in pts], float)
    # Solve #2: ANALYTIC kernel, SAME container, NO UtiDelAll -> the cross-solve K cache MUST rebuild
    rad.SolverConfig(moment_analytic_kernel=True)
    rad.Solve(cont, 1e-8, 3000, 2)
    Ba = np.asarray([rad.Fld(cont, "b", p) for p in pts], float)
    rel = np.linalg.norm(Ba - Bg) / max(np.linalg.norm(Bg), 1e-30)
    assert rel > 1e-9, (f"cross-solve analytic flip IGNORED (rel {rel:.2e}) -- stale moment K reused; "
                        f"kernel flag missing from cross-solve cache key / lifecycle hook")
    assert rel < 5e-3, f"cross-solve analytic vs gauss not same physics (rel {rel:.2e})"
