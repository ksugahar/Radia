"""Golden lock for the opt-in analytic closed-form moment kernel (moment_analytic_kernel, 2026-06-26).

The multipole-moment surface-charge kernel computes, per source face, the demag field H and field-gradient
gH at a target centroid.  Default: a 64pt (8x8) Gauss bilinear-quad quadrature.  Opt-in
(rad.SolverConfig(moment_analytic_kernel=True)): each face is fan-triangulated and integrated with the
CLOSED FORM (FieldGradFromChargedTriangleLocal) -- H = van Oosterom-Strackee, gH = its Mathematica-verified
symbolic gradient (the quadrupole field-gradient).  ~64x fewer kernel evals/face; EXACT for planar faces.

WIRING (2026-06-26): the analytic path is wired into the matrix-BUILD paths -- method 0 (dense
BuildCentroidFieldGrad) and method 2 (HACApK MomentSystemEntry), both delegating to the single analytic-
capable CentroidFieldGradFromFace.  The method-1 matrix-free matvec + method-2 block-Jacobi precond still
use the precomputed 64 Gauss samples (they would need the face corners cached) -- a documented follow-up.

These tests lock that (a) the analytic path is actually LIVE (changes the result vs Gauss, not a no-op),
(b) it reproduces the Gauss demag physics, (c) it produces the correct cube demag, and (d) the flag
round-trips + defaults off + does not leak.  Derivation + verification (gH rel ~ Gauss accuracy, symmetric
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
    rad.SolverConfig(moment_analytic_kernel=False)   # never leak the global kernel flag to other goldens
    rad.set_demag_backend("auto"); rad.UtiDelAll()


def _iron_block_extB(analytic, method, n=2, mu_r=200.0, L=0.01):
    """Solve an n x n x n iron-hex block in a uniform Hz with the chosen face kernel + method; return ext B."""
    rad.UtiDelAll(); rad.set_demag_backend("yano")
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
    rad.UtiDelAll(); rad.set_demag_backend("yano")
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


def test_analytic_kernel_config_roundtrip_and_default_off():
    """The moment_analytic_kernel flag round-trips through SolverConfig/GetSolverConfig and defaults OFF."""
    assert rad.GetSolverConfig()["moment_analytic_kernel"] is False     # default off
    rad.SolverConfig(moment_analytic_kernel=True)
    assert rad.GetSolverConfig()["moment_analytic_kernel"] is True
    rad.SolverConfig(moment_analytic_kernel=False)
    assert rad.GetSolverConfig()["moment_analytic_kernel"] is False
