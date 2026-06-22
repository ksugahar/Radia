"""Golden lock for the moment-yano upgrade (Steps 3-4, 2026-06-22): the parameter-free MOMENT formula
(BuildMomentSystemCore) is the DEFAULT for pure 6-DOF hexahedral soft-iron demag, solved by the direct
dense solver, with a method dispatch in SolveGen:

  - default        : rad.SolverConfig()["yano_moment"] is True.
  - method 0 (LU)  : moment, physical (cube demag ~1/3 -> M_z ~ 3*H0).
  - method 1 (BiCG): reroutes to the dense moment LU -> bit-identical M to method 0.
  - method 2 (HACApK): the scalable moment H-matrix + block-Jacobi BiCGSTAB (Phase-2 Inc 3) -- solves the
    moment system (no longer Error204) and == method 0; LINEAR and NONLINEAR (per-element chi, Inc 4), with
    O(N log N) storage (no dense interaction/base/system matrix, Inc 4 -- see bench_moment_storage_scaling.py).
  - IMA (image=)   : moment-capable (BuildCentroidFieldGrad adds the mirror images) -> reproduces explicit full.
  - opt-out        : yano_moment=False -> EIEM2 collocation (close to moment externally).

These lock the Step-3 default flip + the Step-4 dispatch + the Phase-2 H-matrix path so a future change cannot
silently break them.  Self-contained (mesh-less ObjHexahedron + MatLin/MatSatIsoTab), no NGSolve, fast.
"""
import numpy as np
import pytest
import radia as rad

MU0 = 4e-7 * np.pi
H0 = 1000.0


@pytest.fixture(autouse=True)
def _clean():
    rad.UtiDelAll(); rad.set_demag_backend("auto"); rad.SolverConfig(yano_moment=True)
    yield
    rad.SolverConfig(yano_moment=True); rad.set_demag_backend("auto"); rad.UtiDelAll()


def _cube_Mz(method, moment, image=None, L=0.01, center=(0.0, 0.0, 0.0)):
    rad.UtiDelAll(); rad.set_demag_backend("yano"); rad.SolverConfig(yano_moment=bool(moment))
    cx, cy, cz = center
    v = [[cx - L, cy - L, cz - L], [cx + L, cy - L, cz - L], [cx + L, cy + L, cz - L], [cx - L, cy + L, cz - L],
         [cx - L, cy - L, cz + L], [cx + L, cy - L, cz + L], [cx + L, cy + L, cz + L], [cx - L, cy + L, cz + L]]
    h = rad.ObjHexahedron(v, [0, 0, 0]); rad.MatApl(h, rad.MatLin(1000.0))
    cont = rad.ObjCnt([h, rad.ObjBckg(lambda p: [0.0, 0.0, MU0 * H0])])
    if image is None:
        rad.Solve(cont, 1e-6, 1000, method)
    else:
        rad.Solve(cont, 1e-6, 1000, method, image=image)
    return np.asarray(rad.ObjM(h)["magnetization"], float)


def test_moment_is_default():
    assert rad.GetSolverConfig().get("yano_moment") is True


def test_moment_cube_physical():
    """method 0 moment: a cube in uniform Hz magnetizes with demag ~1/3 -> M_z ~ 3*H0, transverse ~ 0."""
    M = _cube_Mz(0, moment=True)
    assert 2.0 * H0 < M[2] < 4.0 * H0, f"moment cube M_z={M[2]:.1f} not ~3*H0"
    assert abs(M[0]) < 0.05 * abs(M[2]) and abs(M[1]) < 0.05 * abs(M[2])


def test_method1_bicgstab_reroutes_to_moment_lu():
    """method 1 (BiCGSTAB) is rerouted to the dense moment LU -> bit-identical to method 0."""
    M0 = _cube_Mz(0, moment=True)
    M1 = _cube_Mz(1, moment=True)
    assert np.linalg.norm(M1 - M0) <= 1e-9 * max(np.linalg.norm(M0), 1.0), f"M1={M1} != M0={M0}"


def test_method2_hacapk_solves_via_hmatrix():
    """method 2 (HACApK H-matrix + block-Jacobi BiCGSTAB, Phase-2 Increment 3) now SOLVES the moment system
    (no longer raises Error204) and == method 0 (dense LU).  Single cube: the 6x6 block-Jacobi is the exact
    local inverse so BiCGSTAB converges immediately."""
    M0 = _cube_Mz(0, moment=True)
    M2 = _cube_Mz(2, moment=True)
    assert np.all(np.isfinite(M2)) and np.linalg.norm(M2) > 1e-6
    rel = np.linalg.norm(M2 - M0) / max(np.linalg.norm(M0), 1e-30)
    assert rel < 1e-3, f"method2 H-BiCGSTAB M={M2} != method0 M={M0} (rel {rel:.2e})"


def test_method2_hacapk_multihex_external_field():
    """method 2 H-matrix BiCGSTAB == method 0 dense LU on a multi-hex block, compared by the EXTERNAL field
    (formulation/solver-tolerance-independent observable; internal M is BiCGSTAB-tol-limited).  Exercises the
    real off-diagonal H-matvec (not a single 6x6 block)."""
    MU0 = 4e-7 * np.pi; mu_r = 200.0; L = 0.01

    def solve_extB(method):
        rad.UtiDelAll(); rad.set_demag_backend("yano"); rad.SolverConfig(yano_moment=True, bicgstab_tol=1e-9)
        objs = []
        for iz in range(2):
            for ix in range(3):
                for iy in range(3):
                    x0, y0, z0 = ix * L, iy * L, iz * L
                    v = [[x0, y0, z0], [x0 + L, y0, z0], [x0 + L, y0 + L, z0], [x0, y0 + L, z0],
                         [x0, y0, z0 + L], [x0 + L, y0, z0 + L], [x0 + L, y0 + L, z0 + L], [x0, y0 + L, z0 + L]]
                    h = rad.ObjHexahedron(v, [0, 0, 0]); rad.MatApl(h, rad.MatLin(mu_r)); objs.append(h)
        cont = rad.ObjCnt(objs + [rad.ObjBckg(lambda p: [0.0, 0.0, MU0 * 1e3])])
        rad.Solve(cont, 1e-8, 3000, method)
        B = np.asarray([rad.Fld(cont, "b", p) for p in ([0.05, 0.01, 0.01], [0.0, 0.05, 0.02], [0.02, 0.02, 0.06])], float)
        rad.UtiDelAll()
        return B

    B0 = solve_extB(0)
    B2 = solve_extB(2)
    rel = np.linalg.norm(B2 - B0) / max(np.linalg.norm(B0), 1e-30)
    assert rel < 1e-5, f"method2 H-BiCGSTAB external B != method0 (rel {rel:.2e})"


def test_method2_nonlinear_matches_method0():
    """Phase-2 Increment-4: the NONLINEAR moment Picard loop solves via method 2 (HACApK H-matrix +
    PER-ELEMENT-chi block-Jacobi BiCGSTAB) to the SAME saturated state as method 0 (dense moment LU).  A
    compact block of nonlinear soft iron (MatSatIsoTab) is driven past the knee along its long (y) axis;
    both methods run the SAME Picard outer loop (chi(H) recomputed per element each iteration), differing
    only in the linear solver, so the EXTERNAL B must match tightly AND both must take the same (>1) number
    of nonlinear iterations.  (Rigorous saturation sweep: examples/vim/verify_moment_nonlinear.py; storage
    scaling: examples/vim/bench_moment_storage_scaling.py.)"""
    MU0 = 4e-7 * np.pi; L = 0.01
    BH = [[0.0, 0.0], [200.0, 0.75], [500.0, 1.30], [1200.0, 1.70],
          [4000.0, 1.95], [20000.0, 2.08], [100000.0, 2.15]]
    Msat = 2.15 / MU0

    def solve(method):
        rad.UtiDelAll(); rad.set_demag_backend("yano"); rad.SolverConfig(yano_moment=True, bicgstab_tol=1e-9)
        objs = []
        for iz in range(2):
            for ix in range(3):
                for iy in range(4):                                  # 3x4x2 block, long along y -> moderate demag -> saturates
                    x0, y0, z0 = ix * L, iy * L, iz * L
                    v = [[x0, y0, z0], [x0 + L, y0, z0], [x0 + L, y0 + L, z0], [x0, y0 + L, z0],
                         [x0, y0, z0 + L], [x0 + L, y0, z0 + L], [x0 + L, y0 + L, z0 + L], [x0, y0 + L, z0 + L]]
                    h = rad.ObjHexahedron(v, [0, 0, 0]); rad.MatApl(h, rad.MatSatIsoTab(BH)); objs.append(h)
        cont = rad.ObjCnt(objs + [rad.ObjBckg(lambda p: [0.0, MU0 * 6.0e5, 0.0])])
        nit = rad.Solve(cont, 1e-6, 500, method)
        iters = int(round(nit[-1])) if isinstance(nit, (list, tuple)) else int(nit)
        M = np.asarray([m[1] for m in rad.ObjM(rad.ObjCnt(objs))], float)
        B = np.asarray([rad.Fld(cont, "b", p) for p in ([0.05, 0.02, 0.06], [0.0, 0.04, 0.1], [0.02, 0.0, -0.03])], float)
        rad.UtiDelAll()
        return M, B, iters

    M0, B0, n0 = solve(0)
    M2, B2, n2 = solve(2)
    Mmax = np.max(np.linalg.norm(M0, axis=1))
    assert n0 > 1 and abs(n2 - n0) <= 2, f"not a genuine nonlinear iteration (iters m0={n0} m2={n2})"
    assert Mmax < Msat, f"unphysical: |M|max={Mmax:.3e} exceeds Msat={Msat:.3e}"
    relB = np.linalg.norm(B2 - B0) / max(np.linalg.norm(B0), 1e-30)
    relM = np.linalg.norm(M2 - M0) / max(np.linalg.norm(M0), 1e-30)
    assert relB < 1e-4, f"nonlinear method2 external B != method0 (rel {relB:.2e})"
    assert relM < 5e-3, f"nonlinear method2 M != method0 (rel {relM:.2e})"


# NOTE: moment-vs-EIEM2 agreement is NOT tested on a single cube -- the solved M is an INTERNAL,
# formulation-dependent quantity (CLAUDE.md "do not compare MSC internal fields"), and a single cube is
# the coarsest discretization where the two formulations legitimately differ most (~18%); they converge as
# the mesh refines (test_linear_cube_parity: moment is 0.077% from the HDiv reference on 512 hexes, vs
# EIEM2's 0.94% -- moment is the MORE accurate formula, and that accuracy is locked there, not here).


def _ima_boxes_half():
    """A 2x2 layer of hexes entirely in z>0 (disjoint from its z<0 mirror -> no boundary elements on z=0)."""
    out = []
    for ix in range(2):
        for iy in range(2):
            x0, y0 = -0.02 + ix * 0.02, -0.02 + iy * 0.02
            out.append([[x0, y0, 0.006], [x0 + 0.02, y0, 0.006], [x0 + 0.02, y0 + 0.02, 0.006], [x0, y0 + 0.02, 0.006],
                        [x0, y0, 0.026], [x0 + 0.02, y0, 0.026], [x0 + 0.02, y0 + 0.02, 0.026], [x0, y0 + 0.02, 0.026]])
    return out


def _ima_solve(boxes, Happ, image):
    rad.UtiDelAll(); rad.set_demag_backend("yano"); rad.SolverConfig(yano_moment=True)
    objs = [rad.ObjHexahedron(b, [0, 0, 0]) for b in boxes]
    for h in objs:
        rad.MatApl(h, rad.MatLin(200.0))
    cont = rad.ObjCnt(objs + [rad.ObjBckg(lambda p: [MU0 * Happ[0], MU0 * Happ[1], MU0 * Happ[2]])])
    if image is None:
        rad.Solve(cont, 1e-8, 500, 0)
    else:
        rad.Solve(cont, 1e-8, 500, 0, image=image)
    return np.asarray([m[1] for m in rad.ObjM(rad.ObjCnt(objs))], float)


@pytest.mark.parametrize("Happ,image", [([1000.0, 0.0, 0.0], "+z"),    # parallel to z=0 -> symmetric
                                        ([0.0, 0.0, 1000.0], "-z")])   # perpendicular -> antisymmetric
def test_ima_image_uses_moment(Happ, image):
    """image= (IMA) is now moment-capable: BuildCentroidFieldGrad adds the mirror images, so a HALF model
    solved with image= reproduces the EXPLICIT FULL model (half + hand-mirrored z<0 copy) to machine
    precision -- same moment formulation, IMA is just the computational shortcut for the mirror."""
    half = _ima_boxes_half()
    full = half + [[[p[0], p[1], -p[2]] for p in b] for b in half]
    M_ref = _ima_solve(full, Happ, None)[:len(half)]      # explicit full, no image, z>0 elements
    M_ima = _ima_solve(half, Happ, image)                 # half + image
    rel = np.linalg.norm(M_ima - M_ref) / max(np.linalg.norm(M_ref), 1e-30)
    assert rel < 1e-6, f"moment IMA {image} != explicit full (rel {rel:.2e})"


def test_moment_entry_reproduces_system():
    """Phase-2 Increment-1: the on-demand moment H-matrix entry (MomentSystemDenseRaw, built entry-by-entry
    via MomentSystemEntry) reproduces the moment system: (1) re-normalizing A_raw's rows == the normalized
    BuildMomentSystem A (machine precision); (2) the UN-normalized A_raw solves to the SAME magnetization
    (the row 2-norm is a diagonal scaling -> direct solve invariant -- the premise the H-LU path rests on)."""
    rad.UtiDelAll(); rad.set_demag_backend("yano"); rad.SolverConfig(yano_moment=True)
    mu_r = 50.0; chi = mu_r - 1.0; L = 0.01
    objs = []                                          # 2x2x1 grid of hexes (mutual + local entries to test)
    for ix in range(2):
        for iy in range(2):
            x0, y0 = ix * L, iy * L
            v = [[x0, y0, 0.0], [x0 + L, y0, 0.0], [x0 + L, y0 + L, 0.0], [x0, y0 + L, 0.0],
                 [x0, y0, L], [x0 + L, y0, L], [x0 + L, y0 + L, L], [x0, y0 + L, L]]
            h = rad.ObjHexahedron(v, [0, 0, 0]); rad.MatApl(h, rad.MatLin(mu_r)); objs.append(h)
    handle = rad.BuildMatrix(rad.ObjCnt(objs))
    A_norm, rhs_norm, dof = rad.BuildMomentSystem(handle, chi, 0.0, 0.0, 1.0e3)
    A_norm = np.asarray(A_norm, float); rhs_norm = np.asarray(rhs_norm, float)
    A_raw, _ = rad.MomentSystemDenseRaw(handle, chi)
    A_raw = np.asarray(A_raw, float)

    rownorm = np.linalg.norm(A_raw, axis=1)
    A_renorm = A_raw / np.where(rownorm > 1e-300, rownorm, 1.0)[:, None]
    assert np.max(np.abs(A_renorm - A_norm)) < 1e-9, "on-demand entry != BuildMomentSystem (renormalized)"

    x_norm = np.linalg.solve(A_norm, rhs_norm)
    x_raw = np.linalg.solve(A_raw, rhs_norm * rownorm)
    rel = np.linalg.norm(x_raw - x_norm) / max(np.linalg.norm(x_norm), 1e-30)
    assert rel < 1e-9, f"un-normalized A_raw solves to a different x (rel {rel:.2e})"
    rad.UtiDelAll()


def test_moment_hmatrix_matvec_equals_dense():
    """Phase-2 Increment-2: the moment system A_raw built as a HACApK H-matrix (RadHACApKMomentSystem)
    reproduces the dense A_raw matvec.  MomentHMatrixProbe builds the H-matrix + compares H-matvec(x) to
    dense A_raw @ x.  (Compression-at-scale -- n_lowrank/compression growing with N -- is exercised in
    C:/temp/verify_moment_hmatrix.py on the larger C-yoke; here we lock matvec correctness.)"""
    rad.UtiDelAll(); rad.set_demag_backend("yano"); rad.SolverConfig(yano_moment=True)
    mu_r = 200.0; chi = mu_r - 1.0; L = 0.01
    objs = []                                          # 3x3x2 grid of hexes (some off-diagonal structure)
    for iz in range(2):
        for ix in range(3):
            for iy in range(3):
                x0, y0, z0 = ix * L, iy * L, iz * L
                v = [[x0, y0, z0], [x0 + L, y0, z0], [x0 + L, y0 + L, z0], [x0, y0 + L, z0],
                     [x0, y0, z0 + L], [x0 + L, y0, z0 + L], [x0 + L, y0 + L, z0 + L], [x0, y0 + L, z0 + L]]
                h = rad.ObjHexahedron(v, [0, 0, 0]); rad.MatApl(h, rad.MatLin(mu_r)); objs.append(h)
    handle = rad.BuildMatrix(rad.ObjCnt(objs))
    d = rad.MomentHMatrixProbe(handle, chi, 1e-6, 32, 2.0)
    rad.UtiDelAll()
    if not d["ok"]:
        pytest.skip("HACApK not available in this build")
    assert d["ndof"] == 6 * len(objs)
    assert d["matvec_relerr"] < 1e-6, f"moment H-matvec != dense A_raw (relerr {d['matvec_relerr']:.2e})"


def test_moment_nonlinear_picard_matches_linear_in_linear_regime():
    """The moment LU path drives a NONLINEAR material (MatSatIsoTab) through the Picard loop, reading chi(H)
    from ctx.CurrentChiArray each step.  Locked robustly WITHOUT saturation-extrapolation tuning: at a field
    weak enough that the operating point stays on the curve's FIRST (linear) segment, MatSatIsoTab(curve)
    must give the same M as MatLin(mu_r = curve's initial slope).  This exercises the nonlinear Picard
    plumbing of the moment branch.  Saturation itself is cross-checked in examples/vim/verify_moment_nonlinear.py
    (C-yoke driven to 94% of Msat, moment vs EIEM2 external B within 2.5e-4)."""
    L, H_app = 0.01, 1.0                                   # tiny field -> internal H stays in segment 1 (H<200)
    mu_r0 = 0.75 / (MU0 * 200.0)                           # initial slope of the B-H curve = relative permeability
    bh = [[0.0, 0.0], [200.0, 0.75], [500.0, 1.30], [1200.0, 1.70], [4000.0, 1.95]]
    v = [[-L, -L, -L], [L, -L, -L], [L, L, -L], [-L, L, -L], [-L, -L, L], [L, -L, L], [L, L, L], [-L, L, L]]

    def solve(make_mat):
        rad.UtiDelAll(); rad.set_demag_backend("yano"); rad.SolverConfig(yano_moment=True)
        h = rad.ObjHexahedron(v, [0, 0, 0]); rad.MatApl(h, make_mat())   # build material AFTER UtiDelAll
        rad.Solve(rad.ObjCnt([h, rad.ObjBckg(lambda p: [0.0, 0.0, MU0 * H_app])]), 1e-8, 500, 0)
        return rad.ObjM(h)["magnetization"][2]

    mz_nl = solve(lambda: rad.MatSatIsoTab(bh))
    mz_lin = solve(lambda: rad.MatLin(mu_r0))
    assert np.isfinite(mz_nl) and mz_nl > 0
    rel = abs(mz_nl - mz_lin) / max(abs(mz_lin), 1e-30)
    assert rel < 1e-3, f"nonlinear Picard M_z={mz_nl:.4e} != linear M_z={mz_lin:.4e} in the linear regime (rel {rel:.2e})"
