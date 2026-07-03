"""Golden lock for multipole-moment MMM (Steps 3-4 + Phase 2-3, 2026-06-22): the parameter-free MOMENT formula
(BuildMomentSystemCore) is the SOLE surface-charge soft-iron demag (hex 6-DOF + wedge/pyramid 5-DOF).  The EIEM2
collocation kernel and its old EIEM2/moment opt-out were removed (Phase 3b).  Method dispatch in SolveGen:

  - method 0 (LU)  : moment, physical (cube demag ~1/3 -> M_z ~ 3*H0).
  - method 1 (BiCG): matrix-free moment BiCGSTAB + element-block Jacobi for pure hex and mixed 5/6 DOF.
  - method 2 (HACApK): the scalable moment H-matrix + block-Jacobi BiCGSTAB (Phase-2 Inc 3) -- solves the
    moment system and == method 0; LINEAR and NONLINEAR (per-element chi, Inc 4), with O(N log N) storage
    (no dense interaction/base/system matrix, Inc 4 -- see bench_moment_storage_scaling.py).
  - IMA (image=)   : moment-capable (BuildCentroidFieldGrad adds the mirror images) -> reproduces explicit full.
  - wedge (5-DOF)  : moment 5-row axial-quad (Phase 3a); hex+wedge mixed via the dense moment path.

These lock the Step-3 default flip + the Step-4 dispatch + the Phase-2 H-matrix path + the Phase-3 wedge
extension so a future change cannot silently break them.  Self-contained (mesh-less ObjHexahedron / ObjWedge
+ MatLin/MatSatIsoTab), no NGSolve, fast.
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
    rad.set_demag_backend("auto"); rad.UtiDelAll()


def _cube_Mz(method, image=None, L=0.01, center=(0.0, 0.0, 0.0)):
    rad.UtiDelAll(); rad.set_demag_backend("collocation_mmmm")
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


def test_moment_cube_physical():
    """method 0 moment: a cube in uniform Hz magnetizes with demag ~1/3 -> M_z ~ 3*H0, transverse ~ 0."""
    M = _cube_Mz(0)
    assert 2.0 * H0 < M[2] < 4.0 * H0, f"moment cube M_z={M[2]:.1f} not ~3*H0"
    assert abs(M[0]) < 0.05 * abs(M[2]) and abs(M[1]) < 0.05 * abs(M[2])


def test_method1_bicgstab_matrix_free_moment():
    """method 1 (BiCGSTAB) solves the pure-hex moment system with matrix-free matvec + block Jacobi."""
    M0 = _cube_Mz(0)
    M1 = _cube_Mz(1)
    rel = np.linalg.norm(M1 - M0) / max(np.linalg.norm(M0), 1.0)
    assert rel <= 1e-8, f"M1={M1} != M0={M0} (rel {rel:.2e})"


def test_method1_matrix_free_multihex_external_field():
    """method 1 matrix-free matvec exercises off-diagonal moment blocks and matches method 0 externally."""
    MU0 = 4e-7 * np.pi; mu_r = 200.0; L = 0.01

    def solve_extB(method):
        rad.UtiDelAll(); rad.set_demag_backend("collocation_mmmm"); rad.SolverConfig(bicgstab_tol=1e-10)
        objs = []
        for ix in range(2):
            for iy in range(2):
                x0, y0, z0 = ix * L, iy * L, 0.0
                v = [[x0, y0, z0], [x0 + L, y0, z0], [x0 + L, y0 + L, z0], [x0, y0 + L, z0],
                     [x0, y0, z0 + L], [x0 + L, y0, z0 + L], [x0 + L, y0 + L, z0 + L], [x0, y0 + L, z0 + L]]
                h = rad.ObjHexahedron(v, [0, 0, 0]); rad.MatApl(h, rad.MatLin(mu_r)); objs.append(h)
        cont = rad.ObjCnt(objs + [rad.ObjBckg(lambda p: [0.0, 0.0, MU0 * 1e3])])
        rad.Solve(cont, 1e-8, 2000, method)
        B = np.asarray([rad.Fld(cont, "b", p) for p in ([0.04, 0.01, 0.01], [0.0, 0.04, 0.02], [0.015, 0.015, 0.05])], float)
        rad.UtiDelAll()
        return B

    B0 = solve_extB(0)
    B1 = solve_extB(1)
    rel = np.linalg.norm(B1 - B0) / max(np.linalg.norm(B0), 1e-30)
    assert rel < 1e-5, f"method1 matrix-free external B != method0 (rel {rel:.2e})"


def test_method1_parallel_hotpath_12hex_matches_method0():
    """method 1 with nHex>8 (2x2x3=12 hex) exercises the PARALLEL MomentKernelMatVec6x6 branch.

    The pure-hex method-1 matvec runs serially only for nHex<=8 and otherwise parallelizes over rows
    (rad_interaction.cpp MomentKernelMatVec6x6).  The 4-hex test above only reaches the serial path; this
    12-hex case forces the parallel hot path, which must still match method 0 (dense LU) externally.
    """
    MU0 = 4e-7 * np.pi; mu_r = 200.0; L = 0.01

    def solve_extB(method):
        rad.UtiDelAll(); rad.set_demag_backend("collocation_mmmm"); rad.SolverConfig(bicgstab_tol=1e-10)
        objs = []
        for ix in range(2):
            for iy in range(2):
                for iz in range(3):
                    x0, y0, z0 = ix * L, iy * L, iz * L
                    v = [[x0, y0, z0], [x0 + L, y0, z0], [x0 + L, y0 + L, z0], [x0, y0 + L, z0],
                         [x0, y0, z0 + L], [x0 + L, y0, z0 + L], [x0 + L, y0 + L, z0 + L], [x0, y0 + L, z0 + L]]
                    h = rad.ObjHexahedron(v, [0, 0, 0]); rad.MatApl(h, rad.MatLin(mu_r)); objs.append(h)
        assert len(objs) == 12  # > 8 -> MomentKernelMatVec6x6 takes the ParallelFor branch
        cont = rad.ObjCnt(objs + [rad.ObjBckg(lambda p: [0.0, 0.0, MU0 * 1e3])])
        rad.Solve(cont, 1e-8, 4000, method)
        B = np.asarray([rad.Fld(cont, "b", p) for p in
                        ([0.06, 0.01, 0.015], [0.0, 0.06, 0.015], [0.01, 0.01, 0.07])], float)
        rad.UtiDelAll()
        return B

    B0 = solve_extB(0)
    B1 = solve_extB(1)
    rel = np.linalg.norm(B1 - B0) / max(np.linalg.norm(B0), 1e-30)
    assert rel < 1e-5, f"method1 (12 hex, parallel hot path) external B != method0 (rel {rel:.2e})"


def test_method2_hacapk_solves_via_hmatrix():
    """method 2 (HACApK H-matrix + block-Jacobi BiCGSTAB, Phase-2 Increment 3) now SOLVES the moment system
    (no longer raises Error204) and == method 0 (dense LU).  Single cube: the 6x6 block-Jacobi is the exact
    local inverse so BiCGSTAB converges immediately."""
    M0 = _cube_Mz(0)
    M2 = _cube_Mz(2)
    assert np.all(np.isfinite(M2)) and np.linalg.norm(M2) > 1e-6
    rel = np.linalg.norm(M2 - M0) / max(np.linalg.norm(M0), 1e-30)
    assert rel < 1e-3, f"method2 H-BiCGSTAB M={M2} != method0 M={M0} (rel {rel:.2e})"


def test_method2_gmres_comparison_path():
    """The explicit moment_krylov='gmres' comparison path uses the same moment H-matvec + block Jacobi and
    matches the production BiCGSTAB path on the one-hex lock case.  Bad solver names fail loud in SolverConfig
    instead of silently switching to BiCGSTAB."""
    with pytest.raises(ValueError):
        rad.SolverConfig(moment_krylov="not-a-solver")
    try:
        rad.SolverConfig(moment_krylov="gmres", moment_gmres_restart=8)
        M0 = _cube_Mz(0)
        M2 = _cube_Mz(2)
        rel = np.linalg.norm(M2 - M0) / max(np.linalg.norm(M0), 1e-30)
        assert rel < 1e-3, f"method2 H-GMRES M={M2} != method0 M={M0} (rel {rel:.2e})"
    finally:
        rad.SolverConfig(moment_krylov="bicgstab", moment_anderson_depth=0)


def test_removed_inexact_bicgstab_knob_fails_loud():
    """The failed inexact-BiCGSTAB branch is documented in MEMORY, but removed from SolverConfig."""
    with pytest.raises(ValueError):
        rad.SolverConfig(moment_inexact_bicgstab=False)


def test_removed_two_level_preconditioner_knob_fails_loud():
    """The failed two-level coarse-space branch is documented in MEMORY, but removed from SolverConfig."""
    with pytest.raises(ValueError):
        rad.SolverConfig(moment_two_level_precond=True)


def test_method2_hacapk_multihex_external_field():
    """method 2 H-matrix BiCGSTAB == method 0 dense LU on a multi-hex block, compared by the EXTERNAL field
    (formulation/solver-tolerance-independent observable; internal M is BiCGSTAB-tol-limited).  Exercises the
    real off-diagonal H-matvec (not a single 6x6 block)."""
    MU0 = 4e-7 * np.pi; mu_r = 200.0; L = 0.01

    def solve_extB(method):
        rad.UtiDelAll(); rad.set_demag_backend("collocation_mmmm"); rad.SolverConfig(bicgstab_tol=1e-9)
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
    of nonlinear iterations.  (Rigorous saturation sweep: validation_test/feec/vim_legacy/verify_moment_nonlinear.py; storage
    scaling: validation_test/feec/vim_legacy/bench_moment_storage_scaling.py.)"""
    MU0 = 4e-7 * np.pi; L = 0.01
    BH = [[0.0, 0.0], [200.0, 0.75], [500.0, 1.30], [1200.0, 1.70],
          [4000.0, 1.95], [20000.0, 2.08], [100000.0, 2.15]]
    Msat = 2.15 / MU0

    def solve(method):
        rad.UtiDelAll(); rad.set_demag_backend("collocation_mmmm"); rad.SolverConfig(bicgstab_tol=1e-9)
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


def test_wedge_moment_matches_hex_externally():
    """Phase-3a: multipole-moment MMM extended to 5-face WEDGE elements (3 dipole + 1 monopole + 1 AXIAL quad about the
    prism axis = 5 rows).  A unit cube as ONE hex vs the SAME cube tiled by TWO triangular-prism wedges must
    give the same EXTERNAL field (same geometry + uniform applied field).  The AXIAL quad (3*(d.a)^2-|d|^2,
    a = prism axis) is REQUIRED: the naive antisymmetric dx^2-dy^2 leaves a symmetric near-null sigma mode
    that blows up (|M| ~ 1e9).  Locks method 0 (LU) and method 1 matrix-free variable-DOF moment path for wedges."""
    MU0 = 4e-7 * np.pi; L = 0.02; Happ = 1000.0
    probes = [[0.05, 0.0, 0.0], [0.0, 0.05, 0.01], [0.01, 0.01, 0.05], [-0.04, 0.02, 0.03]]
    hexv = [[0, 0, 0], [L, 0, 0], [L, L, 0], [0, L, 0], [0, 0, L], [L, 0, L], [L, L, L], [0, L, L]]
    wA = [[0, 0, 0], [L, 0, 0], [0, L, 0], [0, 0, L], [L, 0, L], [0, L, L]]      # lower-left triangle prism
    wB = [[L, 0, 0], [L, L, 0], [0, L, 0], [L, 0, L], [L, L, L], [0, L, L]]      # upper-right triangle prism

    def solve(build, method):
        rad.UtiDelAll(); rad.set_demag_backend("collocation_mmmm"); rad.SolverConfig(bicgstab_tol=1e-9)
        objs = build()
        for o in objs:
            rad.MatApl(o, rad.MatLin(1000.0))
        cont = rad.ObjCnt(objs + [rad.ObjBckg(lambda p: [0.0, 0.0, MU0 * Happ])])
        rad.Solve(cont, 1e-8, 2000, method)
        M = np.asarray([rad.ObjM(o)["magnetization"] for o in objs], float)
        B = np.asarray([rad.Fld(cont, "b", p) for p in probes], float)
        rad.UtiDelAll()
        return M, B

    _, B_hex = solve(lambda: [rad.ObjHexahedron(hexv, [0, 0, 0])], 0)
    for method in (0, 1):
        Mw, Bw = solve(lambda: [rad.ObjWedge(wA, [0, 0, 0]), rad.ObjWedge(wB, [0, 0, 0])], method)
        assert np.all(np.isfinite(Mw)) and np.max(np.abs(Mw)) < 1e6, \
            f"wedge sigma blew up (near-null mode unpinned): |M|max={np.max(np.abs(Mw)):.2e}"
        Mzw = float(np.mean(Mw[:, 2]))
        assert 2.0 * Happ < Mzw < 4.0 * Happ, f"wedge M_z={Mzw:.1f} not ~3*Happ (cube demag ~1/3)"
        assert np.all(np.abs(Mw[:, :2]) < 0.1 * abs(Mzw)), "wedge transverse M not small (symmetric mode leak)"
        relB = np.linalg.norm(Bw - B_hex) / max(np.linalg.norm(B_hex), 1e-30)
        assert relB < 0.02, f"wedge method{method} external B != hex (rel {relB:.2e})"


def test_mixed_hex_wedge_moment():
    """Phase-3a: a MIXED container (hex 6-DOF + wedge 5-DOF; same 5-DOF family as pyramid) solves via the moment dense path with variable
    DOF offsets.  A 2-cube stack with the top cube tiled by 2 wedges == both cubes as hexes, externally."""
    MU0 = 4e-7 * np.pi; L = 0.02; Happ = 1000.0
    probes = [[0.06, 0.01, 0.02], [0.0, 0.06, 0.03], [0.01, 0.01, 0.08], [-0.05, 0.02, 0.04]]

    def hx(z0):
        return [[0, 0, z0], [L, 0, z0], [L, L, z0], [0, L, z0], [0, 0, z0+L], [L, 0, z0+L], [L, L, z0+L], [0, L, z0+L]]

    def wg(z0):
        return [[[0, 0, z0], [L, 0, z0], [0, L, z0], [0, 0, z0+L], [L, 0, z0+L], [0, L, z0+L]],
                [[L, 0, z0], [L, L, z0], [0, L, z0], [L, 0, z0+L], [L, L, z0+L], [0, L, z0+L]]]

    def solve(mixed, method):
        rad.UtiDelAll(); rad.set_demag_backend("collocation_mmmm"); rad.SolverConfig(bicgstab_tol=1e-9)
        objs = [rad.ObjHexahedron(hx(0.0), [0, 0, 0])]
        objs += ([rad.ObjWedge(w, [0, 0, 0]) for w in wg(L)] if mixed else [rad.ObjHexahedron(hx(L), [0, 0, 0])])
        for o in objs:
            rad.MatApl(o, rad.MatLin(500.0))
        cont = rad.ObjCnt(objs + [rad.ObjBckg(lambda p: [0.0, 0.0, MU0 * Happ])])
        rad.Solve(cont, 1e-8, 3000, method)
        B = np.asarray([rad.Fld(cont, "b", p) for p in probes], float)
        rad.UtiDelAll()
        return B

    B_hex = solve(False, 0)
    for method in (0, 1):
        B_mix = solve(True, method)
        relB = np.linalg.norm(B_mix - B_hex) / max(np.linalg.norm(B_hex), 1e-30)
        assert np.all(np.isfinite(B_mix)) and relB < 0.01, f"mixed hex+wedge method{method} external B != all-hex (rel {relB:.2e})"


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
    rad.UtiDelAll(); rad.set_demag_backend("collocation_mmmm")
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
    rad.UtiDelAll(); rad.set_demag_backend("collocation_mmmm")
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


def test_moment_nonlinear_picard_matches_linear_in_linear_regime():
    """The moment LU path drives a NONLINEAR material (MatSatIsoTab) through the Picard loop, reading chi(H)
    from ctx.CurrentChiArray each step.  Locked robustly WITHOUT saturation-extrapolation tuning: at a field
    weak enough that the operating point stays on the curve's FIRST (linear) segment, MatSatIsoTab(curve)
    must give the same M as MatLin(mu_r = curve's initial slope).  This exercises the nonlinear Picard
    plumbing of the moment branch.  Saturation itself is cross-checked in validation_test/feec/vim_legacy/verify_moment_nonlinear.py
    (C-yoke driven to 94% of Msat, moment vs EIEM2 external B within 2.5e-4)."""
    L, H_app = 0.01, 1.0                                   # tiny field -> internal H stays in segment 1 (H<200)
    mu_r0 = 0.75 / (MU0 * 200.0)                           # initial slope of the B-H curve = relative permeability
    bh = [[0.0, 0.0], [200.0, 0.75], [500.0, 1.30], [1200.0, 1.70], [4000.0, 1.95]]
    v = [[-L, -L, -L], [L, -L, -L], [L, L, -L], [-L, L, -L], [-L, -L, L], [L, -L, L], [L, L, L], [-L, L, L]]

    def solve(make_mat):
        rad.UtiDelAll(); rad.set_demag_backend("collocation_mmmm")
        h = rad.ObjHexahedron(v, [0, 0, 0]); rad.MatApl(h, make_mat())   # build material AFTER UtiDelAll
        rad.Solve(rad.ObjCnt([h, rad.ObjBckg(lambda p: [0.0, 0.0, MU0 * H_app])]), 1e-8, 500, 0)
        return rad.ObjM(h)["magnetization"][2]

    mz_nl = solve(lambda: rad.MatSatIsoTab(bh))
    mz_lin = solve(lambda: rad.MatLin(mu_r0))
    assert np.isfinite(mz_nl) and mz_nl > 0
    rel = abs(mz_nl - mz_lin) / max(abs(mz_lin), 1e-30)
    assert rel < 1e-3, f"nonlinear Picard M_z={mz_nl:.4e} != linear M_z={mz_lin:.4e} in the linear regime (rel {rel:.2e})"
