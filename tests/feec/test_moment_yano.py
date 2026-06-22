"""Golden lock for the moment-yano upgrade (Steps 3-4, 2026-06-22): the parameter-free MOMENT formula
(BuildMomentSystemCore) is the DEFAULT for pure 6-DOF hexahedral soft-iron demag, solved by the direct
dense solver, with a method dispatch in SolveGen:

  - default        : rad.SolverConfig()["yano_moment"] is True.
  - method 0 (LU)  : moment, physical (cube demag ~1/3 -> M_z ~ 3*H0).
  - method 1 (BiCG): reroutes to the dense moment LU -> bit-identical M to method 0.
  - method 2 (HACApK): fails loud (Error204) with the EIEM2 opt-out (yano_moment=False).
  - IMA (image=)   : NOT moment-eligible -> EIEM2 path runs (finite), no Error204.
  - opt-out        : yano_moment=False -> EIEM2 collocation (close to moment externally).

These lock the Step-3 default flip + the Step-4 dispatch so a future change cannot silently break them.
Self-contained (mesh-less ObjHexahedron + MatLin), no NGSolve dependency, fast.
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


def test_method2_hacapk_fails_loud():
    """method 2 (HACApK) has no moment path yet -> raise with the EIEM2 opt-out hint (Error204)."""
    with pytest.raises(RuntimeError) as ei:
        _cube_Mz(2, moment=True)
    assert "yano_moment=False" in str(ei.value)


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
