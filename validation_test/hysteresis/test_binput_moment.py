#!/usr/bin/env python
"""
Golden validation for the B-input PLAY-model hysteresis solver on MMMM moment
elements (tetrahedron 4-DOF, hexahedron 6-DOF; wedge/pyramid 5-DOF share the
same dof>=4 code path).

Background
----------
After the tet/RecMag -> MSC unification there are no 3-DOF soft-iron elements
left, so B-input hysteresis (rad.SolverConfig(b_input_newton=True)) can no
longer reach a soft-iron body through the dense 3-DOF AutoRelax_BInput_Newton
path.  This solver routes the all-moment B-input case through the normal MMMM
moment Picard loop with a B-INPUT chi update:

    B = mu0*(H + M)                         (element flux density from solved state)
    chi = |M_play| / |H_play|  with  H_play = Forward(B),  M_play = B/mu0 - H_play

Per-element play state is saved at the start of the solve and restored before
each material evaluation so each element keeps its own play trajectory across
Picard iterations; states are committed once at the end of the converged solve.

THE physics-consistency gate
----------------------------
A wrong hysteresis solver produces SILENTLY WRONG numbers.  The mandatory gate
verifies the converged state actually satisfies the play relation H = Forward(B):

  (A) For a single hexahedron (axis-aligned, so the centroid field IS the
      effective internal H), assert ||H - Forward(B)|| / ||H|| < 1e-3 directly,
      where Forward(B) is reconstructed from the COMMITTED play state via the
      exact identity  Forward(B) = nu_rev*B + H_irr(B)  (MatHysGetNuRev +
      MatHysIrreversible).
  (A') For BOTH hex and tet, cross-validate the B-input solve against the
      trusted H-input hysteresis solve (ComputeChiFromH path) element-by-element.
      Both solve the same nonlinear constitutive system at the same applied
      field from the virgin state, so they MUST converge to the same M.  This is
      the geometry-independent all-element-types consistency gate (it does not
      depend on per-element field-evaluation quality, which is poor for skewed
      MSC tets).
  (B) Drive a cyclic applied field and assert the (H_z, B_z) loop has nonzero
      area (hysteretic, not single-valued) with branch separation; save to JSON.
  (C) tet(4) and hex(6) both solve without Error205 and pass (A)/(A').
  (D) With b_input OFF and a linear MatLin material, the same geometry solves
      (no regression of the linear / H-input moment path).
"""

import sys
import json
import time
from pathlib import Path
from datetime import datetime, timezone

import numpy as np
import pytest

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT / "src"))

import radia as rad  # noqa: E402

MU_0 = 4.0 * np.pi * 1e-7
HERE = Path(__file__).resolve().parent


# --------------------------------------------------------------------------
# Synthetic play model: K=3, eta=[0,0.5,1], piecewise-linear shape functions.
# f_0 (reversible, eta=0) sets the small-signal mu_r; f_1,f_2 (negative, the
# Hane-Sugahara Potter sign convention) add the irreversible / hysteretic part.
# --------------------------------------------------------------------------
def _synthetic_play(mu_r=200.0, irr1=0.30, irr2=0.15, eta=(0.0, 0.3, 0.7),
                    rmax=2.0, n=41):
    r = np.linspace(0.0, rmax, n)
    a0 = 1.0 / (MU_0 * mu_r)            # reversible slope dH/dB ~ a0 -> mu_r
    f0 = a0 * r
    f1 = -irr1 * a0 * r
    f2 = -irr2 * a0 * r
    eta = np.asarray(eta, dtype=float)
    tables = [(r.tolist(), f0.tolist()),
              (r.tolist(), f1.tolist()),
              (r.tolist(), f2.tolist())]
    return len(eta), eta, tables


_CUBE_V = None


def _cube_verts(s=0.02):
    rr = s / 2.0
    return [[-rr, -rr, -rr], [rr, -rr, -rr], [rr, rr, -rr], [-rr, rr, -rr],
            [-rr, -rr, rr], [rr, -rr, rr], [rr, rr, rr], [-rr, rr, rr]]


_TET_IDX = [[0, 1, 3, 4], [1, 2, 3, 6], [1, 4, 5, 6], [3, 4, 6, 7], [1, 3, 4, 6]]


def _build_hex(mat_factory):
    """Single hexahedron cube with a per-call material instance."""
    h = rad.ObjHexahedron(_cube_verts(), [0, 0, 0])
    rad.MatApl(h, mat_factory())
    return h, [h]


def _build_tet(mat_factory):
    """5-tetrahedron decomposition of the same cube; each tet gets its own
    duplicated material instance (MatApl duplicates hysteresis materials)."""
    v = _cube_verts()
    elems = []
    for idx in _TET_IDX:
        e = rad.ObjTetrahedron([v[i] for i in idx], [0, 0, 0])
        rad.MatApl(e, mat_factory())
        elems.append(e)
    return rad.ObjCnt(elems), elems


def _forward_from_committed(mat, nu_rev, B):
    """Forward(B) = H from the material's COMMITTED play state, via the exact
    identity Forward(B) = nu_rev*B + H_irr(B) (does not depend on the solver)."""
    B = np.asarray(B, dtype=float)
    H_irr = np.asarray(rad.MatHysIrreversible(mat, B.tolist())).ravel()[:3]
    return H_irr + nu_rev * B


def _objm_list(obj):
    """Return [(center(3,), M(3,)), ...] for a single object or container."""
    d = rad.ObjM(obj)
    if isinstance(d, dict):
        return [(np.asarray(d['center'], float), np.asarray(d['magnetization'], float))]
    return [(np.asarray(c, float), np.asarray(m, float)) for (c, m) in d]


# ==========================================================================
# (A) Physics consistency on a single hexahedron: H == Forward(B) exactly.
# ==========================================================================
def test_A_hex_play_relation_satisfied():
    rad.UtiDelAll()
    K, eta, tables = _synthetic_play()
    h = rad.ObjHexahedron(_cube_verts(), [0, 0, 0])
    mat = rad.MatPlayHysteresis(K, eta, tables)
    rad.MatApl(h, mat)
    # Each element keeps its OWN duplicated material; recover the duplicated
    # handle's nu_rev (identical model -> identical nu_rev as the template).
    nu_rev = rad.MatHysGetNuRev(mat)

    B_EXT = 0.3
    cont = rad.ObjCnt([h, rad.ObjBckg(lambda p: [0, 0, B_EXT])])
    rad.SolverConfig(b_input_newton=True)
    res = rad.Solve(cont, 1e-7, 500, 0)
    rad.SolverConfig(b_input_newton=False)
    n_iter = int(res[3])
    assert n_iter < 500, f"B-input hex solve did not converge (niter={n_iter})"

    (center, M), = _objm_list(h)
    # Effective internal H at the (axis-aligned) hex centroid is reliable here.
    H = np.asarray(rad.Fld(cont, 'h', center.tolist())).ravel()[:3]
    B = MU_0 * (H + M)
    H_fwd = _forward_from_committed(mat, nu_rev, B)
    resid = np.linalg.norm(H - H_fwd) / max(np.linalg.norm(H), 1e-30)
    print(f"[A hex] niter={n_iter} |M|={np.linalg.norm(M):.0f} |B|={np.linalg.norm(B):.4f} "
          f"resid(||H-Forward(B)||/||H||)={resid:.3e}")
    assert resid < 1e-3, f"hex play relation violated: resid={resid:.3e}"
    rad.UtiDelAll()


# ==========================================================================
# (A') / (C) B-input == trusted H-input, element-by-element, hex AND tet.
# ==========================================================================
def _solve_get_M(build, b_input, B_EXT=0.3):
    rad.UtiDelAll()
    K, eta, tables = _synthetic_play()
    obj, elems = build(lambda: rad.MatPlayHysteresis(K, eta, tables))
    cont = rad.ObjCnt([obj, rad.ObjBckg(lambda p: [0, 0, B_EXT])])
    if b_input:
        rad.SolverConfig(b_input_newton=True)
    res = rad.Solve(cont, 1e-7, 500, 0)
    if b_input:
        rad.SolverConfig(b_input_newton=False)
    Ms = np.array([m for (_, m) in _objm_list(obj)])
    n_iter = int(res[3])
    rad.UtiDelAll()
    return Ms, n_iter


@pytest.mark.parametrize("name,build,dof", [
    ("hex", _build_hex, 6),
    ("tet", _build_tet, 4),
])
def test_Aprime_binput_matches_hinput(name, build, dof):
    Mb, nb = _solve_get_M(build, b_input=True)
    Mh, nh = _solve_get_M(build, b_input=False)
    assert nb < 500 and nh < 500, f"{name}: non-convergence (b={nb}, h={nh})"
    rel = np.linalg.norm(Mb - Mh) / max(np.linalg.norm(Mh), 1e-30)
    print(f"[A' {name}] dof={dof} b_input niter={nb} h_input niter={nh} "
          f"per-elem rel diff={rel:.3e}")
    # B-input moment Picard must reproduce the trusted H-input moment solve.
    assert rel < 5e-3, f"{name}: B-input vs H-input disagree (rel={rel:.3e})"


# ==========================================================================
# (C) tet(4) AND hex(6) both solve without Error205.
# ==========================================================================
@pytest.mark.parametrize("name,build", [("hex", _build_hex), ("tet", _build_tet)])
def test_C_no_error205(name, build):
    rad.UtiDelAll()
    K, eta, tables = _synthetic_play()
    obj, elems = build(lambda: rad.MatPlayHysteresis(K, eta, tables))
    cont = rad.ObjCnt([obj, rad.ObjBckg(lambda p: [0, 0, 0.2])])
    rad.SolverConfig(b_input_newton=True)
    # rad.Solve raises on Error205; reaching a finite iteration count == no error.
    res = rad.Solve(cont, 1e-5, 500, 0)
    rad.SolverConfig(b_input_newton=False)
    n_iter = int(res[3])
    Ms = np.array([m for (_, m) in _objm_list(obj)])
    assert n_iter > 0 and n_iter < 500, f"{name}: bad niter={n_iter}"
    assert np.all(np.isfinite(Ms)) and np.linalg.norm(Ms) > 0, f"{name}: non-physical M"
    print(f"[C {name}] solved b_input, niter={n_iter}, |M|={np.linalg.norm(Ms):.0f} (no Error205)")
    rad.UtiDelAll()


# ==========================================================================
# (B) Hysteresis loop: nonzero area + branch separation; saved to JSON.
# ==========================================================================
def test_B_hysteresis_loop():
    rad.UtiDelAll()
    K, eta, tables = _synthetic_play(mu_r=200.0, irr1=0.45, irr2=0.25,
                                     eta=(0.0, 0.25, 0.6))
    h = rad.ObjHexahedron(_cube_verts(), [0, 0, 0])
    mat = rad.MatPlayHysteresis(K, eta, tables)
    rad.MatApl(h, mat)
    rad.SolverConfig(b_input_newton=True)

    Bmax = 0.8
    seq = (list(np.linspace(0, Bmax, 6))
           + list(np.linspace(Bmax, -Bmax, 11))[1:]
           + list(np.linspace(-Bmax, Bmax, 11))[1:])
    loop = []
    center = [0.0, 0.0, 0.0]
    for Bext in seq:
        cont = rad.ObjCnt([h, rad.ObjBckg(lambda p, b=Bext: [0, 0, b])])
        rad.Solve(cont, 1e-6, 500, 0)
        (_, M), = _objm_list(h)
        H = np.asarray(rad.Fld(cont, 'h', center)).ravel()[:3]
        Bz = MU_0 * (H[2] + M[2])
        loop.append((float(H[2]), float(Bz), float(Bext)))
    rad.SolverConfig(b_input_newton=False)

    arr = np.array(loop)
    Hz, Bz = arr[:, 0], arr[:, 1]
    # Shoelace area of the (H,B) loop.
    area = 0.5 * abs(np.sum(Hz * np.roll(Bz, -1) - np.roll(Hz, -1) * Bz))

    # Branch separation: at a mid-field |H|, ascending and descending B differ.
    # Compare the +Bmax->-Bmax (descending) leg vs the -Bmax->+Bmax (ascending)
    # leg at a common positive H bin.
    asc = arr[16:]              # last leg: -Bmax -> +Bmax (ascending)
    desc = arr[5:16]            # +Bmax -> -Bmax (descending)
    # interpolate B at H = +2000 A/m on each branch where it exists
    def B_at_H(branch, Htgt):
        Hh, Bb = branch[:, 0], branch[:, 1]
        order = np.argsort(Hh)
        return np.interp(Htgt, Hh[order], Bb[order])
    Htest = 2000.0
    B_asc = B_at_H(asc, Htest)
    B_desc = B_at_H(desc, Htest)
    branch_gap = abs(B_desc - B_asc)

    print(f"[B] loop area={area:.1f} (H*B)  branch gap at H={Htest:.0f}: "
          f"B_desc={B_desc*1e3:.1f}mT B_asc={B_asc*1e3:.1f}mT gap={branch_gap*1e3:.1f}mT")

    assert area > 1.0, f"loop area must be nonzero (hysteretic): area={area:.3e}"
    assert branch_gap > 1e-3, (
        f"ascending/descending branches must differ (multi-valued): "
        f"gap={branch_gap:.3e} T")

    # remanence proxy: largest |B| reached on the saturated corner; coercivity
    # proxy: |H| span. (Loop is B_ext-controlled so it passes through origin at
    # the reversal; the open branches above are the hysteresis evidence.)
    out = {
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "radia_version": getattr(rad, "__version__", "unknown"),
        "python_version": sys.version.split()[0],
        "model": {"K": K, "eta": eta.tolist(), "mu_r": 200.0,
                  "irr1": 0.45, "irr2": 0.25},
        "geometry": "single_hexahedron_2cm",
        "B_ext_max": Bmax,
        "loop_Hz_Bz_Bext": loop,
        "loop_area_H_times_B": float(area),
        "branch_gap_T_at_H2000": float(branch_gap),
        "B_max_T": float(np.max(np.abs(Bz))),
        "H_max_Am": float(np.max(np.abs(Hz))),
    }
    out_path = HERE / "binput_moment_loop.json"
    out_path.write_text(json.dumps(out, indent=2), encoding="utf-8")
    print(f"[B] loop saved to {out_path}")
    rad.UtiDelAll()


# ==========================================================================
# (D) Regression: b_input OFF + linear MatLin still solves the same geometry.
# ==========================================================================
@pytest.mark.parametrize("name,build", [("hex", _build_hex), ("tet", _build_tet)])
def test_D_linear_no_binput_regression(name, build):
    rad.UtiDelAll()
    obj, elems = build(lambda: rad.MatLin(200.0))
    cont = rad.ObjCnt([obj, rad.ObjBckg(lambda p: [0, 0, 0.3])])
    # b_input is OFF (default); the moment LINEAR/H-input path must be intact.
    res = rad.Solve(cont, 1e-6, 500, 0)
    n_iter = int(res[3])
    Ms = np.array([m for (_, m) in _objm_list(obj)])
    assert n_iter > 0 and n_iter < 500, f"{name}: linear solve niter={n_iter}"
    assert np.all(np.isfinite(Ms)) and np.linalg.norm(Ms) > 0, f"{name}: bad linear M"
    print(f"[D {name}] linear (b_input OFF) solved niter={n_iter}, |M|={np.linalg.norm(Ms):.0f}")
    rad.UtiDelAll()


if __name__ == "__main__":
    t0 = time.time()
    test_A_hex_play_relation_satisfied()
    for nm, bld, dof in [("hex", _build_hex, 6), ("tet", _build_tet, 4)]:
        test_Aprime_binput_matches_hinput(nm, bld, dof)
    for nm, bld in [("hex", _build_hex), ("tet", _build_tet)]:
        test_C_no_error205(nm, bld)
    test_B_hysteresis_loop()
    for nm, bld in [("hex", _build_hex), ("tet", _build_tet)]:
        test_D_linear_no_binput_regression(nm, bld)
    print(f"\nALL B-input moment validation passed in {time.time()-t0:.1f}s")
