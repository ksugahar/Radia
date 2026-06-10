"""Wire extraction from a foliated volume current  (Stage B / B2).

Turns the foliated-Clebsch solenoid's continuous solved lambda (a VOLUME H1
scalar, examples/feec_vim/foliated_clebsch_solenoid.py) into WINDABLE DISCRETE
EQUAL-CURRENT WIRES, and validates that their straight-segment Biot-Savart
reproduces the target Bz.

Pipeline (per the Stage-B design):
  [GATE 1]  helicity_relative ~ 0  (Clebsch-feasible; helicity_diagnostic.py)
  [STAGE 2] foliated lambda, mu = r FIXED  (linear ACA+TSVD; foliated #3)
  [STAGE 3] per mu-leaf r_k: sample lambda(phi,z), marching-squares contour at
            GLOBAL equal-Delta-lambda -> (phi,z) wires; each carries
            I = Delta-lambda * Delta-mu  (signed flux = Jacobian, verified).
  [STAGE 4] Biot-Savart of the wires vs B0, cross-checked against Radia
            (rad.ObjFlmCur + rad.Fld) -- the two-codebase agreement invariant.

For the solenoid lambda is monotone in z (J_phi = d lambda/dz) so the contours
are RINGS (the b1=1 net-azimuthal-current case); chaining them is the already-
solved Gz "smooth helix".  We validate the EXTRACTED EQUAL-CURRENT RINGS here;
single-stroke chaining is the downstream step.

Run:  python examples/feec_vim/foliated_solenoid_wires.py
"""
import os
import sys
import math

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "../../src/radia"))

from ngsolve import H1, GridFunction, TaskManager, CF, x, y, sqrt
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from foliated_clebsch_solenoid import (
    build_mesh, assemble_response, R0, R1, ZL, B0, MU0,
)
from helicity_diagnostic import helicity_relative
from stream_function import aca_tsvd, pseudo_inverse_solve

PI = math.pi


def solve_lambda(maxh=0.03):
    mesh = build_mesh(maxh)
    fes = H1(mesh, order=2)
    zt = np.linspace(-0.6 * ZL, 0.6 * ZL, 9)
    targets = [(0.0, 0.0, float(zz)) for zz in zt]
    A = assemble_response(mesh, fes, targets)
    res = aca_tsvd(len(targets), fes.ndof, lambda i, j: A[i, j],
                   modes=len(targets), kmax=len(targets), aca_eps=1e-10, method=3)
    lam = pseudo_inverse_solve(res, np.full(len(targets), B0))
    gf = GridFunction(fes)
    gf.vec.FV().NumPy()[:] = lam
    return mesh, gf, targets


def _sample_cylinder(mesh, gf, r_k, nphi, nz):
    phis = np.linspace(0, 2 * PI, nphi, endpoint=True)
    zs = np.linspace(-0.95 * ZL, 0.95 * ZL, nz)
    LAM = np.zeros((nz, nphi))
    for j, ph in enumerate(phis):
        cx, sy = r_k * math.cos(ph), r_k * math.sin(ph)
        for i, zz in enumerate(zs):
            LAM[i, j] = gf(mesh(cx, sy, zz))
    return phis, zs, LAM


def extract_wires(mesh, gf, radii, dmu, dlam, nphi=49, nz=41):
    """Contour lambda on each cylinder at equal Delta-lambda -> 3D wires.

    Returns (wires, currents): wires = list of (n,3) arrays (closed loops),
    currents = list of I = dlam * dmu (equal-current by construction).
    """
    wires, currents = [], []
    # global level set (same Delta-lambda on every leaf = equal current/wire)
    lo, hi = +1e30, -1e30
    samples = []
    for r_k in radii:
        phis, zs, LAM = _sample_cylinder(mesh, gf, r_k, nphi, nz)
        samples.append((r_k, phis, zs, LAM))
        lo, hi = min(lo, LAM.min()), max(hi, LAM.max())
    levels = np.arange(math.ceil(lo / dlam) * dlam, hi, dlam)
    for (r_k, phis, zs, LAM) in samples:
        PHI, Z = np.meshgrid(phis, zs)
        cs = plt.contour(PHI, Z, LAM, levels=levels)
        for segs in cs.allsegs:
            for seg in segs:
                if len(seg) < 2:
                    continue
                ph, zz = seg[:, 0], seg[:, 1]
                pts = np.column_stack([r_k * np.cos(ph), r_k * np.sin(ph), zz])
                # close a ring that spans (nearly) the full azimuth
                if (ph.max() - ph.min()) > 1.8 * PI:
                    pts = np.vstack([pts, pts[0]])
                wires.append(pts)
                currents.append(dlam * dmu)
        plt.close("all")
    return wires, currents


def _bz_segments(wires, currents, obs):
    """Bz at obs points from straight current segments (vectorised over segs).

    Finite-segment law: B = (mu0 I/4pi) (a x b)(|a|+|b|) /
        (|a||b|(|a||b| + a.b)),  a=P0-P1, b=P0-P2.
    """
    out = []
    for p in obs:
        P0 = np.asarray(p)
        Bz = 0.0
        for W, I in zip(wires, currents):
            P1, P2 = W[:-1], W[1:]
            a, b = P0 - P1, P0 - P2
            am = np.linalg.norm(a, axis=1)
            bm = np.linalg.norm(b, axis=1)
            cr = np.cross(a, b)
            denom = am * bm * (am * bm + np.einsum("ij,ij->i", a, b))
            good = denom > 1e-30
            f = np.zeros_like(am)
            f[good] = (am[good] + bm[good]) / denom[good]
            Bz += (MU0 * I / (4 * PI)) * np.sum(f * cr[:, 2])
        out.append(Bz)
    return np.array(out)


def _bz_radia(wires, currents, obs):
    """Cross-check: Bz from Radia rad.ObjFlmCur + rad.Fld (independent codebase)."""
    import radia as rad
    rad.UtiDelAll()
    objs = []
    for W, I in zip(wires, currents):
        pts = [[float(x), float(y), float(z)] for (x, y, z) in W]
        objs.append(rad.ObjFlmCur(pts, float(I)))
    cnt = rad.ObjCnt(objs)
    out = [float(rad.Fld(cnt, "b", list(p))[2]) for p in obs]
    rad.UtiDelAll()
    return np.array(out)


def main():
    print("=== wire extraction from a foliated volume current (Stage B / B2) ===")
    with TaskManager():
        mesh, gf, targets = solve_lambda()

        # Gate 1: helicity (must be ~0 for clean Clebsch/level-set extraction).
        # T = lambda * grad(mu),  mu = r  ->  grad(mu) = (x/r, y/r, 0).
        r = sqrt(x ** 2 + y ** 2)
        T_cf = gf * CF((x / r, y / r, 0.0))
        h_rel = helicity_relative(mesh, T_cf)

        # lambda range -> choose Delta-lambda for ~10 rings per cylinder
        K = 3
        radii = np.linspace(R0 + (R1 - R0) / (2 * K), R1 - (R1 - R0) / (2 * K), K)
        dmu = (R1 - R0) / K
        lam_vals = gf.vec.FV().NumPy()
        lam_span = lam_vals.max() - lam_vals.min()
        dlam = lam_span / 14.0

        wires, currents = extract_wires(mesh, gf, radii, dmu, dlam)
        Bz_wire = _bz_segments(wires, currents, targets)
        Bz_radia = _bz_radia(wires, currents, targets)

    field_res = np.linalg.norm(Bz_wire - B0) / (abs(B0) * math.sqrt(len(targets)))
    twocode = np.linalg.norm(Bz_wire - Bz_radia) / (np.linalg.norm(Bz_wire) + 1e-30)
    I_arr = np.array(currents)
    equal_cur = (I_arr.max() - I_arr.min()) / (I_arr.mean() + 1e-30)

    print(f"   gate-1 helicity H_rel      = {h_rel:.2e}   (clean extraction needs ~0)")
    print(f"   n_wires                    = {len(wires)}  (I = dlam*dmu = {currents[0]:.3e} A each)")
    print(f"   equal-current spread       = {equal_cur:.2e}   (wires equal by construction)")
    print(f"   field residual vs B0       = {field_res:.2e}")
    print(f"   two-codebase (mine-Radia)  = {twocode:.2e}   (straight-seg vs rad.ObjFlmCur+Fld)")
    print(f"   mean Bz_wire               = {Bz_wire.mean():.4e} T   (target {B0:.1e})")

    ok = (h_rel < 5e-2 and field_res < 8e-2 and twocode < 1e-3 and equal_cur < 1e-9
          and len(wires) >= 3 * K)
    print(f"\n=== PASS: {ok} ===")
    return {"h_rel": float(h_rel), "n_wires": len(wires), "field_res": float(field_res),
            "twocode": float(twocode), "equal_cur": float(equal_cur)}, ok


if __name__ == "__main__":
    _, ok = main()
    raise SystemExit(0 if ok else 1)
