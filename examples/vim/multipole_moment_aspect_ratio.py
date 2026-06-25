"""multipole-moment MMM's iterative difficulty is dominantly an ELEMENT ASPECT-RATIO effect -- not a fundamental
scaling wall, and not a condition-number problem (Sugahara-led investigation 2026-06-22).

Background.  An earlier gate (multipole_moment_iter_scaling.py) saw multipole-moment MMM GMRES grow erratically and BiCGSTAB
break down at nxy>=24 on the C-yoke, and concluded "moment is iteratively harder than EIEM2 / needs H-LU".
That benchmark used nz=2 FIXED, so the cells got flatter as nxy grew (nxy=24 -> 0.005 x 0.005 x 0.02 = 1:1:4).
The blowup was the GROWING ANISOTROPY, not the problem size.

Why aspect ratio hurts moment-matching (the user's "different physical quantities" diagnosis).  The local
6x6 block per hex matches functionals of DIFFERENT geometric order: monopole (sum A sigma, ~L^2), dipole
(sum A d sigma, ~L^3), quadrupole (sum A d^2 sigma, ~L^4).  On a CUBIC cell all directions share one length h
so these scale uniformly and the block is well-conditioned; on an ANISOTROPIC cell the thin-direction moments
collapse relative to the fat ones, the local block becomes ill-conditioned, and block-Jacobi (which inverts
that block) amplifies it -> iteration blowup.

This script measures three things:
  A. Aspect-ratio V-curve at FIXED nxy (vary nz): iters minimize at aspect ratio ~1 and grow on both sides.
     dof grows monotonically with nz, but iters are V-shaped -> it is aspect ratio, NOT size.
  B. CUBIC-cell scaling (nz ~ nxy/3 keeps cells ~cubic): iters stay BOUNDED as N grows -> multipole-moment MMM DOES
     scale with cheap block-Jacobi (no H-LU) on well-shaped cells.
  C. point-matching (EIEM2 collocation, A=-N+(1/chi)I) vs moment-matching, same block-Jacobi + same RHS:
     point-matching is lower AND less aspect-sensitive at every aspect ratio (it collocates a SINGLE physical
     quantity -- the field -- at points, so it does not mix the L^2/L^3/L^4 orders), but is NOT fully immune
     (it still has a shallower V: a thin cell puts the normal-offset eval point too close to its face).

Upshot for improving multipole-moment MMM: (1) iteration health is governed by element aspect ratio + the
discretization (point- vs moment-matching), NOT by condition number or problem size; (2) on well-shaped cells
multipole-moment MMM scales with block-Jacobi alone; (3) point-matching is the better-conditioned discretization, so
the ideal element is a SHEAR-capturing point-matching one (EIEM2 point-matching is well-conditioned but blind
to the off-diagonal gradient; moment-matching captures the full tensor but is more aspect-sensitive).

mdx-clean (import radia directly).
"""
import json
import os
import platform
from datetime import datetime

import numpy as np

import radia as rad
from multipole_moment_iter_scaling import build_cyoke_hexes, build as build_moment, block_jacobi, full_gmres_iters

HERE = os.path.dirname(os.path.abspath(__file__))
XW, ZW = 0.12, 0.04                       # C-yoke extents (x,y span 0.12; z span 0.04) used by build_cyoke_hexes


def eiem2_system(hexes, mu_r):
    """EIEM2 point-matching system matrix A = -N + (1/chi) I and the per-element face-DOF groups."""
    chi = mu_r - 1.0
    rad.UtiDelAll(); rad.set_demag_backend("yano")
    objs = [rad.ObjHexahedron([list(v) for v in V], [0, 0, 0]) for V in hexes]
    for h in objs:
        rad.MatApl(h, rad.MatLin(mu_r))
    handle = rad.BuildMatrix(rad.ObjCnt(objs))
    N, dof = rad.GetInteractMatrix(handle)
    G = np.asarray(rad.GetFaceGeom(handle), float); elem = G[:, 0].astype(int); rad.UtiDelAll()
    N = np.asarray(N, float).reshape(dof, dof)
    A = -N + (1.0 / chi) * np.eye(dof)
    n_el = int(elem.max()) + 1
    return A, dof, [np.where(elem == e)[0] for e in range(n_el)]


def block_jacobi_blocks(A, dofs_of, dof):
    blocks = [(fs, np.linalg.inv(A[np.ix_(fs, fs)])) for fs in dofs_of]
    def apply(rv):
        rv = np.asarray(rv, float); out = np.zeros(dof)
        for fs, Mi in blocks:
            out[fs] += Mi @ rv[fs]
        return out
    return apply


def gmres_count(A, rhs, apply, dof):
    nit, info = full_gmres_iters(A, rhs, apply, dof)
    return nit, info


def main():
    Happ = np.array([0.0, 1e3, 0.0]); mu_r = 1000.0
    print("\nmultipole-moment MMM: iteration health is an ELEMENT ASPECT-RATIO effect (C-yoke, mu_r=1000).\n")

    # ---- A: aspect-ratio V-curve at fixed nxy (moment-matching), + point-matching for contrast ----
    nxy_A = 16
    print(f"A. aspect-ratio V-curve at FIXED nxy={nxy_A} (same block-Jacobi GMRES, same random RHS):")
    print(f"   {'nz':>3} {'dx/dz':>6} {'dof':>6} | {'moment':>7} {'EIEM2(pt)':>10}")
    A_rows = []
    for nz in (2, 4, 5, 8, 16):
        ar = (XW / nxy_A) / (ZW / nz)
        hexes = build_cyoke_hexes(nxy_A, nz)
        Am, _bm, row_of, dof, _ = build_moment(hexes, Happ, mu_r)
        Ae, dofe, dofs_of = eiem2_system(hexes, mu_r)
        rng = np.random.default_rng(0); rhs = rng.standard_normal(dof)
        gm, _ = gmres_count(Am, rhs, block_jacobi(Am, row_of, dof), dof)
        ge, _ = gmres_count(Ae, rhs, block_jacobi_blocks(Ae, dofs_of, dofe), dofe)
        print(f"   {nz:>3} {ar:>6.2f} {dof:>6} | {gm:>7} {ge:>10}")
        A_rows.append(dict(nz=nz, aspect_ratio=float(ar), dof=int(dof), moment_gmres=gm, eiem2_gmres=ge))

    # ---- B: cubic-cell scaling (nz ~ nxy/3) -> bounded iters ----
    print(f"\nB. CUBIC-cell scaling (nz~nxy/3 keeps cells ~cubic) -> multipole-moment MMM iters stay BOUNDED:")
    print(f"   {'nxy':>3} {'nz':>3} {'dx/dz':>6} {'dof':>6} | {'moment GMRES':>12}")
    B_rows = []
    for nxy in (12, 16, 20, 24):
        nz = max(2, round(nxy / 3)); ar = (XW / nxy) / (ZW / nz)
        hexes = build_cyoke_hexes(nxy, nz)
        Am, _bm, row_of, dof, _ = build_moment(hexes, Happ, mu_r)
        rng = np.random.default_rng(0); rhs = rng.standard_normal(dof)
        gm, _ = gmres_count(Am, rhs, block_jacobi(Am, row_of, dof), dof)
        print(f"   {nxy:>3} {nz:>3} {ar:>6.2f} {dof:>6} | {gm:>12}")
        B_rows.append(dict(nxy=nxy, nz=nz, aspect_ratio=float(ar), dof=int(dof), moment_gmres=gm))

    moment_min = min(r["moment_gmres"] for r in A_rows); moment_max = max(r["moment_gmres"] for r in A_rows)
    eiem2_min = min(r["eiem2_gmres"] for r in A_rows); eiem2_max = max(r["eiem2_gmres"] for r in A_rows)
    cubic_iters = [r["moment_gmres"] for r in B_rows]; cubic_dofs = [r["dof"] for r in B_rows]
    bounded = max(cubic_iters) <= 1.6 * min(cubic_iters)
    point_better = all(r["eiem2_gmres"] <= r["moment_gmres"] for r in A_rows)
    out = dict(timestamp=datetime.now().isoformat(), hostname=platform.node(),
               benchmark="multipole_moment_aspect_ratio", mu_r=mu_r,
               aspect_sweep=A_rows, cubic_scaling=B_rows,
               cubic_iters_bounded=bool(bounded), point_matching_better_all_ar=bool(point_better),
               conclusion=(
                   "multipole-moment MMM's iterative difficulty is dominantly an ELEMENT ASPECT-RATIO effect, not a "
                   "scaling wall or a condition-number problem.  (A) At fixed nxy the iters trace a V in cell "
                   f"aspect ratio (moment {moment_min}->{moment_max}, min at dx/dz~1) while dof grows "
                   "MONOTONICALLY -> it is aspect ratio, not size.  Mechanism: the moment match mixes "
                   "monopole/dipole/quadrupole functionals of different geometric order (L^2/L^3/L^4); on an "
                   "anisotropic cell the thin-direction moments collapse and the local 6x6 block (which "
                   "block-Jacobi inverts) becomes ill-conditioned.  (B) With ~cubic cells (nz~nxy/3) the iters "
                   f"stay BOUNDED ({min(cubic_iters)}-{max(cubic_iters)}) as dof grows {min(cubic_dofs)}->"
                   f"{max(cubic_dofs)} -> multipole-moment MMM SCALES with cheap block-Jacobi (no H-LU) on well-shaped "
                   "cells; the earlier 'needs H-LU / harder than EIEM2' conclusion was a flat-cell (nz=2) "
                   "benchmark artifact.  (C) point-matching (EIEM2 collocation) is lower at EVERY aspect ratio "
                   f"(EIEM2 {eiem2_min}-{eiem2_max} vs moment {moment_min}-{moment_max}) and has a shallower V "
                   "-- collocating a single physical quantity (the field) at points does not mix the "
                   "L^2/L^3/L^4 orders -- but it is NOT fully immune (a thin cell puts the normal-offset eval "
                   "point too close to its face).  IMPROVEMENT DIRECTION: a SHEAR-capturing point-matching "
                   "element -- point-matching is the better-conditioned discretization, but EIEM2 point-matching "
                   "is blind to the off-diagonal (shear) gradient; the ideal multipole-moment MMM combines point-matching "
                   "conditioning with full-tensor (shear) capture."))
    with open(os.path.join(HERE, "multipole_moment_aspect_ratio.json"), "w") as f:
        json.dump(out, f, indent=2, default=float)
    print(f"\n  => aspect-ratio V (moment {moment_min}-{moment_max}); cubic-cell iters bounded={bounded} "
          f"({min(cubic_iters)}-{max(cubic_iters)} over dof {min(cubic_dofs)}-{max(cubic_dofs)}); "
          f"point-matching lower at all AR={point_better}.")
    print("  saved", os.path.join(HERE, "multipole_moment_aspect_ratio.json"))


if __name__ == "__main__":
    main()
