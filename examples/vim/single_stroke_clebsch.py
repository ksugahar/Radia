"""Rigorous single-stroke via the Clebsch / cohomology-secular term.

The single-stroke (一筆書き) problem: a stream-function design gives a family of
iso-contours {psi = k*dpsi} -- SEPARATE closed turns.  To wind them as ONE
continuous wire you must connect them.  The single-stroke-chain skill does this
HEURISTICALLY (ad-hoc 'rung' connectors, field impact estimated + minimised by
RMS).  This puts it on a RIGOROUS Clebsch foundation:

  Make the current potential MULTIVALUED:  Psi = psi + (dpsi/2pi) * theta.
  The SINGLE level-set {Psi = c} is then ONE connected curve -- a continuous
  HELIX threading every turn.  The connection is NOT an ad-hoc rung: it is the
  SECULAR (cohomology b1=1) term -- a net axial 'lead' current dpsi = I -- and
  the single-stroke wire is EXACTLY the level-set of a well-defined potential.

What is rigorous (vs the heuristic):
  1. TOPOLOGY: single-valued psi -> N disconnected contours; Psi = psi + secular
     -> exactly 1 connected component (the single stroke).  By construction.
  2. The connection = the COHOMOLOGY secular term (net axial current I), the same
     object as cohomology_net_current.py -- not an ad-hoc path.
  3. The Clebsch level-set DISTRIBUTES that unavoidable connection current
     axisymmetrically (the helix), giving far LESS near-field stray than the
     heuristic LOCALISED rung -- demonstrated here (rungs >> helix near the join).

Honest scope: for a solenoid the rigorous level-set is the ordinary helix (a
clean check of the framework).  For a general winding pattern the level-set
extraction still defines the wire rigorously, but CHOOSING the multivalued
potential whose single level-set covers an arbitrary contour family is the
routing problem (F1-like: rigorous parameterisation, non-convex routing).

Verified with Radia Biot-Savart (rad.ObjFlmCur).  Run:
    python examples/vim/single_stroke_clebsch.py
"""
from __future__ import annotations

import math
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "src",
                                "radia"))

R = 0.05        # winding radius [m]
L = 0.20        # length [m]
I = 100.0       # current per turn [A]
MU0 = 4e-7 * math.pi


def _rings(rad, n):
    """n separate closed rings = the ideal multi-turn design (N components)."""
    o = []
    for z in np.linspace(-L / 2, L / 2, n):
        th = np.linspace(0, 2 * math.pi, 72, endpoint=True)
        o.append(rad.ObjFlmCur(
            [[R * math.cos(a), R * math.sin(a), float(z)] for a in th], I))
    return rad.ObjCnt(o)


def _helix(rad, n, ns=72):
    """ONE continuous helix = the Clebsch single-stroke {Psi=c} (1 component)."""
    th = np.linspace(0, 2 * math.pi * n, n * ns + 1)
    z = -L / 2 + (th / (2 * math.pi * n)) * L
    return rad.ObjFlmCur(
        [[R * math.cos(a), R * math.sin(a), float(zz)]
         for a, zz in zip(th, z)], I)


def _rings_rungs(rad, n):
    """n rings + (n-1) straight axial rungs at theta=0 = the HEURISTIC stroke."""
    zs = np.linspace(-L / 2, L / 2, n)
    o = []
    for z in zs:
        th = np.linspace(0, 2 * math.pi, 72, endpoint=True)
        o.append(rad.ObjFlmCur(
            [[R * math.cos(a), R * math.sin(a), float(z)] for a in th], I))
    for k in range(n - 1):
        o.append(rad.ObjFlmCur(
            [[R, 0.0, float(zs[k])], [R, 0.0, float(zs[k + 1])]], I))
    return rad.ObjCnt(o)


def _fld(rad, obj, pts):
    return np.array([rad.Fld(obj, "b", list(p)) for p in pts])


def main(n_turns=20):
    import radia as rad

    print("=== rigorous single-stroke via the Clebsch/cohomology-secular term ===")
    n = int(n_turns)

    # on-axis Bz: all three are valid solenoids (same N*I ampere-turns)
    axis = np.array([[0.0, 0.0, z] for z in np.linspace(-0.2 * L, 0.2 * L, 5)])
    # NEAR the connection (theta~0, rho just outside R): the rung localises here
    near = []
    for z in np.linspace(-0.3 * L, 0.3 * L, 3):
        for a in np.linspace(-0.3, 0.3, 7):
            near.append([1.3 * R * math.cos(a), 1.3 * R * math.sin(a), float(z)])
    near = np.array(near)

    rad.UtiDelAll()
    rings, helix, rungs = _rings(rad, n), _helix(rad, n), _rings_rungs(rad, n)
    Bz_rings = _fld(rad, rings, axis)[:, 2]
    Bz_helix = _fld(rad, helix, axis)[:, 2]
    Bz_rungs = _fld(rad, rungs, axis)[:, 2]
    Bn_rings = _fld(rad, rings, near)
    Bn_helix = _fld(rad, helix, near)
    Bn_rungs = _fld(rad, rungs, near)
    rad.UtiDelAll()

    # valid coils: on-axis Bz of the single-stroke wires ~ the ideal rings
    bz_ref = float(np.mean(Bz_rings))
    helix_bz_err = float(np.max(np.abs(Bz_helix - Bz_rings)) / abs(bz_ref))
    rungs_bz_err = float(np.max(np.abs(Bz_rungs - Bz_rings)) / abs(bz_ref))

    # near-connection STRAY (deviation from the ideal rings): helix << rungs
    stray_helix = float(np.linalg.norm(Bn_helix - Bn_rings))
    stray_rungs = float(np.linalg.norm(Bn_rungs - Bn_rings))
    rungs_over_helix = stray_rungs / (stray_helix + 1e-30)

    # the connection = net axial current I (the secular term): one wire makes one
    # net axial traverse -> net axial current = I, INDEPENDENT of n (the secular
    # coefficient dpsi = I).  (analytic statement; the helix carries it
    # axisymmetrically, the rung localises it at theta=0.)
    print(f"turns N={n}, current I={I} A  (N*I = {n*I:.0f} ampere-turns)")
    print(f"on-axis Bz = {bz_ref*1e3:.3f} mT (ideal rings); "
          f"single-stroke Bz error: helix {helix_bz_err:.2e}, rungs {rungs_bz_err:.2e}")
    print(f"near-connection stray (vs ideal rings): "
          f"helix={stray_helix:.3e}  rungs={stray_rungs:.3e}  "
          f"-> rungs/helix = {rungs_over_helix:.2f}")
    print("  => the Clebsch level-set (helix) distributes the unavoidable")
    print("     connection current (net axial I = the secular term) axisymmetrically;")
    print("     the heuristic rung localises it -> larger near-field stray.")

    result = {
        "n_turns": n, "ampere_turns": n * I, "bz_axis_T": bz_ref,
        "helix_bz_err": helix_bz_err, "rungs_bz_err": rungs_bz_err,
        "stray_helix": stray_helix, "stray_rungs": stray_rungs,
        "rungs_over_helix": rungs_over_helix,
        "net_axial_current_A": I,          # the secular coefficient dpsi
    }
    ok = (helix_bz_err < 0.05            # the single-stroke helix IS a valid coil
          and rungs_over_helix > 1.5)    # heuristic rung worse near the join
    print(f"\n=== PASS: {ok} ===")
    return result, ok


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-turns", type=int, default=20)
    a = ap.parse_args()
    _, ok = main(n_turns=a.n_turns)
    raise SystemExit(0 if ok else 1)
