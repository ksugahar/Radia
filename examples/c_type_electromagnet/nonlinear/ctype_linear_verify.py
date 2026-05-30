#!/usr/bin/env python
"""LINEAR foundation check for charge-neutral loop deflation (do this BEFORE nonlinear).

Linear (uniform mu_r) is the clean regime: span(L) is an exact A-invariant
eigenspace, so a correct deflation MUST leave the field unchanged everywhere
while only changing the (field-free) loop content of the magnetization. We verify
on the C-type (quarter + IMA +x-z), mu_r = 1e4:

  1. FIELD-EXACT: max vector difference ||B_off - B_on|| over a cloud of field
     points (gap + outside), relative to max|B_off|. Must be ~ machine/solver tol.
  2. SPEEDUP: linear-iteration ratio off/on (the unit-norm charge-neutral basis
     should recover the cube-like acceleration that the raw 1/area basis lost).
  3. LOOPS ARE FIELD-FREE (answers "do loops carry accuracy?"): the magnetization
     M changes by a LOT (loop content removed) yet B is unchanged everywhere ->
     the removed content produces no field -> it cannot carry accuracy.

Compares OFF vs AUTO (charge-neutral, alpha<=0) and a manual alpha for reference.

Usage: python ctype_linear_verify.py [density]   (default 6x6x6)
"""
import os
import sys
import time

import numpy as np

work_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, work_dir)
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(work_dir))), 'src'))

import radia as rad
from deflation_benchmark import load_geometry, build_model

IMAGE = '+x-z'
MU_R = 1e4


def field_points():
    """A cloud of evaluation points: gap region near origin + a few outside."""
    pts = []
    for x in (-0.02, 0.0, 0.02):
        for y in (-0.02, 0.0, 0.02):
            for z in (0.0, 0.01, 0.02):
                pts.append([x, y, z])
    pts += [[0, 0, 0.05], [0.05, 0, 0], [0, 0.1, 0], [0.1, 0.1, 0.05]]
    return pts


def solve_get(nodes, elements, deflate_alpha, pts):
    if deflate_alpha is None:
        rad.SetDeflateNullspace(False, 0.0)
    else:
        rad.SetDeflateNullspace(True, deflate_alpha)
    rad.SolverConfig(bicgstab_tol=1e-8, hacapk_eps=1e-4, hacapk_leaf=10, hacapk_eta=2.0)
    model, yoke, ne = build_model(nodes, elements, ('linear', MU_R))
    t0 = time.perf_counter()
    rad.Solve(model, 0.0001, 100, 2, image=IMAGE)
    dt = time.perf_counter() - t0
    st = rad.GetSolveStats()
    B = np.array([rad.Fld(model, 'b', p) for p in pts])           # (Npts, 3)
    M = np.array([m[1] for m in rad.ObjM(yoke)])                  # (ne, 3)
    return {'dt': dt, 'it': st.get('linear_iterations', -1),
            'a': st.get('deflation_alpha', float('nan')),
            'cyc': st.get('deflation_cycles', -1), 'B': B, 'M': M}


def main():
    dens = sys.argv[1] if len(sys.argv) > 1 else '6x6x6'
    nodes, elements = load_geometry(dens)
    pts = field_points()
    print(f"C-type {dens} ({len(elements)*6} DOF) LINEAR mu_r={MU_R:.0e}, image={IMAGE}")
    print(f"{len(pts)} field points (gap + outside), BiCGSTAB tol=1e-8\n")

    off = solve_get(nodes, elements, None, pts)
    auto = solve_get(nodes, elements, -1.0, pts)

    Bref = np.max(np.linalg.norm(off['B'], axis=1))
    def field_dev(r):
        return np.max(np.linalg.norm(r['B'] - off['B'], axis=1)) / Bref
    def mag_dev(r):
        mref = np.max(np.linalg.norm(off['M'], axis=1))
        return np.max(np.linalg.norm(r['M'] - off['M'], axis=1)) / mref

    print(f"{'case':<14} {'it':>7} {'t(s)':>7} {'a_used':>10} {'spdup':>6} "
          f"{'max dB/B':>10} {'max dM/M':>10}")
    print(f"{'off':<14} {off['it']:>7} {off['dt']:>7.1f} {0.0:>10.3g} {1.0:>6.2f} "
          f"{0.0:>10.2e} {0.0:>10.2e}")
    print(f"{'auto':<14} {auto['it']:>7} {auto['dt']:>7.1f} {auto['a']:>10.3g} "
          f"{off['it']/max(auto['it'],1):>6.2f} {field_dev(auto):>10.2e} {mag_dev(auto):>10.2e}")

    # per-point breakdown: where is the field changing? (gap vs near-iron)
    print("\nPer-point |B_off| and |dB| (auto):")
    print(f"  {'point (m)':<24} {'|B_off|(mT)':>12} {'|dB|(mT)':>10} {'dB/|Bloc|':>10}")
    dB = np.linalg.norm(auto['B'] - off['B'], axis=1) * 1000
    Bloc = np.linalg.norm(off['B'], axis=1) * 1000
    order = np.argsort(-dB)
    for idx in order[:8]:
        p = pts[idx]
        rel = dB[idx] / Bloc[idx] if Bloc[idx] > 1e-9 else float('inf')
        print(f"  ({p[0]:+.2f},{p[1]:+.2f},{p[2]:+.2f})        "
              f"{Bloc[idx]:>12.2f} {dB[idx]:>10.3f} {rel:>10.3f}")
    print("\nInterpretation:")
    print(f"  - cycles installed: {auto['cyc']}")
    print(f"  - FIELD-EXACT if max dB/B << 1e-3 (loops carry no field).")
    print(f"  - LOOPS REMOVED if max dM/M is sizable (magnetization changed a lot).")
    print(f"  - Both together => loops are a field-free gauge freedom: removing")
    print(f"    them changes M but NOT B, so they cannot 'carry accuracy'.")
    rad.SetDeflateNullspace(False, 0.0); rad.UtiDelAll()


if __name__ == '__main__':
    main()
