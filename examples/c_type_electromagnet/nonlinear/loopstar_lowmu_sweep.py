#!/usr/bin/env python
"""Demonstrate the ALPHA-FREE loop-star gauge (rad.SetLoopStarGauge) is EFFECTIVE
on the FULL C-type magnet (no IMA) across mu_r, anchored at LOW mu_r.

Physics: loop-star removes the cycle-space eigenspace of N, which sits at
eigenvalue 1/(mu_r-1). At LOW mu_r that is O(1) (not the conditioning
bottleneck) -> plain BiCGSTAB already converges, loop-star is field-exact with
~parity iterations (establishes CORRECTNESS). As mu_r rises, 1/(mu_r-1)->0 makes
that eigenspace near-singular -> plain BiCGSTAB iterations GROW while loop-star
stays low (the EFFECT). This sweep shows both ends.

Metrics per (mu_r, loop-star on/off): BiCGSTAB linear_iterations + gap-field Bz
(field-exact check: loop-star Bz must match off Bz).

Usage: python loopstar_lowmu_sweep.py [density]   (default 6x6x6)
"""
import os, sys, time, json
import numpy as np

work_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, work_dir)
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(work_dir))), 'src'))

import radia as rad
from deflation_benchmark import load_geometry, build_model


def solve(nodes, elements, mu_r, loopstar):
    rad.SetDeflateNullspace(False, 0.0)
    rad.SetLoopStarGauge(bool(loopstar))
    rad.SolverConfig(bicgstab_tol=1e-8, hacapk_eps=1e-4, hacapk_leaf=10, hacapk_eta=2.0)
    model, yoke, ne = build_model(nodes, elements, ('linear', mu_r))
    t0 = time.perf_counter()
    rad.Solve(model, 0.0001, 100, 2)          # full model, NO image (no IMA)
    dt = time.perf_counter() - t0
    st = rad.GetSolveStats()
    B = np.array(rad.Fld(model, 'b', [0, 0, 0]))
    return {'ne': ne, 'it': int(st.get('linear_iterations', -1)), 'dt': dt, 'Bz': float(B[2])}


def main():
    dens = sys.argv[1] if len(sys.argv) > 1 else '6x6x6'
    nodes, elements = load_geometry(dens)
    mus = [2.0, 10.0, 100.0, 1000.0]
    print(f"C-type FULL model {dens} ({len(elements)*6} DOF), BiCGSTAB tol=1e-8, method=2 (HACApK)")
    print(f"\n{'mu_r':>7} | {'it_off':>7} {'it_loopstar':>11} {'it_ratio':>8} | "
          f"{'Bz_off[mT]':>11} {'Bz_ls[mT]':>11} {'dBz/Bz':>9}")
    print("-"*78)
    rows = []
    for mu in mus:
        off = solve(nodes, elements, mu, False)
        ls  = solve(nodes, elements, mu, True)
        dBz = abs(ls['Bz'] - off['Bz']) / (abs(off['Bz']) + 1e-30)
        ratio = off['it'] / ls['it'] if ls['it'] > 0 else 0.0
        print(f"{mu:>7.0f} | {off['it']:>7} {ls['it']:>11} {ratio:>7.2f}x | "
              f"{off['Bz']*1000:>11.3f} {ls['Bz']*1000:>11.3f} {dBz:>9.2e}")
        rows.append({'mu_r': mu, 'it_off': off['it'], 'it_loopstar': ls['it'],
                     'Bz_off': off['Bz'], 'Bz_loopstar': ls['Bz'], 'dBz_rel': dBz})
    out = os.path.join(work_dir, f"loopstar_lowmu_sweep_{dens}.json")
    with open(out, 'w') as f:
        json.dump({'density': dens, 'dof': len(elements)*6, 'rows': rows}, f, indent=2)
    print(f"\nsaved -> {out}")
    print("Field-exact if dBz/Bz << 1e-3 at every mu_r; loop-star EFFECT = it_ratio rising with mu_r.")
    rad.SetLoopStarGauge(False); rad.UtiDelAll()


if __name__ == '__main__':
    main()
