#!/usr/bin/env python
"""LINEAR loop-deflation scaling study on the C-type (solidify before nonlinear).

Establishes that, in the linear regime where span(L) is an exact A-invariant
eigenspace, loop deflation is BEAM-FIELD-EXACT and the speedup SCALES with both
mesh density and permeability (higher mu_r = worse loop conditioning = bigger
win). Beam-field fidelity is measured by the B VECTOR at the gap centre (origin),
the accelerator-relevant quantity -- NOT a max over all points (points inside
the iron see the loops' large local field, which is irrelevant to the beam and
was a misleading metric earlier).

For each (density, mu_r): solve OFF vs AUTO (charge..no: +/-1 null basis,
alpha auto-scaled), recording BiCGSTAB iters, wall time, auto-alpha, and the
relative B-vector change at the origin.

PASS per row: dB/B(origin) << 1e-3 (field-exact beam) AND it_off/it_on > 1.

Usage: python ctype_linear_scaling.py [--with-20]
"""
import argparse
import json
import os
import platform
import sys
import time
from datetime import datetime

import numpy as np

work_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, work_dir)
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(work_dir))), 'src'))

import radia as rad
from deflation_benchmark import load_geometry, build_model

IMAGE = '+x-z'


def solve_origin(nodes, elements, mu_r, deflate):
    rad.SetDeflateNullspace(bool(deflate), -1.0 if deflate else 0.0)
    rad.SolverConfig(bicgstab_tol=1e-8, hacapk_eps=1e-4, hacapk_leaf=10, hacapk_eta=2.0)
    model, yoke, ne = build_model(nodes, elements, ('linear', mu_r))
    t0 = time.perf_counter()
    rad.Solve(model, 0.0001, 100, 2, image=IMAGE)
    dt = time.perf_counter() - t0
    st = rad.GetSolveStats()
    B = np.array(rad.Fld(model, 'b', [0, 0, 0]))   # beam centre vector
    return {'ne': ne, 'it': st.get('linear_iterations', -1),
            'a': st.get('deflation_alpha', float('nan')),
            'cyc': st.get('deflation_cycles', -1), 'dt': dt, 'B': B}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--with-20', action='store_true', help='also 20x20x20 (165600 DOF, slow)')
    ap.add_argument('--mu-rs', type=float, nargs='*', default=[1e3, 1e4, 1e5])
    args = ap.parse_args()
    densities = ['6x6x6', '10x10x10'] + (['20x20x20'] if args.with_20 else [])

    print("C-type LINEAR loop-deflation scaling (full model is QUARTER+IMA +x-z, BiCGSTAB tol=1e-8)")
    print(f"{'dens':<9} {'mu_r':>6} {'dof':>7} | {'it_off':>7} {'it_on':>6} {'spd':>5} "
          f"{'a_used':>9} | {'|B|orig':>9} {'dB/B(0)':>9} | {'t_off':>7} {'t_on':>6}")
    results = []
    for dens in densities:
        nodes, elements = load_geometry(dens)
        for mu_r in args.mu_rs:
            off = solve_origin(nodes, elements, mu_r, False)
            on = solve_origin(nodes, elements, mu_r, True)
            spd = off['it'] / max(on['it'], 1)
            Bmag = float(np.linalg.norm(off['B'])) * 1000
            dBrel = float(np.linalg.norm(on['B'] - off['B']) / (np.linalg.norm(off['B']) + 1e-30))
            results.append({'density': dens, 'mu_r': mu_r, 'dof': off['ne']*6,
                            'it_off': off['it'], 'it_on': on['it'], 'speedup': spd,
                            'alpha': on['a'], 'cyc': on['cyc'],
                            'Bz_off_mT': float(off['B'][2])*1000, 'Bz_on_mT': float(on['B'][2])*1000,
                            'dB_rel_origin': dBrel, 't_off': off['dt'], 't_on': on['dt']})
            print(f"{dens:<9} {mu_r:>6.0e} {off['ne']*6:>7} | {off['it']:>7} {on['it']:>6} "
                  f"{spd:>5.2f} {on['a']:>9.3g} | {Bmag:>9.1f} {dBrel:>9.2e} | "
                  f"{off['dt']:>7.1f} {on['dt']:>6.1f}")

    rad.SetDeflateNullspace(False, 0.0); rad.UtiDelAll()
    out = os.path.join(work_dir, "linear_scaling_results.json")
    with open(out, 'w', encoding='utf-8') as f:
        json.dump({'timestamp': datetime.now().isoformat(), 'hostname': platform.node(),
                   'benchmark': 'ctype_linear_deflation_scaling',
                   'problem': {'image': IMAGE, 'bicgstab_tol': 1e-8}, 'results': results}, f,
                  indent=2, ensure_ascii=False)
    print(f"\nResults -> {out}")
    print("PASS per row: dB/B(0) << 1e-3 (beam field-exact) AND spd > 1.")


if __name__ == '__main__':
    main()
