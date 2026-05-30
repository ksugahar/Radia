#!/usr/bin/env python
"""Why does loop-deflation HURT on the C-type while it helps on a regular cube?

Hypothesis: the fixed alpha=1.0 over-shifts on the C-type's irregular (chamfered,
non-uniform connectivity) mesh -- the RAW overcomplete cycle basis L has an
uneven L L^T spectrum, so alpha=1 throws the most-overlapped loop directions far
beyond the physical spectrum, worsening the condition number.

Clean test: the quarter operator WITHOUT IMA is a valid linear system whose loop
null space IS exactly what BuildLoopBasis spans. Sweep alpha on the fast linear
mu_r=1e4 6x6x6 case and watch the BiCGSTAB iteration count. If some alpha beats
the alpha=0 baseline, deflation CAN help with tuning (-> need auto-alpha). If
nothing beats it, the raw-L shift is fundamentally mis-conditioned here and we
need normalized/projection deflation.
"""
import os
import sys
import time

import numpy as np

work_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, work_dir)
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(work_dir))), 'src'))

import radia as rad
from deflation_benchmark import load_geometry, create_racetrack_coil, SCALE, COIL_AT


def build_linear(nodes, elements, mu_r):
    rad.UtiDelAll()
    mat = rad.MatLin(float(mu_r))
    hexes = []
    for nid in elements:
        verts = [[nodes[n][0]*SCALE, nodes[n][1]*SCALE, nodes[n][2]*SCALE] for n in nid]
        h = rad.ObjHexahedron(verts, [0, 0, 0]); rad.MatApl(h, mat); hexes.append(h)
    return rad.ObjCnt([rad.ObjCnt(hexes), create_racetrack_coil(COIL_AT)])


def main():
    mu_r = 1e4
    dens = sys.argv[1] if len(sys.argv) > 1 else '6x6x6'
    nodes, elements = load_geometry(dens)
    alphas = [0.0, 1e-4, 3e-4, 1e-3, 3e-3, 1e-2, 3e-2, 0.1, 0.3, 1.0, 3.0]
    print(f"C-type {dens} ({len(elements)*6} DOF), linear mu_r={mu_r:.0e}, NO IMA, BiCGSTAB tol=1e-8")
    print(f"{'alpha':>8} {'lin_iter':>8} {'cyc':>5} {'t(s)':>6} {'Bz(mT)':>9}")
    base = None
    for a in alphas:
        rad.SetDeflateNullspace(a > 0.0, a)
        rad.SolverConfig(bicgstab_tol=1e-8, hacapk_eps=1e-4, hacapk_leaf=10, hacapk_eta=2.0)
        model = build_linear(nodes, elements, mu_r)
        t0 = time.perf_counter()
        rad.Solve(model, 0.0001, 100, 2)
        dt = time.perf_counter() - t0
        st = rad.GetSolveStats()
        it = st.get('linear_iterations', -1)
        Bz = float(np.array(rad.Fld(model, 'b', [0, 0, 0]))[2]) * 1000.0
        if base is None:
            base = it
        tag = '' if a == 0 else (f'  {it/base:.2f}x base' if base else '')
        print(f"{a:>8.4g} {it:>8} {st.get('deflation_cycles',-1):>5} {dt:>6.1f} {Bz:>9.2f}{tag}")
    rad.SetDeflateNullspace(False, 0.0)
    rad.UtiDelAll()


if __name__ == '__main__':
    main()
