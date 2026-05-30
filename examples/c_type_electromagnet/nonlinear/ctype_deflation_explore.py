#!/usr/bin/env python
"""C-type loop-deflation explorer: material x image(IMA) x alpha.

Findings driving this script:
- alpha=1 (the cube value) HURTS the C-type; the optimal alpha (~0.03 at mu=1e4)
  must lift the loop cluster at 1/chi to the physical eigenvalue scale.
- The paper solves the QUARTER mesh with IMA (+x-z); the no-IMA solve gives only
  1/4 of the field (B_z ~ -119 vs ELF -960 mT). So we must test WITH image to be
  physically correct -- and check that deflation (loop basis built from the
  quarter adjacency) stays correct + helpful under IMA (the IMA operator's null
  space may differ from the bare quarter cycle space).

Usage:
  python ctype_deflation_explore.py --density 6x6x6 --material linear --mu-r 1e4 \
         --image +x-z --alphas 0 0.01 0.03 0.1
  python ctype_deflation_explore.py --density 6x6x6 --material nonlinear \
         --image +x-z --alphas 0 0.03 0.1
"""
import argparse
import os
import sys
import time

import numpy as np

work_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, work_dir)
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(work_dir))), 'src'))

import radia as rad
from deflation_benchmark import load_geometry, build_model, load_elf_bz


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--density', default='6x6x6')
    ap.add_argument('--material', choices=['linear', 'nonlinear'], default='linear')
    ap.add_argument('--mu-r', type=float, default=1e4)
    ap.add_argument('--image', default='', help="IMA spec, e.g. +x-z (empty = full model, no IMA)")
    ap.add_argument('--alphas', type=float, nargs='*', default=[0.0, 0.01, 0.03, 0.1])
    ap.add_argument('--max-nl', type=int, default=100)
    ap.add_argument('--precision', type=float, default=1e-4, help='nonlinear stop (paper uses 1e-3)')
    ap.add_argument('--relax', type=float, default=0.0, help='under-relaxation (0=full step)')
    args = ap.parse_args()

    material = ('nonlinear',) if args.material == 'nonlinear' else ('linear', args.mu_r)
    mlabel = 'nonlinear' if args.material == 'nonlinear' else f'linear mu={args.mu_r:.0e}'
    nodes, elements = load_geometry(args.density)
    elf_bz = load_elf_bz(args.density)

    print(f"C-type {args.density} ({len(elements)*6} DOF), {mlabel}, "
          f"image='{args.image or 'NONE'}', BiCGSTAB tol=1e-8")
    if elf_bz is not None:
        print(f"ELF reference B_z = {elf_bz*1000:.1f} mT")
    print("(alpha: 0 = OFF baseline; <0 = AUTO-scaled; >0 = manual override)")
    print(f"{'alpha_in':>8} {'a_used':>8} {'lin_iter':>8} {'cyc':>5} {'t(s)':>7} {'Bz(mT)':>10} {'note':>14}")

    base = None
    for a in args.alphas:
        rad.SetDeflateNullspace(a != 0.0, a)        # a<0 -> enable + auto; a>0 -> manual
        rad.SolverConfig(bicgstab_tol=1e-8, hacapk_eps=1e-4, hacapk_leaf=10, hacapk_eta=2.0,
                         relax_param=args.relax)
        model, yoke, ne = build_model(nodes, elements, material)
        t0 = time.perf_counter()
        try:
            rad.Solve(model, args.precision, args.max_nl, 2, image=args.image)
            err = ''
        except RuntimeError as e:
            err = ' DID-NOT-CONVERGE'
        dt = time.perf_counter() - t0
        st = rad.GetSolveStats()
        it = st.get('linear_iterations', -1)
        a_used = st.get('deflation_alpha', float('nan'))
        Bz = float(np.array(rad.Fld(model, 'b', [0, 0, 0]))[2]) * 1000.0
        if base is None and a == 0.0:
            base = it
        note = ''
        if base and a != 0.0:
            note = f'{it/base:.2f}x base'
        if a < 0:
            note = ('AUTO ' + note).strip()
        print(f"{a:>8.4g} {a_used:>8.4g} {it:>8} {st.get('deflation_cycles',-1):>5} "
              f"{dt:>7.1f} {Bz:>10.2f} {note:>14}{err}")

    rad.SetDeflateNullspace(False, 0.0)
    rad.UtiDelAll()


if __name__ == '__main__':
    main()
