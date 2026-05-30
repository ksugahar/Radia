#!/usr/bin/env python
"""End-to-end test of the ALPHA-FREE loop-star gauge (rad.SetLoopStarGauge) on the
C-type magnet, LINEAR high mu_r where the loop ill-conditioning is worst.

Loop-star support status (2026-05-30):
  - FULL model (no image=) and SYMMETRIC IMA planes ('+x'/'+y'/'+z'): field-exact,
    alpha-free. The standard star basis spans the complement of the loop null space.
  - ANTISYMMETRIC IMA planes ('-x'/'-y'/'-z', e.g. the '-z' in '+x-z'): SUPPORTED.
    The antisymmetric image folds the operator so ker(N_ima) gains nP-1 extra
    "plane-coupled" modes (nP = #element-faces on the plane) beyond the topological
    loops. ComputePlaneSlabAddedModes finds them SVD-FREE (Gram eigendecomposition
    G=N_sl^T N_sl + a basis-independent plane-restriction Gram, NO LAPACKE_dgesdd),
    and SolveLoopStar deflates them. Validated field-exact on small converging cubes
    (3.4e-10 '-z' 3^3, 3.1e-10 '+x-z' 3^3, 1.2e-9 '-z' 4^3; C:/temp/antisym_cube_svdfree.py).

NOTE on high-mu_r C-type convergence: the reduced-space DIAGONAL-Jacobi preconditioner
is weak for mu_r=1e5 compact iron, so BOTH '+x' and '+x-z' loop-star may not fully
converge here (cap near max_iter) -- the field is then only approximately exact
(~1e-3). This is a PRECONDITIONER limitation (next: block-Jacobi on the reduced
operator), NOT a gauge/deflation error: on a converging model the gauge is field-exact
to ~1e-9 (see antisym_cube_svdfree.py). Use a lower mu_r or the full model for a
cleanly converging loop-star demo.

This script exercises BOTH: '+x' (symmetric) and '+x-z' (antisymmetric, SVD-free
plane-slab deflation).

Usage: python ctype_loopstar_test.py [density]   (default 6x6x6)
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

MU_R = 1e5


def solve(nodes, elements, loopstar, image):
    rad.SetDeflateNullspace(False, 0.0)
    rad.SetLoopStarGauge(bool(loopstar))
    rad.SolverConfig(bicgstab_tol=1e-8, hacapk_eps=1e-4, hacapk_leaf=10, hacapk_eta=2.0)
    model, yoke, ne = build_model(nodes, elements, ('linear', MU_R))
    t0 = time.perf_counter()
    rad.Solve(model, 0.0001, 100, 2, image=image)
    dt = time.perf_counter() - t0
    st = rad.GetSolveStats()
    B = np.array(rad.Fld(model, 'b', [0, 0, 0]))
    return {'ne': ne, 'it': st.get('linear_iterations', -1), 'dt': dt, 'B': B}


def main():
    dens = sys.argv[1] if len(sys.argv) > 1 else '6x6x6'
    nodes, elements = load_geometry(dens)
    print(f"C-type {dens} ({len(elements)*6} DOF) LINEAR mu_r={MU_R:.0e}, BiCGSTAB tol=1e-8")

    # --- SYMMETRIC IMA '+x': loop-star must be field-exact vs OFF ---
    print("\n[1] SYMMETRIC IMA image='+x' (loop-star supported, must be field-exact)")
    off = solve(nodes, elements, False, '+x')
    ls = solve(nodes, elements, True, '+x')
    dBz = abs(ls['B'][2] - off['B'][2]) / (abs(off['B'][2]) + 1e-30)
    print(f"  off       lin_iter={off['it']:>6}  Bz={off['B'][2]*1000:>11.3f} mT")
    print(f"  loop-star lin_iter={ls['it']:>6}  Bz={ls['B'][2]*1000:>11.3f} mT")
    print(f"  dBz/Bz = {dBz:.2e}  -> PASS if << 1e-3 (field-exact, alpha-free)")

    # --- ANTISYMMETRIC IMA '+x-z': now SUPPORTED via the plane-slab added-mode
    # deflation (ComputePlaneSlabAddedModes); must be field-exact vs OFF ---
    print("\n[2] ANTISYMMETRIC IMA image='+x-z' (plane-slab deflation, must be field-exact)")
    off2 = solve(nodes, elements, False, '+x-z')
    ls2 = solve(nodes, elements, True, '+x-z')
    dBz2 = abs(ls2['B'][2] - off2['B'][2]) / (abs(off2['B'][2]) + 1e-30)
    print(f"  off       lin_iter={off2['it']:>6}  Bz={off2['B'][2]*1000:>11.3f} mT")
    print(f"  loop-star lin_iter={ls2['it']:>6}  Bz={ls2['B'][2]*1000:>11.3f} mT")
    print(f"  dBz/Bz = {dBz2:.2e}  -> PASS if << 1e-3 (antisym IMA, field-exact)")

    rad.SetLoopStarGauge(False)
    rad.UtiDelAll()


if __name__ == '__main__':
    main()
