#!/usr/bin/env python
"""Resolve the nonlinear deflation B_z DRIFT so the field is accelerator-grade.

Context: on the nonlinear C-type, loop deflation is 4-5x faster but B_z drifts
~3% (off -963.7 -> alpha=0.03 -934.2). Mechanism (derived): with non-uniform chi
(partial saturation) D=diag(1/chi) is not scalar, so span(L)=ker N is NOT an
A-invariant eigenspace; the symmetric shift A+alpha L L^T then leaks into the
PHYSICAL solution, biasing the Picard fixed point by O(alpha). (Linear/uniform chi
shows no drift -- confirmed 0.01%.) A 3% field error is unacceptable for
accelerator magnets, so we need field-exact + fast.

This script measures B_z and cost for:
  (1) off                : true reference (slow)
  (2) auto 1e-3          : deflation, loose tol  -> expect drift
  (3) auto 1e-4          : deflation, tight tol  -> if drift persists => H2 (bias),
                           if it recovers => H1 (convergence path)
  (4) two-phase          : auto(1e-3) then OFF(1e-3) warm-start polish. The
                           undeflated Picard fixed point is the TRUE physics, so
                           the polish should pull B_z back to the off value while
                           the bulk speed came from phase 1.

PASS = a scenario whose final B_z matches (1) within ~0.1% AND is much faster than (1).

Usage: python ctype_drift_fix.py [density]   (default 6x6x6)
"""
import os
import sys
import time

import numpy as np

work_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, work_dir)
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(work_dir))), 'src'))

import radia as rad
from deflation_benchmark import load_geometry, build_model, load_elf_bz

IMAGE = '+x-z'
MAXNL = 200
DEFL_ALPHA = 0.03   # known-good manual shift (auto-alpha tuning is separate)


def solve(model, deflate_alpha, tol):
    """One Solve phase on an existing model (M persists for warm-start)."""
    if deflate_alpha is None:
        rad.SetDeflateNullspace(False, 0.0)
    else:
        rad.SetDeflateNullspace(True, deflate_alpha)
    rad.SolverConfig(bicgstab_tol=1e-8, hacapk_eps=1e-4, hacapk_leaf=10, hacapk_eta=2.0)
    t0 = time.perf_counter()
    conv = True
    try:
        rad.Solve(model, tol, MAXNL, 2, image=IMAGE)
    except RuntimeError:
        conv = False
    dt = time.perf_counter() - t0
    st = rad.GetSolveStats()
    return dt, st.get('linear_iterations', -1), st.get('deflation_alpha', float('nan')), conv


def main():
    dens = sys.argv[1] if len(sys.argv) > 1 else '6x6x6'
    nodes, elements = load_geometry(dens)
    elf_bz = load_elf_bz(dens)
    mat = ('nonlinear',)
    print(f"C-type {dens} ({len(elements)*6} DOF) nonlinear, image={IMAGE}")
    print(f"ELF B_z = {elf_bz*1000:.2f} mT\n")

    rows = []

    # (1) off reference
    m, _, _ = build_model(nodes, elements, mat)
    dt, it, _, ok = solve(m, None, 1e-3)
    bz_ref = float(np.array(rad.Fld(m, 'b', [0, 0, 0]))[2]) * 1000
    rows.append(('off 1e-3 (ref)', it, dt, bz_ref, 0.0, ok))

    # (2) auto, loose
    m, _, _ = build_model(nodes, elements, mat)
    dt, it, a, ok = solve(m, DEFL_ALPHA, 1e-3)
    bz = float(np.array(rad.Fld(m, 'b', [0, 0, 0]))[2]) * 1000
    rows.append((f'defl 1e-3 (a={a:.3g})', it, dt, bz, 100*(bz-bz_ref)/bz_ref, ok))

    # (3) auto, tight
    m, _, _ = build_model(nodes, elements, mat)
    dt, it, a, ok = solve(m, DEFL_ALPHA, 1e-4)
    bz = float(np.array(rad.Fld(m, 'b', [0, 0, 0]))[2]) * 1000
    rows.append((f'defl 1e-4 (a={a:.3g})', it, dt, bz, 100*(bz-bz_ref)/bz_ref, ok))

    # (4) two-phase: auto(1e-3) then OFF(1e-3) warm polish
    m, _, _ = build_model(nodes, elements, mat)
    dt1, it1, a, ok1 = solve(m, -1.0, 1e-3)       # phase 1 (deflated, fast)
    bz_p1 = float(np.array(rad.Fld(m, 'b', [0, 0, 0]))[2]) * 1000
    dt2, it2, _, ok2 = solve(m, None, 1e-3)        # phase 2 (off, warm polish)
    bz = float(np.array(rad.Fld(m, 'b', [0, 0, 0]))[2]) * 1000
    rows.append((f'two-phase (p1 bz={bz_p1:.1f})', it1+it2, dt1+dt2, bz,
                 100*(bz-bz_ref)/bz_ref, ok1 and ok2))

    rad.SetDeflateNullspace(False, 0.0); rad.UtiDelAll()

    print(f"{'scenario':<26} {'lin_it':>7} {'t(s)':>7} {'Bz(mT)':>10} {'dBz%':>8} {'conv':>5}")
    for name, it, dt, bz, d, ok in rows:
        print(f"{name:<26} {it:>7} {dt:>7.1f} {bz:>10.2f} {d:>+8.3f} {'y' if ok else 'NO':>5}")
    print("\nPASS = |dBz%| < 0.1 AND much faster than the off reference.")


if __name__ == '__main__':
    main()
