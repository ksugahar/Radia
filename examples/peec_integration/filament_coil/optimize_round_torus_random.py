"""optimize_round_torus_random.py

Baseline random-search optimization over round-wire torus IH coils.

Design space:
  r_wire   in [1, 7]   mm
  R_major  in [28, 55] mm  (outside the 25mm-radius workpiece)
  arc_deg  in [90, 359]

Objective: maximize P_total in the steel cylinder workpiece at 7 kHz
with total coil current I = 1 A.

This is the baseline ClaudeProposer will be compared against: for a
fair head-to-head, run this first to record best-of-N random.
"""

import os
import sys
import time
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.abspath(os.path.join(HERE, '..', '..', '..', 'src'))
SRC_RADIA = os.path.join(SRC, 'radia')
PANELS = os.path.join(SRC, 'radia', 'panels')
for p in (SRC, SRC_RADIA, PANELS, HERE):
    if p not in sys.path:
        sys.path.insert(0, p)

from ngsolve import Mesh

from radia.ih_pipeline import IHWorkpieceContext
from radia.ih_optimize import IHOptimizer, RandomProposer

from radia.panels.surface_mesh_extract import _extract_surface_mesh_filtered


MM = 1e-3
VOL_FILE = os.path.join(SRC, 'radia', 'panels', 'samples',
                        'ih_bem_sample.vol')


def main(n_iter: int = 20, seed: int = 0):
    print(f"Loading workpiece...")
    mesh_full = Mesh(VOL_FILE)
    wp = _extract_surface_mesh_filtered(mesh_full, keep_label='workpiece')
    print(f"  {wp.nv} vertices")

    print("Building IHWorkpieceContext (BEM assembly, ~11s one-time)...")
    t0 = time.perf_counter()
    ctx = IHWorkpieceContext(wp, frequency=7e3, sigma=2e6, mu_r=100.0)
    print(f"  done in {(time.perf_counter() - t0)*1e3:.0f} ms\n")

    proposer = RandomProposer(
        ranges={
            'r_wire_mm': (1.0, 7.0),
            'R_major_mm': (28.0, 55.0),
            'arc_deg': (90.0, 359.0),
        },
        profile_type='circle',
        current=1.0,
        seed=seed,
    )
    opt = IHOptimizer(ctx, proposer, nw=5, nh=8)

    print(f"Random search: n_iter = {n_iter}, seed = {seed}\n")
    opt.run(n_iter=n_iter)

    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)
    print(opt.summary())


if __name__ == '__main__':
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument('--n-iter', type=int, default=20)
    ap.add_argument('--seed', type=int, default=0)
    args = ap.parse_args()
    main(n_iter=args.n_iter, seed=args.seed)
