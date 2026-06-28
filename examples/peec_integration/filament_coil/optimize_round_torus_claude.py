"""optimize_round_torus_claude.py

LLM-driven (Claude Opus 4.6) coil optimization over round-wire torus
IH coils, head-to-head against the random-search baseline.

Uses the ``claude`` CLI backend -- subscription auth, no
ANTHROPIC_API_KEY required. See src/radia/ih_claude_proposer.py.

Usage:
    python optimize_round_torus_claude.py --n-iter 10

Same design space, objective, and workpiece as
optimize_round_torus_random.py, so the results are directly
comparable.
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
from radia.ih_claude_proposer import ClaudeProposer

from radia.panels.surface_mesh_extract import _extract_surface_mesh_filtered


MM = 1e-3
VOL_FILE = os.path.join(SRC, 'radia', 'panels', 'samples',
                        'ih_bem_sample.vol')


RANGES = {
    'r_wire_mm': (1.0, 7.0),
    'R_major_mm': (28.0, 55.0),
    'arc_deg': (90.0, 359.0),
}

WORKPIECE_DESC = (
    "Solid steel cylinder, outer radius 25 mm, height 25 mm, "
    "conductivity sigma = 2e6 S/m, relative permeability mu_r = 100. "
    "Centered at the origin with its axis along +Z. The coil torus must "
    "surround this cylinder (R_major > 25 mm + r_wire).")


def main(n_iter: int = 10, seed: int = 0):
    print("Loading workpiece...")
    wp = _extract_surface_mesh_filtered(Mesh(VOL_FILE),
                                         keep_label='workpiece')
    print(f"  {wp.nv} vertices")

    print("Building IHWorkpieceContext (BEM assembly, ~11 s)...")
    t0 = time.perf_counter()
    ctx = IHWorkpieceContext(wp, frequency=7e3, sigma=2e6, mu_r=100.0)
    print(f"  done in {(time.perf_counter() - t0)*1e3:.0f} ms")

    fallback = RandomProposer(RANGES, profile_type='circle',
                               current=1.0, seed=seed)
    proposer = ClaudeProposer(
        ranges=RANGES,
        workpiece_desc=WORKPIECE_DESC,
        frequency_hz=7e3,
        fallback_proposer=fallback,
        verbose=True,
    )
    opt = IHOptimizer(ctx, proposer, nw=5, nh=8)

    print(f"\nClaude-driven search: n_iter = {n_iter}\n")
    opt.run(n_iter=n_iter)

    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)
    print(opt.summary())
    print("\nClaude telemetry:")
    for k, v in proposer.telemetry().items():
        print(f"  {k}: {v}")


if __name__ == '__main__':
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument('--n-iter', type=int, default=10)
    ap.add_argument('--seed', type=int, default=0,
                     help='Seed for the RandomProposer fallback.')
    args = ap.parse_args()
    main(n_iter=args.n_iter, seed=args.seed)
