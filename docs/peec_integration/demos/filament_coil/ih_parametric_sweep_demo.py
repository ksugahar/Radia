"""ih_parametric_sweep_demo.py

Demonstrates the per-coil cost win of caching the BEM-SIBC solver via
IHWorkpieceContext. Sweeps a set of round-wire torus coils over a
grid of (wire_radius, n_turns_arc) values against one cached workpiece
and reports per-coil time + H_t_rms / P_total.

Baseline (no cache, rebuilding SIBC each coil): ~14 s / coil on the
original 166-node reference workpiece.
With cache (this script): ~3 s / coil after a single ~11 s one-time
assembly.
"""

import argparse
from pathlib import Path

import time
import numpy as np

from ngsolve import Mesh, BND

from radia.coil_builder import CoilBuilder
from radia.coil_profile import CircleProfile
from radia.ih_pipeline import IHWorkpieceContext

from radia.panels.surface_mesh_extract import _extract_surface_mesh_filtered

MM = 1e-3

R_MAJOR = 30 * MM
FREQUENCY = 7e3
COIL_CURRENT = 1.0
WP_SIGMA = 2e6
WP_MU_R = 100.0


def build_round_torus(wire_r_mm, arc_deg):
    return (CoilBuilder(current=COIL_CURRENT)
            .set_start([R_MAJOR, 0.0, 0.0])
            .set_profile(CircleProfile(r=wire_r_mm * MM))
            .add_arc(radius=R_MAJOR, arc_angle=arc_deg))


def main():
    parser = argparse.ArgumentParser(
        description="Sweep CAD-first coils against a checked workpiece mesh."
    )
    parser.add_argument("vol", type=Path, help="checked workpiece .vol mesh")
    args = parser.parse_args()
    if not args.vol.is_file():
        parser.error(f"workpiece mesh does not exist: {args.vol}")

    print("Loading workpiece mesh...")
    mesh_full = Mesh(str(args.vol))
    wp = _extract_surface_mesh_filtered(mesh_full, keep_label='workpiece')
    print(f"  {wp.nv} vertices, {wp.GetNE(BND)} surface elements")

    print("\nBuilding IHWorkpieceContext (one-time BEM assembly)...")
    t0 = time.perf_counter()
    ctx = IHWorkpieceContext(wp, frequency=FREQUENCY,
                             sigma=WP_SIGMA, mu_r=WP_MU_R)
    t_ctx = time.perf_counter() - t0
    print(f"  SIBC assembly = {t_ctx*1e3:.0f} ms (ndof = {ctx.solver.ndof})")

    # Parametric grid
    wire_r_mm = [2.0, 2.5, 3.0, 3.5, 4.0]
    arc_deg = [355.0, 300.0]
    combos = [(r, a) for r in wire_r_mm for a in arc_deg]
    print(f"\nSweeping {len(combos)} coil designs (cached SIBC)...\n")

    header = ("  {r:>8}  {a:>8}  {nfil:>5}  {Hrms:>10}  "
              "{P:>12}  {tper:>10}")
    print(header.format(r='r[mm]', a='arc[deg]', nfil='N_fil',
                        Hrms='H_t_rms', P='P_total', tper='t_coil'))

    t_sweep_start = time.perf_counter()
    results = []
    for (r_mm, a_deg) in combos:
        coil = build_round_torus(r_mm, a_deg)
        res = ctx.evaluate(coil, nw=5, nh=8)
        results.append((r_mm, a_deg, res))
        print(header.format(
            r=f"{r_mm:.1f}", a=f"{a_deg:.0f}",
            nfil=f"{res['n_fil']}",
            Hrms=f"{res['H_t_rms']:.4f}",
            P=f"{res['P_total']:.3e}",
            tper=f"{res['timings']['total_ms']:.0f} ms"))
    t_sweep = time.perf_counter() - t_sweep_start

    # Aggregate timing
    per_coil = [r[2]['timings']['total_ms'] for r in results]
    mean_per_coil = np.mean(per_coil)
    total_without_cache = len(combos) * (t_ctx * 1e3 + mean_per_coil)
    total_with_cache = t_ctx * 1e3 + sum(per_coil)

    print("\n" + "=" * 70)
    print("Timing summary")
    print("=" * 70)
    print(f"  one-time BEM assembly       : {t_ctx*1e3:8.0f} ms")
    print(f"  mean per-coil (cached)      : {mean_per_coil:8.0f} ms")
    print(f"  sweep total ({len(combos):2d} coils)        : "
          f"{t_sweep*1e3:8.0f} ms")
    print(f"  projected without cache     : "
          f"{total_without_cache:8.0f} ms  "
          f"({total_without_cache/total_with_cache:.1f}x slower)")


if __name__ == '__main__':
    main()
