"""ih_cad_first_demo.py

End-to-end CAD-first IH (induction heating) demo.

Demonstrates the full workflow:
  (1) Design the coil in build123d / CoilBuilder (CAD-first)
  (2) PEEC filament solve on the coil -- NO MESH GENERATION
  (3) phi_inc from PEEC filament currents on workpiece surface
  (4) BEM-SIBC on workpiece -> H_t_rms, P_total, Q_total

Two coil designs are compared on the SAME steel-cylinder workpiece:
  A) Circular wire torus (CircleProfile, r = 3 mm) -- simple round wire
  B) Loft arc: rectangular 5x5 mm straight leg into a round wire 180
     degree return -- showcases Stage 3c LoftArcSegment

Both are at 7 kHz, I_total = 1 A, and run against the same checked
workpiece .vol supplied on the command line.

The BEM-SIBC solver is CACHED across both coils via
IHWorkpieceContext -- see per-coil timing breakdown in the output.
"""

import argparse
from pathlib import Path

import numpy as np

from ngsolve import Mesh, BND

from radia.coil_builder import CoilBuilder
from radia.coil_profile import RectProfile, CircleProfile
from radia.ih_pipeline import IHWorkpieceContext

from radia.panels.surface_mesh_extract import _extract_surface_mesh_filtered

MM = 1e-3

R_MAJOR = 30 * MM
COIL_CURRENT = 1.0
FREQUENCY = 7e3

WP_SIGMA = 2e6
WP_MU_R = 100.0


def load_workpiece(vol_file: Path):
    mesh_full = Mesh(str(vol_file))
    wp = _extract_surface_mesh_filtered(mesh_full, keep_label='workpiece')
    return wp


def build_coil_A():
    """Round wire torus, R_major = 30 mm, r = 3 mm, 355 deg arc."""
    return (CoilBuilder(current=COIL_CURRENT)
            .set_start([R_MAJOR, 0.0, 0.0])
            .set_profile(CircleProfile(r=3 * MM))
            .add_arc(radius=R_MAJOR, arc_angle=355.0))


def build_coil_B():
    """Rect 5x5 mm -> round 3 mm wire LoftArc 355 deg."""
    return (CoilBuilder(current=COIL_CURRENT)
            .set_start([R_MAJOR, 0.0, 0.0])
            .set_profile(RectProfile(5 * MM, 5 * MM))
            .add_loft_arc(profile_end=CircleProfile(r=3 * MM),
                          radius=R_MAJOR, arc_angle=355.0, n_sub=40))


def report(tag, coil, result):
    print(f"\n{'=' * 70}")
    print(f"Coil: {tag}")
    print(f"  segments: {len(coil.segments)}")
    for i, s in enumerate(coil.segments):
        print(f"    [{i}] {type(s).__name__}  "
              f"w={s.width*1e3:.2f} h={s.height*1e3:.2f} mm")

    I = result['I_peec']
    mag = np.abs(I)
    ph = np.angle(I, deg=True)
    print(f"  filaments N = {result['n_fil']}  n_loop = {result['n_loop']}")
    print(f"  |I_k| range [{mag.min():.4f}, {mag.max():.4f}]  "
          f"max/min = {mag.max()/mag.min():.2f}")
    print(f"  arg range [{ph.min():.1f}, {ph.max():.1f}] deg")

    t = result['timings']
    print(f"  timing (per-coil, SIBC solver cached):")
    print(f"    to_filaments         : {t['to_filaments_ms']:7.1f} ms")
    print(f"    build_bundle_solver  : {t['build_bundle_ms']:7.1f} ms")
    print(f"    compute_branch_currents: {t['branch_currents_ms']:7.1f} ms")
    print(f"    compute_phi_inc      : {t['phi_inc_ms']:7.1f} ms")
    print(f"    BEM-SIBC solve       : {t['sibc_solve_ms']:7.1f} ms")
    print(f"    TOTAL                : {t['total_ms']:7.1f} ms")

    print(f"  workpiece response (I = {COIL_CURRENT} A, f = {FREQUENCY:g} Hz):")
    print(f"    H_t_rms = {result['H_t_rms']:.4f}  A/m")
    print(f"    P_total = {result['P_total']:.4e}  W")
    print(f"    Q_total = {result['Q_total']:.4e}  VAr")


def main():
    parser = argparse.ArgumentParser(
        description="Compare two CAD-first coils against a checked workpiece mesh."
    )
    parser.add_argument("vol", type=Path, help="checked workpiece .vol mesh")
    args = parser.parse_args()
    if not args.vol.is_file():
        parser.error(f"workpiece mesh does not exist: {args.vol}")

    print("Loading workpiece mesh...")
    wp = load_workpiece(args.vol)
    print(f"  {wp.nv} vertices, {wp.GetNE(BND)} surface elements")

    print("\nBuilding IHWorkpieceContext (one-time BEM assembly)...")
    ctx = IHWorkpieceContext(wp, frequency=FREQUENCY,
                             sigma=WP_SIGMA, mu_r=WP_MU_R)
    print(f"  SIBC assembly = {ctx.t_build*1e3:.0f} ms  "
          f"(ndof = {ctx.solver.ndof})")
    print(f"  delta = {ctx.delta*1e3:.3f} mm  "
          f"Z_s = {ctx.Z_s:.3e} ohm")

    coil_A = build_coil_A()
    resA = ctx.evaluate(coil_A, nw=5, nh=8)
    report("A: round wire torus (r=3mm)", coil_A, resA)

    coil_B = build_coil_B()
    resB = ctx.evaluate(coil_B, nw=5, nh=5)
    report("B: rect->round loft arc (5x5 mm -> r=3 mm)", coil_B, resB)

    print("\n" + "=" * 70)
    print("Summary (IH workpiece response at 7 kHz, I = 1 A)")
    print("=" * 70)
    fmt = "  {tag:50s}  {H:>10}  {P:>12}  {t:>10}"
    print(fmt.format(tag='coil', H='H_t_rms', P='P_total', t='t_per_coil'))
    for tag, r in [("A: round wire torus (r=3mm)", resA),
                   ("B: rect->round loft arc", resB)]:
        print(fmt.format(
            tag=tag, H=f"{r['H_t_rms']:.4f}",
            P=f"{r['P_total']:.3e}",
            t=f"{r['timings']['total_ms']:.0f} ms"))

    print(f"\nOne-time BEM assembly cost: {ctx.t_build*1e3:.0f} ms")
    print("Per-coil cost above excludes BEM assembly; for sweeps the "
          "per-coil cost is what dominates at scale.")


if __name__ == '__main__':
    main()
