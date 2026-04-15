"""Matched-mesh 2D axisym SIBC vs 3D FEM-SIBC comparison.

Principle: with gap_deg=0 (full 360 torus), identical mesh sizes, identical
scalar Z_s, and identical Curve order, 2D axisym and 3D FEM-SIBC MUST give
the same L and P. Any residual gap is a formulation bug, not physics.

Knobs matched:
  - gap_deg = 0              (3D only; 2D is inherently axisymmetric)
  - maxh_wp_bnd              (WP surface/boundary mesh)
  - maxh_air                 (air domain)
  - maxh_coil                (coil volume)
  - curve_order              (Curve(3) both)
  - scalar Z_s               (cylindrical Bessel, ESIM)

Usage:
    python compare_2d_vs_3d_matched_mesh.py --freq 7000
    python compare_2d_vs_3d_matched_mesh.py --freq 50000 --maxh-wp-bnd 0.001
"""

import argparse
import math
import os
import sys
import time

sys.path.insert(0, os.path.join(os.path.dirname(__file__),
                                '../cubit_panels/inductance'))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../src/radia'))

SIGMA_MAP = {'copper': 5.8e7, 'aluminum': 3.5e7}
MU_0 = 4e-7 * math.pi


def skin_depth(sigma, mu_r, freq):
    return math.sqrt(2.0 / (2 * math.pi * freq * mu_r * MU_0 * sigma))


def run_2d_sibc(R_coil, a_coil, R_wp, H_wp, sigma, freq, mu_r,
                maxh_air, maxh_coil, maxh_wp_bnd, order, a_kelvin, z_offset):
    from fem_esim_kelvin import solve_fem_sibc
    return solve_fem_sibc(
        R_coil=R_coil, a_coil=a_coil, R_wp=R_wp, H_wp=H_wp,
        sigma=sigma, frequency=freq, mu_r=mu_r,
        a_kelvin=a_kelvin, z_offset=z_offset,
        maxh=maxh_air, order=order,
        max_iter=8, tol=1e-4,
        coil_shape='circular',
        maxh_wp_bnd=maxh_wp_bnd,
        maxh_coil=maxh_coil,
    )


def run_3d_sibc(R_coil, a_coil, R_wp, H_wp, sigma, freq, material,
                gap_deg, maxh_air, maxh_coil, maxh_wp_bnd,
                R_air, R_kelvin, order):
    from fem_esim_3d import run
    return run(
        R_coil=R_coil, a_coil=a_coil, gap_deg=gap_deg,
        R_wp=R_wp, H_wp=H_wp,
        sigma=sigma, frequency=freq, material=material,
        R_air=R_air, maxh_air=maxh_air,
        maxh_coil=maxh_coil,
        maxh_wp_bnd=maxh_wp_bnd,
        R_kelvin=R_kelvin, order=order,
        esim_geometry='cylinder', approach='hole',
    )


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--freq", type=float, default=7000)
    ap.add_argument("--material", default="copper")
    ap.add_argument("--R-coil", type=float, default=0.030)
    ap.add_argument("--a-coil", type=float, default=0.003)
    ap.add_argument("--R-wp", type=float, default=0.025)
    ap.add_argument("--H-wp", type=float, default=0.025)
    ap.add_argument("--gap-deg", type=float, default=0.0,
                    help="3D gap angle (0 = full axisymmetric)")
    ap.add_argument("--maxh-air", type=float, default=0.012,
                    help="Air mesh size [m] (both 2D and 3D)")
    ap.add_argument("--maxh-wp-bnd", type=float, default=0.002,
                    help="WP boundary/surface mesh [m]")
    ap.add_argument("--maxh-coil", type=float, default=0.0005,
                    help="Coil volume mesh [m] (default 0.5mm = 12 elem/6mm dia)")
    ap.add_argument("--order", type=int, default=1,
                    help="HCurl order for 3D, H1 order is fixed at 3 for 2D")
    ap.add_argument("--R-air", type=float, default=0.060,
                    help="3D air sphere radius [m]")
    ap.add_argument("--R-kelvin", type=float, default=0.120,
                    help="3D Kelvin exterior domain outer radius [m]")
    ap.add_argument("--skip-2d", action="store_true")
    ap.add_argument("--skip-3d", action="store_true")
    args = ap.parse_args()

    sigma = SIGMA_MAP[args.material]
    mu_r = 1.0
    delta = skin_depth(sigma, mu_r, args.freq)

    print()
    print("=" * 92)
    print(f"  2D axisym SIBC vs 3D FEM-SIBC   --   MATCHED MESH")
    print("=" * 92)
    print(f"  {args.material} @ {args.freq} Hz")
    print(f"  delta = {delta*1e3:.3f} mm, R/delta = {args.R_wp/delta:.1f}")
    print(f"  Coil: R={args.R_coil*1e3:.1f}mm a={args.a_coil*1e3:.1f}mm "
          f"gap={args.gap_deg}deg (3D)")
    print(f"  WP:   R={args.R_wp*1e3:.1f}mm H={args.H_wp*1e3:.1f}mm")
    print(f"  Mesh: maxh_air={args.maxh_air*1e3:.1f}mm  "
          f"maxh_coil={args.maxh_coil*1e3:.2f}mm  "
          f"maxh_wp_bnd={args.maxh_wp_bnd*1e3:.2f}mm")
    print(f"  Curve(3) both.  Z_s: cylindrical Bessel (ESIM).")
    print("=" * 92)

    # For 2D: Kelvin radius ~3*R_coil, z_offset = 2.5*a_kelvin
    a_kelvin_2d = max(args.R_coil * 4, args.H_wp / 2 + args.R_coil * 2, 0.120)
    z_offset_2d = a_kelvin_2d * 2.5

    out = {}
    if not args.skip_2d:
        print("\n[a] 2D axisym SIBC (circular coil, WP = hole, Robin BC)")
        t0 = time.perf_counter()
        try:
            r2d = run_2d_sibc(
                R_coil=args.R_coil, a_coil=args.a_coil,
                R_wp=args.R_wp, H_wp=args.H_wp,
                sigma=sigma, freq=args.freq, mu_r=mu_r,
                maxh_air=args.maxh_air, maxh_coil=args.maxh_coil,
                maxh_wp_bnd=args.maxh_wp_bnd, order=3,
                a_kelvin=a_kelvin_2d, z_offset=z_offset_2d,
            )
            dt = time.perf_counter() - t0
            print(f"  wall: {dt:.1f}s")
            out['2d'] = r2d
        except Exception as e:
            import traceback; traceback.print_exc()
            out['2d'] = None

    if not args.skip_3d:
        print(f"\n[b] 3D FEM-SIBC (full 360 torus, gap_deg={args.gap_deg})")
        t0 = time.perf_counter()
        try:
            r3d = run_3d_sibc(
                R_coil=args.R_coil, a_coil=args.a_coil,
                R_wp=args.R_wp, H_wp=args.H_wp,
                sigma=sigma, freq=args.freq, material=args.material,
                gap_deg=args.gap_deg,
                maxh_air=args.maxh_air, maxh_coil=args.maxh_coil,
                maxh_wp_bnd=args.maxh_wp_bnd,
                R_air=args.R_air, R_kelvin=args.R_kelvin, order=args.order,
            )
            dt = time.perf_counter() - t0
            print(f"  wall: {dt:.1f}s")
            out['3d'] = r3d
        except Exception as e:
            import traceback; traceback.print_exc()
            out['3d'] = None

    print()
    print("-" * 80)
    print("  SUMMARY (matched mesh)")
    print("-" * 80)
    print(f"  {'solver':>14s}  {'L [nH]':>10s}  {'P [W]':>12s}  "
          f"{'ndof':>8s}  {'t_solve':>8s}")
    for key, label in [('2d', '2D SIBC'), ('3d', '3D FEM-SIBC')]:
        r = out.get(key)
        if r is None:
            print(f"  {label:>14s}  {'--':>10s}  {'--':>12s}  "
                  f"{'--':>8s}  {'--':>8s}")
            continue
        L = r['L'] * 1e9
        P = r.get('P_total', 0.0)
        ndof = r.get('ndof', -1)
        t_solve = r.get('t_solve', 0.0)
        print(f"  {label:>14s}  {L:10.4f}  {P:12.4e}  "
              f"{ndof:>8d}  {t_solve:>7.1f}s")

    if out.get('2d') and out.get('3d'):
        L2 = out['2d']['L'] * 1e9
        L3 = out['3d']['L'] * 1e9
        P2 = out['2d'].get('P_total', 0.0)
        P3 = out['3d'].get('P_total', 0.0)
        print()
        print(f"  dL (3D - 2D)/2D = {(L3-L2)/L2*100:+.3f}%")
        print(f"  dP (3D - 2D)/2D = {(P3-P2)/P2*100:+.3f}%")


if __name__ == "__main__":
    main()
