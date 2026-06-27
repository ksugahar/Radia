"""Compare 2D axisymmetric SIBC vs 3D FEM-SIBC on the same problem.

Geometry:
  Coil:      ring at r=30 mm, cross-section radius a=3 mm
             (2D: square cross-section 2a x 2a; 3D: circular torus)
  Workpiece: cylinder R=25 mm, H=25 mm (axisymmetric, no gap)
  3D gap:    default 5 deg (same as ih_fem_kelvin_sample.jou)

Solvers:
  (a) 2D full-res  : fem_esim_kelvin.solve_fem_full
                     workpiece meshed, direct eddy currents
                     -- REFERENCE (validated <0.6% vs Dodd-Deeds)
  (b) 2D SIBC      : fem_esim_kelvin.solve_fem_sibc (Karl iteration)
                     workpiece = hole, Robin BC with cylindrical Bessel Z_s
                     -- for linear Cu/Al this converges in ~1 iter to Bessel
  (c) 3D FEM-SIBC  : fem_esim_3d.run (hole approach + Kelvin + HCurl)
                     gapped torus coil + cylinder hole + Robin BC
                     -- the solver under test

All three use Kelvin open-boundary (no Dirichlet truncation).

The 2D side is known to match full-res within 1% (see
eddy_current_analytical_validation/README.md). Any 3D-vs-2D discrepancy
is attributable to:
  - The 5-deg gap in 3D (removes ~1.4% of L)
  - 3D mesh coarseness relative to the skin depth
  - The Robin BC sign/scaling in HCurl (potential formulation bug)
  - High-order curving on the torus/sphere

Usage:
    python compare_sibc_2d_vs_3d.py                      # Cu @ 7 kHz
    python compare_sibc_2d_vs_3d.py --material aluminum --freq 1000
    python compare_sibc_2d_vs_3d.py --sweep              # Cu/Al x 1k,7k,50k Hz
    python compare_sibc_2d_vs_3d.py --skip-full          # skip expensive full-res
"""

import argparse
import math
import os
import sys
import time

sys.path.insert(0, os.path.join(os.path.dirname(__file__),
                                '../cubit_panels/inductance'))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../src/radia'))

SIGMA_MAP = {'copper': 5.8e7, 'aluminum': 3.5e7, 'steel': 2e6}
MU_R_MAP = {'copper': 1.0, 'aluminum': 1.0, 'steel': None}  # steel uses BH
MU_0 = 4e-7 * math.pi


def skin_depth(sigma, mu_r, freq):
    mu_r_eff = mu_r if mu_r else 100.0
    return math.sqrt(2.0 / (2 * math.pi * freq * mu_r_eff * MU_0 * sigma))


def run_2d_full(R_coil, a_coil, R_wp, H_wp, sigma, freq,
                maxh=0.004, order=3, a_kelvin=0.12, z_offset=0.40,
                coil_shape='circular'):
    """2D axisym full-res (workpiece meshed).

    Uses reference_2d_axisym.solve_2d_axisym which auto-refines WP
    mesh to delta/5 for accurate skin-effect resolution. The
    fem_esim_kelvin.solve_fem_full alternative uses a coarser WP
    mesh (maxh/3) that under-resolves skin depth at high frequency,
    causing ~8% L under-estimation.
    """
    from em_material import EMMaterial
    from reference_2d_axisym import solve_2d_axisym
    # Map sigma -> EMMaterial (only sigma and mu_r matter here)
    # reference_2d_axisym takes a mat object; fake one with our sigma.
    mu_r = 1.0  # linear materials (Cu/Al) only for now
    mat = EMMaterial(name='custom', sigma=sigma, mu_r=mu_r)
    return solve_2d_axisym(r_coil=R_coil, a_coil=a_coil,
                           r_wp=R_wp, h_wp=H_wp,
                           mat=mat, frequency=freq,
                           I_total=1.0, maxh_air=0.010,
                           kelvin=True, order=order)


def run_2d_sibc(R_coil, a_coil, R_wp, H_wp, sigma, freq, mu_r=1.0,
                maxh=0.004, order=3, a_kelvin=0.12, z_offset=0.40,
                max_iter=8, tol=1e-4, coil_shape='circular'):
    """2D axisym SIBC (Karl iteration, WP = hole + Robin).

    coil_shape='circular' matches the 3D circular torus cross-section.
    """
    from fem_esim_kelvin import solve_fem_sibc
    return solve_fem_sibc(R_coil, a_coil, R_wp, H_wp, sigma, freq,
                          bh_curve=None, mu_r=mu_r,
                          a_kelvin=a_kelvin, z_offset=z_offset,
                          maxh=maxh, order=order,
                          max_iter=max_iter, tol=tol,
                          coil_shape=coil_shape)


def run_3d_sibc(R_coil, a_coil, R_wp, H_wp, sigma, freq, material,
                gap_deg=5.0, R_air=0.060, R_kelvin=0.120,
                maxh_air=0.012, maxh_coil=0.005, order=1,
                esim_geometry='cylinder'):
    """3D FEM-SIBC (gapped torus, circular cross-section, WP hole + Kelvin)."""
    from fem_esim_3d import run
    return run(R_coil=R_coil, a_coil=a_coil, gap_deg=gap_deg,
               R_wp=R_wp, H_wp=H_wp,
               sigma=sigma, frequency=freq, material=material,
               R_air=R_air, maxh_air=maxh_air, maxh_coil=maxh_coil,
               R_kelvin=R_kelvin, order=order,
               esim_geometry=esim_geometry, approach='hole')


def print_result(tag, r):
    """Uniform one-line summary for a solver result."""
    if r is None:
        print(f"  {tag:>16s}: FAILED")
        return
    L_nH = r['L'] * 1e9
    P = r.get('P_total', float('nan'))
    H_t = r.get('H_t_rms', r.get('H_t_avg', float('nan')))
    ndof = r.get('ndof', -1)
    t = r.get('t_solve', 0)
    print(f"  {tag:>16s}: L={L_nH:8.3f} nH  P={P:10.3e} W  "
          f"H_t={H_t:8.2f} A/m  ndof={ndof:>7d}  t_solve={t:5.1f}s")


def compare_one(material, freq, R_coil, a_coil, R_wp, H_wp,
                gap_deg, maxh_air_3d, order_3d, skip_full=False,
                skip_3d=False):
    sigma = SIGMA_MAP[material]
    mu_r = MU_R_MAP[material]
    if mu_r is None:
        # Steel: 2D SIBC needs nonlinear ESIM too (and currently
        # diverges badly against analytical) — skip it for first pass.
        print(f"\n=== {material} @ {freq:g} Hz  (steel: nonlinear, SKIPPED) ===")
        return None
    delta = skin_depth(sigma, mu_r, freq)
    print()
    print("=" * 80)
    print(f"  Case: material={material}  freq={freq:g} Hz  "
          f"delta={delta*1e3:.3f} mm  R/delta={R_wp/delta:.1f}")
    print("=" * 80)

    out = {}

    # (a) 2D full-resolution reference
    if not skip_full:
        print("\n[a] 2D axisym FULL-RES (square coil, workpiece meshed)")
        try:
            t0 = time.perf_counter()
            r_full = run_2d_full(R_coil, a_coil, R_wp, H_wp, sigma, freq)
            print(f"  wall time: {time.perf_counter()-t0:.1f}s")
            out['full'] = r_full
            print_result("2D full-res", r_full)
        except Exception as e:
            import traceback
            traceback.print_exc()
            out['full'] = None
    else:
        out['full'] = None

    # (b) 2D SIBC (Karl; converges in ~1 iter for linear mu_r=1)
    print("\n[b] 2D axisym SIBC (square coil, Karl iter, Bessel Z_s)")
    try:
        t0 = time.perf_counter()
        r_sibc2d = run_2d_sibc(R_coil, a_coil, R_wp, H_wp, sigma, freq,
                               mu_r=mu_r)
        print(f"  wall time: {time.perf_counter()-t0:.1f}s")
        out['sibc2d'] = r_sibc2d
        print_result("2D SIBC", r_sibc2d)
    except Exception as e:
        import traceback
        traceback.print_exc()
        out['sibc2d'] = None

    # (c) 3D FEM-SIBC
    print(f"\n[c] 3D FEM-SIBC (circular torus, gap={gap_deg} deg, "
          f"maxh_air={maxh_air_3d*1e3:.1f}mm, order={order_3d})")
    try:
        t0 = time.perf_counter()
        r_3d = run_3d_sibc(R_coil, a_coil, R_wp, H_wp, sigma, freq, material,
                           gap_deg=gap_deg, maxh_air=maxh_air_3d,
                           order=order_3d)
        print(f"  wall time: {time.perf_counter()-t0:.1f}s")
        out['sibc3d'] = r_3d
        print_result("3D FEM-SIBC", r_3d)
    except Exception as e:
        import traceback
        traceback.print_exc()
        out['sibc3d'] = None

    # Comparison table
    print()
    print("-" * 80)
    print(f"  SUMMARY  ({material} @ {freq:g} Hz)")
    print("-" * 80)

    ref = out.get('full') or out.get('sibc2d')
    if ref is None:
        print("  No reference available.")
        return out
    ref_name = "2D full-res" if out.get('full') else "2D SIBC"
    L_ref = ref['L'] * 1e9
    P_ref = ref.get('P_total', 0.0)

    header = f"  {'solver':>16s}  {'L [nH]':>9s}  {'dL %':>7s}  " \
             f"{'P [W]':>12s}  {'dP %':>7s}"
    print(header)
    for key, label in [('full', '2D full-res'),
                       ('sibc2d', '2D SIBC'),
                       ('sibc3d', '3D FEM-SIBC')]:
        r = out.get(key)
        if r is None:
            print(f"  {label:>16s}  {'--':>9s}  {'--':>7s}  "
                  f"{'--':>12s}  {'--':>7s}")
            continue
        L = r['L'] * 1e9
        P = r.get('P_total', 0.0)
        dL = (L - L_ref) / L_ref * 100 if L_ref else 0
        dP = (P - P_ref) / P_ref * 100 if P_ref else 0
        print(f"  {label:>16s}  {L:9.3f}  {dL:+7.2f}  "
              f"{P:12.4e}  {dP:+7.2f}")
    print(f"  (errors are vs {ref_name})")
    return out


def main():
    ap = argparse.ArgumentParser(
        description="Compare 2D axisym SIBC vs 3D FEM-SIBC")
    ap.add_argument("--material", default="copper",
                    choices=list(SIGMA_MAP.keys()))
    ap.add_argument("--freq", type=float, default=7000)
    ap.add_argument("--sweep", action="store_true",
                    help="Cu/Al x 1k,7k,50k Hz")
    ap.add_argument("--skip-full", action="store_true",
                    help="Skip 2D full-res (slow, already validated)")
    ap.add_argument("--R-coil", type=float, default=0.030)
    ap.add_argument("--a-coil", type=float, default=0.003)
    ap.add_argument("--R-wp", type=float, default=0.025)
    ap.add_argument("--H-wp", type=float, default=0.025)
    ap.add_argument("--gap-deg", type=float, default=5.0,
                    help="3D gap angle [deg]")
    ap.add_argument("--maxh-air-3d", type=float, default=0.012,
                    help="3D air mesh size [m]")
    ap.add_argument("--order-3d", type=int, default=1,
                    help="3D HCurl polynomial order")
    args = ap.parse_args()

    common = dict(R_coil=args.R_coil, a_coil=args.a_coil,
                  R_wp=args.R_wp, H_wp=args.H_wp,
                  gap_deg=args.gap_deg, maxh_air_3d=args.maxh_air_3d,
                  order_3d=args.order_3d, skip_full=args.skip_full)

    print()
    print("=" * 80)
    print("  2D axisym SIBC vs 3D FEM-SIBC -- cross-validation")
    print(f"  Coil:      R_coil={args.R_coil*1e3:.1f}mm, "
          f"a_coil={args.a_coil*1e3:.1f}mm, gap_deg={args.gap_deg}")
    print(f"  Workpiece: R_wp={args.R_wp*1e3:.1f}mm, H_wp={args.H_wp*1e3:.1f}mm")
    print(f"  3D mesh:   maxh_air={args.maxh_air_3d*1e3:.1f}mm, "
          f"order={args.order_3d}")
    print("=" * 80)

    if args.sweep:
        for material in ['copper', 'aluminum']:
            for freq in [1000, 7000, 50000]:
                compare_one(material, freq, **common)
    else:
        compare_one(args.material, args.freq, **common)


if __name__ == "__main__":
    main()
