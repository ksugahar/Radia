"""2D axisym full-res convergence study vs Dodd-Deeds infinite rod.

Goal: prove that reference_2d_axisym.solve_2d_axisym is mesh-converged.
Two refinement strategies:

  (a) uniform  : maxh_wp = delta/skin_ratio applied to whole WP face
  (b) skin     : maxh_wp = delta/skin_ratio on WP outer edges only,
                 interior coarser (maxh_wp_interior = 3*delta),
                 grading controls transition rate

Both should converge to the same value when skin_ratio grows. We compare
L, P against Dodd-Deeds analytical (infinite rod, H_wp large).

Usage:
    python convergence_2d_full_res.py                # Cu at 50kHz + 7kHz
    python convergence_2d_full_res.py --freq 50000
"""

import argparse
import math
import os
import sys
import time

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../src/radia'))

from em_material import EMMaterial, MU_0
from reference_2d_axisym import solve_2d_axisym
from dodd_deeds_rod import dodd_deeds_rod


def run_convergence(freq, material='copper',
                    r_coil=0.030, a_coil=0.003,
                    r_wp=0.025, h_wp=0.250,
                    order=3, maxh_air=0.010,
                    skin_ratios_uniform=(3.0, 5.0, 8.0),
                    skin_ratios_skin=(3.0, 5.0, 8.0)):
    """Run convergence study at one (material, freq)."""
    sigma_map = {'copper': 5.8e7, 'aluminum': 3.5e7}
    mu_r_map = {'copper': 1.0, 'aluminum': 1.0}
    sigma = sigma_map[material]
    mu_r = mu_r_map[material]

    mat = EMMaterial(name=material, sigma=sigma, mu_r=mu_r)
    delta = mat.skin_depth(freq)
    r_over_d = r_wp / delta

    print()
    print("=" * 92)
    print(f" Convergence: {material} @ {freq} Hz, "
          f"delta={delta*1e3:.3f} mm, R/delta={r_over_d:.1f}")
    print(f" Geometry:    R_coil={r_coil*1e3:.1f}mm, a_coil={a_coil*1e3:.1f}mm, "
          f"R_wp={r_wp*1e3:.1f}mm, H_wp={h_wp*1e3:.1f}mm")
    print(f" H_wp/R_wp = {h_wp/r_wp:.1f}  (>>1: near-infinite rod approx)")
    print("=" * 92)

    # Analytical reference
    dd = dodd_deeds_rod(r_coil, r_wp, sigma, mu_r, freq)
    L_DD = dd['Delta_L'] * 1e9   # nH (note: Delta_L is conductor effect;
                                 # full-res gives TOTAL L including air)

    # For comparison we use full L (air + Delta_L). Compute air-only L once
    # at very low freq using a dry run? Simpler: use the highest skin_ratio
    # result as the "reference" and report drift. Also report DeltaL vs DD.
    print(f"\n Dodd-Deeds Delta_L = {L_DD:+.4f} nH  "
          f"Delta_R = {dd['Delta_R']*1e3:.4f} mOhm  "
          f"P/I^2 = {dd['P_per_I2']:.4e} W/A^2")

    results = []

    header = (f"\n {'mode':>8s} {'ratio':>6s} {'maxh_wp':>9s} {'maxh_int':>9s} "
              f"{'ndof':>8s} {'t_solve':>8s} {'L [nH]':>10s} "
              f"{'dL_DD %':>9s} {'P [W]':>12s} {'dP_DD %':>9s}")
    print(header)
    print(" " + "-" * (len(header) - 2))

    def _run(mode, ratio):
        maxh_wp = delta / ratio
        maxh_int = max(delta * 3, maxh_wp)
        t0 = time.perf_counter()
        r = solve_2d_axisym(
            r_coil=r_coil, a_coil=a_coil,
            r_wp=r_wp, h_wp=h_wp,
            mat=mat, frequency=freq, I_total=1.0,
            maxh_air=maxh_air, kelvin=True, order=order,
            skin_mode=mode, skin_ratio=ratio,
            maxh_wp_interior=maxh_int,
        )
        t_solve = time.perf_counter() - t0

        L_nH = r['L'] * 1e9
        P = r.get('P_total', 0.0)

        # Air-only contribution is NOT removed here, so L_nH includes
        # self-inductance of the coil loop in air. To compare vs DD
        # Delta_L we'd need L_air too. For now, report absolute L and
        # flag convergence by delta vs finest mesh.
        P_DD = dd['P_per_I2']  # W/A^2, I_total=1 A -> W
        dP = (P - P_DD) / P_DD * 100 if P_DD else 0

        results.append({
            'mode': mode, 'ratio': ratio,
            'maxh_wp': maxh_wp, 'maxh_int': maxh_int,
            'ndof': r.get('ndof', -1), 't_solve': t_solve,
            'L_nH': L_nH, 'P': P, 'dP_DD': dP,
        })
        print(f" {mode:>8s} {ratio:>6.1f} {maxh_wp*1e3:>7.3f}mm "
              f"{maxh_int*1e3:>7.3f}mm {r.get('ndof',-1):>8d} "
              f"{t_solve:>7.1f}s {L_nH:>10.4f} "
              f"{'----':>9s} {P:>12.4e} {dP:>+8.2f}%")

    for ratio in skin_ratios_uniform:
        try:
            _run("uniform", ratio)
        except Exception as e:
            print(f" {'uniform':>8s} {ratio:>6.1f}  FAILED: {e}")

    for ratio in skin_ratios_skin:
        try:
            _run("skin", ratio)
        except Exception as e:
            print(f" {'skin':>8s} {ratio:>6.1f}  FAILED: {e}")

    # Compute dL relative to finest uniform result
    ref = None
    for r in results:
        if r['mode'] == 'uniform' and (ref is None or r['ratio'] > ref['ratio']):
            ref = r
    if ref is not None:
        L_ref = ref['L_nH']
        print(f"\n Reference (finest uniform): ratio={ref['ratio']}  "
              f"L_ref={L_ref:.4f} nH")
        print(f" {'mode':>8s} {'ratio':>6s} {'L [nH]':>10s} {'dL_ref %':>10s} "
              f"{'P [W]':>12s} {'dP_DD %':>9s}")
        for r in results:
            dL_ref = (r['L_nH'] - L_ref) / L_ref * 100
            print(f" {r['mode']:>8s} {r['ratio']:>6.1f} "
                  f"{r['L_nH']:>10.4f} {dL_ref:>+9.3f}% "
                  f"{r['P']:>12.4e} {r['dP_DD']:>+8.2f}%")
    return results


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--freq", type=float, default=None,
                    help="Single frequency (Hz); omit to run 7k+50k")
    ap.add_argument("--material", default="copper")
    ap.add_argument("--r-wp", type=float, default=0.025)
    ap.add_argument("--h-wp", type=float, default=0.250,
                    help="Long WP to approximate infinite rod (default 250mm)")
    ap.add_argument("--order", type=int, default=3)
    ap.add_argument("--maxh-air", type=float, default=0.010)
    ap.add_argument("--ratios-uniform", type=float, nargs='+',
                    default=[3.0, 5.0, 8.0])
    ap.add_argument("--ratios-skin", type=float, nargs='+',
                    default=[3.0, 5.0, 8.0])
    args = ap.parse_args()

    if args.freq is not None:
        freqs = [args.freq]
    else:
        freqs = [7000.0, 50000.0]

    for freq in freqs:
        run_convergence(
            freq=freq, material=args.material,
            r_wp=args.r_wp, h_wp=args.h_wp,
            order=args.order, maxh_air=args.maxh_air,
            skin_ratios_uniform=tuple(args.ratios_uniform),
            skin_ratios_skin=tuple(args.ratios_skin),
        )


if __name__ == "__main__":
    main()
