"""2D axisym reference vs 3D PEEC+BEM vs 3D FEM coilmesh.

Reference geometry (matches src/radia/panels/samples/ih_closed_torus.jou):
  Coil:      closed torus, R=30mm, a=3mm, Cu sigma=5.8e7
  Workpiece: cylinder R=25mm H=25mm, Cu sigma=5.8e7, mu_r=1
  Frequency: 7 kHz
  BC:        open (Kelvin for 2D axisym, Kelvin for 3D FEM)

Outputs:
  | Method              | L [nH]    | P_total [W]  |
  | 2D axisym full-eddy | reference | reference    |
  | 3D PEEC+BEM         | --        | test (5%)    |
  | 3D FEM coilmesh     | test (1%) | test         |

BEM only reports P (L is mesh-sensitive in 1-way forward — see memory
project peec_bem_production_decision_2026_04_18).  FEM reports both.
"""
from __future__ import annotations

import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
SRC_RADIA = os.path.join(ROOT, "src", "radia")
PANELS = os.path.join(SRC_RADIA, "panels")
BEM_REF = os.path.join(ROOT, "examples", "induction_heating",
                        "bem_reference")
for p in (SRC_RADIA, PANELS, BEM_REF, HERE):
    if p not in sys.path:
        sys.path.insert(0, p)

import radia  # DLL setup
from em_material import EMMaterial


def run_2d_axisym(freq, coil_R=0.030, coil_a=0.003,
                   wp_R=0.025, wp_H=0.025):
    from reference_2d_axisym import solve_2d_axisym
    mat = EMMaterial.from_name("copper")
    r = solve_2d_axisym(
        r_coil=coil_R, a_coil=coil_a, r_wp=wp_R, h_wp=wp_H,
        mat=mat, frequency=freq, I_total=1.0,
        kelvin=True, order=2, sigma_coil=5.8e7, mu_r_coil=1.0)
    return {"L_nH": r["L"] * 1e9, "P_W": r["P_total"],
            "H_t_Am": r["H_t_rms"], "ndof": r["ndof"], "ne": r["ne"]}


def run_peec_bem_1way(vol, step, freq, coil_a=0.003):
    from calc_peec_bem import solve_peec_bem_forward
    r = solve_peec_bem_forward(
        peec_step=step, peec_nwinc=3, peec_nhinc=3,
        frequency=freq, current=1.0, coil_sigma=5.8e7,
        vol=vol, wp_label="sibc",
        sigma=5.8e7, half_thickness=wp_half_thickness(0.025, 0.025),
        mu_r=1.0, impedance_model="sibc", h1_order=1)
    return {"L_nH": None,  # not evaluated for BEM
             "P_W": r["P_wp_W"],
             "H_t_Am": r["H_t_rms_Am"],
             "L_coil_peec_nH": r["L_coil_nH"],
             "bem_ndof": r["bem_ndof"]}


def wp_half_thickness(R, H):
    # SIBC reference length: smallest half-dimension, so the Dowell
    # formula reduces to the semi-infinite slab limit for the cylinder
    # wall (where the tangential H is largest).
    return min(R, H / 2)


def run_fem_coilmesh(vol, freq):
    from calc_fem_coilmesh import solve_fem_coilmesh
    r = solve_fem_coilmesh(
        vol=vol, frequency=freq, I_target=1.0,
        coil_sigma=5.8e7, coil_mu_r=1.0,
        wp_sigma=5.8e7, wp_mu_r=1.0,
        half_thickness=wp_half_thickness(0.025, 0.025),
        fes_order=1, solver="pardiso",
        sibc_bnd="sibc", source_bnd="source", sink_bnd="sink",
        coil_mat="coil")
    return {"L_nH": r["L_total_nH"], "P_W": r["P_total_W"],
             "H_t_Am": r["H_t_rms_wp_Am"],
             "L_vol_nH": r["L_vol_nH"],
             "L_skin_wp_nH": r["L_skin_wp_nH"],
             "P_coil_W": r["P_coil_W"], "P_surf_W": r["P_wp_W"],
             "ndof": r["ndof"], "ne": r["ne"]}


def main():
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--freq", type=float, default=7000.0)
    # The two 3D methods need DIFFERENT coil topologies:
    #   - FEM coilmesh imposes uniform J_phi over the coil volume.
    #     This is consistent only for a CLOSED axisymmetric torus.
    #   - BEM (scalar BIE SIBC) requires a single-valued magnetic
    #     scalar potential, which only exists when the coil has a GAP.
    # We therefore run each method on its own sample and use the 2D
    # axisym (which is geometrically closed) as the common reference.
    # The 3D BEM case has a ~5 deg gap (1.4% of circumference) — an
    # O(1.4%) asymmetry vs the 2D closed reference.
    ap.add_argument("--vol-fem", default=os.path.join(
        SRC_RADIA, "panels", "samples", "ih_fem_kelvin_skin_fine.vol"),
        help="FEM A-V .vol (gapped torus, source/sink, skin-resolved "
             "coil mesh h~0.16mm)")
    ap.add_argument("--vol-bem", default=os.path.join(
        SRC_RADIA, "panels", "samples", "ih_peec_bem_coarse.vol"),
        help="BEM .vol (gapped torus, wp as sibc sideset, coarse "
             "mesh for fast assembly)")
    ap.add_argument("--step-bem", default=os.path.join(
        SRC_RADIA, "panels", "samples",
        "ih_peec_bem_coarse_coil.step"),
        help="BEM STEP (gapped torus coil, matches vol-bem)")
    ap.add_argument("--skip-bem", action="store_true")
    ap.add_argument("--skip-fem", action="store_true")
    ap.add_argument("--skip-2d", action="store_true")
    args = ap.parse_args()

    print("=" * 72)
    print("  2D axisym vs 3D PEEC+BEM (P) vs 3D FEM coilmesh (L+P)")
    print(f"  Freq={args.freq:.0f} Hz, Cu wp (R=25mm H=25mm), "
          f"Cu coil (R=30mm a=3mm)")
    print(f"  FEM .vol : {os.path.relpath(args.vol_fem, ROOT)} (closed)")
    print(f"  BEM .vol : {os.path.relpath(args.vol_bem, ROOT)} (gap=5 deg)")
    print(f"  BEM .step: {os.path.relpath(args.step_bem, ROOT)}")
    print("=" * 72)

    rows = []
    if not args.skip_2d:
        print()
        print("[1/3] 2D axisym full-eddy (reference_2d_axisym.py)")
        r = run_2d_axisym(args.freq)
        rows.append(("2D axisym full-eddy", r))
        print(f"      L={r['L_nH']:.2f} nH  P={r['P_W']:.3e} W  "
              f"(ndof={r['ndof']}, ne={r['ne']})")

    if not args.skip_bem:
        print()
        print("[2/3] 3D PEEC+BEM 1-way  (gapped sample, calc_peec_bem.py)")
        r = run_peec_bem_1way(args.vol_bem, args.step_bem, args.freq)
        rows.append(("3D PEEC+BEM (P only)", r))
        print(f"      L={r['L_nH']}  P={r['P_W']:.3e} W  "
              f"(BEM ndof={r['bem_ndof']}, L_coil_peec="
              f"{r['L_coil_peec_nH']:.2f} nH for ref)")

    if not args.skip_fem:
        print()
        print("[3/3] 3D FEM A-V (gapped sample w/ source/sink, "
              "calc_fem_coilmesh.py)")
        r = run_fem_coilmesh(args.vol_fem, args.freq)
        rows.append(("3D FEM A-V (L+P)", r))
        print(f"      L={r['L_nH']:.2f} nH "
              f"(vol={r['L_vol_nH']:.2f} + skin_wp={r['L_skin_wp_nH']:+.2f})")
        print(f"      P={r['P_W']:.3e} W "
              f"(coil={r['P_coil_W']:.3e} + wp_surf={r['P_surf_W']:.3e})")
        print(f"      ndof={r['ndof']}, ne={r['ne']}")

    # Summary
    print()
    print("-" * 72)
    print(f"  {'Method':<24s} {'L [nH]':>10s} {'P [W]':>12s} "
           f"{'H_t [A/m]':>12s}")
    print("-" * 72)
    L_ref = None
    P_ref = None
    for name, r in rows:
        if "2D axisym" in name:
            L_ref = r["L_nH"]
            P_ref = r["P_W"]
        L_str = f"{r['L_nH']:.2f}" if r["L_nH"] is not None else "--"
        print(f"  {name:<24s} {L_str:>10s} {r['P_W']:>12.3e} "
              f"{r['H_t_Am']:>12.2f}")

    # Error table
    if L_ref is not None:
        print()
        print(f"  Error vs 2D ref ({L_ref:.2f} nH, {P_ref:.3e} W):")
        for name, r in rows:
            if "2D axisym" in name:
                continue
            dL = (f"{(r['L_nH'] - L_ref) / L_ref * 100:+.2f}%"
                   if r["L_nH"] is not None else "--")
            dP = (f"{(r['P_W'] - P_ref) / P_ref * 100:+.2f}%"
                   if r["P_W"] is not None else "--")
            print(f"    {name:<24s} dL={dL:>8s}  dP={dP:>8s}")


if __name__ == "__main__":
    main()
