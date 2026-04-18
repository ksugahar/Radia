"""compare_fem_vs_peec_skin.py

2-way IH inductance + heating comparison, same mesh and same SIBC workpiece:

  (1) A-V: full 3D FEM, coil volume meshed with coil_sigma
      - Coil skin/proximity eddies captured by the FEM itself
      - L_av = 2 * int nu |curl A|^2 / I^2
      - P_av = 0.5 * Re(Z_s) * H_t_rms^2 * A_wp
      - Heavy (ndof ~ 280k for this sample), but 3D reference truth

  (2) PEEC + back-reaction: PEEC circuit (L_peec) + FEM scattered A_r
      - PEEC filaments from STEP give Biot-Savart excitation A_s
      - FEM solves scattered A_r with SIBC Robin on workpiece surface
      - L_total = L_peec + Delta_L (line-integral of A_r along filaments)
      - P_wp    = 0.5 * Delta_R * |I|^2

Both paths use SIBC for the workpiece.  3D without a 2D-axisym shortcut
needs SIBC to avoid meshing the skin layer.
"""
import argparse
import math
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.abspath(os.path.join(HERE, '..', '..', '..', 'src'))
for p in (SRC, os.path.join(SRC, 'radia'), os.path.join(SRC, 'radia', 'panels'), HERE):
    if p not in sys.path:
        sys.path.insert(0, p)

SAMPLES = os.path.join(SRC, 'radia', 'panels', 'samples')

MU_0 = 4 * math.pi * 1e-7
NU_0 = 1.0 / MU_0


def fmt_pct(val, ref):
    if ref == 0 or ref is None:
        return "n/a"
    return f"{(val - ref) / abs(ref) * 100:+.2f}%"


def solve_peec_back_reaction(vol_file, step, frequency, mat,
                             half_thickness, nwinc, nhinc,
                             peec_sigma=5.8e7, I_total=1.0):
    """Thin wrapper around calc_fem_kelvin.solve_fem_biot_savart.

    Kept for the 3-way compare below and for backward API compat with any
    downstream callers.
    """
    from calc_fem_kelvin import solve_fem_biot_savart
    r = solve_fem_biot_savart(
        vol_file=vol_file, frequency=frequency, mat=mat,
        I_total=I_total, half_thickness=half_thickness,
        peec_step=step, peec_sigma=peec_sigma,
        peec_nwinc=nwinc, peec_nhinc=nhinc)
    return {
        'L_peec': r['L_peec'],
        'R_peec': r['R_peec'],
        'Delta_L': r['Delta_L'],
        'L_total': r['L'],
        'Delta_R': r['Delta_R'],
        'P_wp': r['P_total'],
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--vol", default=os.path.join(SAMPLES, "ih_fem_kelvin_skin.vol"))
    ap.add_argument("--step", default=os.path.join(SAMPLES, "ih_fem_kelvin_skin_coil.step"))
    ap.add_argument("--frequency", type=float, default=7000)
    ap.add_argument("--peec-nwinc", type=int, default=3)
    ap.add_argument("--peec-nhinc", type=int, default=3)
    ap.add_argument("--material", default="steel")
    ap.add_argument("--half-thickness", type=float, default=0.0125)
    args = ap.parse_args()

    from calc_fem_kelvin import solve_fem
    from calc_common import EMMaterial
    mat = EMMaterial.from_name(args.material)
    f = args.frequency

    runs = {}

    # --- Run 1: A-V (full FEM with coil_sigma + SIBC) ---
    print("=" * 60)
    print(f"Run 1: A-V (coil volume meshed + coil_sigma=5.8e7, SIBC wp)")
    print("=" * 60)
    t0 = time.perf_counter()
    runs["av"] = solve_fem(vol_file=args.vol, frequency=f, mat=mat,
                           impedance_model="sibc",
                           half_thickness=args.half_thickness,
                           coil_sigma=5.8e7, solver="pardiso")
    runs["t_av"] = time.perf_counter() - t0

    # --- Run 2: PEEC + back-reaction (scattered A_r via SIBC Robin) ---
    print("=" * 60)
    print(f"Run 2: PEEC nwinc={args.peec_nwinc} + back-reaction (scattered A_r)")
    print("=" * 60)
    t0 = time.perf_counter()
    runs["peec"] = solve_peec_back_reaction(
        vol_file=args.vol, step=args.step, frequency=f, mat=mat,
        half_thickness=args.half_thickness,
        nwinc=args.peec_nwinc, nhinc=args.peec_nhinc)
    runs["t_peec"] = time.perf_counter() - t0

    # --- Run 3: PEEC+BR with nwinc=1 centerline (no skin in PEEC) ---
    # Single-centerline filament: L_peec is the DC Neumann self-L
    # (~97 nH), FEM captures workpiece back-reaction only (no coil
    # eddy). Contrast with Run 2 (nwinc=3, PEEC captures skin) to
    # isolate the nwinc-driven skin correction. Both rows are T0-free.
    print("=" * 60)
    print(f"Run 3: PEEC nwinc=1 centerline + back-reaction (no skin in PEEC)")
    print("=" * 60)
    t0 = time.perf_counter()
    runs["avbs"] = solve_peec_back_reaction(
        vol_file=args.vol, step=args.step, frequency=f, mat=mat,
        half_thickness=args.half_thickness,
        nwinc=1, nhinc=1)
    runs["t_avbs"] = time.perf_counter() - t0

    # --- Results ---
    av = runs["av"]
    pe = runs["peec"]
    bs = runs["avbs"]

    L_av = av.get("L", 0) or 0
    L_pe = pe["L_total"]
    L_bs = bs["L_total"]
    P_av = av.get("P_total", 0) or 0
    P_pe = pe["P_wp"]
    P_bs = bs["P_wp"]

    print()
    print("=" * 60)
    print("COMPARISON")
    print("=" * 60)
    print(f"{'':18} {'A-V T0':>14} {'PEEC+BR nw3':>14} {'PEEC+BR nw1':>14}")
    print("-" * 74)
    print(f"{'L [nH]':18} {L_av*1e9:14.3f} {L_pe*1e9:14.3f} {L_bs*1e9:14.3f}")
    print(f"{'P [W]':18} {P_av:14.4e} {P_pe:14.4e} {P_bs:14.4e}")
    print()
    print(f"  Row 2 (nwinc=3, PEEC captures skin) decomposition:")
    print(f"    L_peec (circuit): {pe['L_peec']*1e9:.3f} nH  (HF skin-corrected)")
    print(f"    Delta_L (FEM BR): {pe['Delta_L']*1e9:+.3f} nH  "
          f"{'(Lenz, PASS)' if pe['Delta_L'] < 0 else '(positive - flux conc.?)'}")
    print(f"    L_total         : {pe['L_total']*1e9:.3f} nH")
    print()
    print(f"  Row 3 (nwinc=1 centerline) decomposition:")
    print(f"    L_peec (DC self): {bs['L_peec']*1e9:.3f} nH  (no skin correction)")
    print(f"    Delta_L (FEM BR): {bs['Delta_L']*1e9:+.3f} nH")
    print(f"    L_total         : {bs['L_total']*1e9:.3f} nH")
    print()
    print(f"  Times: A-V(T0) {runs['t_av']:.1f}s, PEEC+BR nw3 {runs['t_peec']:.1f}s, "
          f"PEEC+BR nw1 {runs['t_avbs']:.1f}s")

    if "error" in av:
        print(f"  A-V ERROR: {av['error']}")


if __name__ == "__main__":
    main()
