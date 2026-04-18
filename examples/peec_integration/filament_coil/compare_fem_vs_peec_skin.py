"""compare_fem_vs_peec_skin.py

2-way IH inductance + heating comparison, same mesh and same SIBC workpiece.

Both runs use the PEEC filament Biot-Savart excitation + FEM scattered
A_r with SIBC Robin on the workpiece surface (solve_fem_biot_savart in
calc_fem_kelvin.py). They differ only in the PEEC filament-bundle
discretization:

  (1) PEEC nwinc=3 : 3x3 cross-section subdivision, captures coil skin
      effect through the PEEC circuit (L_peec from Z_port).
  (2) PEEC nwinc=1 : single centerline filament, no skin in PEEC;
      L_peec is the DC Neumann self-inductance.

Both paths return L_total = L_peec + Delta_L, where Delta_L is the line
integral of A_r along the filaments (workpiece back-reaction). Comparing
the two isolates the skin-effect contribution in L_peec from the
workpiece back-reaction in Delta_L.

The legacy A-V T0 row was removed 2026-04-18 when T0 / source-sink
Laplace was retired (unreliable on gapped geometries: 1/r cusps at the
gap corners inflated the FEM self-energy by +30 nH at 7 kHz on this
sample). The PEEC+BR path is now the trustworthy reference.
"""
import argparse
import math
import os
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.abspath(os.path.join(HERE, '..', '..', '..', 'src'))
for p in (SRC, os.path.join(SRC, 'radia'), os.path.join(SRC, 'radia', 'panels'), HERE):
    if p not in sys.path:
        sys.path.insert(0, p)

SAMPLES = os.path.join(SRC, 'radia', 'panels', 'samples')


def solve_peec_back_reaction(vol_file, step, frequency, mat,
                             half_thickness, nwinc, nhinc,
                             peec_sigma=5.8e7, I_total=1.0):
    """Thin wrapper around calc_fem_kelvin.solve_fem_biot_savart."""
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

    from calc_common import EMMaterial
    mat = EMMaterial.from_name(args.material)
    f = args.frequency

    runs = {}

    print("=" * 60)
    print(f"Run A: PEEC nwinc={args.peec_nwinc} + back-reaction (scattered A_r)")
    print("=" * 60)
    t0 = time.perf_counter()
    runs["peec"] = solve_peec_back_reaction(
        vol_file=args.vol, step=args.step, frequency=f, mat=mat,
        half_thickness=args.half_thickness,
        nwinc=args.peec_nwinc, nhinc=args.peec_nhinc)
    runs["t_peec"] = time.perf_counter() - t0

    print("=" * 60)
    print(f"Run B: PEEC nwinc=1 centerline + back-reaction (no skin in PEEC)")
    print("=" * 60)
    t0 = time.perf_counter()
    runs["center"] = solve_peec_back_reaction(
        vol_file=args.vol, step=args.step, frequency=f, mat=mat,
        half_thickness=args.half_thickness,
        nwinc=1, nhinc=1)
    runs["t_center"] = time.perf_counter() - t0

    pe = runs["peec"]
    ce = runs["center"]

    print()
    print("=" * 60)
    print("COMPARISON")
    print("=" * 60)
    print(f"{'':18} {'PEEC+BR nw3':>14} {'PEEC+BR nw1':>14}")
    print("-" * 60)
    print(f"{'L [nH]':18} {pe['L_total']*1e9:14.3f} {ce['L_total']*1e9:14.3f}")
    print(f"{'P [W]':18} {pe['P_wp']:14.4e} {ce['P_wp']:14.4e}")
    print()
    print(f"  Run A (nwinc={args.peec_nwinc}, PEEC captures skin):")
    print(f"    L_peec (circuit): {pe['L_peec']*1e9:.3f} nH  (HF skin-corrected)")
    print(f"    Delta_L (FEM BR): {pe['Delta_L']*1e9:+.3f} nH  "
          f"{'(Lenz, PASS)' if pe['Delta_L'] < 0 else '(positive - flux conc.?)'}")
    print(f"    L_total         : {pe['L_total']*1e9:.3f} nH")
    print()
    print(f"  Run B (nwinc=1 centerline):")
    print(f"    L_peec (DC self): {ce['L_peec']*1e9:.3f} nH  (no skin correction)")
    print(f"    Delta_L (FEM BR): {ce['Delta_L']*1e9:+.3f} nH")
    print(f"    L_total         : {ce['L_total']*1e9:.3f} nH")
    print()
    print(f"  Times: nw{args.peec_nwinc} {runs['t_peec']:.1f}s, "
          f"nw1 {runs['t_center']:.1f}s")


if __name__ == "__main__":
    main()
