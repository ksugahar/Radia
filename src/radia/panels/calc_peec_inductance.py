"""calc_peec_inductance.py -- PEEC coil inductance from STEP (Layer 4, no GUI).

Layer 4 subprocess calc for the IH panel "PEEC inductance (coil only, STEP)"
method.  No workpiece, no BEM, no FEM mesh.  Just STEP -> filaments -> L, R.

Pipeline:
  1. STEP -> filaments (CoilBuilder, perimeter-only placement)
  2. PEEC Loop-bundle solve at given frequency
  3. Report L_coil, R_coil (vacuum)

Filaments are placed on the cross-section PERIMETER only (thin-skin limit,
d / delta >= 3 — the typical IH operating regime).  Use n_peri filaments
spread around the arc-length perimeter of each cross-section.

Usage (from the IH panel):
    python calc_peec_inductance.py --peec-step coil.step \\
        --peec-n-peri 16 \\
        --frequency 7000 --current 1.0 --coil-sigma 5.8e7

Output: JSON to stdout (calc_main contract).
"""
from __future__ import annotations

import argparse
import math
import os
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
SRC_RADIA = os.path.abspath(os.path.join(HERE, ".."))
for p in (SRC_RADIA, HERE):
    if p not in sys.path:
        sys.path.insert(0, p)

from calc_common import calc_main, progress


def solve_peec_inductance(peec_step, n_peri,
                          frequency, current, coil_sigma):
    """STEP -> perimeter filaments -> PEEC Loop-bundle -> L_coil, R_coil (vacuum)."""
    # Trigger Radia's MKL DLL path setup before peec_matrices loads.
    import radia  # noqa: F401

    from coil_from_cad import filaments_from_step
    from peec_bundle import (build_loop_bundle_impedance,
                              solve_loop_bundle)

    omega = 2 * math.pi * frequency

    progress("PEEC", f"STEP -> perimeter filaments (n_peri={n_peri})")
    t0 = time.perf_counter()
    topo = filaments_from_step(peec_step, sigma=coil_sigma,
                                n_peri=n_peri,
                                use_coil_builder=True)
    paths = topo["filament_paths"]
    seg_of_fil = topo["seg_of_filament"]
    solver = topo["solver"]
    t_topo = time.perf_counter() - t0
    progress("PEEC", f"{len(paths)} filaments, {t_topo:.1f}s")

    progress("PEEC", f"Loop-bundle solve @ {frequency:.0f} Hz")
    t0 = time.perf_counter()
    R_f, L_f = build_loop_bundle_impedance(solver, seg_of_fil)
    I_fil, V_port = solve_loop_bundle(R_f, L_f, frequency, I_port=current)
    t_peec = time.perf_counter() - t0
    Z_coil = V_port / current
    L_coil = Z_coil.imag / omega if omega > 0 else 0.0
    R_coil = Z_coil.real
    progress("PEEC", f"L_coil={L_coil*1e9:.3f} nH, "
                      f"R_coil={R_coil*1e3:.4f} mOhm ({t_peec:.1f}s)")

    return {
        "status": "ok",
        "method": "PEEC inductance (coil only, STEP)",
        "placement": "perimeter",
        "n_peri": int(n_peri),
        "frequency_hz": float(frequency),
        "current_A": float(current),
        "n_filaments": int(len(paths)),
        "coil_sigma_Sm": float(coil_sigma),
        "L_coil_H": float(L_coil),
        "L_coil_nH": float(L_coil * 1e9),
        "R_coil_Ohm": float(R_coil),
        "R_coil_mOhm": float(R_coil * 1e3),
        "impedance_real": float(Z_coil.real),
        "impedance_imag": float(Z_coil.imag),
        "abs_Z": float(abs(Z_coil)),
        "inductance_H": float(L_coil),
        "t_topology_s": float(t_topo),
        "t_peec_solve_s": float(t_peec),
        "t_solve_s": float(t_topo + t_peec),
    }


def main():
    parser = argparse.ArgumentParser(
        description="PEEC coil inductance from STEP (vacuum, no workpiece)")
    parser.add_argument("--peec-step", required=True,
                        help="STEP file for PEEC coil")
    parser.add_argument("--peec-n-peri", type=int, default=16,
                        help="Number of filaments on the cross-section "
                             "perimeter (thin-skin regime; d/delta>=3)")
    parser.add_argument("--frequency", type=float, required=True,
                        help="Frequency [Hz]")
    parser.add_argument("--current", type=float, default=1.0,
                        help="Port current [A]")
    parser.add_argument("--coil-sigma", type=float, default=5.8e7,
                        help="Coil conductivity [S/m]")

    def run(args):
        return solve_peec_inductance(
            peec_step=args.peec_step,
            n_peri=args.peec_n_peri,
            frequency=args.frequency,
            current=args.current,
            coil_sigma=args.coil_sigma,
        )

    calc_main(run, parser)


if __name__ == "__main__":
    main()
