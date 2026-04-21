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


def solve_peec_inductance(peec_input, n_peri,
                          frequency, current, coil_sigma):
    """Build PEEC topology from STEP or JOU input, then solve for L, R.

    Routing:
      - `.jou`            -> `coil_from_jou.filaments_from_jou` (explicit
                              centerline from `move Surface` commands).
                              Use for multi-turn / 3turnCoil-class coils
                              where STEP walker heuristics fail.
      - `.step` / `.stp`  -> `coil_from_cad.filaments_from_step` (walker
                              centerline extraction).  Use for clean
                              single-loop torus / helix coils.
      - other             -> error.
    """
    # Trigger Radia's MKL DLL path setup before peec_matrices loads.
    import radia  # noqa: F401

    from peec_bundle import (build_loop_bundle_impedance,
                              solve_loop_bundle)

    omega = 2 * math.pi * frequency

    ext = os.path.splitext(peec_input)[1].lower()

    # Auto-prefer sibling .jou (exact stem match, case-insensitive).
    # When user happens to have both foo.step and foo.jou in the same
    # directory — the common Cubit export pattern, since the panel's
    # ensure_jou_path() saves the .jou before every STEP export — we
    # use the .jou for correct L on multi-turn lofts.  The STEP
    # longest-edge path currently mis-estimates cross-section area on
    # tight-pancake multi-turn lofts (Kubota's 3turncoil.stp: STEP
    # path L=4.8 nH WRONG vs sibling 3turnCoil.jou L=426 nH correct).
    # Fixing that is a cross-section-geometry bug (see TODO below);
    # meanwhile sibling-match gives the user the right answer when
    # they have the right file.
    if ext in (".step", ".stp"):
        peec_dir = os.path.dirname(peec_input) or "."
        base = os.path.splitext(os.path.basename(peec_input))[0]
        try:
            entries = os.listdir(peec_dir)
        except OSError:
            entries = []
        jou_sibling = None
        for name in entries:
            stem, e = os.path.splitext(name)
            if e.lower() == ".jou" and stem.lower() == base.lower():
                jou_sibling = os.path.join(peec_dir, name)
                break
        if jou_sibling is not None:
            progress("PEEC", f"found sibling .jou, preferring it: "
                              f"{os.path.basename(jou_sibling)}")
            peec_input = jou_sibling
            ext = ".jou"

    if ext == ".jou":
        from coil_from_jou import filaments_from_jou
        source_kind = "JOU"
        progress("PEEC", f"JOU -> explicit centerline (n_peri={n_peri})")
        t0 = time.perf_counter()
        topo = filaments_from_jou(peec_input, sigma=coil_sigma,
                                    n_peri=n_peri)
    elif ext in (".step", ".stp"):
        from coil_from_cad import filaments_from_step
        source_kind = "STEP"
        progress("PEEC", f"STEP -> perimeter filaments (n_peri={n_peri})")
        t0 = time.perf_counter()
        topo = filaments_from_step(peec_input, sigma=coil_sigma,
                                    n_peri=n_peri,
                                    use_coil_builder=True)
    else:
        raise ValueError(
            f"Unsupported input extension {ext!r}. "
            f"Expected .jou (Cubit journal with explicit centerline) or "
            f".step / .stp (geometry file for walker extraction).")

    paths = topo["filament_paths"]
    seg_of_fil = topo["seg_of_filament"]
    solver = topo["solver"]
    t_topo = time.perf_counter() - t0
    progress("PEEC", f"{len(paths)} filaments via {source_kind}, {t_topo:.1f}s")

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
        "input_kind": source_kind,
        "input_path": peec_input,
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
                        help="Input file for PEEC coil: .step / .stp "
                             "(walker extraction) or .jou (explicit "
                             "centerline from Cubit journal). The arg "
                             "name stays --peec-step for backward "
                             "compatibility; extension chooses the path.")
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
            peec_input=args.peec_step,
            n_peri=args.peec_n_peri,
            frequency=args.frequency,
            current=args.current,
            coil_sigma=args.coil_sigma,
        )

    calc_main(run, parser)


if __name__ == "__main__":
    main()
