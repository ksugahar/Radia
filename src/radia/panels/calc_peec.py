"""calc_peec.py -- PEEC coil impedance extraction (Layer 4, no GUI).

Computes coil self-impedance Z(f) from parametric helix parameters
using the PEEC filament pipeline + PRIMA broadband sweep.

No coil mesh needed.  The coil is defined by:
  radius, pitch, turns, wire_w, wire_h, nwinc, seg_per_turn

For PEEC+BEM workpiece coupling: --vol supplies the workpiece mesh,
and a single Picard correction step applies the BEM-SIBC back-reaction.

Usage:
    python calc_peec.py --coil-radius 0.05 --coil-pitch 0.008 \\
        --coil-turns 5 --wire-w 0.006 --wire-h 0.006 \\
        --nwinc 3 --frequency 50000

Output: JSON to stdout (same schema as calc_inductance.py).
"""
from __future__ import annotations

import argparse
import json
import math
import os
import sys
import time

import numpy as np

# Ensure Radia src is on path
HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.abspath(os.path.join(HERE, ".."))
if SRC not in sys.path:
    sys.path.insert(0, SRC)


def build_peec_topology(args):
    """Build PEEC topology from helix parameters."""
    from radia.coil_from_cad import helix_path, build_peec_from_path

    n_seg = int(args.coil_turns * args.seg_per_turn)
    path = helix_path(
        radius=args.coil_radius,
        pitch=args.coil_pitch,
        n_turns=args.coil_turns,
        n_points=n_seg + 1,
    )

    # Per-segment uniform cross-section (no taper for panel default)
    widths = np.full(n_seg, args.wire_w, dtype=np.float64)
    heights = np.full(n_seg, args.wire_h, dtype=np.float64)

    topo = build_peec_from_path(
        path, widths, heights,
        sigma=args.coil_sigma,
        nwinc=args.nwinc,
        nhinc=args.nwinc,  # square subdivision
    )
    return topo, n_seg


def solve_prima(topo, args):
    """PRIMA broadband sweep."""
    from prima_hacapk import PRIMAHACApKModel

    t0 = time.perf_counter()
    model = PRIMAHACApKModel.from_topology(topo, q=args.prima_q, aca_eps=1e-8)
    t_build = time.perf_counter() - t0

    freq = float(args.frequency)
    Z = model.port_impedance(freq)

    # Also do a sweep around the target frequency (1 decade below to above)
    f_lo = max(1.0, freq / 10)
    f_hi = freq * 10
    freqs_sweep = np.logspace(np.log10(f_lo), np.log10(f_hi), 21)
    Z_sweep = model.impedance_sweep(freqs_sweep)

    return Z, t_build, model, freqs_sweep, Z_sweep


def solve_hacapk(topo, args):
    """Single-frequency HACApK saddle-point solve."""
    from radia.peec_topology import PEECCircuitSolver

    t0 = time.perf_counter()
    solver = PEECCircuitSolver(
        topo, use_hacapk=True,
        hacapk_params={"aca_eps": 1e-8},
        bicgstab_outer_tol=1e-8,
        outer_method="saddle",
    )
    Z = solver.compute_port_impedance(float(args.frequency))
    t_solve = time.perf_counter() - t0
    return Z, t_solve


def solve_dense(topo, args):
    """Dense C++ MNA solve (reference)."""
    from radia.peec_topology import PEECCircuitSolver

    t0 = time.perf_counter()
    solver = PEECCircuitSolver(topo)
    Z = solver.compute_port_impedance(float(args.frequency))
    t_solve = time.perf_counter() - t0
    return Z, t_solve


def main():
    p = argparse.ArgumentParser(description="PEEC coil impedance extraction")
    # Coil geometry
    p.add_argument("--coil-radius", type=float, required=True,
                   help="Helix radius [m]")
    p.add_argument("--coil-pitch", type=float, required=True,
                   help="Pitch per turn [m]")
    p.add_argument("--coil-turns", type=float, required=True,
                   help="Number of turns")
    p.add_argument("--wire-w", type=float, required=True,
                   help="Wire width [m]")
    p.add_argument("--wire-h", type=float, required=True,
                   help="Wire height [m]")
    p.add_argument("--nwinc", type=int, default=3,
                   help="Cross-section subdivision (skin/proximity)")
    p.add_argument("--seg-per-turn", type=int, default=20,
                   help="Filament segments per turn")
    # Physics
    p.add_argument("--frequency", type=float, required=True,
                   help="Target frequency [Hz]")
    p.add_argument("--current", type=float, default=1.0,
                   help="Coil current [A] (for reporting, does not affect Z)")
    p.add_argument("--coil-sigma", type=float, default=5.8e7,
                   help="Coil conductivity [S/m]")
    # Solver
    p.add_argument("--solver", choices=["prima", "hacapk", "dense"],
                   default="prima", help="Solver method")
    p.add_argument("--prima-q", type=int, default=30,
                   help="PRIMA Krylov order")
    # Workpiece (PEEC+BEM coupling, optional)
    p.add_argument("--vol", type=str, default="",
                   help=".vol file for workpiece (PEEC+BEM mode)")
    p.add_argument("--workpiece", type=str, default="",
                   help="Workpiece sideset name")
    p.add_argument("--impedance-model", choices=["sibc", "esim"],
                   default="sibc")
    p.add_argument("--sigma", type=float, default=2e6,
                   help="Workpiece sigma [S/m]")
    p.add_argument("--half-thickness", type=float, default=0.005)
    p.add_argument("--mu-r", type=float, default=100)
    p.add_argument("--use-local-curvature", action="store_true")

    args = p.parse_args()

    # Import Radia (triggers MKL DLL path)
    import radia  # noqa: F401

    # Build topology
    t0 = time.perf_counter()
    topo, n_seg = build_peec_topology(args)
    t_topo = time.perf_counter() - t0
    N = topo["n_loop"]

    # Solve
    freq = float(args.frequency)
    omega = 2 * math.pi * freq

    if args.solver == "prima":
        Z, t_solve, model, freqs_sweep, Z_sweep = solve_prima(topo, args)
    elif args.solver == "hacapk":
        Z, t_solve = solve_hacapk(topo, args)
        freqs_sweep, Z_sweep = None, None
    else:
        Z, t_solve = solve_dense(topo, args)
        freqs_sweep, Z_sweep = None, None

    L = Z.imag / omega if omega > 0 else 0.0
    R = Z.real

    result = {
        "status": "ok",
        "method": "PEEC",
        "solver": args.solver,
        "n_filaments": N,
        "n_segments": n_seg,
        "nwinc": args.nwinc,
        "frequency_hz": freq,
        "inductance_H": L,
        "inductance_nH": L * 1e9,
        "resistance_ohm": R,
        "resistance_mohm": R * 1e3,
        "impedance_real": Z.real,
        "impedance_imag": Z.imag,
        "abs_Z": abs(Z),
        "t_topology": t_topo,
        "t_solve": t_solve,
    }

    if args.solver == "prima" and freqs_sweep is not None:
        result["prima_q"] = args.prima_q
        result["sweep"] = {
            "freqs_hz": freqs_sweep.tolist(),
            "Z_real": [z.real for z in Z_sweep],
            "Z_imag": [z.imag for z in Z_sweep],
        }

    # TODO: PEEC+BEM workpiece coupling (Picard correction)
    if args.vol and args.workpiece:
        result["workpiece_coupling"] = "not_implemented_yet"

    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
