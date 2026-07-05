"""
PCB impedance extraction via PEEC (Partial Element Equivalent Circuit).

Called as subprocess from Cubit panel (PCBPEECDialog):
    python calc_pcb_peec.py --inp circuit.inp --output result.json

Pipeline:
  FastHenry .inp -> FastHenryParser -> PEECBuilder -> PEECCircuitSolver
  -> frequency sweep -> Z(f), L(f), R(f)

Supports:
  - conductor PEEC extraction; magnetic-material coupling is application-specific
  - Multi-filament (nwinc/nhinc) for skin/proximity effect
  - Multi-port Z-parameter extraction
  - SPICE netlist output (PRIMA Lanczos MOR)

IMPORTANT: NGSolve import not needed for basic PEEC (only for BEM shield).
"""

import argparse
import math
import os
import sys
import time

import numpy as np

_this_dir = os.path.dirname(os.path.abspath(__file__))
if _this_dir not in sys.path:
    sys.path.insert(0, _this_dir)

from calc_common import setup_paths, progress, calc_main


def _log(msg):
    progress("PEEC", msg)


def solve_peec(inp_file="", inp_text="",
               freq_min=1e3, freq_max=1e9, n_freq=50,
               solver_method=0, solver_prec=1e-4, solver_maxiter=1000,
               spice_output=""):
    """PEEC solver: FastHenry .inp -> impedance extraction.

    Args:
        inp_file: Path to FastHenry .inp file
        inp_text: Inline .inp text (alternative to inp_file)
        freq_min: Minimum frequency [Hz]
        freq_max: Maximum frequency [Hz]
        n_freq: Number of frequency points (log-spaced)
        solver_method: retained CLI slot; magnetic-material coupling is handled
            by HDiv-VIM / reduced-FEM workflows, not this PEEC panel.
        solver_prec: Solver precision
        solver_maxiter: Max solver iterations
        spice_output: Optional SPICE netlist output path

    Returns:
        dict with freqs, Z, L, R, topology info, etc.
    """
    setup_paths()
    from fasthenry_parser import FastHenryParser
    from radia.peec_topology import PEECCircuitSolver

    t_total_start = time.perf_counter()

    # ============================================================
    # Step 1: Parse input
    # ============================================================
    _log("PARSE:loading .inp")
    parser = FastHenryParser()

    if inp_file and os.path.exists(inp_file):
        parser.parse_file(inp_file)
    elif inp_text:
        parser.parse_string(inp_text)
    else:
        return {"error": "No .inp file or text provided"}

    summary = parser.get_summary()
    _log(f"PARSE:{summary.get('n_nodes', 0)} nodes, "
         f"{summary.get('n_segments', 0)} segments, "
         f"{summary.get('n_ports', 0)} ports")

    n_ports = summary.get('n_ports', 0)
    if n_ports == 0:
        return {"error": "No .external ports defined in .inp file"}

    # Use .freq from file if available
    freq_spec = parser.get_frequencies()
    if freq_spec is not None and len(freq_spec) > 0:
        freqs = freq_spec
        _log(f"FREQ:from .inp: {len(freqs)} points, "
             f"{freqs[0]:.1e}-{freqs[-1]:.1e} Hz")
    else:
        freqs = np.logspace(np.log10(freq_min), np.log10(freq_max), n_freq)
        _log(f"FREQ:generated: {n_freq} points, {freq_min:.1e}-{freq_max:.1e} Hz")

    # ============================================================
    # Step 2: Build topology
    # ============================================================
    _log("BUILD:topology")
    t0 = time.perf_counter()

    has_magnetic = bool(parser.magnetic_blocks)

    if has_magnetic:
        return {
            "error": (
                ".magnetic blocks used the retired PEEC magnetic-material "
                "coupling path. Use this panel for conductor/shield PEEC and "
                "HDiv-VIM / reduced FEM for magnetic cores."
            ),
            "has_magnetic": True,
        }

    else:
        # Standard PEEC (no magnetic coupling)
        builder = parser.to_peec_builder()
        topo = builder.build_topology()
        n_loops = topo['n_loop']
        _log(f"BUILD:{n_loops} loops, {topo['n_nodes']} nodes")

        solver = PEECCircuitSolver(topo)
        t_build = time.perf_counter() - t0
        _log(f"BUILD:done in {t_build:.2f}s")

        # Frequency sweep
        _log(f"SOLVE:sweep {len(freqs)} frequencies")
        t0 = time.perf_counter()

        if n_ports == 1:
            Z_port = solver.frequency_sweep(freqs)
            R_arr = np.real(Z_port)
            L_arr = np.imag(Z_port) / (2 * np.pi * np.maximum(freqs, 1e-30))
        else:
            # Multi-port
            sweep = solver.frequency_sweep_multiport(freqs)
            Z_port = sweep.get('Z_matrix', np.zeros((len(freqs), n_ports, n_ports),
                                                      dtype=complex))
            # Extract diagonal for summary
            R_arr = np.real(Z_port[:, 0, 0]) if Z_port.ndim == 3 else np.real(Z_port)
            L_arr = (np.imag(Z_port[:, 0, 0]) / (2 * np.pi * np.maximum(freqs, 1e-30))
                     if Z_port.ndim == 3 else
                     np.imag(Z_port) / (2 * np.pi * np.maximum(freqs, 1e-30)))

        t_solve = time.perf_counter() - t0

    _log(f"SOLVE:done in {t_solve:.2f}s")

    # DC values
    L_dc = float(L_arr[0]) if len(L_arr) > 0 else 0
    R_dc = float(R_arr[0]) if len(R_arr) > 0 else 0
    _log(f"RESULT:L_DC={L_dc*1e9:.4f} nH, R_DC={R_dc*1e3:.4f} mOhm")

    # ============================================================
    # Step 3: Optional SPICE extraction
    # ============================================================
    spice_file = ""
    if spice_output:
        try:
            from lanczos_reduction import SPICEExtractionConfig, PRIMASchurExtractor
            config = SPICEExtractionConfig(
                n_lanczos_loop=min(20, n_loops),
                port_indices=list(range(n_ports)),
                f_min=float(freqs[0]),
                f_max=float(freqs[-1]),
                n_freq=min(50, len(freqs)),
            )
            builder = parser.to_peec_builder()
            topo = builder.build_topology(include_star=True)

            extractor = PRIMASchurExtractor(config)
            R_mat = np.diag(topo['R']) if topo['R'].ndim == 1 else topo['R']
            P_mat = topo.get('P', np.zeros((1, 1)))
            K_LS = topo.get('M_LS', np.zeros((len(topo['R']), max(1, P_mat.shape[0]))))
            n_L = len(topo['R'])

            mor_result = extractor.extract(
                L=topo['L'], R=R_mat, P=P_mat,
                Z_MM=np.zeros((1, 1)), Z_DD=np.zeros((1, 1)),
                Z_CC=np.zeros((1, 1)),
                K_LM=np.zeros((n_L, 1)), K_LD=np.zeros((n_L, 1)),
                K_LC=np.zeros((n_L, 1)), K_LS=K_LS,
                frequencies=freqs,
            )

            netlist = mor_result.get('netlist', '')
            if netlist:
                with open(spice_output, 'w') as f:
                    f.write(netlist)
                spice_file = spice_output
                _log(f"SPICE:{spice_output}")
        except Exception as e:
            _log(f"SPICE:extraction failed: {e}")

    t_total = time.perf_counter() - t_total_start

    # Build result (JSON-serializable)
    result_dict = {
        "freqs": freqs.tolist(),
        "Z_real": np.real(Z_port).tolist() if np.ndim(Z_port) <= 1 else None,
        "Z_imag": np.imag(Z_port).tolist() if np.ndim(Z_port) <= 1 else None,
        "R": R_arr.tolist(),
        "L": L_arr.tolist(),
        "L_dc_nH": L_dc * 1e9,
        "R_dc_mOhm": R_dc * 1e3,
        "n_nodes": summary.get('n_nodes', 0),
        "n_segments": summary.get('n_segments', 0),
        "n_loops": n_loops,
        "n_ports": n_ports,
        "has_magnetic": has_magnetic,
        "t_solve": t_solve,
        "t_total": t_total,
        "spice_file": spice_file,
    }

    return result_dict


def build_argparser():
    """argparse factory shared by main() and the panel generator.

    The panel (PCBPanel.__init__) calls this to auto-generate widgets
    via ModePanel.bind_argparser, so adding/removing/renaming an option
    here propagates to the GUI without any panel-side edit.
    """
    parser = argparse.ArgumentParser(
        description="PCB impedance extraction via PEEC")
    parser.add_argument("--inp", default="",
                        help="FastHenry .inp file")
    parser.add_argument("--freq-min", type=float, default=1e3,
                        help="Min frequency [Hz]")
    parser.add_argument("--freq-max", type=float, default=1e9,
                        help="Max frequency [Hz]")
    parser.add_argument("--n-freq", type=int, default=50,
                        help="Number of frequency points")
    parser.add_argument("--solver-method", type=int, default=0,
                        choices=[0, 1, 2],
                        help="Retained option; magnetic-material PEEC coupling is retired")
    parser.add_argument("--spice-output", default="",
                        help="SPICE netlist output path")
    parser.add_argument("--output", default="",
                        help="JSON output file")
    return parser


def main():
    parser = build_argparser()

    def run(args):
        return solve_peec(
            inp_file=args.inp,
            freq_min=args.freq_min,
            freq_max=args.freq_max,
            n_freq=args.n_freq,
            solver_method=args.solver_method,
            spice_output=args.spice_output,
        )

    calc_main(run, parser)


if __name__ == "__main__":
    main()
