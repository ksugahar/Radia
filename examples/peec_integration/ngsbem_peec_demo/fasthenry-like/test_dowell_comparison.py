"""
test_dowell_comparison.py

Compare ngbem PEEC vs FastHenry PEEC with Dowell skin effect
for a 10x10mm copper plate (35um thick, sigma=5.8e7 S/m).

Part of Radia project
"""

import sys
import os
import time
import numpy as np

_script_dir = os.path.dirname(os.path.abspath(__file__))
_src_radia = os.path.normpath(os.path.join(_script_dir, "..", "..", "..", "src", "radia"))
sys.path.insert(0, _src_radia)
sys.path.insert(0, _script_dir)

# Add MKL DLL directory for peec_matrices.pyd
import sys as _sys
_mkl_bin = os.path.join(_sys.prefix, "Library", "bin")
if os.path.isdir(_mkl_bin):
    os.add_dll_directory(_mkl_bin)

MU_0 = 4.0 * np.pi * 1e-7

from ngbem_peec import (
    NGBEMPEECSolver, create_plate_mesh,
    dowell_F_R, make_dowell_Zs_func,
)
from fasthenry_parser import FastHenryParser
from peec_topology import PEECCircuitSolver


def extract_R_L(Z, freqs):
    R = np.real(Z)
    omega = 2.0 * np.pi * freqs
    L = np.zeros_like(freqs)
    mask = omega > 0
    L[mask] = np.imag(Z[mask]) / omega[mask]
    return R, L


def fmt_z(z):
    return "%+.4e%+.4ej" % (z.real, z.imag)


def main():
    print("=" * 70)
    print("ngbem PEEC vs FastHenry PEEC: Dowell Skin Effect Comparison")
    print("=" * 70)
    print()

    width = 0.010
    height = 0.010
    thickness = 35e-6
    sigma = 5.8e7
    maxh = 0.003

    print("Geometry:")
    print("  Plate:     %.1f x %.1f mm" % (width*1e3, height*1e3))
    print("  Thickness: %.0f um" % (thickness*1e6))
    print("  Sigma:     %.1e S/m" % sigma)
    print("  ngbem maxh: %.1f mm" % (maxh*1e3))
    print()

    key_freqs = np.array([1e3, 1e5, 1e6])
    key_labels = ["1 kHz", "100 kHz", "1 MHz"]

    print("-" * 70)
    print("Section 1: ngbem PEEC (BEM Galerkin)")
    print("-" * 70)

    t0_ngbem = time.perf_counter()
    mesh = create_plate_mesh(width, height, maxh, label="conductor")

    ngbem_solver = NGBEMPEECSolver(
        mesh, conductor_label="conductor",
        sigma=sigma, thickness=thickness,
        order=0, intorder=5,
    )

    ngbem_solver.assemble()
    t_assemble_ngbem = time.perf_counter() - t0_ngbem

    print("  n_loop (edges):  %d" % ngbem_solver.n_loop)
    print("  n_star (cells):  %d" % ngbem_solver.n_star)
    print("  Assembly time:   %.3f s" % t_assemble_ngbem)
    print()

    print("  1a) Frequency sweep: DC resistance only (Zs_func=None)")
    t0 = time.perf_counter()
    Z_ngbem_dc = ngbem_solver.solve_frequency(key_freqs, mode="mqs", Zs_func=None)
    t_dc_ngbem = time.perf_counter() - t0
    print("      Sweep time: %.4f s" % t_dc_ngbem)

    print("  1b) Frequency sweep: Dowell Zs_func")
    Zs_func_ngbem = make_dowell_Zs_func(
        ngbem_solver.n_loop, ngbem_solver.R_loop,
        thickness, sigma, mu_r=1.0,
    )
    t0 = time.perf_counter()
    Z_ngbem_dowell = ngbem_solver.solve_frequency(key_freqs, mode="mqs", Zs_func=Zs_func_ngbem)
    t_dowell_ngbem = time.perf_counter() - t0
    print("      Sweep time: %.4f s" % t_dowell_ngbem)

    t_total_ngbem = time.perf_counter() - t0_ngbem
    print("  Total ngbem time: %.3f s" % t_total_ngbem)
    print()

    print("-" * 70)
    print("Section 2: FastHenry PEEC (filament model)")
    print("-" * 70)

    inp_lines = [
        "* 10x10mm copper plate as 4-segment rectangular loop",
        ".Units m",
        ".default sigma=5.8e7 nwinc=1 nhinc=1",
        "",
        "N1 x=0     y=0     z=0",
        "N2 x=0.01  y=0     z=0",
        "N3 x=0.01  y=0.01  z=0",
        "N4 x=0     y=0.01  z=0",
        "",
        "E1 N1 N2 w=3.5e-5 h=0.01",
        "E2 N2 N3 w=3.5e-5 h=0.01",
        "E3 N3 N4 w=3.5e-5 h=0.01",
        "E4 N4 N1 w=3.5e-5 h=0.01",
        "",
        ".external N1 N2",
        ".freq fmin=1e3 fmax=1e6 ndec=3",
        ".end",
    ]
    inp_text = chr(10).join(inp_lines)

    t0_fh = time.perf_counter()

    parser = FastHenryParser()
    parser.parse_string(inp_text)

    summary = parser.get_summary()
    print("  Model: %s" % str(summary))

    builder = parser.to_peec_builder()
    topo = builder.build_topology()

    n_loop_fh = topo["n_loop"]
    R_dc_fh = np.array(topo["R"])
    L_fh = np.array(topo["L"])

    t_build_fh = time.perf_counter() - t0_fh

    print("  n_loop (filaments): %d" % n_loop_fh)
    print("  R_dc per segment: %s" % str(R_dc_fh))
    print("  L matrix diagonal [nH]: %s" % str(np.diag(L_fh)*1e9))
    print("  Build time: %.4f s" % t_build_fh)
    print()

    fh_solver = PEECCircuitSolver(topo)

    print("  2a) Frequency sweep: DC resistance only (Zs_func=None)")
    t0 = time.perf_counter()
    Z_fh_dc = fh_solver.frequency_sweep(key_freqs, Zs_func=None)
    t_dc_fh = time.perf_counter() - t0
    print("      Sweep time: %.6f s" % t_dc_fh)

    print("  2b) Frequency sweep: Dowell Zs_func")

    def fh_dowell_Zs_func(freq):
        if freq <= 0:
            return np.zeros(n_loop_fh, dtype=complex)
        omega = 2.0 * np.pi * freq
        delta = np.sqrt(2.0 / (omega * MU_0 * sigma))
        xi = thickness / (2.0 * delta)
        F_R_val = dowell_F_R(xi)
        return R_dc_fh * (F_R_val - 1.0) + 0j

    t0 = time.perf_counter()
    Z_fh_dowell = fh_solver.frequency_sweep(key_freqs, Zs_func=fh_dowell_Zs_func)
    t_dowell_fh = time.perf_counter() - t0
    print("      Sweep time: %.6f s" % t_dowell_fh)

    t_total_fh = time.perf_counter() - t0_fh
    print("  Total FastHenry time: %.4f s" % t_total_fh)
    print()

    print("-" * 70)
    print("Section 3: Dowell F_R values at key frequencies")
    print("-" * 70)
    print()
    print("  Conductor thickness d = %.0f um" % (thickness*1e6))
    print()
    print("  %12s  %12s  %14s  %10s" % ("Frequency", "delta [um]", "xi=d/(2*delta)", "F_R"))
    print("  %s  %s  %s  %s" % ("-"*12, "-"*12, "-"*14, "-"*10))

    for i, (f_val, label) in enumerate(zip(key_freqs, key_labels)):
        omega = 2.0 * np.pi * f_val
        delta = np.sqrt(2.0 / (omega * MU_0 * sigma))
        xi = thickness / (2.0 * delta)
        F_R_val = dowell_F_R(xi)
        print("  %12s  %12.2f  %14.6f  %10.6f" % (label, delta*1e6, xi, F_R_val))
    print()

    print("-" * 70)
    print("Section 4: Comparison Table")
    print("-" * 70)
    print()

    R_ngbem_dc, L_ngbem_dc = extract_R_L(Z_ngbem_dc, key_freqs)
    R_ngbem_dw, L_ngbem_dw = extract_R_L(Z_ngbem_dowell, key_freqs)
    R_fh_dc_arr, L_fh_dc_arr = extract_R_L(Z_fh_dc, key_freqs)
    R_fh_dw_arr, L_fh_dw_arr = extract_R_L(Z_fh_dowell, key_freqs)

    print("  4a) DC resistance only (no skin effect)")
    print()
    print("  %10s  %16s  %14s  %14s  %12s" % ("Freq", "ngbem R [mOhm]", "ngbem L [nH]", "FH R [mOhm]", "FH L [nH]"))
    print("  %s  %s  %s  %s  %s" % ("-"*10, "-"*16, "-"*14, "-"*14, "-"*12))
    for i, label in enumerate(key_labels):
        print("  %10s  %16.4f  %14.4f  %14.4f  %12.4f" % (
            label, R_ngbem_dc[i]*1e3, L_ngbem_dc[i]*1e9,
            R_fh_dc_arr[i]*1e3, L_fh_dc_arr[i]*1e9))
    print()

    print("  4b) With Dowell skin effect")
    print()
    print("  %10s  %16s  %14s  %14s  %12s" % ("Freq", "ngbem R [mOhm]", "ngbem L [nH]", "FH R [mOhm]", "FH L [nH]"))
    print("  %s  %s  %s  %s  %s" % ("-"*10, "-"*16, "-"*14, "-"*14, "-"*12))
    for i, label in enumerate(key_labels):
        print("  %10s  %16.4f  %14.4f  %14.4f  %12.4f" % (
            label, R_ngbem_dw[i]*1e3, L_ngbem_dw[i]*1e9,
            R_fh_dw_arr[i]*1e3, L_fh_dw_arr[i]*1e9))
    print()

    print("  4c) Full impedance Z = R + jwL")
    print()
    print("  %10s  %-30s  %-30s" % ("Freq", "ngbem Z (DC)", "FH Z (DC)"))
    print("  %s  %s  %s" % ("-"*10, "-"*30, "-"*30))
    for i, label in enumerate(key_labels):
        print("  %10s  %s  %s" % (label, fmt_z(Z_ngbem_dc[i]), fmt_z(Z_fh_dc[i])))
    print()
    print("  %10s  %-30s  %-30s" % ("Freq", "ngbem Z (Dowell)", "FH Z (Dowell)"))
    print("  %s  %s  %s" % ("-"*10, "-"*30, "-"*30))
    for i, label in enumerate(key_labels):
        print("  %10s  %s  %s" % (label, fmt_z(Z_ngbem_dowell[i]), fmt_z(Z_fh_dowell[i])))
    print()

    print("-" * 70)
    print("Section 5: Timing Summary")
    print("-" * 70)
    print()
    print("  ngbem PEEC:")
    print("    Assembly:       %.3f s" % t_assemble_ngbem)
    print("    DC sweep:       %.4f s" % t_dc_ngbem)
    print("    Dowell sweep:   %.4f s" % t_dowell_ngbem)
    print("    Total:          %.3f s" % t_total_ngbem)
    print()
    print("  FastHenry PEEC:")
    print("    Build:          %.4f s" % t_build_fh)
    print("    DC sweep:       %.6f s" % t_dc_fh)
    print("    Dowell sweep:   %.6f s" % t_dowell_fh)
    print("    Total:          %.4f s" % t_total_fh)
    print()
    if t_total_fh > 0:
        print("  Speedup (ngbem / FastHenry): %.1fx" % (t_total_ngbem / t_total_fh))
    print()

    print("-" * 70)
    print("Section 6: Skin Effect R_ac/R_dc Ratio")
    print("-" * 70)
    print()
    print("  %10s  %16s  %14s  %12s" % ("Freq", "ngbem R_ac/R_dc", "FH R_ac/R_dc", "Dowell F_R"))
    print("  %s  %s  %s  %s" % ("-"*10, "-"*16, "-"*14, "-"*12))
    for i, label in enumerate(key_labels):
        r_ngbem_ratio = R_ngbem_dw[i] / R_ngbem_dc[i] if R_ngbem_dc[i] > 0 else 0
        r_fh_ratio = R_fh_dw_arr[i] / R_fh_dc_arr[i] if R_fh_dc_arr[i] > 0 else 0
        omega = 2.0 * np.pi * key_freqs[i]
        delta = np.sqrt(2.0 / (omega * MU_0 * sigma))
        xi = thickness / (2.0 * delta)
        F_R_val = dowell_F_R(xi)
        print("  %10s  %16.6f  %14.6f  %12.6f" % (label, r_ngbem_ratio, r_fh_ratio, F_R_val))
    print()

    print("=" * 70)
    print("Done.")
    print("=" * 70)


if __name__ == "__main__":
    main()
