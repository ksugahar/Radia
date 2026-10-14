"""
FastHenry-like PEEC: Material coupling (ferrite + shield) on a flat plate.

Models a 10x10mm copper plate (35 um thick) as a 4-segment rectangular
loop using Radia's FastHenry parser. Demonstrates:
  1. Air-only baseline (MQS loop impedance)
  2. Ferrite core coupling (analytical image method: Delta_L = L / (mu_r + 1))
  3. Al shield coupling (BEM+SIBC eddy current solve via ShieldBEMSIBC)
  4. Combined ferrite + shield

Shield geometry: 12x12 mm Al plate, 5 mm above conductor, 0.5 mm thick.

Requires: radia (peec_matrices.pyd), fasthenry_parser, peec_topology,
          ngsolve, ngsolve-ngsbem
"""
import sys
import os
import time
import numpy as np

# Path setup
_script_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(_script_dir, '..'))
sys.path.insert(0, os.path.join(_script_dir, '../../../../src/radia'))

# MKL DLL directory for peec_matrices.pyd
_mkl_bin = os.path.join(sys.prefix, "Library", "bin")
if os.path.isdir(_mkl_bin):
    os.add_dll_directory(_mkl_bin)

MU_0 = 4.0 * np.pi * 1e-7

from fasthenry_parser import FastHenryParser
from peec_topology import PEECCircuitSolver


def compute_delta_L(L_air, mu_r):
    """Image method coupling: Delta_L = L_air / (mu_r + 1)."""
    return L_air / (mu_r + 1.0)


def uniform_port_Z(Z_branch, n):
    """Port impedance via uniform excitation: Z = 1 / (e^T Y e)."""
    Y = np.linalg.inv(Z_branch)
    e = np.ones(n) / n
    return 1.0 / (e @ Y @ e)


def main():
    # --- Geometry (same as ngbem/demo_material_coupling.py) ---
    width = 0.01        # 10 mm
    height = 0.01       # 10 mm
    thickness = 35e-6   # 35 um copper
    sigma = 5.8e7       # Cu conductivity [S/m]

    w_mm = width * 1e3
    h_mm = height * 1e3
    t_mm = thickness * 1e3

    print("=== FastHenry PEEC: Material coupling demo ===")
    print(f"Conductor: {w_mm:.0f} x {h_mm:.0f} mm, t = {thickness*1e6:.0f} um Cu")

    # --- FastHenry model: 4-segment rectangular loop ---
    inp_text = f""".Units mm
.default sigma={sigma:.2e}

N1 x=0       y=0       z=0
N2 x={w_mm}  y=0       z=0
N3 x={w_mm}  y={h_mm}  z=0
N4 x=0       y={h_mm}  z=0

E1 N1 N2 w={w_mm} h={t_mm}
E2 N2 N3 w={h_mm} h={t_mm}
E3 N3 N4 w={w_mm} h={t_mm}
E4 N4 N1 w={h_mm} h={t_mm}

.external N1 N1
.freq fmin=100 fmax=100e3 ndec=3
.end"""

    t0 = time.perf_counter()
    parser = FastHenryParser()
    parser.parse_string(inp_text)
    builder = parser.to_peec_builder()
    topo = builder.build_topology()
    t_build = time.perf_counter() - t0

    n_seg = topo['n_loop']
    R_dc = np.array(topo['R'])
    L = np.array(topo['L'])

    print(f"DOFs: {n_seg} segments (filaments)")
    print(f"Build time: {t_build*1e3:.1f} ms")

    solver = PEECCircuitSolver(topo)

    # --- Case 1: Air only ---
    freqs = np.array([100.0, 1e3, 10e3, 100e3])
    Z_air = solver.frequency_sweep(freqs)

    # --- Case 2: +Ferrite (analytical image method) ---
    mu_r = 1000.0
    Delta_L = compute_delta_L(L, mu_r)
    L_with_core = L + Delta_L * (mu_r - 1.0)

    # Create solver with modified L matrix
    topo_core = dict(topo)
    topo_core['L'] = L_with_core
    solver_core = PEECCircuitSolver(topo_core)
    Z_core = solver_core.frequency_sweep(freqs)

    # --- Case 3: +Shield (BEM) ---
    from netgen.occ import Box, Pnt, OCCGeometry
    from ngsolve import Mesh as NGMesh
    from ngsolve import TaskManager
    from ngbem_eddy import ShieldBEMSIBC

    print("\nAssembling shield BEM solver...")
    t0 = time.perf_counter()
    # Shield: 12x12mm Al plate, 5mm above conductor, 0.5mm thick
    shield_plate = Box(Pnt(-0.006, -0.006, 0.005), Pnt(0.006, 0.006, 0.0055))
    shield_plate.solids.name = "conductor"
    shield_plate.faces.name = "surface"
    with TaskManager():
        shield_mesh = NGMesh(OCCGeometry(shield_plate).GenerateMesh(maxh=0.003))

        shield = ShieldBEMSIBC(shield_mesh, sigma=3.7e7)
        shield.assemble(intorder=4)
        t_shield = time.perf_counter() - t0
        print(f"Shield: {shield._loop.n_loops} loops, {shield._loop.n_active} active DOFs")
        print(f"Shield assembly: {t_shield*1e3:.0f} ms")

        # Compute shielded impedance using uniform excitation
        Z_shield = np.zeros(len(freqs), dtype=complex)
        for k, f in enumerate(freqs):
            omega = 2 * np.pi * f
            Delta_Z = shield.compute_impedance_matrix(f, topo)
            Z_branch = np.diag(R_dc.astype(complex)) + 1j * omega * L + Delta_Z
            Z_shield[k] = uniform_port_Z(Z_branch, n_seg)

        # --- Case 4: +Both (ferrite + shield) ---
        Z_both = np.zeros(len(freqs), dtype=complex)
        for k, f in enumerate(freqs):
            omega = 2 * np.pi * f
            Delta_Z = shield.compute_impedance_matrix(f, topo)
            Z_branch = np.diag(R_dc.astype(complex)) + 1j * omega * L_with_core + Delta_Z
            Z_both[k] = uniform_port_Z(Z_branch, n_seg)

        # --- Results ---
        print(f"\n{'Case':>20s}  {'f [Hz]':>10s}  {'L [nH]':>10s}  "
              f"{'R [mOhm]':>10s}  {'dL [nH]':>10s}")
        print("=" * 70)
        for k, f in enumerate(freqs):
            omega = 2 * np.pi * f
            L_a = np.imag(Z_air[k]) / omega * 1e9
            L_c = np.imag(Z_core[k]) / omega * 1e9
            L_s = np.imag(Z_shield[k]) / omega * 1e9
            L_b = np.imag(Z_both[k]) / omega * 1e9
            R_a = np.real(Z_air[k]) * 1e3
            R_c = np.real(Z_core[k]) * 1e3
            R_s = np.real(Z_shield[k]) * 1e3
            R_b = np.real(Z_both[k]) * 1e3
            print(f"{'Air':>20s}  {f:10.0f}  {L_a:10.2f}  {R_a:10.4f}  {'--':>10s}")
            print(f"{'+ Ferrite':>20s}  {f:10.0f}  {L_c:10.2f}  {R_c:10.4f}  {L_c-L_a:+10.2f}")
            print(f"{'+ Shield':>20s}  {f:10.0f}  {L_s:10.2f}  {R_s:10.4f}  {L_s-L_a:+10.2f}")
            print(f"{'+ Both':>20s}  {f:10.0f}  {L_b:10.2f}  {R_b:10.4f}  {L_b-L_a:+10.2f}")
            if k < len(freqs) - 1:
                print("-" * 70)

        # --- Physics checks ---
        L_air_vals = np.imag(Z_air) / (2 * np.pi * freqs) * 1e9
        L_core_vals = np.imag(Z_core) / (2 * np.pi * freqs) * 1e9
        L_shield_vals = np.imag(Z_shield) / (2 * np.pi * freqs) * 1e9
        R_shield_vals = np.real(Z_shield) * 1e3
        R_air_vals = np.real(Z_air) * 1e3

        dL_core = L_core_vals - L_air_vals
        dL_shield = L_shield_vals - L_air_vals
        dR_shield = R_shield_vals - R_air_vals

        print(f"\n--- Physics checks ---")
        print(f"Ferrite dL > 0:  {np.all(dL_core > 0)}  "
              f"(min dL = {np.min(dL_core):+.2f} nH)")
        print(f"Shield  dL < 0:  {np.all(dL_shield < 0)}  "
              f"(min dL = {np.min(dL_shield):+.2f} nH)")
        print(f"Shield  dR > 0:  {np.all(dR_shield > 0)}  "
              f"(min dR = {np.min(dR_shield):+.4f} mOhm)")


if __name__ == '__main__':
    main()
