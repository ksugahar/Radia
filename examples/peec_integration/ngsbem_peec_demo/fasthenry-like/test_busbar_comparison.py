"""
Bus bar Dowell comparison: ngbem PEEC vs FastHenry PEEC

Same geometry for both:
  - Rectangular bus bar: 100 mm x 10 mm x 1 mm copper
  - sigma = 5.8e7 S/m
  - Single port (current in at one end, out at other end)

ngbem PEEC:  plate surface mesh (100x10 mm), thickness=1 mm, Dowell Zs_func
FastHenry:   single filament segment, w=10mm, h=1mm, L=100mm, Dowell Zs_func

Expected DC resistance: R_dc = L/(sigma*w*h) = 0.1/(5.8e7*0.01*0.001) = 1.724 uOhm
"""

import sys
import os
import time
import numpy as np

# Path setup
script_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, script_dir)
sys.path.insert(0, os.path.join(script_dir, '../../../src/radia'))

MU_0 = 4.0 * np.pi * 1e-7

# ============================================================
# Geometry parameters
# ============================================================
length = 0.100    # 100 mm
bar_width = 0.010  # 10 mm
bar_thick = 0.001  # 1 mm
sigma = 5.8e7      # copper S/m

R_dc_analytical = length / (sigma * bar_width * bar_thick)
print(f"Bus bar: {length*1e3:.0f} mm x {bar_width*1e3:.0f} mm x {bar_thick*1e3:.1f} mm Cu")
print(f"R_dc (analytical) = {R_dc_analytical*1e6:.3f} uOhm")
print()

# Frequencies where Dowell matters
freqs_test = np.logspace(1, 7, 60)  # 10 Hz to 10 MHz
freqs_key = np.array([1e3, 10e3, 100e3, 1e6])

# Print Dowell F_R at key frequencies
from ngbem_peec import dowell_F_R
print("Dowell F_R for 1 mm Cu bus bar:")
print(f"  {'Freq':>10s}  {'delta [um]':>10s}  {'xi':>8s}  {'F_R':>8s}")
print("  " + "-" * 42)
for f in freqs_key:
    delta = np.sqrt(2.0 / (2*np.pi*f * MU_0 * sigma))
    xi = bar_thick / (2 * delta)
    F = dowell_F_R(xi)
    print(f"  {f:10.0f}  {delta*1e6:10.0f}  {xi:8.3f}  {F:8.4f}")
print()

# ============================================================
# A) ngbem PEEC (surface mesh BEM)
# ============================================================
print("=" * 70)
print("A) ngbem PEEC (surface mesh, Galerkin BEM)")
print("=" * 70)

from ngbem_peec import (NGBEMPEECSolver, create_plate_mesh,
                         make_dowell_Zs_func)

mesh = create_plate_mesh(length, bar_width, maxh=0.008, label="conductor")

solver_ngbem = NGBEMPEECSolver(mesh, conductor_label="conductor",
                                sigma=sigma, thickness=bar_thick,
                                order=0, intorder=5)
t0 = time.perf_counter()
solver_ngbem.assemble()
t_ngbem_asm = time.perf_counter() - t0

print(f"  DOFs: {solver_ngbem.n_loop} loop, {solver_ngbem.n_star} star")
print(f"  Assembly: {t_ngbem_asm:.3f} s")
print(f"  R_dc (per-edge mean): {np.mean(solver_ngbem.R_loop)*1e6:.3f} uOhm")

# DC only sweep
t0 = time.perf_counter()
Z_ngbem_dc = solver_ngbem.solve_frequency(freqs_test, mode='mqs')
t_ngbem_dc_sweep = time.perf_counter() - t0

# Dowell sweep
Zs_ngbem = make_dowell_Zs_func(
    n_loop=solver_ngbem.n_loop,
    R_dc_per_edge=solver_ngbem.R_loop,
    thickness=bar_thick,
    sigma=sigma,
)
t0 = time.perf_counter()
Z_ngbem_dw = solver_ngbem.solve_frequency(freqs_test, mode='mqs', Zs_func=Zs_ngbem)
t_ngbem_dw_sweep = time.perf_counter() - t0

print(f"  Sweep DC:     {t_ngbem_dc_sweep*1e3:.1f} ms")
print(f"  Sweep Dowell: {t_ngbem_dw_sweep*1e3:.1f} ms")

# ============================================================
# B) FastHenry PEEC (filament, Neumann formula)
# ============================================================
print(f"\n{'=' * 70}")
print("B) FastHenry PEEC (filament, Neumann formula)")
print("=" * 70)

import radia as rad
from fasthenry_parser import FastHenryParser
from peec_topology import PEECCircuitSolver

rad.UtiDelAll()

# Single segment: same bus bar
inp_text = f""".Units mm
N1 x=0   y=0 z=0
N2 x={length*1e3:.1f} y=0 z=0

E1 N1 N2 w={bar_width*1e3:.1f} h={bar_thick*1e3:.1f} sigma={sigma:.2e}

.external N1 N2
.freq fmin=10 fmax=10e6 ndec=10
.end"""

t0 = time.perf_counter()
parser = FastHenryParser()
parser.parse_string(inp_text)
builder = parser.to_peec_builder()
topo = builder.build_topology()
t_fh_build = time.perf_counter() - t0

fh_solver = PEECCircuitSolver(topo)
n_seg = topo['n_loop']
seg_lengths = np.asarray(topo['segment_lengths'])

print(f"  DOFs: {n_seg} segment(s)")
print(f"  Build: {t_fh_build*1e3:.1f} ms")
print(f"  Segment length: {seg_lengths[0]*1e3:.1f} mm")

# DC only sweep
t0 = time.perf_counter()
Z_fh_dc = fh_solver.frequency_sweep(freqs_test)
t_fh_dc_sweep = time.perf_counter() - t0

# Dowell sweep
R_dc_fh = np.array([length / (sigma * bar_width * bar_thick)])

def dowell_Zs_fh(freq):
    if freq <= 0:
        return np.zeros(n_seg, dtype=complex)
    omega = 2.0 * np.pi * freq
    delta = np.sqrt(2.0 / (omega * MU_0 * sigma))
    xi = bar_thick / (2.0 * delta)
    F_R = dowell_F_R(xi)
    return R_dc_fh * (F_R - 1.0) + 0j

t0 = time.perf_counter()
Z_fh_dw = fh_solver.frequency_sweep(freqs_test, Zs_func=dowell_Zs_fh)
t_fh_dw_sweep = time.perf_counter() - t0

print(f"  R_dc: {R_dc_fh[0]*1e6:.3f} uOhm")
print(f"  Sweep DC:     {t_fh_dc_sweep*1e3:.1f} ms")
print(f"  Sweep Dowell: {t_fh_dw_sweep*1e3:.1f} ms")

rad.UtiDelAll()

# ============================================================
# C) Comparison at key frequencies
# ============================================================
print(f"\n{'=' * 70}")
print("C) Comparison at key frequencies")
print("=" * 70)

print(f"\n  {'Freq':>10s}  {'':3s}  {'L_ngbem[nH]':>12s}  {'L_fh[nH]':>10s}  "
      f"{'L_diff%':>8s}  {'R_ngbem[uOhm]':>14s}  {'R_fh[uOhm]':>12s}  {'R_diff%':>8s}")
print("  " + "-" * 90)

for f in freqs_key:
    idx = np.argmin(np.abs(freqs_test - f))

    # DC resistance results
    L_ng_dc = np.imag(Z_ngbem_dc[idx]) / (2*np.pi*freqs_test[idx]) * 1e9
    L_fh_dc = np.imag(Z_fh_dc[idx]) / (2*np.pi*freqs_test[idx]) * 1e9
    R_ng_dc = np.real(Z_ngbem_dc[idx]) * 1e6
    R_fh_dc_v = np.real(Z_fh_dc[idx]) * 1e6

    L_diff = (L_ng_dc - L_fh_dc) / abs(L_fh_dc) * 100 if abs(L_fh_dc) > 1e-15 else 0
    R_diff = (R_ng_dc - R_fh_dc_v) / abs(R_fh_dc_v) * 100 if abs(R_fh_dc_v) > 1e-15 else 0

    print(f"  {f:10.0f}  DC   {L_ng_dc:12.3f}  {L_fh_dc:10.3f}  "
          f"{L_diff:+7.1f}%  {R_ng_dc:14.3f}  {R_fh_dc_v:12.3f}  {R_diff:+7.1f}%")

    # Dowell results
    L_ng_dw = np.imag(Z_ngbem_dw[idx]) / (2*np.pi*freqs_test[idx]) * 1e9
    L_fh_dw = np.imag(Z_fh_dw[idx]) / (2*np.pi*freqs_test[idx]) * 1e9
    R_ng_dw = np.real(Z_ngbem_dw[idx]) * 1e6
    R_fh_dw_v = np.real(Z_fh_dw[idx]) * 1e6

    L_diff_dw = (L_ng_dw - L_fh_dw) / abs(L_fh_dw) * 100 if abs(L_fh_dw) > 1e-15 else 0
    R_diff_dw = (R_ng_dw - R_fh_dw_v) / abs(R_fh_dw_v) * 100 if abs(R_fh_dw_v) > 1e-15 else 0

    print(f"  {'':10s}  Dwl  {L_ng_dw:12.3f}  {L_fh_dw:10.3f}  "
          f"{L_diff_dw:+7.1f}%  {R_ng_dw:14.3f}  {R_fh_dw_v:12.3f}  {R_diff_dw:+7.1f}%")

# ============================================================
# D) Dowell F_R ratio check
# ============================================================
print(f"\n{'=' * 70}")
print("D) Dowell R_ac/R_dc ratio")
print("=" * 70)

print(f"\n  {'Freq':>10s}  {'ngbem':>8s}  {'FastHenry':>10s}  {'Dowell F_R':>10s}")
print("  " + "-" * 45)

for f in freqs_key:
    idx = np.argmin(np.abs(freqs_test - f))
    ratio_ng = np.real(Z_ngbem_dw[idx]) / np.real(Z_ngbem_dc[idx])
    ratio_fh = np.real(Z_fh_dw[idx]) / np.real(Z_fh_dc[idx])
    delta = np.sqrt(2.0 / (2*np.pi*f * MU_0 * sigma))
    xi = bar_thick / (2 * delta)
    F_R = dowell_F_R(xi)
    print(f"  {f:10.0f}  {ratio_ng:8.4f}  {ratio_fh:10.4f}  {F_R:10.4f}")

# ============================================================
# E) Timing summary
# ============================================================
print(f"\n{'=' * 70}")
print("E) Timing summary")
print("=" * 70)

t_ngbem_total = t_ngbem_asm + t_ngbem_dw_sweep
t_fh_total = t_fh_build + t_fh_dw_sweep

print(f"\n  {'':20s}  {'ngbem PEEC':>12s}  {'FastHenry':>12s}  {'Ratio':>8s}")
print("  " + "-" * 58)
print(f"  {'Assembly/Build':20s}  {t_ngbem_asm*1e3:10.1f} ms  {t_fh_build*1e3:10.1f} ms  "
      f"{t_ngbem_asm/t_fh_build:7.0f}x")
print(f"  {'Dowell sweep':20s}  {t_ngbem_dw_sweep*1e3:10.1f} ms  {t_fh_dw_sweep*1e3:10.1f} ms  "
      f"{t_ngbem_dw_sweep/t_fh_dw_sweep:7.1f}x")
print(f"  {'Total':20s}  {t_ngbem_total*1e3:10.1f} ms  {t_fh_total*1e3:10.1f} ms  "
      f"{t_ngbem_total/t_fh_total:7.0f}x")
print(f"  {'DOFs':20s}  {solver_ngbem.n_loop:10d}     {n_seg:10d}")

print("\nDone.")
