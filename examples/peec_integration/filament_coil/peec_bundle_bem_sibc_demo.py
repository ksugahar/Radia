"""peec_bundle_bem_sibc_demo.py

Third reference for the matched-mesh 2D axisym vs 3D FEM-SIBC
inductance / workpiece-heating comparison:

  (a) coil     : filament-PEEC bundle (peec_bundle.build_bundle_solver)
  (b) workpiece: BEM-SIBC scalar-BIE (ScalarBIESIBCSolver) on a surface
                 mesh of a copper cylinder
  (c) back-reaction: wp surface J_s -> vector potential A at the coil
                     filament midpoints -> flux linkage per filament
                     -> port-impedance change Delta_Z, Delta_L, Delta_R

Refactored 2026-04-15 to reuse existing helpers:
  - compute_phi_inc_from_filaments  (bem_sibc_solver)
  - ScalarBIESIBCSolver             (bem_sibc_solver)
  - extract_surface_J_from_phi      (bem_sibc_solver, new)
  - flux_linkage_in_filaments       (bem_sibc_solver, new)
  - _extract_surface_mesh_filtered  (radia.panels.surface_mesh_extract,
                                     to strip an OCC volume mesh to
                                     surface-only)
  - peec_bundle.build_bundle_solver

Usage:
    python peec_bundle_bem_sibc_demo.py [--nw 3] [--nh 3] [--wp-maxh 0.004]

Matched-mesh reference (Cu @ 7 kHz, R/delta = 31.7):

    solver       L [nH]   P [W]
    2D SIBC       57.61   6.63e-05
    3D FEM-SIBC   56.94   5.78e-05
    PEEC + BEM    (this script's output)
"""

from __future__ import annotations

import argparse
import math
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.abspath(os.path.join(HERE, '..', '..', '..', 'src'))
SRC_RADIA = os.path.join(SRC, 'radia')
PANELS = os.path.join(SRC_RADIA, 'panels')
for p in (SRC, SRC_RADIA, PANELS):
    if p not in sys.path:
        sys.path.insert(0, p)

MU_0 = 4e-7 * math.pi


def build_torus_filaments(R, a, current, gap_deg, nw, nh, n_arc,
                          frequency, sigma):
    """CoilBuilder circular torus (z-axis) -> filament paths."""
    from coil_builder import CoilBuilder
    cb = (CoilBuilder(current=current)
          .set_start([R, 0, 0],
                     orientation=np.array([[1, 0, 0],
                                            [0, 1, 0],
                                            [0, 0, 1]]))
          .set_cross_section(width=2 * a, height=2 * a)
          .add_arc(radius=R, arc_angle=360.0 - gap_deg, tilt=0))
    paths, _ = cb.to_filaments(
        nw, nh, frequency=frequency, sigma=sigma, n_arc=n_arc)
    return cb, paths


def build_wp_surface_mesh(R_wp, H_wp, maxh):
    """Copper cylinder -> 2D surface mesh for BEM-SIBC.

    Extracts the boundary of an OCC-meshed volume via the standard
    radia panel helper. `keep_label='wp'` filters out any phantom
    faces (there should be none for a standalone cylinder).
    """
    from netgen.occ import Cylinder, OCCGeometry, Pnt, Dir
    from ngsolve import Mesh
    from radia.panels.surface_mesh_extract import _extract_surface_mesh_filtered

    cyl = Cylinder(Pnt(0, 0, -H_wp / 2), Dir(0, 0, 1), R_wp, H_wp)
    cyl.name = "wp"
    for f in cyl.faces:
        f.name = "wp_surface"
        f.maxh = maxh
    geo = OCCGeometry(cyl)
    with TaskManager():
        vol_mesh = Mesh(geo.GenerateMesh(maxh=maxh))
        return _extract_surface_mesh_filtered(vol_mesh, keep_label="wp")


def solve_peec_bundle(paths, dw, dh, sigma, I_port, frequency):
    """Build PEEC bundle + Loop-form reduction, return per-fil currents.

    For a parallel-filament series-chain bundle the nodal MNA collapses
    exactly to a K x K Loop-form (K = #filaments) via
    build_loop_bundle_impedance. This scales O(N_seg^3) -> O(K^3) with
    zero approximation error (verified in loop_bundle_verify.py).

    Returns (R_f, L_f, I_fil, V_port). V_port is the port voltage when
    driven by the specified I_port current source, so the air-only port
    impedance is V_port / I_port.
    """
    from peec_bundle import (build_bundle_solver,
                              build_loop_bundle_impedance,
                              solve_loop_bundle)
    solver, seg_of_fil, _, _ = build_bundle_solver(paths, dw, dh, sigma)
    R_f, L_f = build_loop_bundle_impedance(solver, seg_of_fil)
    I_fil, V_port = solve_loop_bundle(R_f, L_f, frequency, I_port=I_port)
    return R_f, L_f, I_fil, V_port


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--nw", type=int, default=3)
    ap.add_argument("--nh", type=int, default=3)
    ap.add_argument("--n-arc", type=int, default=30)
    ap.add_argument("--wp-maxh", type=float, default=0.004)
    ap.add_argument("--freq", type=float, default=7000.0)
    args = ap.parse_args()

    R, a = 0.030, 0.003
    R_wp, H_wp = 0.025, 0.025
    sigma_cu = 5.8e7
    sigma_wp = 5.8e7
    mu_r_wp = 1.0
    omega = 2 * math.pi * args.freq
    delta = math.sqrt(2.0 / (omega * mu_r_wp * MU_0 * sigma_wp))

    print("=" * 70)
    print("  PEEC bundle + BEM-SIBC  (matched-mesh 2D vs 3D third reference)")
    print("=" * 70)
    print(f"  Coil : R={R*1e3:.0f}mm  a={a*1e3:.1f}mm  "
          f"nw x nh = {args.nw} x {args.nh}")
    print(f"  WP   : cylinder R={R_wp*1e3:.0f}mm H={H_wp*1e3:.0f}mm  Cu")
    print(f"  Freq : {args.freq:.0f} Hz, delta={delta*1e3:.3f} mm "
          f"(R/delta = {R_wp / delta:.1f})")

    print("\n[1/5] CoilBuilder -> filament paths")
    t0 = time.perf_counter()
    _, paths = build_torus_filaments(
        R, a, 1.0, gap_deg=5.0,
        nw=args.nw, nh=args.nh, n_arc=args.n_arc,
        frequency=args.freq, sigma=sigma_cu)
    dw = (2 * a) / args.nw
    dh = (2 * a) / args.nh
    print(f"  {len(paths)} filaments, sub-cell {dw*1e3:.3f} x {dh*1e3:.3f} mm "
          f"({time.perf_counter() - t0:.1f}s)")

    print("\n[2/5] PEEC bundle + Loop-form reduction -> per-filament currents")
    t0 = time.perf_counter()
    R_f, L_f, I_fil, V_port = solve_peec_bundle(
        paths, dw, dh, sigma_cu, I_port=1.0, frequency=args.freq)
    # With I_port = 1 A, V_port IS the port impedance.
    Z_port = V_port
    L_air = Z_port.imag / omega
    R_coil = Z_port.real
    print(f"  done ({time.perf_counter() - t0:.1f}s, K x K = "
          f"{R_f.shape[0]} x {R_f.shape[0]})")
    print(f"  sum(I_fil) = {np.sum(I_fil):.4f}  "
          f"|I_fil| mean = {np.mean(np.abs(I_fil)):.3e}")
    print(f"  L_air (PEEC, vacuum) = {L_air * 1e9:.3f} nH")
    print(f"  R_coil              = {R_coil * 1e3:.4f} mOhm")

    print("\n[3/5] Workpiece surface mesh + ScalarBIESIBCSolver")
    t0 = time.perf_counter()
    wp_mesh = build_wp_surface_mesh(R_wp, H_wp, args.wp_maxh)
    from ngsolve import BND
    from ngsolve import TaskManager
    print(f"  wp nv={wp_mesh.nv}  ne(BND)={wp_mesh.GetNE(BND)}  "
          f"({time.perf_counter() - t0:.1f}s)")
    from radia.bem_sibc_solver import (ScalarBIESIBCSolver,
                                  compute_phi_inc_from_filaments,
                                  extract_surface_J_from_phi,
                                  flux_linkage_in_filaments)
    t0 = time.perf_counter()
    bem = ScalarBIESIBCSolver(wp_mesh, order=1)
    print(f"  BEM assembly: ndof={bem.ndof} ({bem.t_assembly:.1f}s)")

    print("\n[4/5] phi_inc from filaments -> BIE + SIBC")
    obs = np.array([[wp_mesh.vertices[i].point[j] for j in range(3)]
                    for i in range(wp_mesh.nv)])
    t0 = time.perf_counter()
    phi_inc = compute_phi_inc_from_filaments(obs, paths, I_fil)
    print(f"  phi_inc ({time.perf_counter() - t0:.1f}s)")
    rho_wp = 1.0 / sigma_wp
    Z_s = (1.0 + 1j) * rho_wp / delta
    t0 = time.perf_counter()
    res = bem.solve(phi_inc, Z_s=Z_s, omega=omega)
    print(f"  BIE solve ({time.perf_counter() - t0:.1f}s)")

    A_wp = 2 * math.pi * R_wp * H_wp + 2 * math.pi * R_wp ** 2
    P_total = res['P_density'] * A_wp
    print(f"  H_t_rms  = {res['H_t_rms']:.2f} A/m")
    print(f"  P_total  = {P_total:.4e} W")

    print("\n[5/5] Back-reaction -> Delta_L, Delta_R")
    t0 = time.perf_counter()
    wp_c, wp_a, wp_J = extract_surface_J_from_phi(wp_mesh, res['phi_vec'])
    Phi = flux_linkage_in_filaments(paths, wp_c, wp_a, wp_J)
    Delta_Z = 1j * omega * np.sum(I_fil * Phi) / (1.0 ** 2)
    Delta_L = Delta_Z.imag / omega
    Delta_R = Delta_Z.real
    L_total = L_air + Delta_L
    print(f"  |Phi| min/max/mean = {np.min(np.abs(Phi)):.3e} / "
          f"{np.max(np.abs(Phi)):.3e} / {np.mean(np.abs(Phi)):.3e}")
    print(f"  Delta_L = {Delta_L * 1e9:+.3f} nH  "
          f"(Delta_R = {Delta_R * 1e3:+.4f} mOhm)")
    print(f"  L_total = L_air + Delta_L = {L_total * 1e9:.3f} nH  "
          f"({time.perf_counter() - t0:.1f}s)")

    print()
    print("-" * 70)
    print(f"  {'solver':<14s} {'L [nH]':>10s} {'P [W]':>12s}")
    print(f"  {'2D SIBC':<14s} {57.61:>10.2f} {6.63e-5:>12.2e}")
    print(f"  {'3D FEM-SIBC':<14s} {56.94:>10.2f} {5.78e-5:>12.2e}")
    print(f"  {'PEEC + BEM':<14s} {L_total * 1e9:>10.2f} "
          f"{P_total:>12.2e}")
    print(f"  dL vs 2D = {(L_total * 1e9 - 57.61) / 57.61 * 100:+.2f}%    "
          f"dP vs 2D = {(P_total - 6.63e-5) / 6.63e-5 * 100:+.2f}%")


if __name__ == "__main__":
    main()
