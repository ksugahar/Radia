"""PEEC filament coil + BEM-SIBC workpiece coupled demo.

Phase 1 minimum: a third-reference computation of the matched-mesh
2D vs 3D inductance/heating comparison.

Pipeline
--------
1. Build a slightly-gapped CoilBuilder torus with circular cross-section.
   Gap is needed because PEEC requires a port (open-chain filament).
2. Generate N = nw x nh filament paths via CoilBuilder.to_filaments().
3. Build the PEEC bundle solver (R + jwL with full mutual M).
4. Solve for per-filament currents at the target frequency.
5. Build a workpiece BEM-SIBC mesh (Cu cylinder, surface only).
6. Compute phi_inc on workpiece nodes from the filament currents.
7. Solve scalar BIE + SIBC -> J_s on workpiece surface.
8. Report:
   - P_wp (workpiece eddy heating)
   - L_air (coil air-only inductance from PEEC, comparison to 2D L_volume)

This script does NOT yet do back-reaction (workpiece eddy current does
not feed back into PEEC currents). That is the next step.
"""

from __future__ import annotations

import math
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
SRC_RADIA = os.path.abspath(os.path.join(HERE, '..', '..', '..',
                                         'src', 'radia'))
sys.path.insert(0, SRC_RADIA)


def build_torus_filaments(R, a, current, gap_deg=5.0, nw=9, nh=9, n_arc=60,
                          frequency=0.0, sigma=5.8e7):
    """Build a CoilBuilder torus with tiny gap and emit filament paths."""
    from radia_coil_builder import CoilBuilder
    from coil_profile import CircleProfile

    # Square 2a x 2a bounding the circular wire (used by to_filaments grid).
    # The actual cross-section is circular; to_filaments samples the
    # bounding rectangle, so corner filaments will fall outside the
    # circular wire. For Phase 1 we accept that; a circular profile
    # with a true circular sampling will follow once peec_bundle is
    # extended.
    # Orientation at start:
    #   x_axis row 0 = radial direction (+x, pointing away from torus axis)
    #   y_axis row 1 = motion direction (+y, tangent to circle)
    #   z_axis row 2 = torus axis (+z)
    # This makes add_arc curve in the world xy plane -> torus around z.
    cb = (CoilBuilder(current=current)
          .set_start([R, 0, 0],
                     orientation=np.array([[1, 0, 0],
                                            [0, 1, 0],
                                            [0, 0, 1]]))
          .set_cross_section(width=2 * a, height=2 * a)
          .add_arc(radius=R, arc_angle=360.0 - gap_deg, tilt=0))

    paths, currents_tier_a = cb.to_filaments(
        nw, nh, frequency=frequency, sigma=sigma, n_arc=n_arc)
    return cb, paths, currents_tier_a


def solve_peec_bundle(filament_paths, dw, dh, sigma, current, frequency):
    """Build PEEC bundle solver and return per-filament currents."""
    from peec_bundle import build_bundle_solver

    solver, seg_of_fil, port_p, port_m = build_bundle_solver(
        filament_paths, dw, dh, sigma)

    I_branch = solver.compute_branch_currents(frequency, [current])
    # Each filament is a series chain: average over segs gives that
    # filament's current.
    n_fil = len(seg_of_fil)
    I_fil = np.zeros(n_fil, dtype=complex)
    for k, segs in enumerate(seg_of_fil):
        I_fil[k] = np.mean(I_branch[segs])
    return I_fil, solver


def build_wp_mesh(R_wp, H_wp, maxh):
    """Build a SURFACE mesh of a copper cylinder workpiece (for BEM)."""
    from netgen.occ import Cylinder, OCCGeometry, Pnt, Dir
    from ngsolve import Mesh
    cyl = Cylinder(Pnt(0, 0, -H_wp / 2), Dir(0, 0, 1), R_wp, H_wp)
    cyl.name = "wp"
    for f in cyl.faces:
        f.name = "wp_surface"
        f.maxh = maxh
    geo = OCCGeometry(cyl)
    ngmesh = geo.GenerateMesh(maxh=maxh)
    vol_mesh = Mesh(ngmesh)

    # Extract surface-only (the BIE solver requires it; volume H1 mesh
    # makes the surface mass matrix singular due to interior dofs).
    sys.path.insert(0, os.path.join(SRC_RADIA, 'panels'))
    from calc_heating_bem import _extract_surface_mesh_filtered
    surf_mesh = _extract_surface_mesh_filtered(vol_mesh, keep_label="wp")
    return surf_mesh


def main():
    R, a = 0.030, 0.003     # torus
    R_wp, H_wp = 0.025, 0.025
    sigma_cu = 5.8e7
    sigma_wp = 5.8e7
    mu_r_wp = 1.0
    freq = 7000.0
    omega = 2 * math.pi * freq

    nw, nh = 3, 3            # 81 filaments

    print("=" * 70)
    print("  PEEC bundle + BEM-SIBC coupled demo (Phase 1, no back-reaction)")
    print("=" * 70)
    print(f"  Coil:  R={R*1e3:.0f} mm, a={a*1e3:.1f} mm, "
          f"nw x nh = {nw} x {nh}")
    print(f"  WP:    R={R_wp*1e3:.0f} mm, H={H_wp*1e3:.0f} mm, Cu")
    print(f"  Freq:  {freq} Hz")
    delta = math.sqrt(2.0 / (omega * 4e-7 * math.pi * mu_r_wp * sigma_wp))
    print(f"  delta: {delta*1e3:.3f} mm  (R_wp/delta = {R_wp/delta:.1f})")
    print()

    print("[1/4] Building filament paths via CoilBuilder...")
    t0 = time.perf_counter()
    cb, paths, _ = build_torus_filaments(
        R, a, current=1.0, gap_deg=5.0, nw=nw, nh=nh, n_arc=30,
        frequency=freq, sigma=sigma_cu)
    dw = (2 * a) / nw
    dh = (2 * a) / nh
    print(f"  {len(paths)} filaments, sub-cell {dw*1e3:.3f} x {dh*1e3:.3f} mm "
          f"({time.perf_counter() - t0:.1f}s)")

    print("[2/4] PEEC solve for filament currents...")
    t0 = time.perf_counter()
    I_fil, peec_solver = solve_peec_bundle(
        paths, dw, dh, sigma_cu, current=1.0, frequency=freq)
    print(f"  PEEC solve done ({time.perf_counter() - t0:.1f}s)")
    print(f"  |I_fil| min/max/mean: {np.min(np.abs(I_fil)):.3e} / "
          f"{np.max(np.abs(I_fil)):.3e} / {np.mean(np.abs(I_fil)):.3e}")
    print(f"  sum(I_fil) = {np.sum(I_fil):.4f} (expected ~1.0)")

    print("[3/4] Building WP surface mesh + BEM-SIBC solver...")
    t0 = time.perf_counter()
    wp_mesh = build_wp_mesh(R_wp, H_wp, maxh=0.004)
    print(f"  WP mesh: {wp_mesh.nv} vertices, {wp_mesh.ne} elements "
          f"({time.perf_counter() - t0:.1f}s)")
    from bem_sibc_solver import (ScalarBIESIBCSolver,
                                  compute_phi_inc_from_filaments)
    t0 = time.perf_counter()
    bem = ScalarBIESIBCSolver(wp_mesh, order=1)
    print(f"  BEM assembly: ndof={bem.ndof} ({bem.t_assembly:.1f}s)")

    print("[4/4] phi_inc + BIE solve...")
    obs = np.array([[wp_mesh.vertices[i].point[j] for j in range(3)]
                    for i in range(wp_mesh.nv)])
    t0 = time.perf_counter()
    phi_inc = compute_phi_inc_from_filaments(obs, paths, I_fil)
    print(f"  phi_inc done ({time.perf_counter() - t0:.1f}s)")
    rho = 1.0 / sigma_wp
    Z_s = (1.0 + 1j) * rho / delta
    print(f"  Z_s = {Z_s:.4e}")

    t0 = time.perf_counter()
    res = bem.solve(phi_inc, Z_s=Z_s, omega=omega)
    print(f"  BIE solve done ({time.perf_counter() - t0:.1f}s)")
    print(f"  H_t_rms     = {res.get('H_t_rms', float('nan')):.3f} A/m")
    print(f"  P_density   = {res.get('P_density', float('nan')):.4e} W/m^2")

    A_wp = 2 * math.pi * R_wp * H_wp + 2 * math.pi * R_wp ** 2
    P_total = res.get('P_density', float('nan')) * A_wp
    print(f"  A_wp        = {A_wp:.4e} m^2 (analytical)")
    print(f"  P_total     = {P_total:.4e} W")
    print()
    print("Reference values (matched-mesh test, Cu @ 7 kHz, R/delta=31.7):")
    print("  2D SIBC:  L=57.61 nH, P=6.63e-05 W")
    print("  3D FEM:   L=56.94 nH, P=5.78e-05 W  (-12.81% P)")
    print(f"  PEEC+BEM: P={P_total:.3e} W  (vs 2D dP={(P_total - 6.63e-5)/6.63e-5*100:+.2f}%)")


if __name__ == "__main__":
    main()
