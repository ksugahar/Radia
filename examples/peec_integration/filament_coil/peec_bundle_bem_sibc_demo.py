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


def extract_wp_surface_J(mesh, phi_vec_complex, omega, Z_s):
    """Extract per-element complex surface current density J_s.

    SIBC scalar BIE: phi on surface is magnetic scalar potential. H_t =
    -grad_s(phi); J_s = n x H_t = -n x grad_s(phi) (tangential vector).

    Returns (centroids, areas, J_s) as (M,3), (M,), (M,3) arrays; J_s
    is complex because phi_vec is complex.
    """
    from ngsolve import (H1, GridFunction, Integrate, CF, BND,
                          InnerProduct, grad, specialcf, Cross)
    fes = H1(mesh, order=1)
    gf_re = GridFunction(fes)
    gf_im = GridFunction(fes)
    gf_re.vec.FV().NumPy()[:] = phi_vec_complex.real
    gf_im.vec.FV().NumPy()[:] = phi_vec_complex.imag

    n_cf = specialcf.normal(3)
    # J_s_complex = -n x grad_s(phi) = -n x (grad phi_re + j grad phi_im)
    Jre_cf = -Cross(n_cf, grad(gf_re).Trace())
    Jim_cf = -Cross(n_cf, grad(gf_im).Trace())

    # Per-element averages via element_wise integrate / area
    elem_A = Integrate(CF(1), mesh, VOL_or_BND=BND, element_wise=True)
    Jre = [Integrate(Jre_cf[i], mesh, VOL_or_BND=BND, element_wise=True)
           for i in range(3)]
    Jim = [Integrate(Jim_cf[i], mesh, VOL_or_BND=BND, element_wise=True)
           for i in range(3)]

    centroids, areas, J_s = [], [], []
    for el in mesh.Elements(BND):
        area = abs(elem_A[el.nr])
        if area < 1e-30:
            continue
        verts = [mesh.vertices[v.nr].point for v in el.vertices]
        c = np.mean([(v[0], v[1], v[2]) for v in verts], axis=0)
        Jvec = np.array([(Jre[i][el.nr] + 1j * Jim[i][el.nr]) / area
                         for i in range(3)], dtype=complex)
        centroids.append(c)
        areas.append(area)
        J_s.append(Jvec)
    return np.array(centroids), np.array(areas), np.array(J_s)


def A_back_at_points(obs_points, wp_c, wp_a, wp_J):
    """Vector potential from wp surface currents at obs points.

    A(r) = (mu_0/4pi) sum_j wp_J[j] * wp_a[j] / |r - wp_c[j]|

    Args:
        obs_points: (N,3) float
        wp_c: (M,3) panel centroids
        wp_a: (M,) panel areas
        wp_J: (M,3) complex panel surface current density

    Returns:
        (N,3) complex A_back.
    """
    r = np.asarray(obs_points, float)
    c = np.asarray(wp_c, float)
    a = np.asarray(wp_a, float)
    J = np.asarray(wp_J, complex)

    # Pair-wise distances (N, M)
    dx = r[:, None, 0] - c[None, :, 0]
    dy = r[:, None, 1] - c[None, :, 1]
    dz = r[:, None, 2] - c[None, :, 2]
    dist = np.sqrt(dx * dx + dy * dy + dz * dz)
    dist[dist < 1e-12] = 1e-12
    MU_0 = 4e-7 * math.pi
    weight = (MU_0 / (4 * math.pi)) * (a / dist)   # (N, M)
    A = np.einsum('nm,mc->nc', weight, J)          # (N, 3) complex
    return A


def compute_flux_linkage(filament_paths, wp_c, wp_a, wp_J):
    """Flux Phi_k = integral of A_back.dl along each filament."""
    N = len(filament_paths)
    Phi = np.zeros(N, dtype=complex)
    for k, path in enumerate(filament_paths):
        midpoints = np.array([0.5 * (np.array(p1) + np.array(p2))
                              for p1, p2 in path])
        dls = np.array([np.array(p2) - np.array(p1) for p1, p2 in path])
        A_mid = A_back_at_points(midpoints, wp_c, wp_a, wp_J)
        Phi[k] = np.sum(A_mid * dls)
    return Phi


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
    print("[5/5] Back-reaction: extract J_s, compute flux linkage, DeltaL...")
    t0 = time.perf_counter()
    wp_c, wp_a, wp_J = extract_wp_surface_J(
        wp_mesh, res['phi_vec'], omega, Z_s)
    print(f"  Extracted J_s on {len(wp_c)} wp elements "
          f"({time.perf_counter() - t0:.1f}s)")

    # Flux linkage per filament: Phi_k = integral A_back . dl
    t0 = time.perf_counter()
    Phi = compute_flux_linkage(paths, wp_c, wp_a, wp_J)
    print(f"  Flux linkage per filament "
          f"({time.perf_counter() - t0:.1f}s):")
    print(f"    |Phi| min/max/mean = {np.min(np.abs(Phi)):.3e} / "
          f"{np.max(np.abs(Phi)):.3e} / {np.mean(np.abs(Phi)):.3e}")

    # Port impedance change from wp mutual coupling:
    #   Delta_Z_port = +jw * sum_k (I_k * Phi_wp_k) / I_port^2
    # The +jw sign follows from mutual flux linkage V = jw Phi adding
    # to the port voltage (NOT the -jw EMF notation). Validated by the
    # sign of Delta_L < 0 for Cu (Lenz) and Delta_R > 0 for eddy loss.
    I_port = 1.0
    Sum_IPhi = np.sum(I_fil * Phi)
    Delta_Z = 1j * omega * Sum_IPhi / I_port ** 2
    Delta_L = Delta_Z.imag / omega
    Delta_R = Delta_Z.real
    print(f"  Delta_Z = {Delta_Z:.4e}")
    print(f"  Delta_L = {Delta_L * 1e9:+.3f} nH  (wp-induced inductance change)")
    print(f"  Delta_R = {Delta_R * 1e3:+.4f} mOhm  (wp-induced resistance)")

    # Air-only L from PEEC (port impedance at this frequency)
    Z_port = peec_solver.compute_port_impedance(freq)
    L_air = Z_port.imag / omega
    R_coil = Z_port.real
    L_total = L_air + Delta_L
    print(f"  L_air (PEEC)     = {L_air * 1e9:.3f} nH")
    print(f"  L_total (+ wp)   = {L_total * 1e9:.3f} nH")

    print()
    print("Reference values (matched-mesh test, Cu @ 7 kHz, R/delta=31.7):")
    print(f"  {'solver':<12s} {'L [nH]':>10s} {'P [W]':>12s}")
    print(f"  {'2D SIBC':<12s} {57.61:>10.2f} {6.63e-5:>12.2e}")
    print(f"  {'3D FEM':<12s} {56.94:>10.2f} {5.78e-5:>12.2e}")
    print(f"  {'PEEC+BEM':<12s} {L_total * 1e9:>10.2f} "
          f"{P_total:>12.2e}")
    print(f"  dL vs 2D = {(L_total * 1e9 - 57.61) / 57.61 * 100:+.2f}%")
    print(f"  dP vs 2D = {(P_total - 6.63e-5) / 6.63e-5 * 100:+.2f}%")


if __name__ == "__main__":
    main()
