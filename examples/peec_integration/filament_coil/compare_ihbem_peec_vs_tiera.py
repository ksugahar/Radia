"""compare_ihbem_peec_vs_tiera.py

Phase 2 validation on the ih_bem_sample workpiece: compare the downstream
heating response (H_t_rms, P_total, Q_total) for THREE ways of feeding
phi_inc into the SIBC solver on the same workpiece mesh:

  (a) BEM truth           -- EFIE on coil surface -> surface J -> phi_inc.
                             DC-uniform coil current; workpiece physics only.
  (b) Tier A filaments    -- analytical 1/R_k * exp(-(1+j)d/delta) currents
                             on the CoilBuilder filament discretization.
  (c) PEEC filaments      -- full (R + jwM) MNA solve on the SAME filament
                             discretization (AC-truth for the coil alone),
                             then phi_inc from those complex currents.

Why this matters
----------------
Phase 1 showed Tier A vs BEM-truth differ by ~1.5% on H_t_rms / ~3% on
P for this geometry. Phase 1b showed Tier A vs PEEC on the bare coil
differ by ~40-100% in per-filament currents at 7 kHz. The open question:
does the workpiece "average out" the per-filament differences, or does
switching to PEEC currents close / widen / flip the gap in the workpiece
quantities?

If (c) lands between (a) and (b) or nearly on top of (a), the workpiece
is insensitive to filament-level redistribution -> Tier A is fine for
IH design. If (c) drifts far from (a), we need PEEC currents (or a
coupled coil+workpiece solve) even for the workpiece metrics.

Usage
-----
    python compare_ihbem_peec_vs_tiera.py
"""

import os
import sys
import time
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.abspath(os.path.join(HERE, '..', '..', '..', 'src'))
SRC_RADIA = os.path.join(SRC, 'radia')
PANELS = os.path.join(SRC, 'radia', 'panels')
for p in (SRC, SRC_RADIA, PANELS, HERE):
    if p not in sys.path:
        sys.path.insert(0, p)

from ngsolve import Mesh, Integrate, CF, BND
from ngsolve import TaskManager

from radia.coil_builder import CoilBuilder
from radia.bem_sibc_solver import (ScalarBIESIBCSolver,
                                    compute_phi_inc_from_surface_J,
                                    compute_phi_inc_from_filaments)
from radia.bem_inductance import compute_inductance_source_sink

from calc_heating_bem import _extract_surface_mesh_filtered

from peec_bundle import build_bundle_solver

MU_0 = 4.0 * np.pi * 1e-7
MM = 1e-3

# --- Match ih_bem_sample.jou ---------------------------------------
VOL_FILE = os.path.join(SRC, 'radia', 'panels', 'samples',
                        'ih_bem_sample.vol')

R_MAJOR = 30 * MM
R_MINOR = 3 * MM
GAP_DEG = 5.0
ARC_DEG = 360.0 - GAP_DEG

COIL_CURRENT = 1.0               # A
FREQUENCY = 7e3                  # Hz

WP_SIGMA = 2e6                   # S/m, steel
WP_MU_R = 100.0                  # relative

CU_SIGMA = 5.8e7                 # coil
CU_MU_R = 1.0

TIER_A_SIDE = np.sqrt(np.pi) * R_MINOR
NW = 5
NH = 5
N_ARC = 40


# ------------------------------------------------------------------
# Mesh loading
# ------------------------------------------------------------------

def load_meshes():
    print(f"Loading {VOL_FILE}")
    mesh_full = Mesh(VOL_FILE)
    coil_mesh = _extract_surface_mesh_filtered(mesh_full, keep_label='coil')
    wp_mesh = _extract_surface_mesh_filtered(mesh_full, keep_label='workpiece')
    print(f"  coil surface: nse={coil_mesh.GetNE(BND)}, nv={coil_mesh.nv}")
    print(f"  wp   surface: nse={wp_mesh.GetNE(BND)},   nv={wp_mesh.nv}")
    return coil_mesh, wp_mesh


# ------------------------------------------------------------------
# (a) BEM truth: EFIE coil -> surface J -> phi_inc
# ------------------------------------------------------------------

def phi_inc_from_bem(coil_mesh, wp_nodes, current):
    print("\n[BEM truth] EFIE on coil (source/sink)...")
    t0 = time.perf_counter()
    sol = compute_inductance_source_sink(
        coil_mesh, source_label='source', sink_label='sink', fes_order=0)
    if 'error' in sol:
        raise RuntimeError(f"Coil EFIE failed: {sol['error']}")
    L_coil = sol['L']
    gf_J = sol['gf_J']
    print(f"  L_coil = {L_coil*1e9:.2f} nH, "
          f"t_solve = {sol['t_solve']:.1f} s")

    with TaskManager():
        elem_A = Integrate(CF(1), coil_mesh, VOL_or_BND=BND, element_wise=True)
        elem_Jx = Integrate(gf_J[0], coil_mesh, VOL_or_BND=BND, element_wise=True)
        elem_Jy = Integrate(gf_J[1], coil_mesh, VOL_or_BND=BND, element_wise=True)
        elem_Jz = Integrate(gf_J[2], coil_mesh, VOL_or_BND=BND, element_wise=True)

        centroids, areas, Jvecs = [], [], []
        for el in coil_mesh.Elements(BND):
            a = abs(elem_A[el.nr])
            if a < 1e-30:
                continue
            j = np.array([elem_Jx[el.nr], elem_Jy[el.nr],
                          elem_Jz[el.nr]]) / a
            verts = [coil_mesh.vertices[v.nr].point for v in el.vertices]
            c = np.mean([(v[0], v[1], v[2]) for v in verts], axis=0)
            centroids.append(c)
            areas.append(a)
            Jvecs.append(j)

        centroids = np.asarray(centroids)
        areas = np.asarray(areas)
        Jvecs = np.asarray(Jvecs) * current

        print(f"  {len(centroids)} source panels, "
              f"total area {areas.sum()*1e6:.1f} mm^2")
        print("  integrating phi_inc...")
        t0 = time.perf_counter()
        phi = compute_phi_inc_from_surface_J(
            wp_nodes, centroids, areas, Jvecs, n_quad=20)
        print(f"  phi range [{phi.min():.4f}, {phi.max():.4f}] "
              f"({time.perf_counter()-t0:.1f} s)")
        return phi


# ------------------------------------------------------------------
# Common: build filament discretization via CoilBuilder
# ------------------------------------------------------------------

def build_coil_and_filaments(current, frequency_for_tier_a):
    """Returns (coil, filament_paths, I_tier_a, dw, dh).

    Geometry: circular arc 355 deg, area-equivalent square cross-section.
    filament_paths + I_tier_a are both outputs of CoilBuilder.to_filaments
    at the specified frequency (Tier A currents). The same filament_paths
    are reused as input to the PEEC bundle solver.
    """
    coil = (CoilBuilder(current=current)
            .set_start([R_MAJOR, 0.0, 0.0])
            .set_cross_section(width=TIER_A_SIDE, height=TIER_A_SIDE)
            .add_arc(radius=R_MAJOR, arc_angle=ARC_DEG))
    paths, I_tier_a = coil.to_filaments(
        nw=NW, nh=NH, frequency=frequency_for_tier_a,
        sigma=CU_SIGMA, mu_r=CU_MU_R, n_arc=N_ARC)
    info = coil._last_filament_info
    dw = info['w'] / NW
    dh = info['h'] / NH
    return coil, paths, I_tier_a, dw, dh


# ------------------------------------------------------------------
# (b) Tier A: analytical filament currents -> phi_inc
# ------------------------------------------------------------------

def phi_inc_from_tier_a(wp_nodes, paths, I_tier_a, frequency):
    print(f"\n[Tier A] f={frequency:.0f} Hz, K={len(paths)}")
    I_sum = np.sum(I_tier_a)
    I_mag = np.abs(I_tier_a)
    ph = np.angle(I_tier_a, deg=True) if np.iscomplexobj(I_tier_a) else np.zeros_like(I_mag)
    print(f"  sum(I) = {I_sum}  "
          f"|I_k| in [{I_mag.min():.3f}, {I_mag.max():.3f}]  "
          f"arg in [{ph.min():.1f}, {ph.max():.1f}] deg")
    t0 = time.perf_counter()
    phi = compute_phi_inc_from_filaments(
        wp_nodes, paths, I_tier_a, n_quad=20)
    print(f"  phi Re range [{phi.real.min():.4f}, {phi.real.max():.4f}], "
          f"|Im| max {np.abs(phi.imag).max():.4e} "
          f"({time.perf_counter()-t0:.1f} s)")
    return phi


# ------------------------------------------------------------------
# (c) PEEC: full MNA on filament bundle -> per-filament I_peec -> phi_inc
# ------------------------------------------------------------------

def phi_inc_from_peec(wp_nodes, paths, dw, dh, current, frequency):
    print(f"\n[PEEC bundle] f={frequency:.0f} Hz, K={len(paths)}")
    t0 = time.perf_counter()
    solver, seg_of_fil, port_p, port_m = build_bundle_solver(
        paths, dw, dh, CU_SIGMA)
    print(f"  PEEC n_loop = {solver.n_loop}, "
          f"n_nodes = {solver.n_nodes}  "
          f"(build {time.perf_counter()-t0:.1f} s)")

    t0 = time.perf_counter()
    I_branch = solver.compute_branch_currents(frequency,
                                               [complex(current)])
    print(f"  compute_branch_currents: {time.perf_counter()-t0:.2f} s")

    N = len(paths)
    I_peec = np.zeros(N, dtype=complex)
    for k in range(N):
        I_peec[k] = np.mean(I_branch[seg_of_fil[k]])

    I_sum = np.sum(I_peec)
    I_mag = np.abs(I_peec)
    ph = np.angle(I_peec, deg=True)
    print(f"  sum(I) = {I_sum}  "
          f"|I_k| in [{I_mag.min():.3f}, {I_mag.max():.3f}]  "
          f"arg in [{ph.min():.1f}, {ph.max():.1f}] deg")

    t0 = time.perf_counter()
    phi = compute_phi_inc_from_filaments(
        wp_nodes, paths, I_peec, n_quad=20)
    print(f"  phi Re range [{phi.real.min():.4f}, {phi.real.max():.4f}], "
          f"|Im| max {np.abs(phi.imag).max():.4e} "
          f"({time.perf_counter()-t0:.1f} s)")
    return phi, I_peec


# ------------------------------------------------------------------
# SIBC run
# ------------------------------------------------------------------

def run_sibc(wp_mesh, phi_inc, frequency, sigma, mu_r):
    omega = 2.0 * np.pi * frequency
    mu = MU_0 * mu_r
    rho = 1.0 / sigma
    delta = np.sqrt(2.0 * rho / (omega * mu))
    Z_s = complex(1, 1) * rho / delta

    solver = ScalarBIESIBCSolver(wp_mesh, order=1)
    result = solver.solve(phi_inc, Z_s=Z_s, omega=omega)

    area = result['area']
    H_t_rms = result['H_t_rms']
    P_density = result['P_density']
    P_total = P_density * area
    Q_density = 0.5 * Z_s.imag * H_t_rms ** 2
    Q_total = Q_density * area
    return {
        'Z_s': Z_s,
        'delta_mm': delta * 1e3,
        'H_t_rms': H_t_rms,
        'P_total': P_total,
        'Q_total': Q_total,
    }


# ------------------------------------------------------------------
# Main
# ------------------------------------------------------------------

def main():
    coil_mesh, wp_mesh = load_meshes()

    wp_nodes = np.array(
        [[wp_mesh.vertices[i].point[j] for j in range(3)]
         for i in range(wp_mesh.nv)])

    # Build filament discretization (paths shared between Tier A and PEEC)
    _, paths, I_tier_a, dw, dh = build_coil_and_filaments(
        COIL_CURRENT, FREQUENCY)

    # Three phi_inc
    phi_bem = phi_inc_from_bem(coil_mesh, wp_nodes, COIL_CURRENT)
    phi_ta = phi_inc_from_tier_a(wp_nodes, paths, I_tier_a, FREQUENCY)
    phi_pe, I_peec = phi_inc_from_peec(
        wp_nodes, paths, dw, dh, COIL_CURRENT, FREQUENCY)

    # Per-filament Tier A vs PEEC difference (context)
    I_ta_c = I_tier_a.astype(complex)
    diff_rel = np.max(np.abs(I_peec - I_ta_c)) / np.max(np.abs(I_ta_c))
    print(f"\n[Per-filament] max|I_PEEC - I_TierA| / max|I_TierA| = "
          f"{diff_rel*100:.2f}%")

    # Run SIBC on each phi_inc
    print("\n[SIBC] BEM truth (DC coil)")
    res_bem = run_sibc(wp_mesh, phi_bem, FREQUENCY, WP_SIGMA, WP_MU_R)
    for k, v in res_bem.items():
        print(f"  {k}: {v}")

    print("\n[SIBC] Tier A")
    res_ta = run_sibc(wp_mesh, phi_ta, FREQUENCY, WP_SIGMA, WP_MU_R)
    for k, v in res_ta.items():
        print(f"  {k}: {v}")

    print("\n[SIBC] PEEC")
    res_pe = run_sibc(wp_mesh, phi_pe, FREQUENCY, WP_SIGMA, WP_MU_R)
    for k, v in res_pe.items():
        print(f"  {k}: {v}")

    # Summary table
    print("\n" + "=" * 78)
    print("Summary  (workpiece metrics, I_total = 1 A, f = 7 kHz)")
    print("=" * 78)
    fmt = "  {name:>10}  {bem:>14}  {ta:>14}  {pe:>14}  {ta_rel:>+7.3f}%  {pe_rel:>+7.3f}%"
    print(f"  {'metric':>10}  {'BEM truth':>14}  {'Tier A':>14}  "
          f"{'PEEC':>14}  {'TA-BEM':>8}  {'PE-BEM':>8}")
    print("  " + "-" * 76)
    for name, key in [('H_t_rms', 'H_t_rms'),
                       ('P_total', 'P_total'),
                       ('Q_total', 'Q_total')]:
        b = res_bem[key]
        t = res_ta[key]
        p = res_pe[key]
        print(fmt.format(name=name,
                         bem=f"{b:.4e}", ta=f"{t:.4e}", pe=f"{p:.4e}",
                         ta_rel=(t - b) / b * 100.0,
                         pe_rel=(p - b) / b * 100.0))


if __name__ == '__main__':
    main()
