"""compare_tier_a_vs_ngbem.py

Phase 1 validation of the Tier A analytical filament current model
against the frequency-aware BEM+SIBC pipeline already used in the
Radia-NGSolve Induction Heating panel (ih_bem_sample).

Validation target
-----------------
The scalar BIE + SIBC solver in src/radia/bem_sibc_solver.py accepts an
incident scalar magnetic potential phi_inc on the workpiece surface and
returns the workpiece heating response (H_t_rms, P_density, Q_density)
at the specified frequency. Two ways of feeding phi_inc are compared:

  (a) BEM truth: solve EFIE (LaplaceSL saddle-point with source/sink
      ports) on the coil surface mesh to obtain the DC surface current J,
      then integrate J with a Biot-Savart kernel to phi_inc on the
      workpiece surface nodes. This is exactly the --coil-vol path of
      calc_heating_bem.py.

  (b) Tier A: represent the coil as nw x nh filaments with analytical
      complex per-filament currents (path-length bias + skin-effect
      factor), and integrate the filament bundle to phi_inc on the
      same workpiece surface nodes via compute_phi_inc_from_filaments.

Both phi_inc are fed into the SAME ScalarBIESIBCSolver with the same
Z_s(omega) and compared on P_total, H_t_rms, and pointwise phi_inc.

Geometry (matched to src/radia/panels/samples/ih_bem_sample.jou):
  - Coil: gapped torus, R_major = 30 mm, R_minor = 3 mm, gap = 5 deg.
    (Circular cross-section. Tier A uses an area-equivalent square
    cross-section of side sqrt(pi) * R_minor ~ 5.32 mm.)
  - Workpiece: steel cylinder, R = 25 mm, H = 25 mm.
  - Frequency: 7 kHz.

Hypothesis
----------
Without strong proximity coupling between coil and workpiece, the Tier A
filament bundle should reproduce the BEM-truth phi_inc on the workpiece
to within a few percent, and the downstream heating metrics should
agree similarly. This is the "no-workpiece-proximity" regime.

When the workpiece is close to the coil and induces substantial image
currents in the coil (proximity effect), Tier A is expected to deviate
because it omits off-diagonal filament-to-workpiece coupling. The
present sample has a ~2 mm coil-to-workpiece gap, so some deviation is
already expected. Larger deviations in future workpiece-close cases
motivate moving from Tier A to full PEEC.

Usage
-----
    python compare_tier_a_vs_ngbem.py
"""

import os
import sys
import time
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.abspath(os.path.join(HERE, '..', '..', '..', 'src'))
if SRC not in sys.path:
    sys.path.insert(0, SRC)
PANELS = os.path.join(SRC, 'radia', 'panels')
if PANELS not in sys.path:
    sys.path.insert(0, PANELS)

from ngsolve import Mesh, Integrate, CF, BND
from ngsolve import TaskManager

from radia.coil_builder import CoilBuilder
from radia.bem_sibc_solver import (ScalarBIESIBCSolver,
                                    compute_phi_inc_from_surface_J,
                                    compute_phi_inc_from_filaments)
from radia.bem_inductance import compute_inductance_source_sink

# Mesh extractors from the panel code (reuse, don't reinvent)
from calc_heating_bem import (_extract_surface_mesh_filtered,
                               _extract_surface_mesh)

MU_0 = 4.0 * np.pi * 1e-7
MM = 1e-3


# ------------------------------------------------------------
# Parameters (matched to ih_bem_sample.jou)
# ------------------------------------------------------------
VOL_FILE = os.path.join(SRC, 'radia', 'panels', 'samples',
                        'ih_bem_sample.vol')

R_MAJOR = 30 * MM
R_MINOR = 3 * MM
GAP_DEG = 5.0
ARC_DEG = 360.0 - GAP_DEG

COIL_CURRENT = 1.0                # [A]
FREQUENCY = 7e3                   # [Hz]

# Workpiece (steel) electrical properties
WP_SIGMA = 2e6                    # [S/m]
WP_MU_R = 100.0                   # relative permeability (linear)

# Coil material (Cu) for Tier A skin-effect factor
CU_SIGMA = 5.8e7
CU_MU_R = 1.0

# Tier A filament grid. Cross-section is 5.32 x 5.32 mm (area-equivalent
# to the 3 mm radius circular BEM cross-section).
TIER_A_SIDE = np.sqrt(np.pi) * R_MINOR
NW = 5
NH = 5


# ------------------------------------------------------------
# Extract coil and workpiece surface meshes from ih_bem_sample.vol
# ------------------------------------------------------------

def load_meshes():
    print(f"Loading {VOL_FILE}")
    mesh_full = Mesh(VOL_FILE)
    print(f"  materials:  {mesh_full.GetMaterials()}")
    print(f"  boundaries: {sorted(set(mesh_full.GetBoundaries()))}")

    coil_mesh = _extract_surface_mesh_filtered(mesh_full, keep_label='coil')
    wp_mesh = _extract_surface_mesh_filtered(mesh_full, keep_label='workpiece')

    print(f"  coil surface: nse={coil_mesh.GetNE(BND)}, nv={coil_mesh.nv}")
    print(f"  wp   surface: nse={wp_mesh.GetNE(BND)},   nv={wp_mesh.nv}")
    return coil_mesh, wp_mesh


# ------------------------------------------------------------
# BEM truth: EFIE on coil -> surface J -> phi_inc on workpiece nodes
# ------------------------------------------------------------

def phi_inc_from_bem_coil(coil_mesh, wp_nodes, current):
    print("\n[BEM truth] Solving EFIE on coil (source/sink)...")
    t0 = time.perf_counter()
    sol = compute_inductance_source_sink(
        coil_mesh, source_label='source', sink_label='sink', fes_order=0)
    if 'error' in sol:
        raise RuntimeError(f"Coil EFIE failed: {sol['error']}")
    L_coil = sol['L']
    gf_J = sol['gf_J']
    print(f"  L_coil = {L_coil*1e9:.2f} nH, n_J = {sol['n_J']}, "
          f"residual = {sol['residual']:.2e}, "
          f"t_solve = {sol['t_solve']:.1f} s")

    # Extract per-element J vector and centroid for Biot-Savart integration
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
        Jvecs = np.asarray(Jvecs) * current   # scale to physical current

        print(f"  {len(centroids)} source panels "
              f"(total area = {areas.sum()*1e6:.1f} mm^2)")

        print("  integrating phi_inc from surface J...")
        t0 = time.perf_counter()
        phi = compute_phi_inc_from_surface_J(
            wp_nodes, centroids, areas, Jvecs, n_quad=20)
        print(f"  phi_inc: range [{phi.min():.4f}, {phi.max():.4f}] "
              f"({time.perf_counter()-t0:.1f} s)")
        return phi, L_coil


# ------------------------------------------------------------
# Tier A: CoilBuilder -> to_filaments -> phi_inc
# ------------------------------------------------------------

def build_tier_a_coil(current):
    """Circular arc 355 deg, area-equivalent square cross-section.

    The CoilBuilder arc is laid out in the x-y plane with the coil
    starting at x=R_MAJOR, y=0 moving in +y initially, matching the
    ih_bem_sample torus orientation (axis along z, starting at +x).
    """
    coil = (CoilBuilder(current=current)
            .set_start([R_MAJOR, 0.0, 0.0])
            .set_cross_section(width=TIER_A_SIDE, height=TIER_A_SIDE)
            .add_arc(radius=R_MAJOR, arc_angle=ARC_DEG))
    return coil


def phi_inc_from_tier_a(wp_nodes, current, frequency, sigma, mu_r,
                         nw, nh, n_arc):
    print(f"\n[Tier A] nw={nw} nh={nh} n_arc={n_arc} "
          f"f={frequency:.0f} Hz")
    coil = build_tier_a_coil(current)
    fil_paths, fil_I = coil.to_filaments(
        nw=nw, nh=nh, frequency=frequency, sigma=sigma,
        mu_r=mu_r, n_arc=n_arc)
    I_sum = np.sum(fil_I)
    I_max = np.max(np.abs(fil_I))
    I_min = np.min(np.abs(fil_I))
    ph_max = np.max(np.angle(fil_I, deg=True))
    ph_min = np.min(np.angle(fil_I, deg=True))
    print(f"  filaments: K = {len(fil_I)}, sum(I) = {I_sum} "
          f"(should be {current})")
    print(f"  |I_k|: [{I_min:.3f}, {I_max:.3f}] A   "
          f"arg(I_k): [{ph_min:.1f}, {ph_max:.1f}] deg")

    # Skin depth at the coil conductor
    delta = np.sqrt(2.0 / (2.0 * np.pi * frequency * MU_0 * mu_r * sigma))
    cell = min(TIER_A_SIDE / nw, TIER_A_SIDE / nh)
    print(f"  delta = {delta*1e3:.3f} mm, cell = {cell*1e3:.3f} mm, "
          f"cell/delta = {cell/delta:.2f}")

    print("  integrating phi_inc from filaments...")
    t0 = time.perf_counter()
    phi = compute_phi_inc_from_filaments(
        wp_nodes, fil_paths, fil_I, n_quad=20)
    print(f"  phi_inc: |Re| range [{np.abs(phi.real).min():.4f}, "
          f"{np.abs(phi.real).max():.4f}], |Im| max "
          f"{np.abs(phi.imag).max():.4e}  "
          f"({time.perf_counter()-t0:.1f} s)")
    return phi


# ------------------------------------------------------------
# Run BEM SIBC with given phi_inc and report heating metrics
# ------------------------------------------------------------

def run_sibc(wp_mesh, phi_inc, frequency, sigma, mu_r):
    omega = 2.0 * np.pi * frequency
    mu = MU_0 * mu_r
    rho = 1.0 / sigma
    delta = np.sqrt(2.0 * rho / (omega * mu))
    Z_s = complex(1, 1) * rho / delta   # linear SIBC

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
        'area_m2': area,
        'H_t_rms': H_t_rms,
        'P_total': P_total,
        'Q_total': Q_total,
    }


# ------------------------------------------------------------
# Main
# ------------------------------------------------------------

def main():
    coil_mesh, wp_mesh = load_meshes()

    # Workpiece node coordinates (used as phi_inc observation points)
    wp_nodes = np.array(
        [[wp_mesh.vertices[i].point[j] for j in range(3)]
         for i in range(wp_mesh.nv)])

    phi_bem, L_coil = phi_inc_from_bem_coil(
        coil_mesh, wp_nodes, COIL_CURRENT)
    phi_tierA = phi_inc_from_tier_a(
        wp_nodes, COIL_CURRENT, FREQUENCY, CU_SIGMA, CU_MU_R,
        NW, NH, n_arc=40)

    # phi_inc pointwise comparison. The scalar potential is defined up
    # to an additive constant AND a sign (source/sink direction
    # convention can flip between CoilBuilder and the BEM jou). The
    # physical BEM output depends only on grad(phi_inc), so we detect
    # the sign by fitting ta = s * bem + c and report in that gauge.
    phi_bem_c = phi_bem.astype(complex)
    # Fit (s, c) to minimize ||Re(Tier A) - (s * bem + c)||^2 with s in
    # {-1, +1}.
    ta_re = phi_tierA.real
    best = None
    for s in (+1, -1):
        c = np.mean(ta_re - s * phi_bem)
        resid = ta_re - (s * phi_bem + c)
        err = np.max(np.abs(resid))
        if best is None or err < best[2]:
            best = (s, c, err, resid)
    s, c, max_abs_err, resid = best
    rel = np.abs(resid) / (np.abs(phi_bem) + 1e-30)
    print("\n[phi_inc comparison]")
    print(f"  BEM truth (DC, real): range [{phi_bem.min():.4f}, "
          f"{phi_bem.max():.4f}]")
    print(f"  Tier A (complex):     Re range "
          f"[{phi_tierA.real.min():.4f}, {phi_tierA.real.max():.4f}], "
          f"|Im| max {np.abs(phi_tierA.imag).max():.4e}")
    print(f"  Best-fit affine map:  Re(Tier A) = {s:+d} * BEM + "
          f"{c:+.4f}")
    print(f"    residual: max|diff| = {max_abs_err:.4e}, "
          f"max rel = {np.max(rel)*100:.2f}%, "
          f"mean rel = {np.mean(rel)*100:.2f}%")

    print("\n[BEM SIBC with phi_inc from BEM truth]")
    res_bem = run_sibc(wp_mesh, phi_bem, FREQUENCY, WP_SIGMA, WP_MU_R)
    for k, v in res_bem.items():
        print(f"  {k}: {v}")

    print("\n[BEM SIBC with phi_inc from Tier A]")
    res_ta = run_sibc(wp_mesh, phi_tierA, FREQUENCY, WP_SIGMA, WP_MU_R)
    for k, v in res_ta.items():
        print(f"  {k}: {v}")

    print("\n[Summary: Tier A vs BEM truth]")
    fmt = "  {name:>12}  {bem:>14}  {ta:>14}  {rel:>+8.3f}%"
    for name, keys in [('H_t_rms', 'H_t_rms'),
                        ('P_total', 'P_total'),
                        ('Q_total', 'Q_total')]:
        b = res_bem[keys]
        t = res_ta[keys]
        rel = (t - b) / b * 100.0
        print(fmt.format(name=name, bem=f"{b:.4e}", ta=f"{t:.4e}", rel=rel))


if __name__ == '__main__':
    main()
