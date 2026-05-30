"""
BEM coil + ESIM workpiece: impedance extraction with skin effect.

Physical setup:
  - Torus coil (gapped, OCC) carrying unit current
  - Cylindrical workpiece at torus center

Pipeline:
  1. BEM (source/sink saddle point EFIE) -> coil surface current J
  2. Biot-Savart from J -> H at workpiece surface panels
  3. ESIM cell problem -> surface impedance Z_s(H), power loss per panel
  4. Total: Z_workpiece = sum(Z_s * dA) integrated over workpiece surface

Usage:
    python impedance_esim.py
    python impedance_esim.py --material steel    # nonlinear BH curve
    python impedance_esim.py --material copper    # linear conductor
    python impedance_esim.py --freq 100000        # 100 kHz
"""

import math
import os
import sys
import time
import numpy as np

# Import from src/radia/
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../../src/radia'))
from bem_inductance import compute_inductance_source_sink, MU_0
from esim_cell_problem import ESIMFiniteSlabSolver

from ngsolve import *
from ngsolve import TaskManager


# ============================================================
# Material database
# ============================================================
# Steel BH curve (typical carbon steel, sigma ~ 2e6 S/m)
STEEL_BH = [
    [0.0, 0.0], [50.0, 0.10], [100.0, 0.25], [200.0, 0.55],
    [500.0, 0.95], [1000.0, 1.20], [2000.0, 1.40], [5000.0, 1.55],
    [10000.0, 1.65], [20000.0, 1.75], [50000.0, 1.90], [100000.0, 2.00],
]
STEEL_SIGMA = 2.0e6  # S/m

# Copper (linear, mu_r=1, sigma=5.8e7 S/m)
COPPER_SIGMA = 5.8e7
COPPER_MU_R = 1.0

# Aluminum (linear, mu_r=1, sigma=3.5e7 S/m)
ALUMINUM_SIGMA = 3.5e7
ALUMINUM_MU_R = 1.0


# ============================================================
# Geometry: gapped torus (coil) + cylinder (workpiece)
# ============================================================
def make_gapped_torus_mesh(R, a, gap_deg=5, maxh=None):
    """Create OCC gapped torus surface mesh with source/sink labels.

    Args:
        R: Major radius [m]
        a: Minor radius [m]
        gap_deg: Gap angle [degrees]
        maxh: Max mesh size (default: a/2)

    Returns:
        NGSolve Mesh with boundary labels "conductor", "source", "sink"
    """
    from netgen.occ import (WorkPlane, Axes, Axis, Pnt, Dir,
                             OCCGeometry, Glue)
    from netgen.meshing import MeshingParameters

    if maxh is None:
        maxh = a / 2

    wp = WorkPlane(Axes(p=Pnt(R, 0, 0), n=Dir(0, 1, 0), h=Dir(0, 0, 1)))
    circle = wp.Circle(a).Face()
    torus = circle.Revolve(Axis(Pnt(0, 0, 0), Dir(0, 0, 1)), 360 - gap_deg)

    for f in torus.faces:
        f.name = "conductor"
    faces_sorted = sorted(torus.faces, key=lambda f: f.mass)
    faces_sorted[0].name = "source"
    faces_sorted[1].name = "sink"

    torus_surf = Glue(torus.faces)
    geo = OCCGeometry(torus_surf)
    ngmesh = geo.GenerateMesh(
        mp=MeshingParameters(maxh=maxh, curvaturesafety=2.0,
                             segmentsperedge=2))
    mesh = Mesh(ngmesh)
    with TaskManager():
        mesh.Curve(2)
        return mesh


def make_workpiece_panels(center, radius, height, n_phi=16, n_z=8):
    """Create cylindrical workpiece surface panels.

    Discretizes a cylinder into rectangular panels on the lateral surface
    and triangular/rectangular panels on the top and bottom caps.

    Args:
        center: [x, y, z] center of cylinder [m]
        radius: Cylinder radius [m]
        height: Cylinder height [m]
        n_phi: Number of circumferential divisions
        n_z: Number of axial divisions on lateral surface

    Returns:
        panels: list of dicts, each with:
            'center': [x, y, z] panel center
            'normal': [nx, ny, nz] outward normal
            'area': panel area [m^2]
            'type': 'lateral' or 'cap'
    """
    cx, cy, cz = center
    z_bot = cz - height / 2
    z_top = cz + height / 2
    dphi = 2 * math.pi / n_phi
    dz = height / n_z

    panels = []

    # Lateral surface panels
    for i in range(n_phi):
        phi_mid = (i + 0.5) * dphi
        nx = math.cos(phi_mid)
        ny = math.sin(phi_mid)
        for j in range(n_z):
            z_mid = z_bot + (j + 0.5) * dz
            px = cx + radius * nx
            py = cy + radius * ny
            pz = z_mid
            area = radius * dphi * dz
            panels.append({
                'center': np.array([px, py, pz]),
                'normal': np.array([nx, ny, 0.0]),
                'area': area,
                'type': 'lateral',
            })

    # Top and bottom cap panels (annular sectors)
    n_r_cap = max(2, n_z // 2)
    dr = radius / n_r_cap
    for cap_sign, cap_z, cap_nz in [(+1, z_top, +1.0), (-1, z_bot, -1.0)]:
        for i in range(n_phi):
            phi_mid = (i + 0.5) * dphi
            cos_p = math.cos(phi_mid)
            sin_p = math.sin(phi_mid)
            for k in range(n_r_cap):
                r_mid = (k + 0.5) * dr
                px = cx + r_mid * cos_p
                py = cy + r_mid * sin_p
                pz = cap_z
                # Annular sector area: (r_outer^2 - r_inner^2)/2 * dphi
                r_inner = k * dr
                r_outer = (k + 1) * dr
                area = 0.5 * (r_outer**2 - r_inner**2) * dphi
                panels.append({
                    'center': np.array([px, py, pz]),
                    'normal': np.array([0.0, 0.0, cap_nz]),
                    'area': area,
                    'type': 'cap',
                })

    return panels


# ============================================================
# Biot-Savart: H at workpiece panels from BEM coil current
# ============================================================
def compute_H_at_panels(mesh_coil, gf_J, panels):
    """Compute H-field at workpiece panel centers via Biot-Savart.

    H = (1/4pi) * sum_elem { J x r / |r|^3 * dA }

    Args:
        mesh_coil: NGSolve surface mesh of coil
        gf_J: GridFunction(HDivSurface) with solved current
        panels: list of panel dicts from make_workpiece_panels

    Returns:
        H_panels: (n_panels, 3) array of H vectors [A/m]
    """
    INV_4PI = 1.0 / (4.0 * np.pi)

    # Extract per-element J and geometry from coil mesh
    elem_A = Integrate(CF(1), mesh_coil, VOL_or_BND=BND, element_wise=True)
    elem_Jx = Integrate(gf_J[0], mesh_coil, VOL_or_BND=BND, element_wise=True)
    elem_Jy = Integrate(gf_J[1], mesh_coil, VOL_or_BND=BND, element_wise=True)
    elem_Jz = Integrate(gf_J[2], mesh_coil, VOL_or_BND=BND, element_wise=True)

    centroids, areas, J_vecs = [], [], []
    for el in mesh_coil.Elements(BND):
        area = abs(elem_A[el.nr])
        if area < 1e-30:
            continue
        jvec = np.array([elem_Jx[el.nr], elem_Jy[el.nr], elem_Jz[el.nr]]) / area
        verts = [mesh_coil.vertices[v.nr].point for v in el.vertices]
        c = np.mean([(v[0], v[1], v[2]) for v in verts], axis=0)
        centroids.append(c)
        areas.append(area)
        J_vecs.append(jvec)

    centroids = np.array(centroids)  # (ne, 3)
    areas = np.array(areas)          # (ne,)
    J_vecs = np.array(J_vecs)        # (ne, 3)

    # Observation points
    obs_pts = np.array([p['center'] for p in panels])  # (np, 3)

    # Vectorized Biot-Savart: H = (1/4pi) * sum { J x (r_obs - r_src) / |r|^3 * dA }
    dx = obs_pts[:, None, :] - centroids[None, :, :]   # (np, ne, 3)
    r = np.sqrt(np.maximum(np.sum(dx**2, axis=2), 1e-30))  # (np, ne)
    r3_inv = areas[None, :] / (r ** 3)                 # (np, ne)
    cross = np.cross(J_vecs[None, :, :], dx)            # (np, ne, 3)
    H_panels = INV_4PI * np.sum(cross * r3_inv[:, :, None], axis=1)  # (np, 3)

    return H_panels


def compute_H_tangential(H_panels, panels):
    """Extract tangential H magnitude at each panel.

    H_t = H - (H . n) * n

    Args:
        H_panels: (n_panels, 3) H vectors
        panels: list of panel dicts

    Returns:
        H_t_mag: (n_panels,) tangential H magnitude [A/m]
        H_t_vec: (n_panels, 3) tangential H vectors
    """
    normals = np.array([p['normal'] for p in panels])
    H_n = np.sum(H_panels * normals, axis=1, keepdims=True)
    H_t_vec = H_panels - H_n * normals
    H_t_mag = np.linalg.norm(H_t_vec, axis=1)
    return H_t_mag, H_t_vec


# ============================================================
# ESIM: surface impedance and power loss
# ============================================================
def compute_esim_impedance(panels, H_t_mag, half_thickness, sigma, frequency,
                           bh_curve=None, mu_r=1.0, n_nodes=100):
    """Compute ESIM surface impedance and power loss for each panel.

    Args:
        panels: list of panel dicts
        H_t_mag: (n_panels,) tangential H magnitude [A/m]
        half_thickness: Workpiece half-thickness [m] (for slab model)
        sigma: Conductivity [S/m]
        frequency: Operating frequency [Hz]
        bh_curve: BH curve data (for nonlinear materials)
        mu_r: Relative permeability (for linear materials)
        n_nodes: ESIM mesh nodes

    Returns:
        result: dict with per-panel and total results
    """
    solver = ESIMFiniteSlabSolver(
        half_thickness=half_thickness,
        bh_curve=bh_curve,
        sigma=sigma,
        frequency=frequency,
        mu_r=mu_r if bh_curve is None else None,
        n_nodes=n_nodes,
    )

    n_panels = len(panels)
    Z_panels = np.zeros(n_panels, dtype=complex)
    P_panels = np.zeros(n_panels)   # Active power per panel [W]
    Q_panels = np.zeros(n_panels)   # Reactive power per panel [W]
    R_ratio_panels = np.zeros(n_panels)
    delta_panels = np.zeros(n_panels)

    for i in range(n_panels):
        H0 = max(float(H_t_mag[i]), 1e-3)  # Avoid zero
        sol = solver.solve(H0)
        Z_panels[i] = sol['Z']
        # Power per panel: P' [W/m^2] * area [m^2]
        P_panels[i] = sol['P_prime'] * panels[i]['area']
        Q_panels[i] = sol['Q_prime'] * panels[i]['area']
        R_ratio_panels[i] = sol['R_ratio']
        delta_panels[i] = sol['delta']

    P_total = np.sum(P_panels)
    Q_total = np.sum(Q_panels)
    S_total = P_total + 1j * Q_total

    # Effective impedance seen from coil (for unit current):
    # Z_wp = sum(Z_s * |H_t|^2 * dA) / I^2
    # For unit current I=1, this simplifies to the sum
    areas = np.array([p['area'] for p in panels])
    Z_effective = np.sum(Z_panels * H_t_mag**2 * areas)

    omega = 2 * np.pi * frequency
    L_effective = Z_effective.imag / omega if omega > 0 else 0.0
    R_effective = Z_effective.real

    return {
        'Z_panels': Z_panels,
        'P_panels': P_panels,
        'Q_panels': Q_panels,
        'R_ratio_panels': R_ratio_panels,
        'delta_panels': delta_panels,
        'P_total': P_total,
        'Q_total': Q_total,
        'Z_effective': Z_effective,
        'R_effective': R_effective,
        'L_effective': L_effective,
    }


# ============================================================
# Main
# ============================================================
def run(material='steel', frequency=50000.0, R=0.05, a=0.005,
        wp_radius=0.01, wp_height=0.02, n_phi=24, n_z=12,
        verbose=True):
    """Run coupled BEM + ESIM analysis.

    Args:
        material: 'steel', 'copper', or 'aluminum'
        frequency: Operating frequency [Hz]
        R: Torus major radius [m]
        a: Torus minor radius [m]
        wp_radius: Workpiece cylinder radius [m]
        wp_height: Workpiece cylinder height [m]
        n_phi: Workpiece circumferential panels
        n_z: Workpiece axial panels
        verbose: Print progress

    Returns:
        dict with all results
    """
    # Material selection
    if material == 'steel':
        bh_curve = STEEL_BH
        sigma = STEEL_SIGMA
        mu_r = None
        mat_label = f"Steel (sigma={sigma:.0e} S/m, nonlinear BH)"
    elif material == 'copper':
        bh_curve = None
        sigma = COPPER_SIGMA
        mu_r = COPPER_MU_R
        mat_label = f"Copper (sigma={sigma:.1e} S/m, mu_r={mu_r})"
    elif material == 'aluminum':
        bh_curve = None
        sigma = ALUMINUM_SIGMA
        mu_r = ALUMINUM_MU_R
        mat_label = f"Aluminum (sigma={sigma:.1e} S/m, mu_r={mu_r})"
    else:
        raise ValueError(f"Unknown material: {material}")

    omega = 2 * np.pi * frequency
    rho = 1.0 / sigma
    mu_eff = MU_0 * (mu_r if mu_r else 1.0)
    delta = math.sqrt(2 * rho / (omega * mu_eff))

    if verbose:
        print("=" * 65)
        print("BEM Coil + ESIM Workpiece: Coupled Impedance Extraction")
        print("=" * 65)
        print(f"Coil:      Torus R={R*1e3:.0f} mm, a={a*1e3:.1f} mm")
        print(f"Workpiece: Cylinder r={wp_radius*1e3:.0f} mm, "
              f"h={wp_height*1e3:.0f} mm")
        print(f"Material:  {mat_label}")
        print(f"Frequency: {frequency/1e3:.0f} kHz")
        print(f"Skin depth (linear): {delta*1e3:.3f} mm")
        print(f"Panels:    {n_phi} x {n_z} (lateral) + caps")
        print()

    # --- Step 1: BEM on torus coil ---
    if verbose:
        print("[1/4] Creating coil mesh (OCC gapped torus)...")
    t0 = time.perf_counter()
    mesh_coil = make_gapped_torus_mesh(R, a, gap_deg=5, maxh=a)
    nse = mesh_coil.GetNE(BND)
    if verbose:
        print(f"       {nse} surface elements, "
              f"{mesh_coil.nv} vertices ({time.perf_counter()-t0:.1f}s)")

    if verbose:
        print("[2/4] BEM solve (saddle point EFIE)...")
    t0 = time.perf_counter()
    sol = compute_inductance_source_sink(mesh_coil)
    if 'error' in sol:
        print(f"  ERROR: {sol['error']}")
        return sol
    L_coil = sol['L']
    gf_J = sol['gf_J']
    t_bem = time.perf_counter() - t0
    L_neumann = MU_0 * R * (math.log(8 * R / a) - 2)
    if verbose:
        print(f"       L_coil = {L_coil*1e9:.2f} nH "
              f"(Neumann: {L_neumann*1e9:.2f} nH, "
              f"err={((L_coil-L_neumann)/L_neumann*100):+.1f}%)")
        print(f"       DOFs: {sol['n_J']} J + {sol['n_f']} constraints "
              f"({t_bem:.1f}s)")

    # --- Step 2: Workpiece panels + Biot-Savart H ---
    if verbose:
        print("[3/4] Biot-Savart: H at workpiece panels...")
    t0 = time.perf_counter()
    panels = make_workpiece_panels([0, 0, 0], wp_radius, wp_height,
                                   n_phi=n_phi, n_z=n_z)
    H_panels = compute_H_at_panels(mesh_coil, gf_J, panels)
    H_t_mag, H_t_vec = compute_H_tangential(H_panels, panels)
    t_biot = time.perf_counter() - t0
    if verbose:
        print(f"       {len(panels)} panels, "
              f"|H_t| range: [{H_t_mag.min():.1f}, {H_t_mag.max():.1f}] A/m "
              f"({t_biot:.2f}s)")

    # --- Step 3: ESIM cell problem ---
    if verbose:
        print("[4/4] ESIM surface impedance...")
    t0 = time.perf_counter()
    # Use cylinder radius as half-thickness for slab model
    esim = compute_esim_impedance(
        panels, H_t_mag,
        half_thickness=wp_radius,
        sigma=sigma,
        frequency=frequency,
        bh_curve=bh_curve,
        mu_r=mu_r if mu_r else 1.0,
    )
    t_esim = time.perf_counter() - t0
    if verbose:
        print(f"       Skin depth range: "
              f"[{esim['delta_panels'].min()*1e3:.3f}, "
              f"{esim['delta_panels'].max()*1e3:.3f}] mm")
        print(f"       R_ac/R_dc range: "
              f"[{esim['R_ratio_panels'].min():.2f}, "
              f"{esim['R_ratio_panels'].max():.2f}]")
        print(f"       ({t_esim:.2f}s)")

    # --- Summary ---
    if verbose:
        print()
        print("-" * 65)
        print("Results (unit current I = 1 A)")
        print("-" * 65)
        print(f"Coil self-inductance:   L_coil = {L_coil*1e9:.4f} nH")
        print(f"Coil reactance:         X_coil = {omega*L_coil*1e6:.4f} uOhm")
        print()
        print(f"Workpiece power loss:   P_wp = {esim['P_total']:.6e} W")
        print(f"Workpiece reactive:     Q_wp = {esim['Q_total']:.6e} var")
        print(f"Workpiece R_effective:  R_wp = {esim['R_effective']:.6e} Ohm")
        print(f"Workpiece L_effective:  L_wp = {esim['L_effective']*1e9:.6f} nH")
        print()
        Z_total = (esim['R_effective'] + 1j * omega * L_coil
                   + 1j * esim['Z_effective'].imag)
        print(f"Total impedance:        Z = {Z_total.real:.6e} "
              f"+ j{Z_total.imag:.6e} Ohm")
        print(f"  (at f = {frequency/1e3:.0f} kHz)")
        print("-" * 65)

    return {
        'frequency': frequency,
        'material': material,
        'L_coil': L_coil,
        'L_neumann': L_neumann,
        'n_coil_elements': nse,
        't_bem': t_bem,
        't_biot': t_biot,
        't_esim': t_esim,
        'n_panels': len(panels),
        'H_t_mag': H_t_mag,
        'panels': panels,
        **esim,
    }


def frequency_sweep(material='steel', freqs=None, **kwargs):
    """Run ESIM analysis at multiple frequencies.

    Args:
        material: Material name
        freqs: Array of frequencies [Hz] (default: 1kHz to 1MHz)
        **kwargs: Passed to run()

    Returns:
        results: list of result dicts
    """
    if freqs is None:
        freqs = np.logspace(3, 6, 13)  # 1 kHz to 1 MHz

    print(f"\nFrequency sweep: {len(freqs)} points "
          f"({freqs[0]/1e3:.0f} kHz - {freqs[-1]/1e3:.0f} kHz)")
    print(f"{'f [kHz]':>10s}  {'delta [mm]':>10s}  {'P_wp [W]':>12s}  "
          f"{'R_wp [Ohm]':>12s}  {'L_wp [nH]':>10s}")
    print("-" * 65)

    results = []
    for f in freqs:
        r = run(material=material, frequency=f, verbose=False, **kwargs)
        if 'error' in r:
            continue
        rho = 1.0 / (STEEL_SIGMA if material == 'steel' else
                      COPPER_SIGMA if material == 'copper' else ALUMINUM_SIGMA)
        mu_val = MU_0
        delta = math.sqrt(2 * rho / (2 * math.pi * f * mu_val))
        print(f"{f/1e3:10.1f}  {delta*1e3:10.3f}  {r['P_total']:12.4e}  "
              f"{r['R_effective']:12.4e}  {r['L_effective']*1e9:10.4f}")
        results.append(r)

    return results


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(
        description="BEM coil + ESIM workpiece impedance extraction")
    parser.add_argument("--material", default="steel",
                        choices=["steel", "copper", "aluminum"])
    parser.add_argument("--freq", type=float, default=50000,
                        help="Frequency [Hz]")
    parser.add_argument("--sweep", action="store_true",
                        help="Run frequency sweep")
    parser.add_argument("--R", type=float, default=0.05,
                        help="Torus major radius [m]")
    parser.add_argument("--a", type=float, default=0.005,
                        help="Torus minor radius [m]")
    parser.add_argument("--wp-radius", type=float, default=0.01,
                        help="Workpiece radius [m]")
    parser.add_argument("--wp-height", type=float, default=0.02,
                        help="Workpiece height [m]")
    args = parser.parse_args()

    if args.sweep:
        frequency_sweep(material=args.material,
                        R=args.R, a=args.a,
                        wp_radius=args.wp_radius,
                        wp_height=args.wp_height)
    else:
        run(material=args.material, frequency=args.freq,
            R=args.R, a=args.a,
            wp_radius=args.wp_radius, wp_height=args.wp_height)
