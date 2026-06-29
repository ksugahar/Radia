"""
EFIE-SIBC: 3D BEM (ngsolve.bem) + surface impedance for induction heating.

Two-way coupled eddy current analysis with surface-only mesh.
Captures impedance screening through Z_s in the BEM system.

Formulation (EFIE + SIBC on workpiece surface):
  [Z_s * Mass + jw*mu0 * SL] * J = -jw * rhs(A_inc)

  SL: LaplaceSL on workpiece surface (self-inductance operator)
  Mass: HDivSurface mass matrix
  Z_s: surface impedance from ESIM cell problem (cylinder geometry)
  J: induced surface current on workpiece (unknown)
  A_inc: incident vector potential from coil Biot-Savart

Physics: surface PEEC equation
  [surface impedance + self/mutual inductance] * current = induced EMF

Karl iteration for nonlinear Z_s(H_t):
  1. Solve EFIE+SIBC system for J
  2. H_t = |J| -> area-weighted RMS
  3. Update Z_s from ESIM at H_t_rms
  4. Repeat until Z_s converges

Usage:
    python pmchwt_sibc.py
    python pmchwt_sibc.py --material steel --freq 7000
"""

import math
import os
import sys
import time
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../../src/radia'))

MU_0 = 4e-7 * np.pi


# ============================================================
# Compute incident A and H at workpiece surface from coil BEM solution
# ============================================================
def compute_incident_fields_at_surface(mesh_coil, gf_J, mesh_wp):
    """Compute incident A and H fields at workpiece element centroids via Biot-Savart.

    A(r) = (mu0/4pi) * sum_elem { J_elem / |r - r'| * dA }
    H(r) = (1/4pi) * sum_elem { J_elem x (r - r') / |r - r'|^3 * dA }

    Args:
        mesh_coil: NGSolve surface mesh of coil
        gf_J: GridFunction(HDivSurface) with solved coil current
        mesh_wp: NGSolve surface mesh of workpiece

    Returns:
        A_inc: (n_wp_elem, 3) vector potential [T*m]
        H_inc: (n_wp_elem, 3) magnetic field [A/m]
    """
    from ngsolve import Integrate, CF, BND

    INV_4PI = 1.0 / (4.0 * np.pi)

    # Extract per-element J and geometry from coil mesh
    elem_A = Integrate(CF(1), mesh_coil, VOL_or_BND=BND, element_wise=True)
    elem_Jx = Integrate(gf_J[0], mesh_coil, VOL_or_BND=BND, element_wise=True)
    elem_Jy = Integrate(gf_J[1], mesh_coil, VOL_or_BND=BND, element_wise=True)
    elem_Jz = Integrate(gf_J[2], mesh_coil, VOL_or_BND=BND, element_wise=True)

    src_centroids, src_areas, src_J = [], [], []
    for el in mesh_coil.Elements(BND):
        area = abs(elem_A[el.nr])
        if area < 1e-30:
            continue
        jvec = np.array([elem_Jx[el.nr], elem_Jy[el.nr], elem_Jz[el.nr]]) / area
        verts = [mesh_coil.vertices[v.nr].point for v in el.vertices]
        c = np.mean([(v[0], v[1], v[2]) for v in verts], axis=0)
        src_centroids.append(c)
        src_areas.append(area)
        src_J.append(jvec)

    src_centroids = np.array(src_centroids)  # (ne_coil, 3)
    src_areas = np.array(src_areas)          # (ne_coil,)
    src_J = np.array(src_J)                  # (ne_coil, 3)

    # Workpiece element centroids
    obs_pts = []
    for el in mesh_wp.Elements(BND):
        verts = [mesh_wp.vertices[v.nr].point for v in el.vertices]
        c = np.mean([(v[0], v[1], v[2]) for v in verts], axis=0)
        obs_pts.append(c)
    obs_pts = np.array(obs_pts)  # (ne_wp, 3)

    # Vectorized Biot-Savart
    dx = obs_pts[:, None, :] - src_centroids[None, :, :]   # (nw, nc, 3)
    r = np.sqrt(np.maximum(np.sum(dx**2, axis=2), 1e-30))  # (nw, nc)

    # A = (mu0/4pi) * sum { J / |r| * dA }
    r_inv_area = src_areas[None, :] / r                     # (nw, nc)
    A_inc = (MU_0 * INV_4PI) * np.sum(
        src_J[None, :, :] * r_inv_area[:, :, None], axis=1)  # (nw, 3)

    # H = (1/4pi) * sum { J x dr / |r|^3 * dA }
    r3_inv = src_areas[None, :] / (r ** 3)                  # (nw, nc)
    cross = np.cross(src_J[None, :, :], dx)                  # (nw, nc, 3)
    H_inc = INV_4PI * np.sum(cross * r3_inv[:, :, None], axis=1)  # (nw, 3)

    return A_inc, H_inc


# ============================================================
# EFIE-SIBC solver using ngsolve.bem
# ============================================================
def _extract_dense_matrix(operator, ndof, fes):
    """Extract BEM or BilinearForm operator to dense NumPy matrix."""
    mat = operator if not hasattr(operator, 'mat') else operator.mat
    try:
        rows, cols, vals = mat.COO()
        M = np.zeros((ndof, ndof), dtype=complex)
        for r, c, v in zip(rows, cols, vals):
            M[int(r), int(c)] = v
        return M
    except Exception:
        pass
    # Fallback: column-by-column
    from ngsolve import GridFunction
    M = np.zeros((ndof, ndof), dtype=complex)
    for i in range(ndof):
        ei = GridFunction(fes)
        ei.vec[:] = 0
        ei.vec[i] = 1.0
        res = ei.vec.CreateVector()
        res.data = mat * ei.vec
        M[:, i] = res.FV().NumPy().copy()
    return M


def solve_pmchwt_sibc_3d(mesh_wp, mesh_coil, gf_J_coil,
                          R_wp, sigma, frequency,
                          bh_curve=None, mu_r=1.0,
                          max_iter=15, tol=1e-3, relax=0.5,
                          esim_geometry='cylinder'):
    """3D EFIE-SIBC solve using surface PEEC formulation.

    System: [Z_s*Mass + jw*mu0*SL] * J = -jw * rhs(A_inc)

    Args:
        mesh_wp: NGSolve surface mesh of workpiece
        mesh_coil: NGSolve surface mesh of coil
        gf_J_coil: GridFunction(HDivSurface) coil current (from BEM)
        R_wp: workpiece radius [m] (for ESIM cylinder) or half-thickness (slab)
        sigma: conductivity [S/m]
        frequency: [Hz]
        bh_curve: BH curve for nonlinear material
        mu_r: relative permeability (linear)

    Returns:
        dict with P_total, H_t, Z_s, etc.
    """
    from ngsolve import (HDivSurface, BilinearForm, LinearForm, GridFunction,
                         Integrate, CF, ds, BND, TaskManager, SurfaceL2, div)
    from ngsolve.bem import LaplaceSL
    from esim_cell_problem import ESIMFiniteSlabSolver

    omega = 2 * np.pi * frequency

    # FE space on workpiece surface
    fes_wp = HDivSurface(mesh_wp, order=0)
    fes_l2 = SurfaceL2(mesh_wp, order=0)
    u, v = fes_wp.TnT()
    ndof = fes_wp.ndof
    n_p = fes_l2.ndof
    print(f"  Workpiece DOFs: {ndof} (J) + {n_p} (p)")

    # --- Assemble SL = LaplaceSL on workpiece (self-inductance) ---
    print("  Assembling LaplaceSL on workpiece...")
    t0 = time.perf_counter()
    with TaskManager():
        SL_bf = LaplaceSL(u.Trace() * ds, use_fmm=False) * v.Trace() * ds
    SL_mat = _extract_dense_matrix(SL_bf, ndof, fes_wp)
    t_sl = time.perf_counter() - t0
    print(f"  SL assembled: {ndof}x{ndof} ({t_sl:.1f}s)")

    # --- Mass matrix ---
    mass_bf = BilinearForm(fes_wp)
    mass_bf += u.Trace() * v.Trace() * ds
    mass_bf.Assemble()
    Mass_mat = _extract_dense_matrix(mass_bf, ndof, fes_wp)

    # --- Divergence matrix D: n_p x ndof (enforces div J = 0) ---
    u_J = fes_wp.TrialFunction()
    q = fes_l2.TestFunction()
    bf_D = BilinearForm(trialspace=fes_wp, testspace=fes_l2)
    bf_D += div(u_J.Trace()) * q * ds
    bf_D.Assemble()
    D_mat = np.zeros((n_p, ndof), dtype=float)
    try:
        rows, cols, vals = bf_D.mat.COO()
        for r, c, val in zip(rows, cols, vals):
            D_mat[int(r), int(c)] = float(val)
    except Exception:
        for i in range(ndof):
            ei = GridFunction(fes_wp)
            ei.vec[:] = 0
            ei.vec[i] = 1.0
            res = ei.vec.CreateVector()
            res.data = bf_D.mat * ei.vec
            D_mat[:, i] = res.FV().NumPy()[:n_p].copy()
    rank_D = np.linalg.matrix_rank(D_mat)
    # Closed surface: one redundant constraint (sum div = 0), remove last row
    D_red = D_mat[:rank_D, :]
    n_p_red = rank_D
    print(f"  D (divergence): {n_p}x{ndof}, rank={rank_D} -> {n_p_red} constraints")

    # --- RHS: -jw * A_inc from coil Biot-Savart ---
    print("  Computing A_inc, H_inc via Biot-Savart...")
    A_inc, H_inc = compute_incident_fields_at_surface(mesh_coil, gf_J_coil, mesh_wp)
    H_inc_mag = np.linalg.norm(H_inc, axis=1)
    A_inc_mag = np.linalg.norm(A_inc, axis=1)
    print(f"  |H_inc| range: {H_inc_mag.min():.1f} - {H_inc_mag.max():.1f} A/m")
    print(f"  |A_inc| range: {A_inc_mag.min():.2e} - {A_inc_mag.max():.2e} T*m")

    # Create piecewise constant CF for A_inc
    ne = mesh_wp.GetNE(BND)
    gf_Ax = GridFunction(fes_l2)
    gf_Ay = GridFunction(fes_l2)
    gf_Az = GridFunction(fes_l2)
    gf_Ax.vec.FV().NumPy()[:ne] = A_inc[:, 0]
    gf_Ay.vec.FV().NumPy()[:ne] = A_inc[:, 1]
    gf_Az.vec.FV().NumPy()[:ne] = A_inc[:, 2]
    A_inc_cf = CF((gf_Ax, gf_Ay, gf_Az))

    # Project A_inc onto HDivSurface, multiply by -jw
    lf_rhs = LinearForm(fes_wp)
    lf_rhs += A_inc_cf * v.Trace() * ds
    lf_rhs.Assemble()
    f_rhs = -1j * omega * lf_rhs.vec.FV().NumPy().copy().astype(complex)
    print(f"  |rhs| = {np.linalg.norm(f_rhs):.4e}")

    # --- ESIM solver (per-element) ---
    esim = ESIMFiniteSlabSolver(
        half_thickness=R_wp, bh_curve=bh_curve, sigma=sigma,
        frequency=frequency,
        mu_r=mu_r if bh_curve is None else None,
        n_nodes=200, geometry=esim_geometry)

    # Element areas (constant across iterations)
    elem_areas = Integrate(CF(1), mesh_wp, VOL_or_BND=BND, element_wise=True)
    area_arr = np.array([abs(elem_areas[el.nr]) for el in mesh_wp.Elements(BND)])

    # Initialize per-element Z_s with uniform value
    Z_s_elem = np.full(ne, esim.solve(5.0)['Z'], dtype=complex)

    # Helper: assemble Z_s-weighted mass matrix
    def assemble_Zs_mass(Z_s_arr):
        gf_re = GridFunction(fes_l2)
        gf_im = GridFunction(fes_l2)
        gf_re.vec.FV().NumPy()[:ne] = Z_s_arr.real
        gf_im.vec.FV().NumPy()[:ne] = Z_s_arr.imag
        bf_re = BilinearForm(fes_wp)
        bf_re += CF(gf_re) * u.Trace() * v.Trace() * ds
        bf_re.Assemble()
        bf_im = BilinearForm(fes_wp)
        bf_im += CF(gf_im) * u.Trace() * v.Trace() * ds
        bf_im.Assemble()
        return (_extract_dense_matrix(bf_re, ndof, fes_wp) +
                1j * _extract_dense_matrix(bf_im, ndof, fes_wp))

    # Helper: per-element H_t from complex J_vec
    def compute_per_element_Ht(J_vec):
        gf_re = GridFunction(fes_wp)
        gf_im = GridFunction(fes_wp)
        gf_re.vec.FV().NumPy()[:] = J_vec.real
        gf_im.vec.FV().NumPy()[:] = J_vec.imag
        Jsq_re = Integrate(gf_re * gf_re, mesh_wp, BND, element_wise=True)
        Jsq_im = Integrate(gf_im * gf_im, mesh_wp, BND, element_wise=True)
        H_t = np.zeros(ne)
        for el in mesh_wp.Elements(BND):
            a = area_arr[el.nr]
            if a < 1e-30:
                continue
            H_t[el.nr] = math.sqrt(
                (abs(Jsq_re[el.nr]) + abs(Jsq_im[el.nr])) / a)
        return H_t

    # --- Karl iteration (per-element Z_s) ---
    Z_s_mean = np.mean(np.abs(Z_s_elem))
    print(f"  |Z_s| init (mean) = {Z_s_mean:.4e}")
    print(f"  {'Iter':>4s} {'|Z_s| mean':>12s} {'P [W]':>12s} "
          f"{'H_t_rms':>10s} {'dZ/Z max':>10s}")

    history = []

    # Pre-build saddle point template
    n_total = ndof + n_p_red
    rhs_block = np.zeros(n_total, dtype=complex)
    rhs_block[:ndof] = f_rhs
    # rhs_block[ndof:] = 0 (divergence-free constraint)

    for iteration in range(max_iter):
        # Saddle point: [Z_s*M + jw*mu0*SL,  D^T] [J] = [f_rhs]
        #               [D,                   0  ] [p]   [0    ]
        ZsMass = assemble_Zs_mass(Z_s_elem)
        A_block = np.zeros((n_total, n_total), dtype=complex)
        A_block[:ndof, :ndof] = ZsMass + 1j * omega * MU_0 * SL_mat
        A_block[:ndof, ndof:] = D_red.T
        A_block[ndof:, :ndof] = D_red
        J_p = np.linalg.solve(A_block, rhs_block)
        J_vec = J_p[:ndof]

        # Per-element H_t
        H_t_elem = compute_per_element_Ht(J_vec)
        H_t_rms = math.sqrt(np.sum(H_t_elem**2 * area_arr) /
                            max(np.sum(area_arr), 1e-30))

        # Per-element ESIM: Z_s and P from cell problem
        Z_s_old = Z_s_elem.copy()
        P_elem = np.zeros(ne)
        for e in range(ne):
            H0 = max(float(H_t_elem[e]), 1e-3)
            sol = esim.solve(H0)
            Z_new = sol['Z']
            Z_s_elem[e] = relax * Z_new + (1 - relax) * Z_s_old[e]
            P_elem[e] = sol['P_prime'] * area_arr[e]

        P_total = float(np.sum(P_elem))
        # Q from last iteration's Z_s
        Q_total = float(np.sum(
            [esim.solve(max(float(H_t_elem[e]), 1e-3))['Q_prime'] * area_arr[e]
             for e in range(ne)]))

        # Convergence: max relative change in |Z_s|
        dZ_arr = np.abs(Z_s_elem - Z_s_old) / np.maximum(np.abs(Z_s_old), 1e-30)
        dZ_max = float(np.max(dZ_arr))
        Z_s_mean = float(np.mean(np.abs(Z_s_elem)))

        history.append({'Z_s_mean': Z_s_mean, 'P': P_total,
                        'H_t_rms': float(H_t_rms), 'dZ_max': dZ_max})
        print(f"  {iteration:4d} {Z_s_mean:12.4e} {P_total:12.4e} "
              f"{H_t_rms:10.2f} {dZ_max:10.4e}")

        if dZ_max < tol and iteration > 0:
            print(f"  Converged at iteration {iteration}")
            break

    A_total = float(np.sum(area_arr))
    return {
        'P_total': P_total,
        'Q_total': Q_total,
        'H_t_rms': float(H_t_rms),
        'H_t_elem': H_t_elem,
        'Z_s_elem': Z_s_elem,
        'P_elem': P_elem,
        'J_sol': J_vec,
        'ndof': ndof,
        'ne': ne,
        'area': A_total,
        'iterations': len(history),
        'history': history,
    }


# ============================================================
# Main
# ============================================================
def run(material='copper', frequency=1000, R_coil=0.030, a_coil=0.003,
        wp_radius=0.010, wp_height=0.020, verbose=True,
        esim_geometry='cylinder'):
    """Run EFIE-SIBC with OCC workpiece mesh."""
    from ngsolve import Mesh
    from netgen.occ import (WorkPlane, Axes, Pnt, Dir, Cylinder,
                             OCCGeometry, Glue)
    from bem_inductance import compute_inductance_source_sink

    STEEL_BH = [
        [0, 0], [50, 0.1], [100, 0.25], [200, 0.55],
        [500, 0.95], [1000, 1.2], [2000, 1.4], [5000, 1.55],
        [10000, 1.65], [20000, 1.75], [50000, 1.9], [100000, 2.0],
    ]

    MATERIALS = {
        'steel': (STEEL_BH, 2e6, None),
        'copper': (None, 5.8e7, 1.0),
        'aluminum': (None, 3.5e7, 1.0),
    }
    bh_curve, sigma, mu_r = MATERIALS[material]

    rho = 1.0 / sigma
    mu_eff = MU_0 * (mu_r if mu_r else 1.0)
    delta = math.sqrt(2 * rho / (2 * np.pi * frequency * mu_eff))

    if verbose:
        print("=" * 65)
        print("EFIE-SIBC: 3D BEM + surface impedance")
        print("=" * 65)
        print(f"Material: {material}, f={frequency}Hz, "
              f"delta={delta*1e3:.3f}mm, xi={wp_radius/delta:.1f}")
        print()

    # Step 1: Coil BEM
    if verbose:
        print("[1/3] Coil BEM solve...")
    from impedance_esim import make_gapped_torus_mesh
    mesh_coil = make_gapped_torus_mesh(R_coil, a_coil, gap_deg=5, maxh=a_coil)
    sol_coil = compute_inductance_source_sink(mesh_coil)
    L_coil = sol_coil['L']
    gf_J = sol_coil['gf_J']
    if verbose:
        print(f"  L_coil = {L_coil*1e9:.2f} nH")

    # Step 2: Workpiece surface mesh (OCC cylinder)
    if verbose:
        print("[2/3] Workpiece surface mesh...")
    cyl = Cylinder(Pnt(0, 0, -wp_height/2), Dir(0, 0, 1),
                   wp_radius, wp_height)
    for f in cyl.faces:
        f.name = "wp"
    wp_surf = Glue(cyl.faces)
    geo_wp = OCCGeometry(wp_surf)
    ngmesh_wp = geo_wp.GenerateMesh(maxh=wp_radius / 2)
    mesh_wp = Mesh(ngmesh_wp)
    mesh_wp.Curve(2)
    from ngsolve import BND as BND_
    nse = mesh_wp.GetNE(BND_)
    if verbose:
        print(f"  {nse} surface elements")

    # Step 3: EFIE-SIBC solve
    if verbose:
        print("[3/3] EFIE-SIBC solve...")
    result = solve_pmchwt_sibc_3d(
        mesh_wp, mesh_coil, gf_J,
        R_wp=wp_radius, sigma=sigma, frequency=frequency,
        bh_curve=bh_curve, mu_r=mu_r if mu_r else 1.0,
        esim_geometry=esim_geometry)

    result['L_coil'] = float(L_coil)

    if verbose:
        Z_mean = np.mean(np.abs(result['Z_s_elem']))
        jw_mu0_a = 2 * np.pi * frequency * MU_0 * wp_radius
        print()
        print("=" * 65)
        print(f"  P (EFIE-SIBC)   = {result['P_total']:.6e} W")
        print(f"  Q               = {result['Q_total']:.6e} var")
        print(f"  H_t_rms         = {result['H_t_rms']:.2f} A/m")
        print(f"  |Z_s| mean      = {Z_mean:.4e} Ohm")
        print(f"  L_coil          = {L_coil*1e9:.2f} nH")
        print(f"  Screening: |Z_s|/jw*mu0*a = {Z_mean/jw_mu0_a:.1f} "
              f"({'strong' if Z_mean/jw_mu0_a > 3 else 'weak'})")
        print("=" * 65)

    return result


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="EFIE-SIBC 3D BEM")
    parser.add_argument("--material", default="copper",
                        choices=["copper", "steel", "aluminum"])
    parser.add_argument("--freq", type=float, default=1000)
    parser.add_argument("--geometry", default="local_curvature",
                        choices=["local_curvature", "none"],
                        help="ESIM curvature: local_curvature (Bessel) or none (flat slab)")
    args = parser.parse_args()

    _geom_map = {"local_curvature": "cylinder", "none": "slab"}
    esim_geom = _geom_map.get(args.geometry, "cylinder")
    r = run(material=args.material, frequency=args.freq,
            esim_geometry=esim_geom)

    # Compare with BEM-ESIM one-way
    print()
    print("BEM-ESIM (one-way) for comparison...")
    from impedance_esim import run as bem_run
    b = bem_run(material=args.material, frequency=args.freq,
                R=0.030, a=0.003, wp_radius=0.010, wp_height=0.020,
                n_phi=24, n_z=24, verbose=False)

    dp = (r['P_total'] - b['P_total']) / b['P_total'] * 100
    print(f"{'':>22s} {'EFIE-SIBC':>14s} {'BEM one-way':>14s} {'diff':>8s}")
    print(f"{'P [W]':>22s} {r['P_total']:14.4e} {b['P_total']:14.4e} {dp:+8.1f}%")
    Z_mean = np.mean(np.abs(r['Z_s_elem']))
    jw_mu0_a = 2 * np.pi * args.freq * MU_0 * 0.010
    print(f"\n  Screening physics: |Z_s|={Z_mean:.2e}, "
          f"jw*mu0*a={jw_mu0_a:.2e}, ratio={Z_mean/jw_mu0_a:.1f}")
    if Z_mean / jw_mu0_a > 10:
        print("  -> Z_s >> jw*mu0*a: strong screening, one-way overestimates")
    elif Z_mean / jw_mu0_a < 0.1:
        print("  -> Z_s << jw*mu0*a: weak screening, one-way ~ EFIE-SIBC")
