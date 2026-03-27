"""
PMCHWT-SIBC: 3D BEM (ngsolve.bem) + surface impedance for induction heating.

Two-way coupled eddy current analysis with surface-only mesh.
Captures impedance screening through Z_s in the BEM system.

Formulation (EFIE + SIBC on workpiece surface):
  [jw*mu0*SL_wp + Z_s*Mass_wp] * J_wp = -jw*mu0*SL_coil_to_wp * J_coil

  SL_wp: LaplaceSL self-interaction on workpiece surface
  Mass_wp: HDivSurface mass matrix (identity for tangential J)
  SL_coil_to_wp: Biot-Savart A from coil at workpiece surface
  Z_s: surface impedance from ESIM cell problem (cylinder geometry)
  J_wp: induced surface current on workpiece (unknown)

Karl iteration for nonlinear Z_s(H_t):
  1. Solve BEM system for J_wp
  2. H_t = |J_wp| -> area-weighted RMS
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
# Compute incident A at workpiece from coil BEM solution
# ============================================================
def compute_A_inc_at_panels(mesh_coil, gf_J, obs_pts):
    """Compute incident vector potential A at observation points.

    A(r) = (mu0/4pi) * sum_elem { J_elem * dA_elem / |r - r_elem| }

    This is simpler than Biot-Savart for H (no cross product, just 1/r).

    Args:
        mesh_coil: NGSolve surface mesh of coil
        gf_J: GridFunction(HDivSurface) with solved coil current
        obs_pts: (N, 3) array of observation points

    Returns:
        A_inc: (N, 3) array of vector potential [T*m]
    """
    from ngsolve import Integrate, CF, BND

    INV_4PI = 1.0 / (4.0 * np.pi)

    elem_A = Integrate(CF(1), mesh_coil, VOL_or_BND=BND, element_wise=True)
    elem_Jx = Integrate(gf_J[0], mesh_coil, VOL_or_BND=BND, element_wise=True)
    elem_Jy = Integrate(gf_J[1], mesh_coil, VOL_or_BND=BND, element_wise=True)
    elem_Jz = Integrate(gf_J[2], mesh_coil, VOL_or_BND=BND, element_wise=True)

    centroids, J_dA = [], []
    for el in mesh_coil.Elements(BND):
        area = abs(elem_A[el.nr])
        if area < 1e-30:
            continue
        jvec = np.array([elem_Jx[el.nr], elem_Jy[el.nr], elem_Jz[el.nr]])
        verts = [mesh_coil.vertices[v.nr].point for v in el.vertices]
        c = np.mean([(v[0], v[1], v[2]) for v in verts], axis=0)
        centroids.append(c)
        J_dA.append(jvec)  # J * area (integrated)

    centroids = np.array(centroids)  # (ne, 3)
    J_dA = np.array(J_dA)            # (ne, 3)

    # Vectorized: A = (mu0/4pi) * sum { J*dA / |r - r'| }
    dx = obs_pts[:, None, :] - centroids[None, :, :]   # (np, ne, 3)
    r = np.sqrt(np.maximum(np.sum(dx**2, axis=2), 1e-30))  # (np, ne)
    r_inv = 1.0 / r                                        # (np, ne)

    # A = mu0/(4pi) * sum_j { J_j * dA_j / r_ij }
    A_inc = MU_0 * INV_4PI * np.einsum('ij,jk->ik', r_inv, J_dA)  # (np, 3)

    return A_inc


# ============================================================
# PMCHWT-SIBC solver using ngsolve.bem
# ============================================================
def solve_pmchwt_sibc_3d(mesh_wp, mesh_coil, gf_J_coil,
                          R_wp, sigma, frequency,
                          bh_curve=None, mu_r=1.0,
                          max_iter=15, tol=1e-3, relax=0.5):
    """3D PMCHWT-SIBC solve using ngsolve.bem operators.

    System: [jw*mu0*SL + Z_s*Mass] * J = RHS

    Args:
        mesh_wp: NGSolve surface mesh of workpiece
        mesh_coil: NGSolve surface mesh of coil
        gf_J_coil: GridFunction(HDivSurface) coil current (from BEM)
        R_wp: workpiece radius [m] (for ESIM cylinder)
        sigma: conductivity [S/m]
        frequency: [Hz]
        bh_curve: BH curve for nonlinear material
        mu_r: relative permeability (linear)

    Returns:
        dict with P_total, H_t, Z_s, etc.
    """
    from ngsolve import (HDivSurface, BilinearForm, LinearForm, GridFunction,
                         Integrate, CF, ds, BND, TaskManager)
    from ngsolve.bem import LaplaceSL
    from esim_cell_problem import ESIMFiniteSlabSolver

    omega = 2 * np.pi * frequency

    # FE space on workpiece surface
    fes_wp = HDivSurface(mesh_wp, order=0)
    u, v = fes_wp.TnT()
    ndof = fes_wp.ndof
    print(f"  Workpiece DOFs: {ndof}")

    # --- Assemble SL operator (workpiece self-interaction) ---
    print("  Assembling LaplaceSL on workpiece...")
    t0 = time.perf_counter()
    with TaskManager():
        SL_bf = LaplaceSL(u.Trace() * ds) * v.Trace() * ds
    # Extract dense matrix via COO (fast path)
    try:
        coo = SL_bf.mat.COO()
        rows, cols, vals = coo
        SL_mat = np.zeros((ndof, ndof), dtype=complex)
        for r, c, v in zip(rows, cols, vals):
            SL_mat[int(r), int(c)] = v
    except:
        # Fallback: column-by-column
        SL_mat = np.zeros((ndof, ndof), dtype=complex)
        for i in range(ndof):
            ei = GridFunction(fes_wp)
            ei.vec[:] = 0
            ei.vec[i] = 1.0
            res = ei.vec.CreateVector()
            res.data = SL_bf.mat * ei.vec
            SL_mat[:, i] = res.FV().NumPy().copy()
    t_sl = time.perf_counter() - t0
    print(f"  SL assembled: {ndof}x{ndof} ({t_sl:.1f}s)")

    # --- Mass matrix ---
    mass_bf = BilinearForm(fes_wp)
    mass_bf += u.Trace() * v.Trace() * ds
    mass_bf.Assemble()
    Mass_mat = np.zeros((ndof, ndof), dtype=complex)
    for i in range(ndof):
        ei = GridFunction(fes_wp)
        ei.vec[:] = 0
        ei.vec[i] = 1.0
        res = ei.vec.CreateVector()
        res.data = mass_bf.mat * ei.vec
        Mass_mat[:, i] = res.FV().NumPy().copy()

    # --- RHS: incident A from coil ---
    # Compute A_inc at workpiece element centroids, project onto HDivSurface
    print("  Computing incident A from coil...")
    # Get workpiece element centroids
    wp_centroids = []
    for el in mesh_wp.Elements(BND):
        verts = [mesh_wp.vertices[v.nr].point for v in el.vertices]
        c = np.mean([(v[0], v[1], v[2]) for v in verts], axis=0)
        wp_centroids.append(c)
    wp_centroids = np.array(wp_centroids)

    A_inc = compute_A_inc_at_panels(mesh_coil, gf_J_coil, wp_centroids)

    # Project A_inc onto HDivSurface: f_i = -jw*mu0 * integral A_inc . v_i ds
    # Approximate: use element-wise integration
    f_rhs = np.zeros(ndof, dtype=complex)
    elem_areas = Integrate(CF(1), mesh_wp, VOL_or_BND=BND, element_wise=True)
    for comp in range(3):
        gf_Ac = GridFunction(fes_wp)
        # Set per-element A_inc component
        lf = LinearForm(fes_wp)
        # Use CF with per-element values
        # Simplified: evaluate A_inc at quadrature points
        # For now, use element centroid value * area as approximation
        pass

    # Simpler RHS: direct matrix-vector with A_inc projected per-DOF
    # Actually, for RWG basis, the DOFs are edge-based. A_inc needs proper projection.
    # Simplest: compute f = -jw * Mass * A_inc_projected
    # where A_inc_projected has per-DOF values from A_inc at edge midpoints.

    # For order=0 HDivSurface (RWG), DOFs correspond to mesh edges.
    # Each RWG basis function spans two triangles sharing an edge.
    # The projection: f_i = -jw * integral A_inc . v_i ds

    # Proper way: create a CoefficientFunction from A_inc and assemble LinearForm
    # But A_inc is not a smooth CF... use per-element constant approximation.

    # Per-element A_inc as a piecewise constant vector CF
    ne = mesh_wp.GetNE(BND)
    A_inc_x = np.zeros(ne)
    A_inc_y = np.zeros(ne)
    A_inc_z = np.zeros(ne)
    for el in mesh_wp.Elements(BND):
        A_inc_x[el.nr] = A_inc[el.nr, 0]
        A_inc_y[el.nr] = A_inc[el.nr, 1]
        A_inc_z[el.nr] = A_inc[el.nr, 2]

    from ngsolve import SurfaceL2
    fes_l2 = SurfaceL2(mesh_wp, order=0)
    gf_Ax = GridFunction(fes_l2)
    gf_Ay = GridFunction(fes_l2)
    gf_Az = GridFunction(fes_l2)
    gf_Ax.vec.FV().NumPy()[:ne] = A_inc_x
    gf_Ay.vec.FV().NumPy()[:ne] = A_inc_y
    gf_Az.vec.FV().NumPy()[:ne] = A_inc_z
    A_inc_cf = CF((gf_Ax, gf_Ay, gf_Az))

    # Assemble real part, then multiply by -jw*mu0
    lf_inc_real = LinearForm(fes_wp)
    lf_inc_real += A_inc_cf * v.Trace() * ds
    lf_inc_real.Assemble()
    f_rhs = -1j * omega * MU_0 * lf_inc_real.vec.FV().NumPy().copy().astype(complex)

    # --- ESIM solver ---
    esim = ESIMFiniteSlabSolver(
        half_thickness=R_wp, bh_curve=bh_curve, sigma=sigma,
        frequency=frequency,
        mu_r=mu_r if bh_curve is None else None,
        n_nodes=200, geometry='cylinder')
    Z_s = esim.solve(5.0)['Z']

    # --- Karl iteration ---
    print(f"  Z_s init = {Z_s:.4e}")
    print(f"  {'Iter':>4s} {'|Z_s|':>12s} {'P [W]':>12s} "
          f"{'H_t_rms':>10s} {'dZ/Z':>10s}")

    history = []
    J_sol = None

    for iteration in range(max_iter):
        # System: [jw*mu0*SL + Z_s*Mass] * J = f_rhs
        A_sys = 1j * omega * MU_0 * SL_mat + Z_s * Mass_mat
        J_vec = np.linalg.solve(A_sys, f_rhs)

        # Evaluate |J| per element for H_t
        gf_J_wp = GridFunction(fes_wp)
        gf_J_wp.vec.FV().NumPy()[:] = J_vec.real
        # Compute |J|^2 per element
        J_sq = Integrate(gf_J_wp * gf_J_wp, mesh_wp, BND, element_wise=True)
        areas = Integrate(CF(1), mesh_wp, BND, element_wise=True)

        # For complex J: need |J|^2 = J_real^2 + J_imag^2
        gf_J_re = GridFunction(fes_wp)
        gf_J_im = GridFunction(fes_wp)
        gf_J_re.vec.FV().NumPy()[:] = J_vec.real
        gf_J_im.vec.FV().NumPy()[:] = J_vec.imag
        J_sq_re = Integrate(gf_J_re * gf_J_re, mesh_wp, BND, element_wise=True)
        J_sq_im = Integrate(gf_J_im * gf_J_im, mesh_wp, BND, element_wise=True)

        H_t_sq_total = 0.0
        A_total = 0.0
        for el in mesh_wp.Elements(BND):
            a = abs(areas[el.nr])
            if a < 1e-30:
                continue
            h2 = (abs(J_sq_re[el.nr]) + abs(J_sq_im[el.nr])) / a
            H_t_sq_total += h2 * a
            A_total += a

        H_t_rms = math.sqrt(max(H_t_sq_total / max(A_total, 1e-30), 0))

        # Power
        P = 0.5 * Z_s.real * H_t_sq_total

        # Update Z_s
        Z_s_old = Z_s
        sol_new = esim.solve(max(float(H_t_rms), 1e-3))
        Z_s = relax * sol_new['Z'] + (1 - relax) * Z_s_old
        dZ = abs(Z_s - Z_s_old) / max(abs(Z_s_old), 1e-30)

        history.append({'Z_s': complex(Z_s), 'P': float(P),
                        'H_t_rms': float(H_t_rms), 'dZ': float(dZ)})
        print(f"  {iteration:4d} {abs(Z_s):12.4e} {P:12.4e} "
              f"{H_t_rms:10.2f} {dZ:10.4e}")

        if dZ < tol and iteration > 0:
            print(f"  Converged at iteration {iteration}")
            break

    J_sol = J_vec

    return {
        'P_total': float(P),
        'H_t_rms': float(H_t_rms),
        'Z_s': complex(Z_s),
        'J_sol': J_sol,
        'ndof': ndof,
        'iterations': len(history),
        'history': history,
    }


# ============================================================
# Main
# ============================================================
def run(material='copper', frequency=1000, R_coil=0.030, a_coil=0.003,
        wp_radius=0.010, wp_height=0.020, verbose=True):
    """Run PMCHWT-SIBC with OCC workpiece mesh."""
    from ngsolve import Mesh
    from netgen.occ import (WorkPlane, Axes, Pnt, Dir, Cylinder,
                             OCCGeometry, Glue)
    from bem_inductance import compute_inductance_source_sink

    STEEL_BH = [
        [0, 0], [50, 0.1], [100, 0.25], [200, 0.55],
        [500, 0.95], [1000, 1.2], [2000, 1.4], [5000, 1.55],
        [10000, 1.65], [20000, 1.75], [50000, 1.9], [100000, 2.0],
    ]

    if material == 'steel':
        bh_curve, sigma, mu_r = STEEL_BH, 2e6, None
    else:
        bh_curve, sigma, mu_r = None, 5.8e7, 1.0

    rho = 1.0 / sigma
    mu_eff = MU_0 * (mu_r if mu_r else 1.0)
    delta = math.sqrt(2 * rho / (2 * np.pi * frequency * mu_eff))

    if verbose:
        print("=" * 65)
        print("PMCHWT-SIBC: 3D BEM + surface impedance")
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

    # Step 3: PMCHWT-SIBC solve
    if verbose:
        print("[3/3] PMCHWT-SIBC solve...")
    result = solve_pmchwt_sibc_3d(
        mesh_wp, mesh_coil, gf_J,
        R_wp=wp_radius, sigma=sigma, frequency=frequency,
        bh_curve=bh_curve, mu_r=mu_r if mu_r else 1.0)

    result['L_coil'] = float(L_coil)

    if verbose:
        print()
        print("=" * 65)
        print(f"  P (PMCHWT-SIBC) = {result['P_total']:.6e} W")
        print(f"  L_coil          = {L_coil*1e9:.2f} nH")
        print("=" * 65)

    return result


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="PMCHWT-SIBC 3D BEM")
    parser.add_argument("--material", default="copper",
                        choices=["copper", "steel", "aluminum"])
    parser.add_argument("--freq", type=float, default=1000)
    args = parser.parse_args()

    r = run(material=args.material, frequency=args.freq)

    # Compare with BEM-ESIM one-way
    print()
    print("BEM-ESIM (one-way) for comparison...")
    from impedance_esim import run as bem_run
    b = bem_run(material=args.material, frequency=args.freq,
                R=0.030, a=0.003, wp_radius=0.010, wp_height=0.020,
                n_phi=24, n_z=24, verbose=False)

    dp = (r['P_total'] - b['P_total']) / b['P_total'] * 100
    print(f"{'':>22s} {'PMCHWT-SIBC':>14s} {'BEM one-way':>14s} {'diff':>8s}")
    print(f"{'P [W]':>22s} {r['P_total']:14.4e} {b['P_total']:14.4e} {dp:+8.1f}%")
